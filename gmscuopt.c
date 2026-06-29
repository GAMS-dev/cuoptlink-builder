#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>
#include "gmomcc.h"
#include "gevmcc.h"
#include "optcc.h"
#include <cuopt/linear_programming/cuopt_c.h>

int 
printOut (gevHandle_t gev, char *fmt, ...)
{
  va_list argp;
  int rc = 0;
  char msg[256];

  va_start (argp, fmt);
  rc = vsnprintf(msg, sizeof(msg), fmt, argp);
  va_end(argp);
  gevLogStatPChar(gev, msg);
  return rc;
}

static char fln_mip_trace[256];
static char mip_trace_id[32] = "";
static FILE *fp_mip_trace = NULL;
static int mip_trace_seq = 0;

static int mip_trace_open(const char *fname, const char *solverID, const int optFileNum, const char *inputName);
static int mip_trace_close();
static int mip_trace_line(char seriesID, double node, int giveint, double seconds, double bestint, double bestbnd);

typedef struct sl_state_s
{
  gevHandle_t gev;
  double tstart;
  int nvars;
} sl_state_t;

static void mip_get_solution_cb(const cuopt_float_t *solution, const cuopt_float_t *objective_value,
                                const cuopt_float_t *solution_bound, void *user_data);

int main(int argc, char *argv[])
{
  gmoHandle_t gmo=NULL;
  gevHandle_t gev=NULL;
  optHandle_t opt=NULL;
  cuopt_int_t status;
  char msg[256], filename[256];

  if (!gevCreate(&gev,msg,sizeof(msg))) {
    printf("Error creating GEV: %s\n", msg);
    goto GAMSDONE;
  }
  if (!gmoCreate(&gmo,msg,sizeof(msg))) {
    printf("Error creating GMO: %s\n", msg);
    goto GAMSDONE;
  }
  if (!optCreate(&opt,msg,sizeof(msg))) {
    printf("Error creating OPT: %s\n", msg);
    goto GAMSDONE;
  }

  status = gevInitEnvironmentLegacy(gev,argv[1]);
  if (status) {
    printf("Could not initialize GEV: %d\n", status);
    goto GAMSDONE;
  }

#if defined(CUOPT_VERSION) && defined(CUOPT_HASH)
  printOut(gev, "GAMS/cuOpt link was built against cuOpt version: %s, git hash: %s\n", CUOPT_VERSION, CUOPT_HASH);
#endif

  status = gmoRegisterEnvironment(gmo, gev, msg);
  if (status) {
    printOut(gev, "Could not register GEV: %s\n", msg);
    goto GAMSDONE;
  }

  status = gmoLoadDataLegacy(gmo, msg);
  if (status) {
    printOut(gev, "Could not register GEV: %s\n", msg);
    goto GAMSDONE;
  }

  gevGetStrOpt(gev, gevNameSysDir, filename);
  strcat(filename, "optcuopt.def");
  if (optReadDefinition(opt, filename)) {
    for (int i=1; i<=optMessageCount(opt); i++) {
      int msg_type=0;
      optGetMessage(opt, i, msg, &msg_type);
      printOut(gev, "%s\n",msg);
    }
    optClearMessages(opt);
    goto GAMSDONE;
  }
  for (int i=1; i<=optMessageCount(opt); i++) {
    int msg_type=0;
    optGetMessage(opt, i, msg, &msg_type);
    printOut(gev, "%s\n",msg);
  }
  optClearMessages(opt);

  if (gmoOptFile(gmo)) {
    gevStatCon(gev);
    optEchoSet(opt, 1);
    optReadParameterFile(opt, gmoNameOptFile(gmo, msg));
    for (int i=1; i<=optMessageCount(opt); i++) {
      int msg_type=0;
      optGetMessage(opt, i, msg, &msg_type);
      if (msg_type<=optMsgFileLeave || msg_type==optMsgUserError) {
        printOut(gev, "%s\n",msg);
      }
    }
    optClearMessages(opt);
    optEchoSet(opt, 0);
  }

  gmoObjStyleSet(gmo, gmoObjType_Fun);
  gmoIndexBaseSet(gmo, 0);
  gmoPinfSet(gmo, CUOPT_INFINITY);
  gmoMinfSet(gmo, -CUOPT_INFINITY);
  gmoSetNRowPerm(gmo);
  gmoSolveStatSet(gmo, gmoSolveStat_Capability);
  gmoModelStatSet(gmo, gmoModelStat_NoSolutionReturned);

  /* Activate Q-mode for QCP and MIQCP models */
  if (gmoModelType(gmo) == gmoProc_qcp || gmoModelType(gmo) == gmoProc_miqcp)
    gmoUseQSet(gmo, 1);

  cuOptOptimizationProblem problem = NULL;
  cuOptSolverSettings settings = NULL;
  cuOptSolution solution = NULL;

  cuopt_int_t num_variables = gmoN(gmo);
  cuopt_int_t num_constraints = gmoM(gmo);
  cuopt_int_t nnz = gmoNZ(gmo);
  //int64_t nnz = gmoNZ64(gmo);

  cuopt_int_t* constraint_matrix_row_offsets=NULL;
  cuopt_int_t* constraint_matrix_column_indices=NULL;
  cuopt_float_t* constraint_matrix_coefficent_values=NULL;
  cuopt_float_t* objective_coefficients=NULL;
  cuopt_float_t* rhs=NULL;
  cuopt_float_t* lower_bounds=NULL;
  cuopt_float_t* upper_bounds=NULL;
  char* constraint_sense=NULL;
  char* variable_types=NULL;

  // QCP specific mapping arrays
  int *gams2cuopt_row = NULL;
  cuopt_float_t *orig_rhs = NULL;
  char *orig_sense = NULL;

  int has_integer_vars = 0;
  fln_mip_trace[0] = '\0';
  sl_state_t context;
  int mipstart = 0;

  // Create solver settings
  status = cuOptCreateSolverSettings(&settings);
  if (status != CUOPT_SUCCESS) {
    printOut(gev, "Error creating solver settings: %d\n", status);
    goto DONE;
  }

  // Setup miptrace facility if option is enabled
  if (optGetDefinedStr(opt, "miptrace"))
  {
    if (gmoModelType(gmo) == gmoProc_mip || gmoModelType(gmo) == gmoProc_miqcp)
    {
      optGetStrStr(opt, "miptrace", fln_mip_trace);
      char sval2[256];
      if (mip_trace_open(fln_mip_trace, "cuOpt", gmoOptFile(gmo), gmoNameInput(gmo, sval2)))
      {
        printOut(gev, "Error opening trace file >%s<!", sval2);
        goto DONE;
      }
    }
    else
      printOut(gev, "WARNING: Enabling a MIP trace is only allowed for model type MIP/MIQCP!\n");
  }
  if (fp_mip_trace)
  {
    context.gev = gev;
    context.tstart = gevTimeJNow(gev);
    context.nvars = num_variables;
    mip_trace_line('S', 0, 0, 0, GMS_SV_NA, GMS_SV_NA);
    status = cuOptSetMIPGetSolutionCallback(settings, mip_get_solution_cb, &context);
    if (status != CUOPT_SUCCESS) {
      printOut(gev, "Error setting get-solution callback\n", status);
      goto DONE;
    }
  }

  // Check for MIP start feature
  mipstart = optGetIntStr(opt, "mipstart");
  if (mipstart && (gmoModelType(gmo) != gmoProc_mip && gmoModelType(gmo) != gmoProc_miqcp))
  {
    printOut(gev, "WARNING: Setting a MIP start is only allowed for model type MIP/MIQCP!\n");
    mipstart = 0;
  }

  // Set solver parameters with GAMS options
  if (gevGetIntOpt(gev, gevThreadsRaw) != 0)
  {
    status = cuOptSetIntegerParameter(settings, CUOPT_NUM_CPU_THREADS, gevGetIntOpt(gev, gevThreadsRaw));
    if (status != CUOPT_SUCCESS) {
      printOut(gev, "Error setting number of CPU threads: %d\n", status);
      goto DONE;
    }
  }
  if (gevGetIntOpt(gev, gevIterLim) < ITERLIM_INFINITY) {
    status = cuOptSetIntegerParameter(settings, CUOPT_ITERATION_LIMIT, gevGetIntOpt(gev, gevIterLim));
    if (status != CUOPT_SUCCESS) {
      printOut(gev, "Error setting iteration limit: %d\n", status);
      goto DONE;
    }
  }
  if (gevGetDblOpt(gev, gevResLim) < RESLIM_INFINITY) {
    status = cuOptSetFloatParameter(settings, CUOPT_TIME_LIMIT, gevGetDblOpt(gev, gevResLim));
    if (status != CUOPT_SUCCESS) {
      printOut(gev, "Error setting time limit: %d\n", status);
      goto DONE;
    }
  }
  if (gmoModelType(gmo) == gmoProc_mip || gmoModelType(gmo) == gmoProc_miqcp)
  {
    status = cuOptSetFloatParameter(settings, CUOPT_MIP_ABSOLUTE_GAP, gevGetDblOpt(gev, gevOptCA));
    if (status != CUOPT_SUCCESS) {
      printOut(gev, "Error setting absolute gap: %d\n", status);
      goto DONE;
    }
    status = cuOptSetFloatParameter(settings, CUOPT_MIP_RELATIVE_GAP, gevGetDblOpt(gev, gevOptCR));
    if (status != CUOPT_SUCCESS) {
      printOut(gev, "Error setting relative gap: %d\n", status);
      goto DONE;
    }
  }

  for (int i = 1; i <= optCount(opt); i++) {
    int defined=0, data_type=0, linkopt=0, ival=0, unused=0;
    double dval=0.0;
    char optname[256], sval[256]="";
    optGetInfoNr(opt, i, &defined, &unused, &linkopt, &data_type, &unused, &unused);
    if (data_type == optDataNone || !defined || linkopt) {
      continue;
    }
    optGetValuesNr(opt, i, optname, &ival, &dval, sval);

    if (data_type == optDataInteger) {
      status = cuOptSetIntegerParameter(settings, optname, ival);
      if (status != CUOPT_SUCCESS) {
        printOut(gev, "Error setting integer option >%s<: %d\n", optname, status);
        goto DONE;
      }
    } else if(data_type == optDataDouble) {
      status = cuOptSetFloatParameter(settings, optname, dval);
      if (status != CUOPT_SUCCESS) {
        printOut(gev, "Error setting float option >%s<: %d\n", optname, status);
        goto DONE;
      }
    }
  }

  // Try taking primal or dual (marginal) values from user (for LPs)
  if (gmoModelType(gmo) == gmoProc_lp)
  {
    cuopt_int_t chosen_method;
    status = cuOptGetIntegerParameter(settings, "method", &chosen_method);
    if (status != CUOPT_SUCCESS)
    {
      printOut(gev, "Error querying method option.\n");
      goto DONE;
    }
    // only when some PDLP is used and we have basis
    if ((chosen_method == CUOPT_METHOD_PDLP || chosen_method == CUOPT_METHOD_CONCURRENT) && gmoHaveBasis(gmo))
    {
      int nvars = gmoN(gmo), nconstraints = gmoM(gmo);
      double *lvls = (double *)malloc(sizeof(double) * nvars);
      double *marginals = (double *)malloc(sizeof(double) * nconstraints);
      gmoGetVarL(gmo, lvls);
      gmoGetEquM(gmo, marginals);
#if defined(CUOPT_INSTANTIATE_DOUBLE)
      status = cuOptSetInitialPrimalSolution(settings, lvls, nvars);
      if (status != CUOPT_SUCCESS)
      {
        printOut(gev, "Error setting primal solution for LP.\n");
        goto DONE;
      }
      status = cuOptSetInitialDualSolution(settings, marginals, nconstraints);
      if (status != CUOPT_SUCCESS)
      {
        printOut(gev, "Error setting dual solution for LP.\n");
        goto DONE;
      }
#else
      cuopt_float_t *lvlsf = (cuopt_float_t *)malloc(sizeof(cuopt_float_t) * nvars);
      cuopt_float_t *marginalsf = (cuopt_float_t *)malloc(sizeof(cuopt_float_t) * nconstraints);
      for (int i = 0; i < nvars; i++)
        lvlsf[i] = (cuopt_float_t)lvls[i];
      for (int i = 0; i < nconstraints; i++)
        marginalsf[i] = (cuopt_float_t)marginals[i];
      status = cuOptSetInitialPrimalSolution(settings, lvlsf, nvars);
      if (status != CUOPT_SUCCESS)
      {
        printOut(gev, "Error setting primal solution for LP.\n");
        goto DONE;
      }
      status = cuOptSetInitialDualSolution(settings, marginalsf, nconstraints);
      if (status != CUOPT_SUCCESS)
      {
        printOut(gev, "Error setting dual solution for LP.\n");
        goto DONE;
      }
      free(lvlsf);
      free(marginalsf);
#endif
      free(lvls);
      free(marginals);
      printOut(gev, "Initial primal and dual solutions have been set.\n");
    }
  }

  if (!optGetDefinedStr(opt, "prob_read"))
  {
    int is_qcp = (gmoModelType(gmo) == gmoProc_qcp || gmoModelType(gmo) == gmoProc_miqcp);
    int num_linear_constraints = 0;
    int num_quad_constraints = 0;

    gams2cuopt_row = malloc(num_constraints * sizeof(int));
    orig_rhs = malloc(num_constraints * sizeof(cuopt_float_t));
    orig_sense = malloc(num_constraints * sizeof(char));

    // First pass: Count linear vs quadratic and build index mapping
    for (int i = 0; i < num_constraints; i++)
    {
      int qnz = 0;
      if (is_qcp)
      {
        qnz = gmoGetRowQNZOne(gmo, i);
        if (qnz < 0)
          qnz = 0; // clamp -1 to 0
      }

      if (qnz == 0)
      {
        gams2cuopt_row[i] = num_linear_constraints++;
      }
      else
      {
        num_quad_constraints++;
      }
    }

    // Append quadratic constraint indices sequentially
    int q_idx = num_linear_constraints;
    for (int i = 0; i < num_constraints; i++)
    {
      int qnz = 0;
      if (is_qcp)
      {
        qnz = gmoGetRowQNZOne(gmo, i);
        if (qnz < 0)
          qnz = 0; // clamp -1 to 0
      }
      if (qnz > 0)
      {
        gams2cuopt_row[i] = q_idx++;
      }
    }

    constraint_matrix_row_offsets = malloc((num_linear_constraints + 1) * sizeof(cuopt_int_t));
    constraint_matrix_column_indices = malloc(nnz * sizeof(cuopt_int_t));
    constraint_matrix_coefficent_values = malloc(nnz * sizeof(cuopt_float_t));
    objective_coefficients = malloc((num_variables) * sizeof(cuopt_float_t));
    rhs = malloc((num_linear_constraints) * sizeof(cuopt_float_t));
    lower_bounds = malloc((num_variables) * sizeof(cuopt_float_t));
    upper_bounds = malloc((num_variables) * sizeof(cuopt_float_t));
    constraint_sense = malloc((num_linear_constraints) * sizeof(char));
    variable_types = malloc((num_variables) * sizeof(char));

    if ((constraint_matrix_row_offsets == NULL) ||
        (constraint_matrix_column_indices == NULL) ||
        (constraint_matrix_coefficent_values == NULL) ||
        (objective_coefficients == NULL) ||
        (rhs == NULL) ||
        (lower_bounds == NULL) ||
        (upper_bounds == NULL) ||
        (constraint_sense == NULL) ||
        (variable_types == NULL)) {
      printOut(gev, "Could not allocate arrays\n");
      goto DONE;
    }

    int *temp_equ_types = malloc(num_constraints * sizeof(int));
    status = gmoGetEquType(gmo, temp_equ_types);
    if (status)
    {
      printOut(gev, "gmoGetEquType failed. Status: %d\n", status);
      free(temp_equ_types);
      goto DONE;
    }

    for (int i = 0; i < num_constraints; i++)
    {
      switch (temp_equ_types[i])
      {
      case gmoequ_E:
        orig_sense[i] = CUOPT_EQUAL;
        break;
      case gmoequ_L:
        orig_sense[i] = CUOPT_LESS_THAN;
        break;
      case gmoequ_G:
        orig_sense[i] = CUOPT_GREATER_THAN;
        break;
      default:
        printOut(gev, "Known row type %d\n", temp_equ_types[i]);
      }
    }
    free(temp_equ_types);

    int *temp_var_types = malloc(num_variables * sizeof(int));
    status = gmoGetVarType(gmo, temp_var_types);
    if (status)
    {
      printOut(gev, "gmoGetVarType failed. Status: %d\n", status);
      free(temp_var_types);
      goto DONE;
    }

    for (int j = 0; j < num_variables; j++)
    {
      switch (temp_var_types[j])
      {
      case gmovar_X:
        variable_types[j] = CUOPT_CONTINUOUS;
        break;
      case gmovar_B:
      case gmovar_I:
        variable_types[j] = CUOPT_INTEGER;
        has_integer_vars = 1;
        break;
      case gmovar_SC:
        variable_types[j] = CUOPT_SEMI_CONTINUOUS;
        has_integer_vars = 1;
        break;
      default:
        printOut(gev, "Known variable type %d\n", temp_var_types[j]);
      }
    }
    free(temp_var_types);

    status = gmoGetVarLower(gmo, lower_bounds);
    if (status)
    {
      printOut(gev, "gmoGetVarLower failed. Status: %d\n", status);
      goto DONE;
    }

    status = gmoGetVarUpper(gmo, upper_bounds);
    if (status)
    {
      printOut(gev, "gmoGetVarUpper failed. Status: %d\n", status);
      goto DONE;
    }

    status = gmoGetRhs(gmo, orig_rhs);
    if (status)
    {
      printOut(gev, "gmoGetRhs failed. Status: %d\n", status);
      goto DONE;
    }

    nnz = 0;
    int lin_row = 0;
    for (int i = 0; i < num_constraints; i++)
    {
      int qnz = 0;
      if (is_qcp)
      {
        qnz = gmoGetRowQNZOne(gmo, i);
        if (qnz < 0)
          qnz = 0; // clamp -1 to 0
      }

      // Pack ONLY purely linear constraints
      if (qnz == 0)
      {
        constraint_matrix_row_offsets[lin_row] = nnz;
        int rnz = 0, rnlnz = 0;
        status = gmoGetRowSparse(gmo, i, constraint_matrix_column_indices + nnz, constraint_matrix_coefficent_values + nnz, NULL, &rnz, &rnlnz);
        if (status)
        {
          printOut(gev, "gmoGetRowSparse %d failed. Status: %d\n", i, status);
          goto DONE;
        }

        constraint_sense[lin_row] = orig_sense[i];
        rhs[lin_row] = orig_rhs[i];

        nnz += rnz;
        lin_row++;
      }
    }
    constraint_matrix_row_offsets[num_linear_constraints] = nnz;

    status = gmoGetObjVector(gmo, objective_coefficients, NULL);
    if (status) {
      printOut(gev, "gmoGetObjVector failed. Status: %d\n", status);
      goto DONE;
    }

    status = cuOptCreateProblem(
        num_linear_constraints, // Use mapped linear size
        num_variables,
        (gmoSense(gmo) == gmoObj_Min) ? CUOPT_MINIMIZE : CUOPT_MAXIMIZE,
        gmoObjConst(gmo),
        objective_coefficients,
        constraint_matrix_row_offsets,
        constraint_matrix_column_indices,
        constraint_matrix_coefficent_values,
        constraint_sense,
        rhs,
        lower_bounds,
        upper_bounds,
        variable_types,
        &problem);

    if (status != CUOPT_SUCCESS)
    {
      printOut(gev, "Error creating linear base problem from GAMS model: %d\n", status);
      goto DONE;
    }

    if (is_qcp)
    {
      int obj_qnz = gmoObjQNZ(gmo);
      if (obj_qnz > 0)
      {
        int *temp_q_row = malloc(obj_qnz * sizeof(int));
        int *temp_q_col = malloc(obj_qnz * sizeof(int));
        double *temp_q_coef = malloc(obj_qnz * sizeof(double));

        // Extract the quadratic objective coefficients from GAMS
        gmoGetObjQ(gmo, temp_q_row, temp_q_col, temp_q_coef);

        cuopt_int_t *q_row = malloc(obj_qnz * sizeof(cuopt_int_t));
        cuopt_int_t *q_col = malloc(obj_qnz * sizeof(cuopt_int_t));
        cuopt_float_t *q_coef = malloc(obj_qnz * sizeof(cuopt_float_t));

        for (int k = 0; k < obj_qnz; k++)
        {
          q_row[k] = (cuopt_int_t)temp_q_row[k];
          q_col[k] = (cuopt_int_t)temp_q_col[k];

          // GAMS provides hessian matrix. Diagonal elements must be halved!
          if (q_row[k] == q_col[k])
            q_coef[k] = (cuopt_float_t)(temp_q_coef[k] / 2.0);
          // Non-diagonal elements can be passed through directly to cuOpt
          else
            q_coef[k] = (cuopt_float_t)temp_q_coef[k];
        }

        // Apply the quadratic terms to the objective using the cuOpt API
        status = cuOptSetQuadraticObjective(problem, obj_qnz, q_row, q_col, q_coef);

        free(temp_q_row);
        free(temp_q_col);
        free(temp_q_coef);
        free(q_row);
        free(q_col);
        free(q_coef);

        if (status != CUOPT_SUCCESS)
        {
          printOut(gev, "Error setting quadratic objective: %d\n", status);
          goto DONE;
        }
      }
    }

    // Append the quadratic constraints dynamically
    for (int i = 0; i < num_constraints; i++)
    {
      int qnz = 0;
      if (is_qcp)
      {
        qnz = gmoGetRowQNZOne(gmo, i);
        if (qnz < 0)
          qnz = 0; // clamp -1 to 0
      }

      if (qnz > 0)
      {
        int lin_nz = 0, rnlnz = 0;
        int *temp_lin_cols = malloc(num_variables * sizeof(int));
        double *temp_lin_vals = malloc(num_variables * sizeof(double));
        gmoGetRowSparse(gmo, i, temp_lin_cols, temp_lin_vals, NULL, &lin_nz, &rnlnz);

        cuopt_int_t *lin_cols = malloc(lin_nz * sizeof(cuopt_int_t));
        cuopt_float_t *lin_vals = malloc(lin_nz * sizeof(cuopt_float_t));
        for (int k = 0; k < lin_nz; k++)
        {
          lin_cols[k] = (cuopt_int_t)temp_lin_cols[k];
          lin_vals[k] = (cuopt_float_t)temp_lin_vals[k];
        }

        int *temp_q_row = malloc(qnz * sizeof(int));
        int *temp_q_col = malloc(qnz * sizeof(int));
        double *temp_q_coef = malloc(qnz * sizeof(double));
        gmoGetRowQ(gmo, i, temp_q_row, temp_q_col, temp_q_coef);

        cuopt_int_t *q_row = malloc(qnz * sizeof(cuopt_int_t));
        cuopt_int_t *q_col = malloc(qnz * sizeof(cuopt_int_t));
        cuopt_float_t *q_coef = malloc(qnz * sizeof(cuopt_float_t));
        for (int k = 0; k < qnz; k++)
        {
          q_row[k] = (cuopt_int_t)temp_q_row[k];
          q_col[k] = (cuopt_int_t)temp_q_col[k];

          // GAMS provides hessian matrix. Diagonal elements must be halved!
          if (q_row[k] == q_col[k])
            q_coef[k] = (cuopt_float_t)(temp_q_coef[k] / 2.0);
          // non-diagonal elements can be passed through directly to cuOpt
          else
            q_coef[k] = (cuopt_float_t)temp_q_coef[k];
        }

        status = cuOptAddQuadraticConstraint(
            problem,
            qnz, q_row, q_col, q_coef,
            lin_nz, lin_cols, lin_vals,
            orig_sense[i], orig_rhs[i]);

        free(temp_lin_cols);
        free(temp_lin_vals);
        free(lin_cols);
        free(lin_vals);
        free(temp_q_row);
        free(temp_q_col);
        free(temp_q_coef);
        free(q_row);
        free(q_col);
        free(q_coef);

        if (status != CUOPT_SUCCESS)
        {
          printOut(gev, "Error adding quadratic constraint %d: %d\n", i, status);
          goto DONE;
        }
      }
    }
  }
  else
  {
    status = cuOptReadProblem(optGetStrStr(opt, "prob_read", filename), &problem);
    if (status != CUOPT_SUCCESS) {
      printOut(gev, "Error creating problem from MPS file: %d\n", status);
      goto DONE;
    }
  }

  char logfilename[256];
  gevGetScratchName(gev, "cuopt", logfilename);
  cuOptSetParameter(settings, CUOPT_LOG_FILE, logfilename);
  int loglevel = gevGetIntOpt(gev, gevLogOption);
  if (loglevel != 3) {
    if ((loglevel==0) || (loglevel==2)) {
      cuOptSetIntegerParameter(settings, CUOPT_LOG_TO_CONSOLE, 0);
    } else if (loglevel==4) {
      cuOptSetIntegerParameter(settings, CUOPT_LOG_TO_CONSOLE, 1);
    }
  }

  // Maybe add mip start
  if(mipstart)
  {
    double *initial_levels = malloc(sizeof(double) * gmoN(gmo));
    gmoGetVarL(gmo, initial_levels);
#ifdef CUOPT_INSTANTIATE_DOUBLE
    status = cuOptAddMIPStart(settings, initial_levels, gmoN(gmo));
    if (status != CUOPT_SUCCESS)
    {
      printOut(gev, "Error setting MIP start.\n");
      goto DONE;
    }
#else
    cuopt_float_t *initial_levelsf = malloc(sizeof(cuopt_float_t) * gmoN(gmo));
    for (int i = 0; i < gmoN(gmo); i++)
      initial_levelsf[i] = (cuopt_float_t)initial_levels[i];
    status = cuOptAddMIPStart(settings, initial_levelsf, gmoN(gmo));
    if (status != CUOPT_SUCCESS)
    {
      printOut(gev, "Error setting MIP start.\n");
      goto DONE;
    }
    free(initial_levelsf);
#endif
    free(initial_levels);
    printOut(gev, "MIP start has been set.\n");
  }

  // Solve the problem
  status = cuOptSolve(problem, settings, &solution);
  if (status != CUOPT_SUCCESS) {
    printOut(gev, "Error solving problem: %d\n", status);
    goto DONE;
  }

  if (loglevel == 2) { // sorry no cuOpt log for logOption=4 in the log file
    FILE *cuoptlogfile=fopen(logfilename,"r");
    if (!cuoptlogfile) {
      printOut(gev, "Error opening cuopt log file\n");
      goto DONE;
    }
    char line[256];
    while (fgets(line, sizeof(line), cuoptlogfile) != NULL) {
      printOut(gev, "%s", line);
    }
    fclose(cuoptlogfile);
  }

  // Get solution information
  cuopt_float_t solution_time;
  cuopt_int_t termination_status;
  cuopt_float_t objective_value, solution_bound;

  status = cuOptGetTerminationStatus(solution, &termination_status);
  if (status != CUOPT_SUCCESS) {
    printOut(gev, "Error getting termination status: %d\n", status);
    goto DONE;
  }

  gmoSolveStatSet(gmo, gmoSolveStat_Normal);
  gmoModelStatSet(gmo, gmoModelStat_NoSolutionReturned);
  if (!optGetDefinedStr(opt, "prob_read")) {
    switch (termination_status) {
    case CUOPT_TERMINATION_STATUS_OPTIMAL:
      gmoModelStatSet(gmo, gmoModelStat_OptimalGlobal);
      break;
    case CUOPT_TERMINATION_STATUS_INFEASIBLE:
      gmoModelStatSet(gmo, gmoModelStat_InfeasibleGlobal);
      break;
    case CUOPT_TERMINATION_STATUS_UNBOUNDED:
      gmoModelStatSet(gmo, gmoModelStat_Unbounded);
      break;
    case CUOPT_TERMINATION_STATUS_ITERATION_LIMIT:
      gmoSolveStatSet(gmo, gmoSolveStat_Iteration);
      break;
    case CUOPT_TERMINATION_STATUS_TIME_LIMIT:
      gmoSolveStatSet(gmo, gmoSolveStat_Resource);
      break;
    case CUOPT_TERMINATION_STATUS_NUMERICAL_ERROR:
      gmoSolveStatSet(gmo, gmoSolveStat_SolverErr);
      break;
    case CUOPT_TERMINATION_STATUS_PRIMAL_FEASIBLE:
      if ((gmoModelType(gmo) == gmoProc_mip || gmoModelType(gmo) == gmoProc_miqcp) && has_integer_vars)
      {
        gmoModelStatSet(gmo, gmoModelStat_Integer);
        } else {
        gmoModelStatSet(gmo, gmoModelStat_Feasible);
      }
      break;
    case CUOPT_TERMINATION_STATUS_FEASIBLE_FOUND:
      if ((gmoModelType(gmo) == gmoProc_mip || gmoModelType(gmo) == gmoProc_miqcp) && has_integer_vars)
      {
        gmoModelStatSet(gmo, gmoModelStat_Integer);
        } else {
        gmoModelStatSet(gmo, gmoModelStat_Feasible);
      }
      break;
    case CUOPT_TERMINATION_STATUS_UNBOUNDED_OR_INFEASIBLE:
    default:
      gmoSolveStatSet(gmo, gmoSolveStat_Solver);
    }
  }

  status = cuOptGetSolveTime(solution, &solution_time);
  if (status != CUOPT_SUCCESS) {
    printOut(gev, "Error getting solve time: %d\n", status);
    goto DONE;
  }
  gmoSetHeadnTail(gmo, gmoHresused, solution_time);

  if (gmoModelStat(gmo) != gmoModelStat_NoSolutionReturned) {
    status = cuOptGetObjectiveValue(solution, &objective_value);
    if (status != CUOPT_SUCCESS) {
      printOut(gev, "Error getting objective value: %d\n", status);
      goto DONE;
    }
    gmoSetHeadnTail(gmo, gmoHobjval, objective_value);

    if ((gmoModelType(gmo) == gmoProc_mip || gmoModelType(gmo) == gmoProc_miqcp) && has_integer_vars)
    {
      status = cuOptGetSolutionBound(solution, &solution_bound);
      if (status != CUOPT_SUCCESS) {
        printOut(gev, "Error getting solution bound: %d\n", status);
        goto DONE;
      }
      gmoSetHeadnTail(gmo, gmoTmipbest, solution_bound);
    }

    status = cuOptGetPrimalSolution(solution, objective_coefficients); // reuse n-vector
    if (status != CUOPT_SUCCESS) {
      printOut(gev, "Error getting primal solution: %d\n", status);
      goto DONE;
    }
    gmoSetVarL(gmo, objective_coefficients);

    if(fp_mip_trace)
    {
      double total_elapsed = (gevTimeJNow(gev) - context.tstart) * 3600.0 * 24.0;
      mip_trace_line('E', 0, 1, total_elapsed, objective_value, solution_bound);
    }

    int request_marginals = gevGetIntOpt(gev, gevRequestMarginals);
    cuopt_int_t presolve = 0, dual_postsolve = 0;
    cuOptGetIntegerParameter(settings, "presolve", &presolve);
    cuOptGetIntegerParameter(settings, "dual_postsolve", &dual_postsolve);

    int is_mip_model = gmoModelType(gmo) == gmoProc_mip || has_integer_vars;
    int is_linear_qcp = gmoModelType(gmo) == gmoProc_qcp && !has_integer_vars;

    // Only extract duals, when it's neither MIP nor QCP (and marginals aren't explicitly NOT WANTED with requestMarginals=0)
    if (!request_marginals)
    {
      gmoSetHeadnTail(gmo, gmoHmarginals, 0.0);
    }
    // Linear/quadratic and fully continuous
    else if ((gmoModelType(gmo) == gmoProc_lp || is_linear_qcp) && !is_mip_model)
    {
      if (!presolve || dual_postsolve)
      {
        status = cuOptGetReducedCosts(solution, objective_coefficients); // reuse n-vector
        if (status != CUOPT_SUCCESS)
        {
          printOut(gev, "Error getting reduced cost: %d\n", status);
          goto DONE;
        }
        if (gmoSense(gmo) == gmoObj_Max)
        {
          // patch duals for max problem
          for (int j = 0; j < num_variables; j++)
          {
            objective_coefficients[j] *= -1.0;
          }
        }
        gmoSetVarM(gmo, objective_coefficients);

        // Extract duals using explicit mapping array
        cuopt_float_t *raw_duals = malloc(num_constraints * sizeof(cuopt_float_t));
        status = cuOptGetDualSolution(solution, raw_duals);
        if (status != CUOPT_SUCCESS)
        {
          printOut(gev, "Error getting dual solution: %d\n", status);
          free(raw_duals);
          goto DONE;
        }

        double *final_duals = malloc(num_constraints * sizeof(double));
        for (int i = 0; i < num_constraints; i++)
        {
          if (gams2cuopt_row)
          {
            final_duals[i] = raw_duals[gams2cuopt_row[i]];
          }
          else
          {
            final_duals[i] = raw_duals[i]; // Fallback if directly read from MPS
          }

          if (gmoSense(gmo) == gmoObj_Max)
          {
            final_duals[i] *= -1.0; // patch duals for max problem
          }
        }
        gmoSetEquM(gmo, final_duals);
        free(raw_duals);
        free(final_duals);
      }
      else
      {
        gmoSetHeadnTail(gmo, gmoHmarginals, 0.0);
      }
    }
    // User explicitly forces marginals on MIQCP / MIP / Problems with discrete vars (requestMarginals = 1)
    else if (request_marginals == 1)
    {
      if (gmoModelType(gmo) == gmoProc_miqcp || has_integer_vars)
      {
        printOut(gev, "WARNING: requestMarginals=1 ignored. cuOpt does not currently support dual solutions for MIQCP/discrete models.\n");
      }
      else
      {
        printOut(gev, "WARNING: requestMarginals=1 ignored. cuOpt link does not currently support continuous subproblem solves for MIP marginals.\n");
      }
      gmoSetHeadnTail(gmo, gmoHmarginals, 0.0);
    }
    // Default fallback for MIP / MIQCP when requestMarginals is -1 (inexpensive only)
    else
    {
      gmoSetHeadnTail(gmo, gmoHmarginals, 0.0);
    }

    gmoCompleteSolution(gmo);
  }
  else if (fp_mip_trace) // gmoModelStat(gmo) == gmoModelStat_NoSolutionReturned
  {
    double total_elapsed = (gevTimeJNow(gev) - context.tstart) * 3600.0 * 24.0;
    solution_bound = GMS_SV_NA;
    if ((gmoModelType(gmo) == gmoProc_mip || gmoModelType(gmo) == gmoProc_miqcp) && has_integer_vars)
    {
      status = cuOptGetSolutionBound(solution, &solution_bound);
      if (status != CUOPT_SUCCESS)
        solution_bound = GMS_SV_NA;
    }
    mip_trace_line('E', 0, 0, total_elapsed, GMS_SV_NA, solution_bound);
  }
  status = gmoUnloadSolutionLegacy(gmo);
  if (status) {
    printOut(gev, "Problems unloading solution\n");
  }

DONE:
  cuOptDestroyProblem(&problem);
  cuOptDestroySolverSettings(&settings);
  cuOptDestroySolution(&solution);

  free(constraint_matrix_row_offsets);
  free(constraint_matrix_column_indices);
  free(constraint_matrix_coefficent_values);
  free(objective_coefficients);
  free(rhs);
  free(lower_bounds);
  free(upper_bounds);
  free(constraint_sense);
  free(variable_types);

  // Free dynamically mapped structures
  free(gams2cuopt_row);
  free(orig_rhs);
  free(orig_sense);

  if (fp_mip_trace)
    mip_trace_close();

GAMSDONE:
  gmoFree(&gmo);
  gevFree(&gev);
  optFree(&opt);

  return 0;

} /* main */

int mip_trace_open(const char *fname, const char *solverID, const int optFileNum, const char *inputName)
{
  if (NULL != fp_mip_trace)
    return 1; /* already open: error */

  strcpy(fln_mip_trace, fname);
  fp_mip_trace = fopen(fln_mip_trace, "w");
  if (NULL == fp_mip_trace)
    return 3;

  strncpy(mip_trace_id, solverID, sizeof(mip_trace_id) - 1);
  mip_trace_id[sizeof(mip_trace_id) - 1] = '\0';
  mip_trace_seq = 1;
  fprintf(fp_mip_trace, "* mip_trace_ file %s: ID = %s.%d Instance = %s\n", fln_mip_trace, mip_trace_id, optFileNum, inputName);
  fprintf(fp_mip_trace, "* fields are lineNum, seriesID, node, seconds, bestFound, bestBound\n");
  fflush(fp_mip_trace);
  return 0;
} /* mip_trace_open */

int mip_trace_close()
{
  int rc;
  if (NULL == fp_mip_trace)
    return 2; /* already closed: error */
  fprintf(fp_mip_trace, "* mip_trace_ file %s closed\n", fln_mip_trace);
  rc = fclose(fp_mip_trace);
  fp_mip_trace = NULL;
  return (0 == rc) ? 0 : 1;
} /* mip_trace_close */

#define bnd_na(x) x == GMS_SV_NA || x == HUGE_VAL || x == -HUGE_VAL

int mip_trace_line(char seriesID, double node, int giveint,
                   double seconds, double bestint, double bestbnd)
{
  int rc;

  if (NULL == fp_mip_trace)
    return -1; /* not open: error */

  if (giveint)
  {
    if (bnd_na(bestbnd))
      rc = fprintf(fp_mip_trace, "%d, %c, %g, %.15g, %.15g, na\n", mip_trace_seq,
                   isalnum(seriesID) ? seriesID : 'X',
                   node, seconds, bestint);
    else
      rc = fprintf(fp_mip_trace, "%d, %c, %g, %.15g, %.15g, %.15g\n", mip_trace_seq,
                   isalnum(seriesID) ? seriesID : 'X',
                   node, seconds, bestint, bestbnd);
  }
  else
  {
    if (bnd_na(bestbnd))
      rc = fprintf(fp_mip_trace, "%d, %c, %g, %.15g, na, na\n", mip_trace_seq,
                   isalnum(seriesID) ? seriesID : 'X',
                   node, seconds);
    else
      rc = fprintf(fp_mip_trace, "%d, %c, %g, %.15g, na, %g\n", mip_trace_seq,
                   isalnum(seriesID) ? seriesID : 'X',
                   node, seconds, bestbnd);
  }
  fflush(fp_mip_trace);
  mip_trace_seq++;

  return rc;
} /* mip_trace_line */

static void mip_get_solution_cb(const cuopt_float_t *solution, const cuopt_float_t *objective_value,
                                const cuopt_float_t *solution_bound, void *user_data){
  sl_state_t *state = (sl_state_t *)user_data;
  double elapsed = (gevTimeJNow(state->gev) - state->tstart) * 3600.0 * 24.0;
  double obj = *objective_value;
  double bnd = *solution_bound;
  mip_trace_line('I', 0, 1, elapsed, obj, bnd);
}

#if 0
t program for cuOpt linear programming solver
 */

// Include the cuOpt linear programming solver header
#include <cuopt/linear_programming/cuopt_c.h>
#include <stdio.h>
#include <stdlib.h>

// Convert termination status to string
const char* termination_status_to_string(cuopt_int_t termination_status)
{
  switch (termination_status) {
    case CUOPT_TERMINATION_STATUS_OPTIMAL:
      return "Optimal";
    case CUOPT_TERMINATION_STATUS_INFEASIBLE:
      return "Infeasible";
    case CUOPT_TERMINATION_STATUS_UNBOUNDED:
      return "Unbounded";
    case CUOPT_TERMINATION_STATUS_ITERATION_LIMIT:
      return "Iteration limit";
    case CUOPT_TERMINATION_STATUS_TIME_LIMIT:
      return "Time limit";
    case CUOPT_TERMINATION_STATUS_NUMERICAL_ERROR:
      return "Numerical error";
    case CUOPT_TERMINATION_STATUS_PRIMAL_FEASIBLE:
      return "Primal feasible";
    case CUOPT_TERMINATION_STATUS_FEASIBLE_FOUND:
      return "Feasible found";
    default:
      return "Unknown";
  }
}

// Test simple LP problem
cuopt_int_t test_simple_lp()
{
  cuOptOptimizationProblem problem = NULL;
  cuOptSolverSettings settings = NULL;
  cuOptSolution solution = NULL;

  /* Solve the following LP:
     minimize -0.2*x1 + 0.1*x2
     subject to:
     3.0*x1 + 4.0*x2 <= 5.4
     2.7*x1 + 10.1*x2 <= 4.9
     x1, x2 >= 0
  */

  cuopt_int_t num_variables = 2;
  cuopt_int_t num_constraints = 2;
  cuopt_int_t nnz = 4;

  // CSR format constraint matrix
  // https://docs.nvidia.com/nvpl/latest/sparse/storage_format/sparse_matrix.html#compressed-sparse-row-csr
  // From the constraints:
  // 3.0*x1 + 4.0*x2 <= 5.4
  // 2.7*x1 + 10.1*x2 <= 4.9
  cuopt_int_t row_offsets[] = {0, 2, 4};
  cuopt_int_t column_indices[] = {0, 1, 0, 1};
  cuopt_float_t values[] = {3.0, 4.0, 2.7, 10.1};

  // Objective coefficients
  // From the objective function: minimize -0.2*x1 + 0.1*x2
  // -0.2 is the coefficient of x1
  // 0.1 is the coefficient of x2
  cuopt_float_t objective_coefficients[] = {-0.2, 0.1};

  // Constraint bounds
  // From the constraints:
  // 3.0*x1 + 4.0*x2 <= 5.4
  // 2.7*x1 + 10.1*x2 <= 4.9
  cuopt_float_t constraint_upper_bounds[] = {5.4, 4.9};
  cuopt_float_t constraint_lower_bounds[] = {-CUOPT_INFINITY, -CUOPT_INFINITY};

  // Variable bounds
  // From the constraints:
  // x1, x2 >= 0
  cuopt_float_t var_lower_bounds[] = {0.0, 0.0};
  cuopt_float_t var_upper_bounds[] = {CUOPT_INFINITY, CUOPT_INFINITY};

  // Variable types (continuous)
  // From the constraints:
  // x1, x2 >= 0
  char variable_types[] = {CUOPT_CONTINUOUS, CUOPT_CONTINUOUS};

  cuopt_int_t status;
  cuopt_float_t time;
  cuopt_int_t termination_status;
  cuopt_float_t objective_value;

  printf("Creating and solving simple LP problem...\n");

  // Create the problem
  status = cuOptCreateRangedProblem(num_constraints,
                                   num_variables,
                                   CUOPT_MINIMIZE,  // minimize=False
                                   0.0,            // objective offset
                                   objective_coefficients,
                                   row_offsets,
                                   column_indices,
                                   values,
                                   constraint_lower_bounds,
                                   constraint_upper_bounds,
                                   var_lower_bounds,
                                   var_upper_bounds,
                                   variable_types,
                                   &problem);
  if (status != CUOPT_SUCCESS) {
    printf("Error creating problem: %d\n", status);
    goto DONE;
  }

  // Create solver settings
  status = cuOptCreateSolverSettings(&settings);
  if (status != CUOPT_SUCCESS) {
    printf("Error creating solver settings: %d\n", status);
    goto DONE;
  }

  // Set solver parameters
  status = cuOptSetFloatParameter(settings, CUOPT_ABSOLUTE_PRIMAL_TOLERANCE, 0.0001);
  if (status != CUOPT_SUCCESS) {
    printf("Error setting optimality tolerance: %d\n", status);
    goto DONE;
  }

  // Solve the problem
  status = cuOptSolve(problem, settings, &solution);
  if (status != CUOPT_SUCCESS) {
    printf("Error solving problem: %d\n", status);
    goto DONE;
  }

  // Get solution information
  status = cuOptGetSolveTime(solution, &time);
  if (status != CUOPT_SUCCESS) {
    printf("Error getting solve time: %d\n", status);
    goto DONE;
  }

  status = cuOptGetTerminationStatus(solution, &termination_status);
  if (status != CUOPT_SUCCESS) {
    printf("Error getting termination status: %d\n", status);
    goto DONE;
  }

  status = cuOptGetObjectiveValue(solution, &objective_value);
  if (status != CUOPT_SUCCESS) {
    printf("Error getting objective value: %d\n", status);
    goto DONE;
  }

  // Print results
  printf("\nResults:\n");
  printf("--------\n");
  printf("Termination status: %s (%d)\n", termination_status_to_string(termination_status), termination_status);
  printf("Solve time: %f seconds\n", time);
  printf("Objective value: %f\n", objective_value);

  // Get and print solution variables
  cuopt_float_t* solution_values = (cuopt_float_t*)malloc(num_variables * sizeof(cuopt_float_t));
  status = cuOptGetPrimalSolution(solution, solution_values);
  if (status != CUOPT_SUCCESS) {
    printf("Error getting solution values: %d\n", status);
    free(solution_values);
    goto DONE;
  }

  printf("\nPrimal Solution: Solution variables \n");
  for (cuopt_int_t i = 0; i < num_variables; i++) {
    printf("x%d = %f\n", i + 1, solution_values[i]);
  }
  free(solution_values);

DONE:
  cuOptDestroyProblem(&problem);
  cuOptDestroySolverSettings(&settings);
  cuOptDestroySolution(&solution);

  return status;
}

int main() {
  // Run the test
  cuopt_int_t status = test_simple_lp();

  if (status == CUOPT_SUCCESS) {
    printf("\nTest completed successfully!\n");
    return 0;
  } else {
    printf("\nTest failed with status: %d\n", status);
    return 1;
  }
}
#endif
