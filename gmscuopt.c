#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <strings.h>
#include <errno.h>
#include <ctype.h>
#include <math.h>
#include "gmomcc.h"
#include "gevmcc.h"
#include "optcc.h"
#include <cuopt/mathematical_optimization/cuopt_c.h>

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
  {
    cuopt_int_t vmajor = 0, vminor = 0, vpatch = 0;
    if (cuOptGetVersion(&vmajor, &vminor, &vpatch) == CUOPT_SUCCESS) {
      printOut(gev, "Using cuOpt library version: %d.%02d.%02d\n", vmajor, vminor, vpatch);
#if defined(CUOPT_VERSION)
      int bmajor = -1, bminor = -1;
      if (sscanf(CUOPT_VERSION, "%d.%d", &bmajor, &bminor) == 2 && (bmajor != vmajor || bminor != vminor))
        printOut(gev, "WARNING: cuOpt library version %d.%02d differs from the version %s the link was built against!\n",
                 vmajor, vminor, CUOPT_VERSION);
#endif
    }
  }

  status = gmoRegisterEnvironment(gmo, gev, msg);
  if (status) {
    printOut(gev, "Could not register GEV: %s\n", msg);
    goto GAMSDONE;
  }

  status = gmoLoadDataLegacy(gmo, msg);
  if (status) {
    printOut(gev, "Could not load model data: %s\n", msg);
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

  /* Activate Q-mode for QCP, RMIQCP and MIQCP models */
  int is_qcp = (gmoModelType(gmo) == gmoProc_qcp || gmoModelType(gmo) == gmoProc_rmiqcp ||
                gmoModelType(gmo) == gmoProc_miqcp);
  if (is_qcp)
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
  int *row_qnz = NULL; // number of Q nonzeros per GAMS row (0 for linear rows)
  cuopt_float_t *orig_rhs = NULL;
  char *orig_sense = NULL;

  int has_integer_vars = 0;
  int num_quad_constraints = 0;
  int obj_qnz = 0;
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
      printOut(gev, "WARNING: Enabling a MIP trace is only allowed for model type MIP!\n");
  }
  if (fp_mip_trace)
  {
    context.gev = gev;
    context.tstart = gevTimeJNow(gev);
    context.nvars = num_variables;
    mip_trace_line('S', 0, 0, 0, GMS_SV_NA, GMS_SV_NA);
    status = cuOptSetMIPGetSolutionCallback(settings, mip_get_solution_cb, &context);
    if (status != CUOPT_SUCCESS) {
      printOut(gev, "Error setting get-solution callback: %d\n", status);
      goto DONE;
    }
  }

  // Check for MIP start feature
  mipstart = optGetIntStr(opt, "mipstart");
  if (mipstart && (gmoModelType(gmo) != gmoProc_mip && gmoModelType(gmo) != gmoProc_miqcp))
  {
    printOut(gev, "WARNING: Setting a MIP start is only allowed for model type MIP!\n");
    mipstart = 0;
  }
  if (mipstart && optGetDefinedStr(opt, "prob_read"))
  {
    printOut(gev, "WARNING: Setting a MIP start is not supported together with 'prob_read'!\n");
    mipstart = 0;
  }

  // Set solver parameters with GAMS options
  // gevThreads resolves non-positive GAMS Threads values (e.g. -2 = all but two cores),
  // which cuOpt's num_cpu_threads (>= -1) would reject.
  if (gevGetIntOpt(gev, gevThreadsRaw) != 0)
  {
    status = cuOptSetIntegerParameter(settings, CUOPT_NUM_CPU_THREADS, gevThreads(gev));
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
    // cuOpt only accepts mip_relative_gap in [0, 0.1]
    double optcr = gevGetDblOpt(gev, gevOptCR);
    if (optcr > 0.1)
    {
      printOut(gev, "WARNING: cuOpt supports a relative gap of at most 0.1. Reducing OptCR from %g to 0.1.\n", optcr);
      optcr = 0.1;
    }
    status = cuOptSetFloatParameter(settings, CUOPT_MIP_RELATIVE_GAP, optcr);
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
    } else if(data_type == optDataString) {
      status = cuOptSetParameter(settings, optname, sval);
      if (status != CUOPT_SUCCESS) {
        printOut(gev, "Error setting string option >%s<: %d\n", optname, status);
        goto DONE;
      }
    }
  }

  // Try taking primal or dual (marginal) values from user (for LPs)
  if (gmoModelType(gmo) == gmoProc_lp && !optGetDefinedStr(opt, "prob_read"))
  {
    cuopt_int_t chosen_method;
    status = cuOptGetIntegerParameter(settings, CUOPT_METHOD, &chosen_method);
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
      cuopt_int_t pstatus, dstatus;
#if CUOPT_INSTANTIATE_DOUBLE
      pstatus = cuOptSetInitialPrimalSolution(settings, lvls, nvars);
      dstatus = cuOptSetInitialDualSolution(settings, marginals, nconstraints);
#else
      cuopt_float_t *lvlsf = (cuopt_float_t *)malloc(sizeof(cuopt_float_t) * nvars);
      cuopt_float_t *marginalsf = (cuopt_float_t *)malloc(sizeof(cuopt_float_t) * nconstraints);
      for (int i = 0; i < nvars; i++)
        lvlsf[i] = (cuopt_float_t)lvls[i];
      for (int i = 0; i < nconstraints; i++)
        marginalsf[i] = (cuopt_float_t)marginals[i];
      pstatus = cuOptSetInitialPrimalSolution(settings, lvlsf, nvars);
      dstatus = cuOptSetInitialDualSolution(settings, marginalsf, nconstraints);
      free(lvlsf);
      free(marginalsf);
#endif
      free(lvls);
      free(marginals);
      if (pstatus != CUOPT_SUCCESS || dstatus != CUOPT_SUCCESS)
      {
        printOut(gev, "Error setting %s solution for LP.\n", pstatus != CUOPT_SUCCESS ? "primal" : "dual");
        goto DONE;
      }
      printOut(gev, "Initial primal and dual solutions have been set.\n");
    }
  }

  if (!optGetDefinedStr(opt, "prob_read"))
  {
    int num_linear_constraints = 0;

    // cuOpt has no notion of SOS1/SOS2 sets or semi-integer variables. Reject such
    // models explicitly instead of silently feeding the solver an incomplete/undefined
    // variable type array.
    if (gmoGetVarTypeCnt(gmo, gmovar_S1) || gmoGetVarTypeCnt(gmo, gmovar_S2) ||
        gmoGetVarTypeCnt(gmo, gmovar_SI))
    {
      printOut(gev, "ERROR: cuOpt does not support SOS1, SOS2, or semi-integer variables.\n");
      gmoSolveStatSet(gmo, gmoSolveStat_Capability);
      gmoModelStatSet(gmo, gmoModelStat_NoSolutionReturned);
      goto UNLOAD;
    }

    // cuOpt has no conic (=C=), external (=X=) or logic (=B=) constraints
    if (gmoGetEquTypeCnt(gmo, gmoequ_C) || gmoGetEquTypeCnt(gmo, gmoequ_X) || gmoGetEquTypeCnt(gmo, gmoequ_B))
    {
      printOut(gev, "ERROR: cuOpt does not support conic (=C=), external (=X=), or logic (=B=) equations.\n");
      gmoSolveStatSet(gmo, gmoSolveStat_Capability);
      gmoModelStatSet(gmo, gmoModelStat_NoSolutionReturned);
      goto UNLOAD;
    }

    // cuOpt uses 32-bit indices
    if (gmoNZ64(gmo) > INT32_MAX || (is_qcp && gmoMaxQNZ64(gmo) > INT32_MAX))
    {
      printOut(gev, "ERROR: cuOpt does not support models with more than 2^31 (quadratic) nonzeros.\n");
      gmoSolveStatSet(gmo, gmoSolveStat_Capability);
      gmoModelStatSet(gmo, gmoModelStat_NoSolutionReturned);
      goto UNLOAD;
    }

    gams2cuopt_row = malloc(num_constraints * sizeof(int));
    row_qnz = calloc(num_constraints > 0 ? num_constraints : 1, sizeof(int));
    orig_rhs = malloc(num_constraints * sizeof(cuopt_float_t));
    orig_sense = malloc(num_constraints * sizeof(char));

    // Classify rows (linear, quadratic, general nonlinear) as GMO sees them in Q-mode.
    // cuOpt can only take linear and quadratic rows and objectives.
    if (is_qcp)
    {
      if (gmoGetObjOrder(gmo) == gmoorder_NL)
      {
        printOut(gev, "ERROR: The objective is not quadratic (or the quadratic information could not be extracted).\n");
        gmoSolveStatSet(gmo, gmoSolveStat_Capability);
        gmoModelStatSet(gmo, gmoModelStat_NoSolutionReturned);
        goto UNLOAD;
      }
      for (int i = 0; i < num_constraints; i++)
      {
        int order = gmoGetEquOrderOne(gmo, i);
        if (order == gmoorder_Q)
          row_qnz[i] = gmoGetRowQNZOne(gmo, i);
        if (order == gmoorder_NL || order == gmoorder_ERR || row_qnz[i] < 0)
        {
          printOut(gev, "ERROR: Row %d is not quadratic (or the quadratic information could not be extracted).\n", i + 1);
          gmoSolveStatSet(gmo, gmoSolveStat_Capability);
          gmoModelStatSet(gmo, gmoModelStat_NoSolutionReturned);
          goto UNLOAD;
        }
      }
    }

    // Linear rows first (in GAMS order), then the quadratic rows
    for (int i = 0; i < num_constraints; i++)
      if (row_qnz[i] == 0)
        gams2cuopt_row[i] = num_linear_constraints++;
    int q_idx = num_linear_constraints;
    for (int i = 0; i < num_constraints; i++)
      if (row_qnz[i] > 0)
      {
        gams2cuopt_row[i] = q_idx++;
        num_quad_constraints++;
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
        printOut(gev, "ERROR: Unsupported row type %d for row %d.\n", temp_equ_types[i], i + 1);
        free(temp_equ_types);
        gmoSolveStatSet(gmo, gmoSolveStat_Capability);
        gmoModelStatSet(gmo, gmoModelStat_NoSolutionReturned);
        goto UNLOAD;
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
        // SOS1/SOS2/semi-integer are rejected above; anything else is unexpected.
        printOut(gev, "ERROR: Unsupported variable type %d for column %d.\n", temp_var_types[j], j);
        variable_types[j] = CUOPT_CONTINUOUS; // avoid passing an uninitialized value to cuOpt
        free(temp_var_types);
        goto DONE;
      }
    }
    free(temp_var_types);

    // cuOpt's MIP solver ignores quadratic objective terms and quadratic constraints
    // (or fails in presolve), so discrete quadratic models must be rejected.
    if (is_qcp)
      obj_qnz = gmoObjQMatNZ(gmo);
    if (has_integer_vars && (num_quad_constraints > 0 || obj_qnz > 0))
    {
      printOut(gev, "ERROR: cuOpt does not support quadratic models with discrete variables (MIQCP/MIQP).\n");
      gmoSolveStatSet(gmo, gmoSolveStat_Capability);
      gmoModelStatSet(gmo, gmoModelStat_NoSolutionReturned);
      goto UNLOAD;
    }

    // requestMarginals=2: marginals are demanded, so fail upfront if cuOpt cannot provide them
    if (gevGetIntOpt(gev, gevRequestMarginals) == 2 && (has_integer_vars || num_quad_constraints > 0))
    {
      printOut(gev, "ERROR: Marginals are demanded (requestMarginals=2), but cuOpt does not provide them for %s.\n",
               has_integer_vars ? "models with discrete variables" : "models with quadratic constraints");
      gmoSolveStatSet(gmo, gmoSolveStat_Capability);
      gmoModelStatSet(gmo, gmoModelStat_NoSolutionReturned);
      goto UNLOAD;
    }

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
      int qnz = row_qnz[i];

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
      if (obj_qnz > 0)
      {
        int *temp_q_row = malloc(obj_qnz * sizeof(int));
        int *temp_q_col = malloc(obj_qnz * sizeof(int));
        double *temp_q_coef = malloc(obj_qnz * sizeof(double));

        // Extract the quadratic objective coefficients from GAMS
        gmoGetObjQMat(gmo, temp_q_row, temp_q_col, temp_q_coef);

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
      int qnz = row_qnz[i];

      if (qnz > 0)
      {
        // cuOpt only supports convex quadratic constraints of type <= or >=
        if (orig_sense[i] == CUOPT_EQUAL)
        {
          char rowname[GMS_SSSIZE];
          if (gmoDictionary(gmo))
            gmoGetEquNameOne(gmo, i, rowname);
          else
            snprintf(rowname, sizeof(rowname), "%d", i + 1);
          printOut(gev, "ERROR: cuOpt does not support quadratic equality constraints (row %s).\n", rowname);
          if (!gmoObjReform(gmo))
            printOut(gev, "       Note: the objective variable could not be eliminated (e.g. because it has bounds or "
                     "appears in other equations), so a quadratic objective definition is kept as an equation.\n");
          gmoSolveStatSet(gmo, gmoSolveStat_Capability);
          gmoModelStatSet(gmo, gmoModelStat_NoSolutionReturned);
          goto UNLOAD;
        }

        int lin_nz = 0, rnlnz = 0, gstatus;
        int *temp_lin_cols = malloc(num_variables * sizeof(int));
        double *temp_lin_vals = malloc(num_variables * sizeof(double));
        gstatus = gmoGetRowSparse(gmo, i, temp_lin_cols, temp_lin_vals, NULL, &lin_nz, &rnlnz);
        if (gstatus)
        {
          printOut(gev, "gmoGetRowSparse %d failed. Status: %d\n", i, gstatus);
          free(temp_lin_cols);
          free(temp_lin_vals);
          goto DONE;
        }

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
        gstatus = gmoGetRowQMat(gmo, i, temp_q_row, temp_q_col, temp_q_coef);
        if (gstatus)
        {
          printOut(gev, "gmoGetRowQMat %d failed. Status: %d\n", i, gstatus);
          free(temp_lin_cols);
          free(temp_lin_vals);
          free(lin_cols);
          free(lin_vals);
          free(temp_q_row);
          free(temp_q_col);
          free(temp_q_coef);
          goto DONE;
        }

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
#if CUOPT_INSTANTIATE_DOUBLE
    status = cuOptAddMIPStart(settings, initial_levels, gmoN(gmo));
#else
    cuopt_float_t *initial_levelsf = malloc(sizeof(cuopt_float_t) * gmoN(gmo));
    for (int i = 0; i < gmoN(gmo); i++)
      initial_levelsf[i] = (cuopt_float_t)initial_levels[i];
    status = cuOptAddMIPStart(settings, initial_levelsf, gmoN(gmo));
    free(initial_levelsf);
#endif
    free(initial_levels);
    if (status != CUOPT_SUCCESS)
    {
      printOut(gev, "Error setting MIP start.\n");
      goto DONE;
    }
    printOut(gev, "MIP start has been set.\n");
  }

  context.tstart = gevTimeJNow(gev);
  // Solve the problem
  status = cuOptSolve(problem, settings, &solution);
  if (status != CUOPT_SUCCESS) {
    char errmsg[1024] = "";
    if (solution)
      cuOptGetErrorString(solution, errmsg, sizeof(errmsg));
    printOut(gev, "Error solving problem (cuOpt status %d): %s\n", status, errmsg);
    if (status == CUOPT_VALIDATION_ERROR) {
      gmoSolveStatSet(gmo, gmoSolveStat_Capability);
      gmoModelStatSet(gmo, gmoModelStat_NoSolutionReturned);
    } else {
      gmoSolveStatSet(gmo, gmoSolveStat_SolverErr);
      gmoModelStatSet(gmo, gmoModelStat_ErrorNoSolution);
    }
    goto UNLOAD;
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
  cuopt_float_t objective_value, solution_bound = GMS_SV_NA;

  status = cuOptGetTerminationStatus(solution, &termination_status);
  if (status != CUOPT_SUCCESS) {
    printOut(gev, "Error getting termination status: %d\n", status);    
    goto DONE;
  }

  status = cuOptGetSolveTime(solution, &solution_time);
  if (status != CUOPT_SUCCESS) {
    printOut(gev, "Error getting solve time: %d\n", status);
    goto DONE;
  }
  gmoSetHeadnTail(gmo, gmoHresused, solution_time);

  int is_mip = has_integer_vars && (gmoModelType(gmo) == gmoProc_mip || gmoModelType(gmo) == gmoProc_miqcp);
  int have_solution = 0;
  int limit_point = 0; // continuous model stopped by an iteration/time limit
  gmoSolveStatSet(gmo, gmoSolveStat_Normal);
  gmoModelStatSet(gmo, gmoModelStat_NoSolutionReturned);
  if (!optGetDefinedStr(opt, "prob_read")) {
    switch (termination_status) {
    case CUOPT_TERMINATION_STATUS_OPTIMAL:
      gmoModelStatSet(gmo, gmoModelStat_OptimalGlobal);
      have_solution = 1;
      break;
    // For infeasible/unbounded problems cuOpt returns no (or an empty) primal solution
    case CUOPT_TERMINATION_STATUS_INFEASIBLE:
      gmoModelStatSet(gmo, is_mip ? gmoModelStat_IntegerInfeasible : gmoModelStat_InfeasibleNoSolution);
      break;
    case CUOPT_TERMINATION_STATUS_UNBOUNDED:
      gmoModelStatSet(gmo, gmoModelStat_UnboundedNoSolution);
      break;
    case CUOPT_TERMINATION_STATUS_UNBOUNDED_OR_INFEASIBLE:
      printOut(gev, "cuOpt proved the model to be infeasible or unbounded.\n");
      gmoModelStatSet(gmo, gmoModelStat_InfeasibleNoSolution);
      break;
    case CUOPT_TERMINATION_STATUS_ITERATION_LIMIT:
      gmoSolveStatSet(gmo, gmoSolveStat_Iteration);
      limit_point = 1;
      break;
    case CUOPT_TERMINATION_STATUS_TIME_LIMIT:
      gmoSolveStatSet(gmo, gmoSolveStat_Resource);
      limit_point = 1;
      break;
    case CUOPT_TERMINATION_STATUS_WORK_LIMIT:
      gmoSolveStatSet(gmo, gmoSolveStat_Resource);
      break;
    case CUOPT_TERMINATION_STATUS_NUMERICAL_ERROR:
      gmoSolveStatSet(gmo, gmoSolveStat_SolverErr);
      break;
    case CUOPT_TERMINATION_STATUS_PRIMAL_FEASIBLE:
      gmoModelStatSet(gmo, is_mip ? gmoModelStat_Integer : gmoModelStat_Feasible);
      have_solution = 1;
      break;
    case CUOPT_TERMINATION_STATUS_FEASIBLE_FOUND:
      gmoModelStatSet(gmo, is_mip ? gmoModelStat_Integer : gmoModelStat_Feasible);
      have_solution = 1;
      if (is_mip)
      {
        // cuOpt reports FeasibleFound (not a limit status) when the MIP stops with an incumbent
        // before the gap tolerances are met, so figure out which limit was hit.
        cuopt_float_t time_limit = CUOPT_INFINITY, work_limit = CUOPT_INFINITY;
        cuopt_int_t node_limit = 0;
        cuOptGetFloatParameter(settings, CUOPT_TIME_LIMIT, &time_limit);
        cuOptGetFloatParameter(settings, CUOPT_WORK_LIMIT, &work_limit);
        cuOptGetIntegerParameter(settings, CUOPT_NODE_LIMIT, &node_limit);
        // cuOpt does not expose the work units spent, so a set work limit is assumed to be the cause
        if (solution_time >= 0.99 * time_limit || work_limit < 1e10)
          gmoSolveStatSet(gmo, gmoSolveStat_Resource);
        else if (node_limit < INT32_MAX)
          gmoSolveStatSet(gmo, gmoSolveStat_Iteration);
        else
          gmoSolveStatSet(gmo, gmoSolveStat_Solver);
      }
      break;
    case CUOPT_TERMINATION_STATUS_CONCURRENT_LIMIT:
    default:
      gmoSolveStatSet(gmo, gmoSolveStat_Solver);
    }

    // A MIP stopped by a limit without incumbent has no point. For continuous models PDLP returns
    // its current iterate (finite objective), while dual simplex and barrier return an all-zero
    // placeholder with a NaN objective.
    limit_point = limit_point && !has_integer_vars && num_quad_constraints == 0;
    if (limit_point)
    {
      cuopt_float_t obj = NAN;
      if (cuOptGetObjectiveValue(solution, &obj) != CUOPT_SUCCESS || !isfinite(obj))
        limit_point = 0;
    }
    if (limit_point)
    {
      // Classify the iterate as feasible (7) or intermediate infeasible (6) by checking each row
      // against cuOpt's primal tolerances (PDLP projects the iterate onto the variable bounds).
      cuopt_float_t abs_tol = 1e-4, rel_tol = 1e-4;
      cuOptGetFloatParameter(settings, CUOPT_ABSOLUTE_PRIMAL_TOLERANCE, &abs_tol);
      cuOptGetFloatParameter(settings, CUOPT_RELATIVE_PRIMAL_TOLERANCE, &rel_tol);
      cuopt_float_t *x = malloc(num_variables * sizeof(cuopt_float_t));
      int feasible = 0;
      if (x && cuOptGetPrimalSolution(solution, x) == CUOPT_SUCCESS)
      {
        feasible = 1;
        for (int r = 0; r < num_constraints && feasible; r++)
        {
          double act = 0.0;
          for (int k = constraint_matrix_row_offsets[r]; k < constraint_matrix_row_offsets[r + 1]; k++)
            act += constraint_matrix_coefficent_values[k] * x[constraint_matrix_column_indices[k]];
          double tol = abs_tol + rel_tol * (rhs[r] < 0 ? -rhs[r] : rhs[r]);
          if ((constraint_sense[r] != CUOPT_GREATER_THAN && act - rhs[r] > tol) ||
              (constraint_sense[r] != CUOPT_LESS_THAN && rhs[r] - act > tol))
            feasible = 0;
        }
      }
      free(x);
      gmoModelStatSet(gmo, feasible ? gmoModelStat_Feasible : gmoModelStat_InfeasibleIntermed);
      have_solution = 1;
    }
  }

  if (have_solution) {
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
      cuopt_float_t mip_gap;
      if (cuOptGetMIPGap(solution, &mip_gap) == CUOPT_SUCCESS)
        gmoSetHeadnTail(gmo, gmoTrelgap, mip_gap);
      gmoSetHeadnTail(gmo, gmoTabsgap, objective_value > solution_bound ? objective_value - solution_bound : solution_bound - objective_value);
    }

    status = cuOptGetPrimalSolution(solution, objective_coefficients); // reuse n-vector
    if (status != CUOPT_SUCCESS) {
      printOut(gev, "Error getting primal solution: %d\n", status);
      goto DONE;
    }

    if(fp_mip_trace)
    {
      double total_elapsed = (gevTimeJNow(gev) - context.tstart) * 3600.0 * 24.0;
      mip_trace_line('E', 0, 1, total_elapsed, objective_value, solution_bound);
    }

    int request_marginals = gevGetIntOpt(gev, gevRequestMarginals);
    cuopt_int_t presolve = 0, dual_postsolve = 0;
    status = cuOptGetIntegerParameter(settings, CUOPT_PRESOLVE, &presolve);
    if (status != CUOPT_SUCCESS)
    {
      printOut(gev, "Error querying presolve option.\n");
      goto DONE;
    }
    status = cuOptGetIntegerParameter(settings, CUOPT_DUAL_POSTSOLVE, &dual_postsolve);
    if (status != CUOPT_SUCCESS)
    {
      printOut(gev, "Error querying dual_postsolve option.\n");
      goto DONE;
    }

    // Marginals are only available for continuous models (LP, RMIP, QP) without quadratic
    // constraints: cuOpt returns NaN duals for problems with quadratic constraints. Duals of an
    // unfinished (limit) iterate are not meaningful, and presolve without dual postsolve drops them.
    // cuOpt already returns duals in the GAMS sign convention (also for maximization problems).
    int have_marginals = 0;
    if (request_marginals && !limit_point && !has_integer_vars && num_quad_constraints == 0 &&
        (!presolve || dual_postsolve))
    {
      cuopt_float_t *raw_duals = malloc(num_constraints * sizeof(cuopt_float_t));
      double *final_duals = malloc(num_constraints * sizeof(double));
      status = cuOptGetDualSolution(solution, raw_duals);
      if (status != CUOPT_SUCCESS)
      {
        printOut(gev, "Error getting dual solution: %d\n", status);
        free(raw_duals);
        free(final_duals);
        goto DONE;
      }
      for (int i = 0; i < num_constraints; i++)
        final_duals[i] = gams2cuopt_row ? raw_duals[gams2cuopt_row[i]] : raw_duals[i];

      cuopt_float_t *reduced_costs = malloc(num_variables * sizeof(cuopt_float_t));
      if (obj_qnz == 0)
      {
        // cuOpt 26.08's PDLP returns all-zero reduced costs, so compute them for LPs from the
        // duals as d = c - A^T y (exact for dual simplex and barrier as well). GMO's own
        // computation (gmoSetSolution2) ignores c when the objective variable is eliminated.
        status = gmoGetObjVector(gmo, reduced_costs, NULL);
        for (int r = 0; r < num_constraints && !status; r++)
          for (int k = constraint_matrix_row_offsets[r]; k < constraint_matrix_row_offsets[r + 1]; k++)
            reduced_costs[constraint_matrix_column_indices[k]] -= constraint_matrix_coefficent_values[k] * raw_duals[r];
      }
      else
        status = cuOptGetReducedCosts(solution, reduced_costs);
      if (status)
      {
        printOut(gev, "Error getting reduced costs: %d\n", status);
        free(raw_duals);
        free(final_duals);
        free(reduced_costs);
        goto DONE;
      }

      // Same as gmoSetSolution, but GMO computes the row levels itself
      gmoSetVarL(gmo, objective_coefficients);
      gmoSetEquM(gmo, final_duals);
      gmoSetVarM(gmo, reduced_costs);
      gmoSetSolutionStatus(gmo, NULL, NULL, NULL, NULL);
      gmoCompleteSolution(gmo);
      free(raw_duals);
      free(final_duals);
      free(reduced_costs);
      have_marginals = 1;
    }
    else
    {
      if (request_marginals == 1 || request_marginals == 2)
      {
        if (limit_point)
          printOut(gev, "WARNING: No marginals are returned for a solve stopped by a limit.\n");
        else if (has_integer_vars)
          printOut(gev, "WARNING: cuOpt link does not currently support continuous subproblem solves for MIP marginals.\n");
        else if (num_quad_constraints > 0)
          printOut(gev, "WARNING: cuOpt does not return dual solutions for models with quadratic constraints.\n");
        else
          printOut(gev, "WARNING: No marginals are returned since presolve is enabled without dual postsolve.\n");
      }
      // sets marginals to NA and computes variable/equation statuses from the levels
      gmoSetSolutionPrimal(gmo, objective_coefficients);
    }

    // requestMarginals=2: marginals are demanded, so report a problem if they are missing
    if (request_marginals == 2 && !have_marginals && gmoSolveStat(gmo) == gmoSolveStat_Normal)
      gmoSolveStatSet(gmo, gmoSolveStat_Solver);
  }
  else if (fp_mip_trace) // no solution
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
UNLOAD:
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
  free(row_qnz);
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
