#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <errno.h>
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

int
main (int argc, char *argv[])
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

  cuOptOptimizationProblem problem = NULL;
  cuOptSolverSettings settings = NULL;
  cuOptSolution solution = NULL;

  cuopt_int_t num_variables = gmoN(gmo);
  cuopt_int_t num_constraints = gmoM(gmo);
  cuopt_int_t nnz = gmoNZ(gmo);

  cuopt_int_t* constraint_matrix_row_offsets=NULL;
  cuopt_int_t* constraint_matrix_column_indices=NULL;
  cuopt_float_t* constraint_matrix_coefficent_values=NULL;
  cuopt_float_t* objective_coefficients=NULL;
  cuopt_float_t* rhs=NULL;
  cuopt_float_t* lower_bounds=NULL;
  cuopt_float_t* upper_bounds=NULL;
  char* constraint_sense=NULL;
  char* variable_types=NULL;

  // Create solver settings
  status = cuOptCreateSolverSettings(&settings);
  if (status != CUOPT_SUCCESS) {
    printOut(gev, "Error creating solver settings: %d\n", status);
    goto DONE;
  }

  // Set solver parameters with GAMS options
  if (gevGetIntOpt(gev, gevThreadsRaw) != 0) {
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
  if (gmoModelType(gmo) == gmoProc_mip) {
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
    } else {
      status = cuOptSetFloatParameter(settings, optname, dval);
      if (status != CUOPT_SUCCESS) {
        printOut(gev, "Error setting float option >%s<: %d\n", optname, status);
        goto DONE;
      }
    }
  }
  
  if (!optGetDefinedStr(opt, "prob_read")) {
    constraint_matrix_row_offsets = malloc((num_constraints+1)*sizeof(cuopt_int_t));
    constraint_matrix_column_indices = malloc(nnz*sizeof(cuopt_int_t));
    constraint_matrix_coefficent_values = malloc(nnz*sizeof(cuopt_float_t));
    objective_coefficients = malloc((num_variables)*sizeof(cuopt_float_t));
    rhs = malloc((num_constraints)*sizeof(cuopt_float_t));
    lower_bounds = malloc((num_variables)*sizeof(cuopt_float_t));
    upper_bounds = malloc((num_variables)*sizeof(cuopt_float_t));
    constraint_sense = malloc((num_constraints)*sizeof(char));
    variable_types = malloc((num_variables)*sizeof(char));
  
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
  
    status = gmoGetEquType(gmo, constraint_matrix_row_offsets); // use constraint_matrix_row_offsets as temp array
    if (status) {
      printOut(gev, "gmoGetEquType failed. Status: %d\n", status);
      goto DONE;
    }

    for (int i=0; i<num_constraints; i++) {
      switch (constraint_matrix_row_offsets[i]) {
         case gmoequ_E: constraint_sense[i] = CUOPT_EQUAL; break;
         case gmoequ_L: constraint_sense[i] = CUOPT_LESS_THAN; break;
         case gmoequ_G: constraint_sense[i] = CUOPT_GREATER_THAN; break;
         default: printOut(gev, "Known row type %d\n", constraint_matrix_row_offsets[i]);
      }
    }

    status = gmoGetVarType(gmo, constraint_matrix_column_indices); // use constraint_matrix_column_indices as temp array
    if (status) {
      printOut(gev, "gmoGetVarType failed. Status: %d\n", status);
      goto DONE;
    }

    for (int j=0; j<num_variables; j++) {
      switch (constraint_matrix_column_indices[j]) {
         case gmovar_X: variable_types[j] = CUOPT_CONTINUOUS; break;
         case gmovar_B:
         case gmovar_I: variable_types[j] = CUOPT_INTEGER; break;
         default: printOut(gev, "Known variable type %d\n", constraint_matrix_column_indices[j]);
      }
    }

    status = gmoGetVarLower(gmo, lower_bounds);
    if (status) {
      printOut(gev, "gmoGetVarLower failed. Status: %d\n", status);
      goto DONE;
    }

    status = gmoGetVarUpper(gmo, upper_bounds);
    if (status) {
      printOut(gev, "gmoGetVarUpper failed. Status: %d\n", status);
      goto DONE;
    }

    nnz = 0;
    for (int i=0; i<num_constraints; i++) {
      int rnz=0, rnlnz=0;
      constraint_matrix_row_offsets[i] = nnz;
      status = gmoGetRowSparse(gmo, i, constraint_matrix_column_indices + nnz, constraint_matrix_coefficent_values + nnz, NULL, &rnz, &rnlnz);
      if (status) {
        printOut(gev, "gmoGetRowSparse %d failed. Status: %d\n", i, status);
        goto DONE;
      }
      nnz += rnz;
    }
    constraint_matrix_row_offsets[num_constraints] = nnz;

    status = gmoGetObjVector(gmo, objective_coefficients, NULL);
    if (status) {
      printOut(gev, "gmoGetObjVector failed. Status: %d\n", status);
      goto DONE;
    }

    status = gmoGetRhs(gmo, rhs);
    if (status) {
      printOut(gev, "gmoGetRhs failed. Status: %d\n", status);
      goto DONE;
    }

    status = cuOptCreateProblem(
      num_constraints,
      num_variables,
      (gmoSense(gmo)==gmoObj_Min)?CUOPT_MINIMIZE:CUOPT_MAXIMIZE,
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
      &problem
    );

    if (status != CUOPT_SUCCESS) {
      printOut(gev, "Error creating problem from GAMS model: %d\n", status);
      goto DONE;
    }
  } else {
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
  cuopt_float_t objective_value;  

  status = cuOptGetTerminationStatus(solution, &termination_status);
  if (status != CUOPT_SUCCESS) {
    printOut(gev, "Error getting termination status: %d\n", status);
    goto DONE;
  }

  gmoSolveStatSet(gmo, gmoSolveStat_Normal);
  gmoModelStatSet(gmo, gmoModelStat_NoSolutionReturned);
  if (!optGetDefinedStr(opt, "prob_read")) {
    switch (termination_status) {
      case CUOPT_TERIMINATION_STATUS_OPTIMAL:
        gmoModelStatSet(gmo, gmoModelStat_OptimalGlobal);
        break;
      case CUOPT_TERIMINATION_STATUS_INFEASIBLE:
        gmoModelStatSet(gmo, gmoModelStat_InfeasibleGlobal);
        break;
      case CUOPT_TERIMINATION_STATUS_UNBOUNDED:
        gmoModelStatSet(gmo, gmoModelStat_Unbounded);
        break;
      case CUOPT_TERIMINATION_STATUS_ITERATION_LIMIT:
        gmoSolveStatSet(gmo, gmoSolveStat_Iteration);
        break;
      case CUOPT_TERIMINATION_STATUS_TIME_LIMIT:
        gmoSolveStatSet(gmo, gmoSolveStat_Resource);
        break;
      case CUOPT_TERIMINATION_STATUS_NUMERICAL_ERROR:
        gmoSolveStatSet(gmo, gmoSolveStat_SolverErr);
        break;
      case CUOPT_TERIMINATION_STATUS_PRIMAL_FEASIBLE:
        if (gmoModelType(gmo) == gmoProc_mip) {
          gmoModelStatSet(gmo, gmoModelStat_Integer);
        } else {
          gmoModelStatSet(gmo, gmoModelStat_Feasible);
        }
        break;
      case CUOPT_TERIMINATION_STATUS_FEASIBLE_FOUND:
        if (gmoModelType(gmo) == gmoProc_mip) {
          gmoModelStatSet(gmo, gmoModelStat_Integer);
        } else {
          gmoModelStatSet(gmo, gmoModelStat_Feasible);
        }
        break;
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

    if (gmoModelType(gmo) == gmoProc_mip) {
      status = cuOptGetSolutionBound(solution, &objective_value);
      if (status != CUOPT_SUCCESS) {
        printOut(gev, "Error getting solution bound: %d\n", status);
        goto DONE;
      }
      gmoSetHeadnTail(gmo, gmoTmipbest, objective_value);
    }

    status = cuOptGetPrimalSolution(solution, objective_coefficients); // reuse n-vector
    if (status != CUOPT_SUCCESS) {
      printOut(gev, "Error getting primal solution: %d\n", status);
      goto DONE;
    }
    gmoSetVarL(gmo, objective_coefficients);

    int presolve;
    cuOptGetIntegerParameter(settings, "presolve", &presolve);
    if (gmoModelType(gmo) != gmoProc_mip && !presolve) {
      status = cuOptGetReducedCosts(solution, objective_coefficients); // reuse n-vector
      if (status != CUOPT_SUCCESS) {
        printOut(gev, "Error getting reduced cost: %d\n", status);
        goto DONE;
      }
      if (gmoSense(gmo) == gmoObj_Max) { // patch duals for max problem
        for (int j=0; j<num_variables; j++) {
          objective_coefficients[j] *= -1.0;
        }
      }
      gmoSetVarM(gmo, objective_coefficients);

      status = cuOptGetDualSolution(solution, rhs); // reuse m-vector
      if (status != CUOPT_SUCCESS) {
        printOut(gev, "Error getting reduced cost: %d\n", status);
        goto DONE;
      }
      if (gmoSense(gmo) == gmoObj_Max) { // patch duals for max problem
        for (int i=0; i<num_constraints; i++) {
          rhs[i] *= -1.0;
        }
      }
      gmoSetEquM(gmo, rhs);
    } else {
      gmoSetHeadnTail(gmo, gmoHmarginals, 0.0); // no duals for mip      
    }
    gmoCompleteSolution(gmo);
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

GAMSDONE:
  gmoFree(&gmo);
  gevFree(&gev);
  optFree(&opt);

  return 0;

} /* main */

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
    case CUOPT_TERIMINATION_STATUS_OPTIMAL:
      return "Optimal";
    case CUOPT_TERIMINATION_STATUS_INFEASIBLE:
      return "Infeasible";
    case CUOPT_TERIMINATION_STATUS_UNBOUNDED:
      return "Unbounded";
    case CUOPT_TERIMINATION_STATUS_ITERATION_LIMIT:
      return "Iteration limit";
    case CUOPT_TERIMINATION_STATUS_TIME_LIMIT:
      return "Time limit";
    case CUOPT_TERIMINATION_STATUS_NUMERICAL_ERROR:
      return "Numerical error";
    case CUOPT_TERIMINATION_STATUS_PRIMAL_FEASIBLE:
      return "Primal feasible";
    case CUOPT_TERIMINATION_STATUS_FEASIBLE_FOUND:
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
