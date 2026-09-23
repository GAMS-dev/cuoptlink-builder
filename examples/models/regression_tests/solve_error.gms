$Title A failing cuOpt solve must report the solver failure to GAMS

* Requesting 2 GPUs (concurrent multi-GPU mode) on a single-GPU machine makes cuOptSolve fail
* with a RuntimeError. On machines with 2+ GPUs the solve succeeds and the test is inconclusive.
* (presolve is disabled, since it would solve this tiny LP before the GPU check)
Positive Variable x;
Free Variable z;
Equations obj, c;
obj.. z =e= x;
c..   x =g= 1;
Model m /all/;
option lp = cuopt;

$onEcho > cuopt.opt
num_gpus 2
presolve 0
$offEcho
m.optfile = 1;
Solve m minimizing z using lp;
if(m.modelstat = %modelStat.optimal%,
   put_utility 'log' / 'NOTE: solve succeeded, this machine has more than one GPU - test inconclusive';
else
   abort$(m.solvestat <> %solveStat.solverFailure%) 'unexpected solve status', m.solvestat;
   abort$(m.modelstat <> %modelStat.errorNoSolution%) 'unexpected model status', m.modelstat;
);
