$Title PDLP column reduction and warm starts through LP presolve

* Random feasible LP: min sum(j, c(j)*x(j)) s.t. A x >= b, 0 <= x <= 5
Set i /r1*r3000/, j /c1*c4000/;
Parameter a(i,j), x0(j), b(i), c(j);
option seed = 1;
a(i,j)$(uniform(0,1) < 0.005) = uniform(1,10);
x0(j) = uniform(0,1);
b(i)  = sum(j, a(i,j)*x0(j)) - 1;
c(j)  = uniform(1,10);
Positive Variable x(j);
Free Variable z;
x.up(j) = 5;
Equations obj, con(i);
obj..    z =e= sum(j, c(j)*x(j));
con(i).. sum(j$a(i,j), a(i,j)*x(j)) =g= b(i);
Model m /all/;
option lp = cuopt;

* Column reduction (cuOpt #2067): forced on, automatic for every LP (threshold 0), and off.
* Concurrent mode with large_nnz_threshold 0 also skips barrier and dual simplex.
$onEcho > cuopt.opt
method 1
pdlp_hyper_reduced_matrix 1
$offEcho
$onEcho > cuopt.op2
method 0
large_nnz_threshold 0
$offEcho
$onEcho > cuopt.op3
method 1
large_nnz_threshold -1
pdlp_hyper_reduced_matrix 0
$offEcho

* Warm starts (cuOpt #2024, #2089). The link passes the levels and marginals as initial primal
* and dual solution. PSLP presolve maps them into the presolved space; with Papilo cuOpt skips
* LP presolve instead. Both must still use the warm start.
$onEcho > cuopt.op4
method 1
presolve 2
$offEcho
$onEcho > cuopt.op5
method 1
presolve 1
$offEcho
$onEcho > cuopt.op6
method 1
presolve 2
pdlp_hyper_project_initial_primal 0
$offEcho
$onEcho > cuopt.op7
method 1
presolve 2
pdlp_hyper_project_initial_primal 1
$offEcho

* PDLP only reaches its default tolerance of 1e-4 (relative)
Set run 'option file' /o1*o7/, reduction(run) /o1*o3/, warmstart(run) /o4*o7/;
Scalar zref 'objective of the reference solve', tol / 1e-3 /, k 'option file number';
Parameter iters(run,*) 'PDLP iterations';

$macro checkOptimal(label) \
   abort$(m.solvestat <> %solveStat.normalCompletion%) 'unexpected solve status', label, m.solvestat; \
   abort$(m.modelstat <> %modelStat.optimal%) 'unexpected model status', label, m.modelstat; \
   abort$(abs(z.l - zref) > tol*max(1, abs(zref))) 'wrong objective', label, z.l, zref;

* bratio 1 discards the levels and marginals of the previous solve (cold start)
option bratio = 1;
m.optfile = 3;
Solve m minimizing z using lp;
zref = z.l;
checkOptimal('reference')

loop(run$reduction(run),
   k = ord(run);
   m.optfile = k;
   Solve m minimizing z using lp;
   checkOptimal(k)
   iters(run,'cold') = m.iterusd;
);

* Each warm start (bratio 0 keeps the previous optimal levels and marginals) must need fewer
* iterations than a cold solve with the same options
loop(run$warmstart(run),
   k = ord(run);
   m.optfile = k;
   option bratio = 1;
   Solve m minimizing z using lp;
   checkOptimal(k)
   iters(run,'cold') = m.iterusd;
   option bratio = 0;
   Solve m minimizing z using lp;
   checkOptimal(k)
   iters(run,'warm') = m.iterusd;
   abort$(iters(run,'warm') >= iters(run,'cold')) 'warm start was not used', k, iters;
);
display iters;
