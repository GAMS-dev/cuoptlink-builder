$Title Options and GAMS settings that previously made the link abort

SemiCont Variable s;
Integer Variable n;
Free Variable z;
s.lo = 2; s.up = 8; n.up = 5;
Equations obj, c;
obj.. z =e= s + n;
c..   s + n =g= 1.5;
Model m /all/;
option mip = cuopt;

$onEcho > cuopt.opt
time_limit 100
$offEcho
$onEcho > cuopt.op2
mip_semi_continuous_big_m 1e6
$offEcho

$macro check(label) \
   abort$(m.solvestat <> %solveStat.normalCompletion%) 'unexpected solve status', label, m.solvestat; \
   abort$(m.modelstat <> %modelStat.optimal% and m.modelstat <> %modelStat.integerSolution%) 'unexpected model status', label, m.modelstat; \
   abort$(abs(z.l - 2) > 1e-6) 'wrong objective', label, z.l;

* option time_limit (double option)
m.optfile = 1;
Solve m minimizing z using mip;
check('time_limit')

* option mip_semi_continuous_big_m
m.optfile = 2;
Solve m minimizing z using mip;
check('big_m')

* OptCR beyond cuOpt's maximum of 0.1 (is clamped by the link)
m.optfile = 0;
option optcr = 0.2;
Solve m minimizing z using mip;
check('optcr')

* negative Threads (all but n cores)
option optcr = 1e-4, threads = -2;
Solve m minimizing z using mip;
check('threads')
