$Title RMIQCP must be solved with its quadratic terms (not linearized)

* Reference (CPLEX): z = 16.124515, x = (3.721, 2.481)
Integer Variable x1, x2;
Free Variable z;
x1.up = 10; x2.up = 10;

Equations obj, q;
obj.. z =e= 3*x1 + 2*x2;
q..   sqr(x1) + sqr(x2) =l= 20;

Model m /all/;
option rmiqcp = cuopt;
Solve m maximizing z using rmiqcp;

abort$(m.solvestat <> %solveStat.normalCompletion%) 'unexpected solve status', m.solvestat;
abort$(m.modelstat <> %modelStat.optimal%) 'unexpected model status', m.modelstat;
abort$(abs(z.l - 16.124515) > 1e-3) 'wrong RMIQCP objective', z.l;
