$Title GMO handling: =N= rows, requestMarginals=2, non-reformable objective, MIP tail records

Positive Variable x1, x2;
Free Variable z, w;
Integer Variable n;
n.up = 5;
Equations obj, c1, free1, qobj, link, objn, cn;
obj..   z =e= x1 + 2*x2;
c1..    x1 + x2 =g= 1;
free1.. x1 - x2 =n= 5;
objn..  z =e= n;
cn..    n =g= 1.5;
Model mN   /obj, c1, free1/
      mMIP /objn, cn/;
option lp = cuopt, mip = cuopt, qcp = cuopt;

* =N= rows are skipped by gmoSetNRowPerm and completed by GMO (reference: CPLEX)
Solve mN minimizing z using lp;
abort$(mN.modelstat <> %modelStat.optimal%) '=N= model not optimal', mN.modelstat;
abort$(abs(z.l - 1) > 1e-6 or abs(free1.l + 4) > 1e-6) 'wrong =N= solution', z.l, free1.l;
abort$(abs(c1.m - 1) > 1e-6 or abs(x2.m - 1) > 1e-6) 'wrong marginals', c1.m, x2.m;

* MIP reports the best bound; iterations and nodes are not available from cuOpt
Solve mMIP minimizing z using mip;
abort$(mMIP.modelstat <> %modelStat.optimal%) 'MIP not optimal', mMIP.modelstat;
abort$(abs(mMIP.objEst - 2) > 1e-6) 'wrong best bound', mMIP.objEst;

* requestMarginals=2: MIP marginals are not available -> capability problem before the solve
option requestMarginals = 2;
Solve mMIP minimizing z using mip;
abort$(mMIP.solvestat <> %solveStat.capabilityProblems%) 'requestMarginals=2 MIP: unexpected solve status', mMIP.solvestat;

* requestMarginals=2: LP with presolve but without dual postsolve -> "terminated by solver"
$onEcho > cuopt.opt
presolve 1
dual_postsolve 0
$offEcho
mN.optfile = 1;
Solve mN minimizing z using lp;
abort$(mN.solvestat <> %solveStat.terminatedBySolver%) 'requestMarginals=2 LP: unexpected solve status', mN.solvestat;
abort$(mapVal(c1.m) <> mapVal(na)) 'marginals must be NA', c1.m;
mN.optfile = 0;
option requestMarginals = -1;

* QP whose objective variable cannot be eliminated: quadratic objective equation (CPLEX also refuses)
qobj.. z =e= sqr(x1) + sqr(x2);
link.. w =e= 2*z;
z.up = 100;
Model mQ /qobj, c1, link/;
Solve mQ minimizing z using qcp;
abort$(mQ.solvestat <> %solveStat.capabilityProblems%) 'non-reformable QP: unexpected solve status', mQ.solvestat;
