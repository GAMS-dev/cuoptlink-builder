$Title Solve/model status mapping for infeasible, unbounded, limit and unsupported cases

option lp = cuopt, mip = cuopt, qcp = cuopt;
Free Variable z;

* --- Infeasible LP: model status 19, no solution loaded (levels stay untouched)
Positive Variable x1, x2;
Equations obj1, inf1, inf2;
obj1.. z =e= x1 + x2;
inf1.. x1 + x2 =l= 1;
inf2.. x1 + x2 =g= 2;
Model mInf /obj1, inf1, inf2/;
x1.l = 7;
Solve mInf minimizing z using lp;
abort$(mInf.solvestat <> %solveStat.normalCompletion%) 'infeasible LP: unexpected solve status', mInf.solvestat;
abort$(mInf.modelstat <> %modelStat.infeasibleNoSolution%) 'infeasible LP: unexpected model status', mInf.modelstat;
abort$(x1.l <> 7) 'infeasible LP: levels must not be overwritten', x1.l;

* --- Unbounded LP: model status 18
Equations obj2, lo2;
obj2.. z =e= x1 - x2;
lo2..  x1 - x2 =g= 1;
Model mUnb /obj2, lo2/;
Solve mUnb maximizing z using lp;
abort$(mUnb.solvestat <> %solveStat.normalCompletion%) 'unbounded LP: unexpected solve status', mUnb.solvestat;
abort$(mUnb.modelstat <> %modelStat.unboundedNoSolution% and mUnb.modelstat <> %modelStat.infeasibleNoSolution%)
   'unbounded LP: unexpected model status', mUnb.modelstat;

* --- Infeasible MIP: model status 10
Integer Variable n;
n.up = 10;
Equations obj3, par;
obj3.. z =e= n;
par..  2*n =e= 3;
Model mInfMip /obj3, par/;
Solve mInfMip minimizing z using mip;
abort$(mInfMip.solvestat <> %solveStat.normalCompletion%) 'infeasible MIP: unexpected solve status', mInfMip.solvestat;
abort$(mInfMip.modelstat <> %modelStat.integerInfeasible% and mInfMip.modelstat <> %modelStat.infeasibleNoSolution%)
   'infeasible MIP: unexpected model status', mInfMip.modelstat;

* --- MIP stopped by the time limit with an incumbent: solve status 3, model status 8
* Market split instance (Cornuejols/Dawande): feasible incumbents are trivial, proving
* optimality is very hard for branch-and-bound.
Set i 'rows' /r1*r4/, j 'columns' /c1*c40/;
Parameter a(i,j), d(i);
option seed = 12345;
a(i,j) = uniformInt(0, 99);
d(i) = floor(sum(j, a(i,j)) / 2);
Binary Variable y(j);
Positive Variable sp(i), sm(i);
Equations obj4, split(i);
obj4..     z =e= sum(i, sp(i) + sm(i));
split(i).. sum(j, a(i,j)*y(j)) + sp(i) - sm(i) =e= d(i);
Model mSplit /obj4, split/;
mSplit.reslim = 3;
option optcr = 0;
Solve mSplit minimizing z using mip;
abort$(mSplit.modelstat = %modelStat.optimal%) 'market split unexpectedly solved to optimality; use a harder instance';
abort$(mSplit.solvestat <> %solveStat.resourceInterrupt%) 'time limit: unexpected solve status', mSplit.solvestat;
abort$(mSplit.modelstat <> %modelStat.integerSolution%) 'time limit: unexpected model status', mSplit.modelstat;

* --- Same MIP stopped by the node limit: solve status 2, model status 8
$onEcho > cuopt.opt
node_limit 20
$offEcho
mSplit.reslim = 60;
mSplit.optfile = 1;
Solve mSplit minimizing z using mip;
abort$(mSplit.solvestat <> %solveStat.iterationInterrupt%) 'node limit: unexpected solve status', mSplit.solvestat;
abort$(mSplit.modelstat <> %modelStat.integerSolution%) 'node limit: unexpected model status', mSplit.modelstat;
mSplit.optfile = 0;

* --- Quadratic equality constraint: capability error instead of a generic failure
Variable w1, w2;
Equations obj5, sphere;
obj5..   z =e= w1 + w2;
sphere.. sqr(w1) + sqr(w2) =e= 1;
Model mQeq /obj5, sphere/;
Solve mQeq minimizing z using qcp;
abort$(mQeq.solvestat <> %solveStat.capabilityProblems%) 'quadratic equality: unexpected solve status', mQeq.solvestat;
abort$(mQeq.modelstat <> %modelStat.noSolutionReturned%) 'quadratic equality: unexpected model status', mQeq.modelstat;
