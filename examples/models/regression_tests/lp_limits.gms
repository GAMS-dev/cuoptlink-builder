$Title LPs stopped by an iteration or time limit return the PDLP iterate (if there is one)

* Random feasible LP: max sum(j, c(j)*x(j)) s.t. A x <= b, 0 <= x <= 5
Set i /r1*r3000/, j /c1*c4000/;
Parameter a(i,j), x0(j), b(i), c(j);
option seed = 1;
a(i,j)$(uniform(0,1) < 0.005) = uniform(1,10);
x0(j) = uniform(0,1);
b(i)  = sum(j, a(i,j)*x0(j)) + 1;
c(j)  = uniform(1,10);
Positive Variable x(j);
Free Variable z;
x.up(j) = 5;
Equations obj, con(i);
obj..    z =e= sum(j, c(j)*x(j));
con(i).. sum(j$a(i,j), a(i,j)*x(j)) =l= b(i);
Model m /all/;
option lp = cuopt;

* optfile 1..3 = PDLP, dual simplex, barrier
$onEcho > cuopt.opt
method 1
$offEcho
$onEcho > cuopt.op2
method 2
$offEcho
$onEcho > cuopt.op3
method 3
$offEcho

* PDLP: iteration limit and time limit both return the current iterate
Set lim 'PDLP limit scenarios' /iterlim, reslim/;
Parameter
   expSS(lim) 'expected solve status' / iterlim %solveStat.iterationInterrupt%, reslim %solveStat.resourceInterrupt% /
   report(lim,*);
Scalar maxviol;
m.optfile = 1;
loop(lim,
   m.iterlim = 2e9; m.reslim = 1e10;
   if(sameAs(lim,'iterlim'), m.iterlim = 20; else m.reslim = 0.05;);
   x.l(j) = 0;
   Solve m maximizing z using lp;
   maxviol = smax(i, con.l(i) - b(i));
   report(lim,'solvestat') = m.solvestat;
   report(lim,'modelstat') = m.modelstat;
   report(lim,'z') = z.l;
   report(lim,'maxviol') = maxviol;
   abort$(m.solvestat <> expSS(lim)) 'unexpected solve status', lim, m.solvestat;
   abort$(m.modelstat <> %modelStat.feasibleSolution% and m.modelstat <> %modelStat.intermediateInfeasible%) 'no point returned', lim, m.modelstat;
   abort$(mapVal(z.l) <> 0) 'objective must be a finite number', lim, z.l;
   abort$(mapVal(con.m('r1')) <> mapVal(na)) 'limit point must not have marginals', lim;
   abort$(m.modelstat = %modelStat.feasibleSolution% and maxviol > 1e-4 + 1e-4*smax(i, b(i))) 'point marked feasible but violated', lim, maxviol;
);
display report;

* Dual simplex and barrier have no usable point at a limit
m.reslim = 1e10; m.iterlim = 5;
Scalar k;
for (k = 2 to 3,
   m.optfile = k;
   Solve m maximizing z using lp;
   abort$(m.solvestat <> %solveStat.iterationInterrupt%) 'unexpected solve status', k, m.solvestat;
   abort$(m.modelstat <> %modelStat.noSolutionReturned%) 'placeholder point must not be returned', k, m.modelstat;
);
