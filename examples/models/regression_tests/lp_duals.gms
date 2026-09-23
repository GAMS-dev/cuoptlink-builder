$Title Marginals of min/max LPs must match CPLEX for all cuOpt methods

* Reference (CPLEX):
*   max: c.m = (11, 0, 6), x.m = (-1, 0, -2, 0)
*   min: c.m = (0, 0, 0),  x.m = (4, 1, 5, 3)
Set i 'constraints' /c1*c3/, j 'variables' /x1*x4/;
Table A(i,j)
        x1  x2  x3  x4
    c1   1  -1  -1   3
    c2   5   1   3   8
    c3  -1   2   3  -5 ;
Parameter b(i) /c1 1, c2 55, c3 3/, c(j) /x1 4, x2 1, x3 5, x4 3/;
Parameter
    ymax(i) /c1 11, c2 0, c3 6/,         dmax(j) /x1 -1, x2 0, x3 -2, x4 0/
    ymin(i) /c1 0,  c2 0, c3 0/,         dmin(j) /x1 4,  x2 1, x3 5,  x4 3/;

Positive Variable x(j);
Free Variable z;
Equations obj, con(i);
obj..    z =e= sum(j, c(j)*x(j));
con(i).. sum(j, A(i,j)*x(j)) =l= b(i);
Model m /all/;
option lp = cuopt;

* optfile 1..4 = method 0 (concurrent), 1 (PDLP), 2 (dual simplex), 3 (barrier)
$onEcho > cuopt.opt
method 0
$offEcho
$onEcho > cuopt.op2
method 1
$offEcho
$onEcho > cuopt.op3
method 2
$offEcho
$onEcho > cuopt.op4
method 3
$offEcho

* PDLP only reaches its default tolerance of 1e-4 (relative)
Scalar tol / 1e-2 /, k;
for (k = 1 to 4,
   m.optfile = k;
   Solve m maximizing z using lp;
   abort$(m.modelstat <> %modelStat.optimal%) 'max LP not optimal', k, m.modelstat;
   abort$(abs(z.l - 29) > tol) 'wrong max objective', k, z.l;
   abort$(smax(i, abs(con.m(i) - ymax(i))) > tol) 'wrong max LP duals', k, con.m;
   abort$(smax(j, abs(x.m(j) - dmax(j))) > tol) 'wrong max LP reduced costs', k, x.m;

   Solve m minimizing z using lp;
   abort$(m.modelstat <> %modelStat.optimal%) 'min LP not optimal', k, m.modelstat;
   abort$(smax(i, abs(con.m(i) - ymin(i))) > tol) 'wrong min LP duals', k, con.m;
   abort$(smax(j, abs(x.m(j) - dmin(j))) > tol) 'wrong min LP reduced costs', k, x.m;
);
