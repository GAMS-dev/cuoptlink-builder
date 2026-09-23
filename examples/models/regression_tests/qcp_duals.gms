$Title QP marginals must be correct; QCQP must not return (NaN) marginals

Positive Variable y1, y2;
Free Variable z;
option qcp = cuopt;

* maximize 4 y1 + y2 - 0.5 (y1^2 + y2^2) s.t. y1 + y2 = 1
* Reference (CPLEX): z = 3.5, e1.m = 3, y.m = (0, -2)
Equations qobj, e1;
qobj.. z =e= 4*y1 + y2 - 0.5*(sqr(y1) + sqr(y2));
e1..   y1 + y2 =e= 1;
Model qp /qobj, e1/;
Solve qp maximizing z using qcp;
abort$(qp.modelstat <> %modelStat.optimal%) 'QP not optimal', qp.modelstat;
abort$(abs(z.l - 3.5) > 1e-4) 'wrong QP objective', z.l;
abort$(abs(e1.m - 3) > 1e-4 or abs(y1.m) > 1e-4 or abs(y2.m + 2) > 1e-4) 'wrong QP marginals', e1.m, y1.m, y2.m;

* maximize y1 + y2 s.t. y1^2 + y2^2 <= 2, y1 <= 0.5
* Reference (CPLEX): z = 1.822876, marginals not available (NA)
Equations lobj, ball, lin;
lobj.. z =e= y1 + y2;
ball.. sqr(y1) + sqr(y2) =l= 2;
lin..  y1 =l= 0.5;
Model qc /lobj, ball, lin/;
Solve qc maximizing z using qcp;
abort$(qc.modelstat <> %modelStat.optimal%) 'QCQP not optimal', qc.modelstat;
abort$(abs(z.l - 1.822876) > 1e-4) 'wrong QCQP objective', z.l;
abort$(mapVal(ball.m) <> mapVal(na) or mapVal(lin.m) <> mapVal(na) or mapVal(y1.m) <> mapVal(na)) 'QCQP must return NA (not NaN) marginals', ball.m, lin.m, y1.m;
