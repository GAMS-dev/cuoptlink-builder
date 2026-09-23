$Title Minimal MIQCP model - should be REJECTED (cuOpt does not support MIQCP)

Integer Variable x1, x2;
Free Variable z;

x1.up = 10; x2.up = 10;

Equations obj, q;
obj.. z =e= 3*x1 + 2*x2;
q..   sqr(x1) + sqr(x2) =l= 20;

Model m /all/;
Solve m maximizing z using miqcp;

Display x1.l, x2.l, z.l;
