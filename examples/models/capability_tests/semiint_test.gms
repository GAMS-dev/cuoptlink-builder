$Title Minimal semi-integer model - should be REJECTED by the cuOpt link

Set i /i1*i2/;

SemiInt Variable x(i) 'either 0 or an integer in [3,8]';
Free Variable z;

x.lo(i) = 3;
x.up(i) = 8;

Equations obj;
obj.. z =e= sum(i, x(i));

Model m /all/;
Solve m maximizing z using mip;

Display x.l, z.l;
