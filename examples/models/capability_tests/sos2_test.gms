$Title Minimal SOS2 model - should be REJECTED by the cuOpt link

Set i /i1*i3/;

SOS2 Variable x(i) 'at most two ADJACENT of these may be nonzero';
Free Variable z;

x.up(i) = 10;

Equations obj;
obj.. z =e= sum(i, x(i));

Model m /all/;
Solve m maximizing z using mip;

Display x.l, z.l;
