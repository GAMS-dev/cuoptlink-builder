$Title Small Semi-Continuous Variable Example

Sets
   i 'Products' / A, B /;

Parameters
   profit(i)    'Profit per unit'        / A 5, B 4 /
   min_prod(i)  'Minimum batch if active'/ A 3, B 4 /
   max_prod(i)  'Maximum production limit'/ A 8, B 8 /;

SemiCont Variable x(i) 'Production quantity';
Free Variable     z    'Total profit';

x.lo(i) = min_prod(i);
x.up(i) = max_prod(i);

Equations
   obj       'Objective function: maximize total profit'
   capacity  'Total factory capacity limit';

obj..       z =e= sum(i, profit(i) * x(i));

capacity..  sum(i, x(i)) =l= 10;

Model FactoryModel /all/;
Solve FactoryModel using mip maximizing z;
Display x.l, z.l;