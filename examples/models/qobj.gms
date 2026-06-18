$Title Simple Quadratic Optimization Model (QCP)

Variables
    x1  "Decision variable 1"
    x2  "Decision variable 2"
    z   "Objective function value";

Positive Variable x2;

Equations
    obj   "Objective function to minimize"
    eq1   "Constraint x1 + x2 >= 5";

* Define the equations
obj.. z =e= sqr(x1) + 4*sqr(x2) - 8*x1 - 16*x2;

eq1.. x1 + x2 =g= 5;

* Set specific bounds for the variables
x1.lo = 3;
x1.up = 10;
x2.up = 10;

* Define the model including all equations
Model myModel /all/;

* Solve the model using a Quadratically Constrained Program (QCP) solver
Solve myModel using qcp minimizing z;

* Display the final values of the variables
Display x1.l, x2.l, z.l;