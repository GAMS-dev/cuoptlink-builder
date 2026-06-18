Variables
    x   'decision variable x'
    y   'decision variable y'
    obj 'objective variable'
;

Equations
    objective    'minimize x + y'
    constraint1  'linear lower bound'
    constraint2  'quadratic upper bound'
;

objective..    obj =e= x + y;

constraint1..  x + y =g= -5;

* Note: You can use x**2 or sqr(x); QCP solvers accept both perfectly.
constraint2..  2*sqr(x) + 2*x*y + 2*sqr(y) =l= 6;

Model quadratic_problem /all/;

* Simply change "nlp" to "qcp"
Solve quadratic_problem minimizing obj using qcp;

Display x.l, y.l, obj.l;