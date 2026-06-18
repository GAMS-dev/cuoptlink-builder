$title Simple MIQCP Example

Set 
    i 'Available items or assets' / 1*3 /;

Alias (i,j);

Parameter
    c(i) 'Linear objective coefficients (e.g., expected return)' 
         / 1 5, 2 8, 3 10 /
    Q(i,j) 'Quadratic constraint coefficients (e.g., covariance matrix)' ;

* Initialize Q as an identity matrix for simplicity (diagonal elements = 1)
Q(i,j)$(ord(i) = ord(j)) = 1;

Variables
    x(i) 'Decision variables (quantity of items to select)'
    z    'Objective variable (total return)';

Integer Variable x;

* It is good practice to bound integer variables
x.up(i) = 5;

Equations
    obj_eq   'Objective function to maximize'
    quad_eq  'Quadratic constraint (e.g., risk limit)'
    lin_eq   'Linear constraint (e.g., total budget or item limit)';

* Maximize linear return
obj_eq..  z =e= sum(i, c(i)*x(i));

* Subject to a quadratic constraint
* Notice the multiplication of two variables: x(i) * x(j)
quad_eq.. sum((i,j), x(i)*Q(i,j)*x(j)) =l= 50;

* Subject to a linear constraint
lin_eq..  sum(i, x(i)) =l= 10;

Model simple_miqcp /all/;

* Specify the model type explicitly as MIQCP
Solve simple_miqcp using miqcp maximizing z;

Display x.l, z.l;