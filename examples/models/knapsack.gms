set i /i1*i4/;
parameter
    weights(i) /i1 3, i2 4, i3 8, i4 10/
    utilities(i) /i1 3, i2 4, i3 8, i4 10/;
scalar capacity /3/;
binary variable x(i);
free variable z;
equation obj, cap;
obj .. z =e= sum(i, utilities(i)*x(i));
cap .. sum(i, x(i)*weights(i)) =l= capacity;
x.l(i) = 0;
x.l('i1') = 1;
model m /all/;
*m.optfile=1;
solve m maximizing z using mip;