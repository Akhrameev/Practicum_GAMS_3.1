$title Performance (task 3)
$ontext
09.2014 Akhrameev Pavel
group 513 CMC MSU
Practicum
Task 3
Performance
$offtext

sets
h / 0 * 100 /
scalars
alpha    /20.0/
rho      /30.0/
x_0_0    /10.0/
x_0_1    /100.0/
dx_0_0   /10.0/
dx_0_1   /12.0/
x_T_0    /700.0/
x_T_1    /900.0/
dx_T_0   /150.0/
dx_T_1   /20.0/

zero     /0.0/

deltaH
maxH;
maxH = card(h) - 1;
* range form 0 to 1 devided by number of h dots + 1 (1/max(h))
deltaH = 1/maxH;

variables
x1(h)
x2(h)
x3(h)
x4(h)
u1(h)
u2(h)

diffSum(h);

equations
eq_x1(h)
eq_x2(h)
eq_x3(h)
eq_x4(h)
eq_u(h)
eq_diffSum(h);

eq_x1(h-1).. x1(h) =e= x1(h-1) + deltaH*x3(h-1);
eq_x2(h-1).. x2(h) =e= x2(h-1) + deltaH*x4(h-1);
eq_x3(h-1).. x3(h) =e= x3(h-1) + deltaH*(-alpha * eDist(x3(h-1),x4(h-1)) *
                         x3(h-1) + rho*u1(h-1) );
eq_x4(h-1).. x4(h) =e= x4(h-1) + deltaH*(-alpha * eDist(x3(h-1),x4(h-1)) *
                         x4(h-1) + rho*u2(h-1) );

eq_u(h-1)..  eDist(u1(h-1),u2(h-1)) =l= 1;

eq_diffSum(h-1).. diffSum(h) =e= diffSum(h-1) + sqr(ord(h)) *
         (1 + eDist(x1(h), x1(h-1), x2(h), x2(h-1)) +
              eDist(x3(h), x3(h-1), x4(h), x4(h-1)));

x1.fx(h)$(ord(h) = 1) =  x_0_0;
x2.fx(h)$(ord(h) = 1) =  x_0_1;
x3.fx(h)$(ord(h) = 1) = dx_0_0;
x4.fx(h)$(ord(h) = 1) = dx_0_1;

x1.fx(h)$(ord(h) = card(h)) =  x_T_0;
x2.fx(h)$(ord(h) = card(h)) =  x_T_1;
x3.fx(h)$(ord(h) = card(h)) = dx_T_0;
x4.fx(h)$(ord(h) = card(h)) = dx_T_1;

diffSum.fx(h)$(ord(h) = 1) = zero;

model speed /all/;

solve speed using dnlp minimizing diffSum(h)$(ord(h) = card(h));

Parameter PLOT_1 data for plotter;
PLOT_1("x1",h,"y")=x1.l(h);
PLOT_1("h",h,"x")=ord(h);
$libinclude gnuplotxyz PLOT_1 x y