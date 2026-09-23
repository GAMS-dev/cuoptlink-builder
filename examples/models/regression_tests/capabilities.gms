$Title Unsupported model features must be rejected with a capability error

Set i /i1*i3/;
Free Variable z;
Equations objS1, objS2, objSI, limSI;

* --- SOS1 variables
SOS1 Variable s1(i);
s1.up(i) = 10;
objS1.. z =e= sum(i, s1(i));
Model mS1 /objS1/;

* --- SOS2 variables
SOS2 Variable s2(i);
s2.up(i) = 10;
objS2.. z =e= sum(i, s2(i));
Model mS2 /objS2/;

* --- semi-integer variables
SemiInt Variable si;
si.lo = 2; si.up = 10;
objSI.. z =e= si;
limSI.. si =g= 1;
Model mSI /objSI, limSI/;

option mip = cuopt;
Solve mS1 maximizing z using mip;
abort$(mS1.solvestat <> %solveStat.capabilityProblems%) 'SOS1: unexpected solve status', mS1.solvestat;
abort$(mS1.modelstat <> %modelStat.noSolutionReturned%) 'SOS1: unexpected model status', mS1.modelstat;
Solve mS2 maximizing z using mip;
abort$(mS2.solvestat <> %solveStat.capabilityProblems%) 'SOS2: unexpected solve status', mS2.solvestat;
abort$(mS2.modelstat <> %modelStat.noSolutionReturned%) 'SOS2: unexpected model status', mS2.modelstat;
Solve mSI minimizing z using mip;
abort$(mSI.solvestat <> %solveStat.capabilityProblems%) 'semi-integer: unexpected solve status', mSI.solvestat;
abort$(mSI.modelstat <> %modelStat.noSolutionReturned%) 'semi-integer: unexpected model status', mSI.modelstat;

* --- MIQCP is not in cuOpt's model types (gamsconfig), so GAMS itself must refuse the solver.
* This aborts the GAMS job, hence it is run as a separate job whose return code is checked.
$onEcho > miqcp_sub.gms
Integer Variable x1, x2;
Free Variable z;
x1.up = 10; x2.up = 10;
Equations obj, q;
obj.. z =e= 3*x1 + 2*x2;
q..   sqr(x1) + sqr(x2) =l= 20;
Model m /all/;
Solve m maximizing z using miqcp;
$offEcho
execute '"%gams.sysDir%gams" miqcp_sub.gms miqcp=cuopt lo=2 lf=miqcp_sub.log o=miqcp_sub.lst';
abort$(errorLevel = 0) 'GAMS must refuse cuopt as MIQCP solver';
execute 'grep -q "cannot solve MIQCPs" miqcp_sub.log';
abort$errorLevel 'MIQCP refusal message not found in the log';
