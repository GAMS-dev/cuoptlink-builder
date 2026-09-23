$Title Every cuOpt option in optcuopt.def must be accepted by cuOpt (name, type and default value)

* Write an option file that sets each non-link option (usermap 0) to its default value
$call awk '/^\*/ || /^ / || NF < 4 {next} $2 == "group" {exit} $3 == 0 && $2 != "string" {print $1, $4}' "%gams.sysDir%optcuopt.def" > cuopt.opt
$if errorLevel 1 $abort 'could not generate option file'
* ... plus string options and link options, which must be handled by the link itself
$onEcho > cuopt.op2
solution_file cuopt.sol
user_problem_file cuopt_user.mps
miptrace cuopt.mtr
mipstart 1
$offEcho

SemiCont Variable s;
Integer Variable n;
Positive Variable x;
Free Variable z;
s.lo = 2; s.up = 8; n.up = 5;
Equations obj, c, objlp, clp;
obj..   z =e= s + n;
c..     s + n =g= 1.5;
objlp.. z =e= x;
clp..   x =g= 1;
Model mip /obj, c/, lp /objlp, clp/;
option mip = cuopt, lp = cuopt;

Scalar k;
for (k = 1 to 2,
   mip.optfile = k;
   Solve mip minimizing z using mip;
   abort$(mip.solvestat <> %solveStat.normalCompletion%) 'MIP: unexpected solve status', k, mip.solvestat;
   abort$(abs(z.l - 2) > 1e-6) 'MIP: wrong objective', k, z.l;
   lp.optfile = k;
   Solve lp minimizing z using lp;
   abort$(lp.solvestat <> %solveStat.normalCompletion%) 'LP: unexpected solve status', k, lp.solvestat;
   abort$(abs(z.l - 1) > 1e-6) 'LP: wrong objective', k, z.l;
);
execute 'test -s cuopt.sol';
abort$errorLevel 'solution_file was not written';
execute 'test -s cuopt_user.mps';
abort$errorLevel 'user_problem_file was not written';
execute 'test -s cuopt.mtr';
abort$errorLevel 'miptrace file was not written';
