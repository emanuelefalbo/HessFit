%nprocs=4
%chk=wat_hard.chk
#P Amber Geom=Connectivity  opt=(calcall,tight,maxstep=7,maxcycles=100) Freq

Title Card Required

0 1
O-OW                  5.10901163    2.31104648    0.00000000
H-HW                  6.06901163    2.31104648    0.00000000
H-HW                  4.78855705    3.21598231    0.00000000

1 2 1.0 3 1.0
2
3

