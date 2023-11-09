%nprocs=8
%chk=wat_qm.chk
#P B3LYP/def2TZVPP Geom=Connectivity  opt=(calcall,tight,maxstep=7,maxcycles=100) Freq=hinderedrotor

Title Card Required

0 1
 O                  5.10901163    2.31104648    0.00000000
 H                  6.06901163    2.31104648    0.00000000
 H                  4.78855705    3.21598231    0.00000000

 1 2 1.0 3 1.0
 2
 3

