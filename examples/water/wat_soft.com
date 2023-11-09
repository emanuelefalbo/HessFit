%nprocs=4
%chk=wat_soft.chk
#P Amber=SoftFirst Geom=Connectivity  opt=(calcall,tight,maxstep=7,maxcycles=100) Freq

Title Card Required

0 1
O-OW                  5.10901163    2.31104648    0.00000000
H-HW                  6.06901163    2.31104648    0.00000000
H-HW                  4.78855705    3.21598231    0.00000000

1 2 1.0 3 1.0
2
3
! Amber FF master function
NonBon 3 1 0 0 0.000 0.000 0.500 0.000 0.000 0.833
!
HrmStr1 OW HW  583.7865152030   0.9606501175      
!                                                     
HrmBnd1 HW O HW   40.7119029217 104.9002751576

