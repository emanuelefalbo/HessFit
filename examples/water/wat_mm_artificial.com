%nprocs=1
%chk=wat_mm_artificial.chk
#P Amber=(SoftOnly,Print) Geom=Connectivity  Freq=intmodes

Title Card Required

0 1
O-OW--0.406659                 5.10901163    2.31104648    0.00000000
H-HW-+0.203329                 6.06901163    2.31104648    0.00000000
H-HW-+0.203329                 4.78855705    3.21598231    0.00000000

1 2 1.0 3 1.0
2
3
! Amber FF master function
NonBon 3 1 0 0 0.000 0.000 0.500 0.000 0.000 0.833
!
HrmStr1 OW HW  1.00    0.9685382838
!
HrmBnd1 HW OW HW   1.00 102.7321784514

