%nprocs=2
%chk=amber.chk
#P Amber=(print,hardfirst) Geom=Connectivity   Freq=intmodes

Title Card Required

0 1
O-OW--0.5054                5.10901163    2.31104648    0.00000000
H-HW-+0.2527                6.06901163    2.31104648    0.00000000
H-HW-+0.2527                4.78855705    3.21598231    0.00000000

1 2 1.0 3 1.0
2
3
! Amber FF master function
NonBon 3 1 0 0 0.000 0.000 0.500 0.000 0.000 0.833
!
!HrmStr1 OW HW  1.000   0.9606501175      
!                                                     
!HrmBnd1 HW OW HW   1.000 104.9002751576
!
!VDW OW 1.7683  0.1520
!VDW HW 0.0000  0.0000

