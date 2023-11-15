%mem=1gb
%chk=mmH2O2_2.chk
#p amber=softfirst geom=connectivity nosymm
iop(4/33=3,7/33=1)
freq=intmodes

MM

0 1
O-oh--0.410452   -0.718633164030   -0.118472295063   -0.054617618503
H-ho-0.410452   -1.023538245240    0.665457463989    0.436707052760
O-oh--0.410452    0.718637032315    0.118468958072   -0.054573365001
H-ho-0.410452    1.023507272500   -0.665430766999    0.436820814218

 1 2 1.0 3 1.0
 2
 3 4 1.0
 4

AmbTrs ho oh oh ho 0 0 0 0 0.0 1.0 0.0 0.0 1.0
HrmBnd1 ho oh oh 1.0 100.2486
HrmStr1 ho oh 1.0 0.97412
HrmStr1 oh oh 1.0 1.45667
Nonbon 3 1 0 0 0.0 0.0 0.5 0.0 0.0 -1.2
VDW ho  0.0000  0.0000
VDW oh  1.7210  0.2104

