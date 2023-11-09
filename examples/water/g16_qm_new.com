%nprocs=8
%chk=g16_qm.chk
#P B3LYP/def2TZVPP Geom=Connectivity  opt=(calcall,tight,maxstep=7,maxcycles=100) Freq 

Title Card Required

 0,1
  O                                                 5.109011630000      2.311046480000      0.000000000000
  H                                                 6.069011630000      2.311046480000      0.000000000000
  H                                                 4.788557050000      3.215982310000      0.000000000000

 1 2 1.000 3 1.000
 2 1 1.000
 3 1 1.000
 
