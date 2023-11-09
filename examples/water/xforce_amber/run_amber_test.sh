#!/bin/bash

mpirun -np 4 $AMBERHOME/bin/sander.MPI -O -i 01_Min.in -o 01_Min.out -p wat.prmtop -c wat.inpcrd -r 01_Min.ncrst \
-inf 01_Min.mdinfo  ;

mpirun -np 4 $AMBERHOME/bin/sander.MPI -O -i 02_Heat.in -o 02_Heat.out -p wat.prmtop -c 01_Min.ncrst \
-r 02_Heat.ncrst -x 02_Heat.nc -inf 02_Heat.mdinfo ;
#
mpirun -np 4 $AMBERHOME/bin/sander.MPI -O -i 03_Prod.in -o 03_Prod.out -p wat.prmtop -c 02_Heat.ncrst \
-r 03_Prod.ncrst -x 03_Prod.nc -inf 03_Prod.info & 

exit 
