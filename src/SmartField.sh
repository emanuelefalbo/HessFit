#!/bin/bash

F1=ind_qm.log
F2=ind_qm.fchk
OPT=opt.json
GPATH=$g09root/g09
SM=Smart_harmonic.py
BS=build_4Smart.py
$BS -f1 ${F1} -f2 ${F2} -path ${GPATH}

for f in GauHarm.gjf GauNonBon.gjf
do
#   if [ -f "$f" ]
#   then
#       echo "Existing ${f} file"
#       echo "Using $f file for ${SM}"
#   else
       echo "Executing Gaussian on ${f}"
       ${GPATH}/g09 ${f}
#   fi 
done

for f in GauHarm.chk GauNonBon.chk
do
   echo "Formatchecking $f file"
   $GPATH/formchk -3 $f "${f%.chk}.fchk"
done
#
#$GPATH/g09 GauHarm.gjf
#$GPATH/g09 GauHarm.gjf
echo "Executing ${SM}"
SmartField_harmonic.py $OPT 
${GPATH}/g09 SmartField4gau.gjf
