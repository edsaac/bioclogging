decomposePar
mpirun -np 24 filtrationFix -parallel > log 
reconstructPar
rm -r processor*
foamToVTK
#foamListTimes -rm

