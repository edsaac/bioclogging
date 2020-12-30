decomposePar
mpirun -np 4 filtrationFix -parallel > log 
reconstructPar
rm -r processor*
foamToVTK
#foamListTimes -rm

