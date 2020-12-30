decomposePar
mpirun -np 4 filtrationFix -parallel > log 
reconstructPar
rm -r processor*
foamToVTK
rm -r $(foamListTimes | head -n -1)
