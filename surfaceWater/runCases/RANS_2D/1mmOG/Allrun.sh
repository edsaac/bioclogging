decomposePar
mpirun -np 4 simpleWithG -parallel > log 
reconstructPar
foamToVTK
rm -r processor*
