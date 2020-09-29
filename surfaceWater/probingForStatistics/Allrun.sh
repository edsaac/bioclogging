decomposePar
mpirun -np 4 pimpleWithG -parallel > log
reconstructPar
foamToVTK
rm -r processor*
