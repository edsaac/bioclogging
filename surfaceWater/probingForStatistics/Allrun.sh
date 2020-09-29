decomposePar
mpirun -np 8 pimpleWithG -parallel > log
reconstructPar
foamToVTK
rm -r processor*
