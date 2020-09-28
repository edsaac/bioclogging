decomposePar
mpirun -np 12 pimpleWithG -parallel > log
reconstructPar
foamToVTK
rm -r processor*
