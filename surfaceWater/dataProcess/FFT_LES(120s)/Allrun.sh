decomposePar
mpirun -np 8 pimpleWithG -parallel
reconstructPar
foamToVTK
rm -r processor*
