NPROC=$1
decomposePar
mpirun -np $NPROC filtrationFix -parallel > log 
reconstructPar
rm -r processor*
foamToVTK -noZero
#foamListTimes -rm

