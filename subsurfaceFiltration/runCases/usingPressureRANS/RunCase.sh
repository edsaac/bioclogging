runCase(){
	cd $1
	decomposePar
	mpirun -np 4 filtrationFix -parallel > /dev/null 2>&1
	reconstructPar
	foamToVTK
	rm -r processor*
	cd ..
}

runCase $1
