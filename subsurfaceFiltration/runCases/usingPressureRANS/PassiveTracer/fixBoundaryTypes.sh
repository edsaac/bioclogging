for BOUNDARY in "back" "front"
do

	foamDictionary constant/polyMesh/boundary -entry entry0.$BOUNDARY.type -set empty
	foamDictionary constant/polyMesh/boundary -entry entry0.$BOUNDARY.physicalType -remove

done

for BOUNDARY in "right" "left"
do

	foamDictionary constant/polyMesh/boundary -entry entry0.$BOUNDARY.type -set cyclicAMI
	foamDictionary constant/polyMesh/boundary -entry entry0.$BOUNDARY.physicalType -remove
	foamDictionary constant/polyMesh/boundary -entry entry0.$BOUNDARY.transform -set translational	

	if [ $BOUNDARY = "right" ]
	then
		flipTo="left"
		vector="(-15 0 0)"
	elif [ $BOUNDARY = "left" ]
	then
		flipTo="right"
		vector="(15 0 0)"
	fi

	foamDictionary constant/polyMesh/boundary -entry entry0.$BOUNDARY.neighbourPatch -set $flipTo
	foamDictionary constant/polyMesh/boundary -entry entry0.$BOUNDARY.separationVector -set "$vector"

done