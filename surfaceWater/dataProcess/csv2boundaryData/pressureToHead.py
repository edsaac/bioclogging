'''
HOW TO USE
-----------------
python3 pressureToHead.py
'''

'''
Write here the path to the folder with the files 
extracted from Paraview (a file per time step).
'''
initFolder  = "last2segsAt0.05" 

'''
Put here the column where the field data is located in those files
First column is 1.
'''
colId = int(2)
'''
Write here the end time and time step of the files
at initFolder. A longer endTime will end up in looped 
boundary data
'''
endTime = 0.05 * 45 * 2
deltaTime = 0.05

'''
Parameters to convert from openFoam pressure to hydraulic
head. 
	head = p_rho*toCentimeters/g
'''
g = 9.81 # m/s2
toCentimeters = 100.

'''
Name of the boundary to build
'''
nameOfBoundary = "top"

#####################################################

import numpy as np
import sys
from os import system
from os import listdir
from natsort import realsorted, natsorted

tempFolder1 = "tempFolder1"   #Keep only P
tempFolder2 = "tempFolder2"   #Save heads
tempFolder3 = "tempFolder3"   #Loops data

'''
Converts the files in initFolder to single files of p_rho
and saves them at tempFolder1
# >> Removes header line 
#     ╚ Sort based on X coordinate (Third column)
#        ╚ Keep only the pressure data (Second column)
'''

system('''
initFolder={0}
tempFolder1={1}
columnIdForPressure={2}
rm -rf $tempFolder1; mkdir $tempFolder1
cd $initFolder

for file in *
do
  sed "1d" $file | sort -k3 -g -t"," | cut -d "," -f {2} > ../$tempFolder1/$file 
done
cd ..
'''.format(initFolder,tempFolder1,colId))

# Create time steps
realTimes = np.arange(0,endTime,deltaTime)

# Retrieve list of pressure files
p_rhoFiles = natsorted(listdir(tempFolder1))
system("F={}; rm -rf $F; mkdir $F".format(tempFolder2))


'''
Creates a folder for each time step and writes the pressure
field converted to hydraulic head
'''
i = 0
for f in p_rhoFiles:
		
		f_in = tempFolder1+"/"+str(f)
		timeFolder = round(realTimes[i],2)

		p_rho = np.loadtxt(f_in)
		head = p_rho*toCentimeters/g

		system('''
			mkdir {}/{}
			'''.format(tempFolder2,timeFolder)
			)

		path=str(tempFolder2)+"/"+str(timeFolder)+"/h"

		# This takes care of the nPoints header that each
		# file must hav
		np.savetxt(path,head,\
			header="{}\n(".format(len(head)),\
			footer=")",\
			comments="",\
			fmt='%.8E')
		
		i+=1

'''
If the endTime is greater than the number of files, loops them
'''
listFilesOG = realsorted(listdir(tempFolder2))
nOriginals = len(listFilesOG)

system('''rm -rf {0}; mkdir {0}'''.format(tempFolder3))

while len(realTimes) > i:
	n = i%nOriginals
	timeFolder = round(realTimes[i],2)
	
	copyFromThis = tempFolder2+"/"+listFilesOG[n]
	pathToThis   = str(tempFolder3)+"/"+str(timeFolder)
	
	system('''cp -r {} {}'''.format(copyFromThis,pathToThis))
	i+=1

'''
Assemble into a boundaryData folder
'''

system('''
	rm -rf boundaryData; mkdir boundaryData
	mkdir boundaryData/{0}
	mv {1}/* boundaryData/{0}/
	mv {2}/* boundaryData/{0}/
	'''.format(nameOfBoundary,tempFolder2,tempFolder3)
	)

'''
Extract the list of points coordinates to complete the
boundaryData folder. It choses a random file since XYZ
was sorted before and all files must come from the same nodes.
'''
system('''
	initFolder={}
	tempFolder1={}
	
	cd $initFolder
	
	file=$(ls | shuf -n 1)
	sed -i "1d" $file | sort -k3 -g -t"," | cut -d "," -f 3,4,5 > ../$tempFolder1/points
	
	cd ..

	pointsFile=$tempFolder1/points
	sed -i "s/,/ /g" $pointsFile
	sed -i "s/^/(/g" $pointsFile
	sed -i "s/$/)/g" $pointsFile
	
	WC=$(wc -l $pointsFile | cut -f 1 -d " ")
	sed -i "$ a )" $pointsFile
	sed -i "1 i (" $pointsFile
	sed -i "1 i $WC" $pointsFile
	
	mv $pointsFile boundaryData/top/
	'''.format(initFolder,tempFolder1)
	)

system('''
	rm -rf {}; rm -rf {}; rm -rf {}
	'''.format(tempFolder1,tempFolder2,tempFolder3)
	)