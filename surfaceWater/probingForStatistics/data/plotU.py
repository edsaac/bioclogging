import numpy as np
import matplotlib.pyplot as plt
from pandas import read_csv
from os import system
from os import listdir
import sys

#FOLDER = "./organizedData/lineAtBot/"
#FOLDER = "./organizedData/lineAtTop/"
FOLDER = "./organizedData/boundaryProbes/"

FILELIST = listdir(FOLDER)
FILELIST.sort()
#print(FILELIST)

plt.figure(figsize=(10,10),facecolor="white")
ax1 = plt.subplot(2,1,1)
ax2 = plt.subplot(2,1,2)

for TXTFILE in FILELIST:
	TIMESERIE = read_csv(FOLDER+TXTFILE,delimiter=" ",header=None)
	
	#Normalized Ux
	meanUx = np.mean(TIMESERIE[1])
	ax1.plot(TIMESERIE[0],TIMESERIE[1]-meanUx,color="purple",lw=1,alpha=0.1)
	
	# #Normalized Uy
	# meanUy = np.mean(TIMESERIE[2])
	# ax2.plot(TIMESERIE[0],TIMESERIE[2]-meanUy,color="green",lw=1,alpha=0.02)
	
plt.tight_layout()
plt.show()