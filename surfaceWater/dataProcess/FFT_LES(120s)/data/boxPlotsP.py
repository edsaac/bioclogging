import numpy as np
import matplotlib.pyplot as plt
from pandas import read_csv
from os import system
from os import listdir
import matplotlib.gridspec as gs
import sys

def rollingAverage(X):
	AV = np.zeros_like(X)
	for i in range(len(X)):
		AV[i] = np.mean(X[:i])
	return AV

FOLDER = "./organizedData/"
FILELIST = listdir(FOLDER)
FILELIST.sort()

fig = plt.figure(figsize=(10,4),facecolor="white")
ncols=len(FILELIST)

gs1 = gs.GridSpec(nrows=ncols,ncols=1)

counter = 0
ax0 = plt.subplot(gs1[counter,0])

for TXTFILE in FILELIST:
	ax1 = plt.subplot(gs1[counter,0],sharex=ax0)
	
	TIMESERIE = read_csv(FOLDER+TXTFILE,delimiter=" ",header=None)
	U = TIMESERIE[1]
	T = TIMESERIE[0]
	ax1.boxplot(U,positions=[-1])

#	ax1.plot(T,U,lw=1,color="pink",alpha=0.6)
#	ax1.plot(T,rollingAverage(U),lw=2,color="purple",alpha=1)

	u = U - np.mean(U)
	ututs = np.correlate(u,u,mode="same")
	varu = np.multiply(u,u)
	AUTOCORR = np.divide(ututs,varu)

	#print(AUTOCORR)
	#input("Wait")
	ax1.set_xlim(-2,35)
	plt.setp(ax1.get_xticklabels(), visible=False)
	counter += 1
plt.show()