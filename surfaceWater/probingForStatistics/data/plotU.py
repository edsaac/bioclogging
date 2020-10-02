import numpy as np
import matplotlib.pyplot as plt
from pandas import read_csv
from os import system
from os import listdir
import sys

def rollingAverage(X,k):
	AV = np.zeros_like(X)
	nElem = len(X)
	
	for i in range(nElem):
		if (i > k) * (i+k < nElem): 
			AV[i] = np.mean(X[i-k:i+k+1])
		else:
			AV[i] = None
	return AV

#FOLDER = "./organizedData/lineAtBot/"
#FOLDER = "./organizedData/lineAtTop/"
FOLDER = "./organizedData/"

FILELIST = listdir(FOLDER)
FILELIST.sort()
d = 0.01 # Readings timestep

for TXTFILE in FILELIST:
	plt.figure(figsize=(8,8),facecolor="white")
	ax0 = plt.subplot(3,1,1)
	ax1 = plt.subplot(3,1,2)
	ax2 = plt.subplot(3,1,3)

	TIMESERIE = read_csv(FOLDER+TXTFILE,delimiter=" ",header=None)
	T = np.array(TIMESERIE[0])
	U = np.array(TIMESERIE[1])

	#All data
	timeForSlice = 20
	
	ax0.plot(T,U,color="purple",lw=0.2,alpha=0.2,marker='o',markersize=1)
	ax0.axvline(x=timeForSlice,lw=1,color="gray")
	ax0.set_xlim(T[0],T[-1])
	ax0.set_ylabel("$p$")
	ax1.set_xlabel("Time [s]")
	
	#Slice data4

	U_slc = U[int(timeForSlice/d):]
	T_slc = T[int(timeForSlice/d):]
	U_slc_Mean = np.mean(U_slc)
	uNorm = (U_slc - U_slc_Mean)/U_slc_Mean

	#Frequency of that stuff
	ax1.plot(T_slc,uNorm,color="green",lw=0.2,alpha=0.2,marker='o',markersize=1)
	ax1.boxplot(uNorm,positions=[timeForSlice/2],widths=[10])
	ax1.axvline(x=timeForSlice,lw=1,color="gray")
	ax1.set_xlim(T[0],T[-1])
	ax1.set_ylabel("$p' = \\dfrac{p-\\langle p \\rangle}{\\langle p \\rangle}$")
	plt.setp(ax1.get_xticklabels(), visible=False)
	
	U_fft = np.absolute(np.fft.rfft(uNorm))/len(uNorm)
	Freq  = np.fft.rfftfreq(uNorm.shape[-1],d=d)
	ax2.plot(Freq,U_fft,color="orange",alpha=0.6,lw=4)
	ax2.set_ylabel("Rel. f.")
	
	characFreq = Freq[int(np.argmax(U_fft))]
	kForAverage = int((1./characFreq)/d)
	ax2.axvline(x=characFreq,lw=1,color="gray",label="Max Freq = %.2f $s^{-1}$" %(characFreq))
	ax2.legend()
	ax2.set_xlabel("Frequency [Hz]")

	ax1.plot(T_slc,rollingAverage(uNorm,kForAverage),color="black",alpha=1,lw=1)
	ax0.plot(T,rollingAverage(U,kForAverage),color="black",alpha=1,lw=1)

	# #Normalized Uy
	# meanUy = np.mean(TIMESERIE[2])
	# ax2.plot(TIMESERIE[0],TIMESERIE[2]-meanUy,color="green",lw=1,alpha=0.02)
	
	plt.tight_layout()
	plt.savefig(TXTFILE+".png")