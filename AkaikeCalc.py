import numpy as np
#import matplotlib.pyplot as mpl
from scipy.optimize import curve_fit
from scipy.optimize import least_squares

#separated or no?
sep = 1

from os import path
#file_path = path.relpath("../for_meredith/CH5plus/Data/5Deut_20K10K/5DeutDDPile.txt")
#file_path = path.relpath("../for_meredith/CH5plus/Data/CH5HHPile.txt")
with open("./Piles/3DeutHDPile.txt",'r') as f:
	histPile = f.readlines()

histPile = np.array(histPile) 
histPile = histPile.astype(float)
standard_bins = np.arange(0.5,2.7,0.02)	
totalHIST, bin_edges = np.histogram(histPile,bins=standard_bins)
#totalHIST, bin_edges = np.histogram(histPile, bins=75)
width1 = bin_edges[1]-bin_edges[0]
newNhist = totalHIST/(np.sum(totalHIST)*width1)
bin_edges = np.delete(bin_edges,len(bin_edges)-1)

xdata = bin_edges
ydata = totalHIST

def Gauss1(x,A,c1,sig1):
	return A*np.exp(-((x-c1)**2)/(2*sig1**2))

def Gauss2(x,A,c1,sig1,B,c2,sig2):
	return A*np.exp(-(x-c1)**2/(2*sig1**2))+B*np.exp(-(x-c2)**2/(2*sig2**2))

def Gauss3(x,A,c1,sig1,B,c2,sig2,C,c3,sig3):
	return A*np.exp(-(x-c1)**2/(2*sig1**2))+B*np.exp(-(x-c2)**2/(2*sig2**2))+C*np.exp(-(x-c3)**2/(2*sig3**2))

def Gauss4(x,A,c1,sig1,B,c2,sig2,C,c3,sig3,D,c4,sig4):
	return A*np.exp(-(x-c1)**2/(2*sig1**2))+B*np.exp(-(x-c2)**2/(2*sig2**2))+C*np.exp(-(x-c3)**2/(2*sig3**2))+D*np.exp(-(x-c4)**2/(2*sig4**2))

def Gauss5(x,A,c1,sig1,B,c2,sig2,C,c3,sig3,D,c4,sig4,E,c5,sig5):
	return A*np.exp(-(x-c1)**2/(2*sig1**2))+B*np.exp(-(x-c2)**2/(2*sig2**2))+C*np.exp(-(x-c3)**2/(2*sig3**2))+D*np.exp(-(x-c4)**2/(2*sig4**2))+E*np.exp(-(x-c5)**2/(2*sig5**2))

def Gauss6(x,A,c1,sig1,B,c2,sig2,C,c3,sig3,D,c4,sig4,E,c5,sig5,F,c6,sig6):
	return A*np.exp(-(x-c1)**2/(2*sig1**2))+B*np.exp(-(x-c2)**2/(2*sig2**2))+C*np.exp(-(x-c3)**2/(2*sig3**2))+D*np.exp(-(x-c4)**2/(2*sig4**2))+E*np.exp(-(x-c5)**2/(2*sig5**2))+F*np.exp(-(x-c6)**2/(2*sig6**2))

def Gauss7(x,A,c1,sig1,B,c2,sig2,C,c3,sig3,D,c4,sig4,E,c5,sig5,F,c6,sig6,G,c7,sig7):
	return A*np.exp(-(x-c1)**2/(2*sig1**2))+B*np.exp(-(x-c2)**2/(2*sig2**2))+C*np.exp(-(x-c3)**2/(2*sig3**2))+D*np.exp(-(x-c4)**2/(2*sig4**2))+E*np.exp(-(x-c5)**2/(2*sig5**2))+F*np.exp(-(x-c6)**2/(2*sig6**2))+G*np.exp(-(x-c7)**2/(2*sig7**2))

def AICc(Gauss,inputs,histPile):
	lnL = 0 
	k = len(inputs)
	for HHdist in histPile:
		lnL += np.log(Gauss(HHdist,*inputs))
	AICvalue = 2*k-2*lnL+2*k*(k+1)/(len(histPile)-k-1)
	return AICvalue


#initparam = [0.5,1,0.5]
#lowerbounds = [0,0.5,0]
#upperbounds = [2,2.2,0.5]
#popt, pcov = curve_fit(Gauss1,xdata,ydata, p0=initparam, bounds=(lowerbounds,upperbounds))


#print "NumGauss\tNorm\tMean\tSigma\tArea"
print "NumGauss\tAICc Value"

numGauss = 2
initparam = [0.5,1.8,0.5,0.5,1.8,0.5]
lowerbounds = [0,0.5,0,0,0.5,0]
upperbounds = [2,2.2,0.5,2,2.2,0.5]
popt, pcov = curve_fit(Gauss2,xdata,ydata, p0=initparam, bounds=(lowerbounds,upperbounds))
p1 = [popt[0],popt[1],popt[2]]
p2 = [popt[3],popt[4],popt[5]]
#for x in range(0,len(popt)):
#	if x%3==1:
#		print "%i\t%f\t%f\t%f\t%f" % (numGauss,popt[x-1],popt[x],popt[x+1],popt[x-1]*popt[x+1]*np.sqrt(2*np.pi))
value = AICc(Gauss2,popt,histPile)
print "%i\t%f" % (numGauss,value)

numGauss = 3
initparam = [0.5,1.8,0.5,0.5,1.8,0.5,0.5,1.8,0.5]
lowerbounds = [0,0.5,0,0,0.5,0,0,0.5,0]
upperbounds = [2,2.2,0.5,2,2.2,0.5,2,2.2,0.5]
popt, pcov = curve_fit(Gauss3,xdata,ydata, p0=initparam, bounds=(lowerbounds,upperbounds))
p1 = [popt[0],popt[1],popt[2]]
p2 = [popt[3],popt[4],popt[5]]
p3 = [popt[6],popt[7],popt[8]]
#for x in range(0,len(popt)):
#	if x%3==1:
#		print "%i\t%f\t%f\t%f\t%f" % (numGauss,popt[x-1],popt[x],popt[x+1],popt[x-1]*popt[x+1]*np.sqrt(2*np.pi))
value = AICc(Gauss3,popt,histPile)
print "%i\t%f" % (numGauss,value)

numGauss = 4
initparam = [0.5,1,0.5,0.5,1,0.5,0.5,2,0.5,0.5,2,0.5]
lowerbounds = [0,0.5,0,0,0.5,0,0,0.5,0,0,0.5,0]
upperbounds = [2,2.2,0.5,2,2.2,0.5,2,2.2,0.5,2,2.2,0.5]
popt, pcov = curve_fit(Gauss4,xdata,ydata, p0=initparam, bounds=(lowerbounds,upperbounds))
p1 = [popt[0],popt[1],popt[2]]
p2 = [popt[3],popt[4],popt[5]]
p3 = [popt[6],popt[7],popt[8]]
p4 = [popt[9],popt[10],popt[11]]
#for x in range(0,len(popt)):
#	if x%3==1:
#		print "%i\t%f\t%f\t%f\t%f" % (numGauss,popt[x-1],popt[x],popt[x+1],popt[x-1]*popt[x+1]*np.sqrt(2*np.pi))
value = AICc(Gauss4,popt,histPile)
print "%i\t%f" % (numGauss,value)

numGauss = 5
initparam = [0.5,1,0.5,0.5,1,0.5,0.5,1.7,0.5,0.5,1.9,0.5,0.5,2,0.5]
lowerbounds = [0,0.5,0,0,0.5,0,0,0.5,0,0,0.5,0,0,0.5,0]
upperbounds = [2,2.2,0.5,2,2.2,0.5,2,2.2,0.5,2,2.2,0.5,2,2.2,0.5]
popt, pcov = curve_fit(Gauss5,xdata,ydata, p0=initparam, bounds=(lowerbounds,upperbounds))
p1 = [popt[0],popt[1],popt[2]]
p2 = [popt[3],popt[4],popt[5]]
p3 = [popt[6],popt[7],popt[8]]
p4 = [popt[9],popt[10],popt[11]]
p5 = [popt[12],popt[13],popt[14]]
#for x in range(0,len(popt)):
#	if x%3==1:
#		print "%i\t%f\t%f\t%f\t%f" % (numGauss,popt[x-1],popt[x],popt[x+1],popt[x-1]*popt[x+1]*np.sqrt(2*np.pi))
value = AICc(Gauss5,popt,histPile)
print "%i\t%f" % (numGauss,value)

numGauss = 6
initparam = [0.5,1,0.5,0.5,1.4,0.5,0.5,1.7,0.5,0.5,1.8,0.5,0.5,1.9,0.5,0.5,2,0.5]
lowerbounds = [0,0.5,0,0,0.5,0,0,0.5,0,0,0.5,0,0,0.5,0,0,0.5,0]
upperbounds = [2,2.2,0.5,2,2.2,0.5,2,2.2,0.5,2,2.2,0.5,2,2.2,0.5,2,2.2,0.5]
popt, pcov = curve_fit(Gauss6,xdata,ydata, p0=initparam, bounds=(lowerbounds,upperbounds))
p1 = [popt[0],popt[1],popt[2]]
p2 = [popt[3],popt[4],popt[5]]
p3 = [popt[6],popt[7],popt[8]]
p4 = [popt[9],popt[10],popt[11]]
p5 = [popt[12],popt[13],popt[14]]
p6 = [popt[15],popt[16],popt[17]]
#for x in range(0,len(popt)):
#	if x%3==1:
#		print "%i\t%f\t%f\t%f\t%f" % (numGauss,popt[x-1],popt[x],popt[x+1],popt[x-1]*popt[x+1]*np.sqrt(2*np.pi))
value = AICc(Gauss6,popt,histPile)
print "%i\t%f" % (numGauss,value)

numGauss = 7
initparam = [0.5,1,0.5,0.5,1.4,0.5,0.5,1.7,0.5,0.5,1.8,0.5,0.5,1.9,0.5,0.5,2,0.5,0.5,2,0.5]
lowerbounds = [0,0.5,0,0,0.5,0,0,0.5,0,0,0.5,0,0,0.5,0,0,0.5,0,0,0.5,0]
upperbounds = [2,2.2,0.5,2,2.2,0.5,2,2.2,0.5,2,2.2,0.5,2,2.2,0.5,2,2.2,0.5,2,2.2,0.5]
popt, pcov = curve_fit(Gauss7,xdata,ydata, p0=initparam, bounds=(lowerbounds,upperbounds))
p1 = [popt[0],popt[1],popt[2]]
p2 = [popt[3],popt[4],popt[5]]
p3 = [popt[6],popt[7],popt[8]]
p4 = [popt[9],popt[10],popt[11]]
p5 = [popt[12],popt[13],popt[14]]
p6 = [popt[15],popt[16],popt[17]]
p7 = [popt[18],popt[19],popt[20]]
#for x in range(0,len(popt)):
#	if x%3==1:
#		print "%i\t%f\t%f\t%f\t%f" % (numGauss,popt[x-1],popt[x],popt[x+1],popt[x-1]*popt[x+1]*np.sqrt(2*np.pi))
value = AICc(Gauss7,popt,histPile)
print "%i\t%f" % (numGauss,value)

