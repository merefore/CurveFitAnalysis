import numpy as np
import matplotlib.pyplot as mpl
from scipy.optimize import curve_fit
from scipy.optimize import least_squares
from scipy.special import gamma
from scipy.special import erf

#separated or no?
sep = 1

from os import path
#file_path = path.relpath("../for_meredith/CH5plus/Data/5Deut_20K10K/5DeutDDPile.txt")
#file_path = path.relpath("../for_meredith/CH5plus/Data/CH5HHPile.txt")
with open("./finalhistpiles/2DD.txt",'r') as f:
	histPile = f.readlines()

histPile = np.array(histPile)
print "len histpile %i" % len(histPile)
histPile = histPile.astype(float)
standard_bins = np.arange(0.5,2.7,0.02)
totalHIST, bin_edges = np.histogram(histPile,bins=standard_bins)
#totalHIST, bin_edges = np.histogram(histPile, bins=75)
width1 = bin_edges[1]-bin_edges[0]
newNhist = totalHIST/(np.sum(totalHIST)*width1)
bin_edges = np.delete(bin_edges,len(bin_edges)-1)

xdata = bin_edges
ydata = newNhist

def Gauss1(x,A,c1,sig1,k1):
	return A*np.exp(-((x-c1)**2)/(2*sig1**2))*(1+erf(k1*(x-c1)/2))

def Gauss2(x,A,c1,sig1,k1,B,c2,sig2,k2):
	return A*np.exp(-(x-c1)**2/(2*sig1**2))*(1+erf(k1*(x-c1)/2))+B*np.exp(-(x-c2)**2/(2*sig2**2))*(1+erf(k2*(x-c2)/2))

def Gauss3(x,A,c1,sig1,k1,B,c2,sig2,k2,C,c3,sig3,k3):
	return A*np.exp(-(x-c1)**2/(2*sig1**2))*(1+erf(k1*(x-c1)/2))+B*np.exp(-(x-c2)**2/(2*sig2**2))*(1+erf(k2*(x-c2)/2))+C*np.exp(-(x-c3)**2/(2*sig3**2))*(1+erf(k3*(x-c3)/2))

def Gauss4(x,A,c1,sig1,k1,B,c2,sig2,k2,C,c3,sig3,k3,D,c4,sig4,k4):
	return A*np.exp(-(x-c1)**2/(2*sig1**2))*(1+erf(k1*(x-c1)/2))+B*np.exp(-(x-c2)**2/(2*sig2**2))*(1+erf(k2*(x-c2)/2))+C*np.exp(-(x-c3)**2/(2*sig3**2))*(1+erf(k3*(x-c3)/2))+D*np.exp(-(x-c4)**2/(2*sig4**2))*(1+erf(k4*(x-c4)/2))

def Gauss5(x,A,c1,sig1,k1,B,c2,sig2,k2,C,c3,sig3,k3,D,c4,sig4,k4,E,c5,sig5,k5):
	return A*np.exp(-(x-c1)**2/(2*sig1**2))+B*np.exp(-(x-c2)**2/(2*sig2**2))+C*np.exp(-(x-c3)**2/(2*sig3**2))+D*np.exp(-(x-c4)**2/(2*sig4**2))+E*np.exp(-(x-c5)**2/(2*sig5**2))

def Gauss6(x,A,c1,sig1,k1,B,c2,sig2,k2,C,c3,sig3,k3,D,c4,sig4,k4,E,c5,sig5,k5,F,c6,sig6,k6):
	return A*np.exp(-(x-c1)**2/(2*sig1**2))+B*np.exp(-(x-c2)**2/(2*sig2**2))+C*np.exp(-(x-c3)**2/(2*sig3**2))+D*np.exp(-(x-c4)**2/(2*sig4**2))+E*np.exp(-(x-c5)**2/(2*sig5**2))+F*np.exp(-(x-c6)**2/(2*sig6**2))

def Gauss7(x,A,c1,sig1,k1,B,c2,sig2,k2,C,c3,sig3,k3,D,c4,sig4,k4,E,c5,sig5,k5,F,c6,sig6,k6,G,c7,sig7,k7):
	return A*np.exp(-(x-c1)**2/(2*sig1**2))+B*np.exp(-(x-c2)**2/(2*sig2**2))+C*np.exp(-(x-c3)**2/(2*sig3**2))+D*np.exp(-(x-c4)**2/(2*sig4**2))+E*np.exp(-(x-c5)**2/(2*sig5**2))+F*np.exp(-(x-c6)**2/(2*sig6**2))+G*np.exp(-(x-c7)**2/(2*sig7**2))

def AICc(Gauss,inputs,histPile):
	lnL = 0
	k = len(inputs)
	for HHdist in histPile:
		lnL += np.log(Gauss(HHdist,*inputs))
	AICvalue = 2*k-2*lnL+2*k*(k+1)/(len(histPile)-k-1)
	return AICvalue

def morse(x,A,k,c1,alpha):
	return A*np.sqrt((alpha*(k-1))/gamma(k))*np.exp(-(k/2)*np.exp(-alpha*(x-c1)))*(k*np.exp(-alpha*(x-c1)))**((k-1)/2)
def skew(x,A,c1,sig1,alpha):
	return A*np.exp(-((x-c1)**2)/(2*sig1**2))*(1+erf(alpha*(x-c1)/2))




# initparammorse = [1,5,-1,3]
# initparamskew = [1,1,1,20]
#
#
#
# poptmorse, pcovmorse = curve_fit(morse, xdata,ydata, p0=initparammorse)
# poptskew, pcovskew = curve_fit(skew, xdata,ydata, p0=initparamskew)
# print "Morse: A, k, r0, alpha"
# print "%f\t%f\t%f\t%f" % (poptmorse[0],poptmorse[1],poptmorse[2],poptmorse[3])
# print "Skew: A, mean, sig, alpha"
# print "%f\t%f\t%f\t%f" % (poptskew[0],poptskew[1],poptskew[2],poptskew[3])
#
# mpl.plot(xdata,ydata, 'k', label="data")
# mpl.plot(xdata,morse(xdata,*poptmorse),'r-',label="morse")
# mpl.plot(xdata,skew(xdata,*poptskew),'b-',label="skew")



numGauss = 1
initparam = [0.5,1,0.5,0]
lowerbounds = [0,0.5,0,-1000]
upperbounds = [2,2.2,0.5,1000]
popt, pcov = curve_fit(Gauss1,xdata,ydata, p0=initparam, bounds=(lowerbounds,upperbounds))
mpl.figure()
mpl.plot(xdata,Gauss1(xdata,*popt),'r-',label="fit")
mpl.plot(xdata,ydata, 'b', label="data")
mpl.xticks(np.arange(0.5,3,0.5))
mpl.axis((0.5,3,0,2.0))
for x in range(0,len(popt)):
	if x%4==0:
		print "%i\t%f\t%f\t%f\t%f" % (numGauss,popt[x],popt[x+1],popt[x+2],popt[x+3])
value = AICc(Gauss1,popt,histPile)
print "%i\t%f" % (numGauss,value)
mpl.plot(xdata,Gauss1(xdata,*popt),'g-',label="gaussian")
mpl.xticks(np.arange(0.5,3,0.5))
mpl.axis((0.5,3,0,2.0))

numGauss = 2
initparam = [0.5,1.8,0.5,0,0.5,1.8,0.5,0]
lowerbounds = [0,0.5,0,-1000,0,0.5,0,-1000]
upperbounds = [2,2.2,0.5,1000,2,2.2,0.5,1000]
popt, pcov = curve_fit(Gauss2,xdata,ydata, p0=initparam, bounds=(lowerbounds,upperbounds))
p1 = [popt[0],popt[1],popt[2],popt[3]]
p2 = [popt[4],popt[5],popt[6],popt[7]]
mpl.figure()
mpl.plot(xdata,Gauss2(xdata,*popt),'r-',label="fit")
mpl.plot(xdata,ydata, 'b', label="data")
mpl.xticks(np.arange(0.5,3,0.5))
mpl.axis((0.5,3,0,2.0))
if sep ==1:
	mpl.plot(xdata,Gauss1(xdata,*p1),'k',label="a(x)")
	mpl.plot(xdata,Gauss1(xdata,*p2),'k',label="b(x)")
for x in range(0,len(popt)):
	if x%4==0:
		print "%i\t%f\t%f\t%f\t%f" % (numGauss,popt[x],popt[x+1],popt[x+2],popt[x+3])
value = AICc(Gauss2,popt,histPile)
print "%i\t%f" % (numGauss,value)

numGauss = 3
initparam = [0.5,1.8,0.5,0,0.5,1.8,0.5,0,0.5,1.8,0.5,0]
lowerbounds = [0,0.5,0,-1000,0,0.5,0,-1000,0,0.5,0,-1000]
upperbounds = [2,2.2,0.5,1000,2,2.2,0.5,1000,2,2.2,0.5,1000]
popt, pcov = curve_fit(Gauss3,xdata,ydata, p0=initparam, bounds=(lowerbounds,upperbounds))
p1 = [popt[0],popt[1],popt[2],popt[3]]
p2 = [popt[4],popt[5],popt[6],popt[7]]
p3 = [popt[8],popt[9],popt[10],popt[11]]
mpl.figure()
mpl.legend()
mpl.plot(xdata,Gauss3(xdata,*popt),'r-',label="fit")
mpl.plot(xdata,ydata, 'b', label="data")
mpl.xticks(np.arange(0.5,3,0.5))
mpl.axis((0.5,3,0,2.0))
if sep ==1:
	mpl.plot(xdata,Gauss1(xdata,*p1),'k',label="a(x)")
	mpl.plot(xdata,Gauss1(xdata,*p2),'k',label="b(x)")
	mpl.plot(xdata,Gauss1(xdata,*p3),'k',label="c(x)")
for x in range(0,len(popt)):
	if x%4==0:
		print "%i\t%f\t%f\t%f\t%f" % (numGauss,popt[x],popt[x+1],popt[x+2],popt[x+3])
value = AICc(Gauss3,popt,histPile)
print "%i\t%f" % (numGauss,value)

numGauss = 4
initparam = [0.5,1.8,0.5,0,0.5,1.8,0.5,0,0.5,1.8,0.5,0,0.5,1.8,0.5,0]
lowerbounds = [0,0.5,0,-1000,0,0.5,0,-1000,0,0.5,0,-1000,0,0.5,0,-1000]
upperbounds = [2,2.2,0.5,1000,2,2.2,0.5,1000,2,2.2,0.5,1000,2,2.2,0.5,1000]
popt, pcov = curve_fit(Gauss4,xdata,ydata, p0=initparam, bounds=(lowerbounds,upperbounds))
p1 = [popt[0],popt[1],popt[2],popt[3]]
p2 = [popt[4],popt[5],popt[6],popt[7]]
p3 = [popt[8],popt[9],popt[10],popt[11]]
p4 = [popt[12],popt[13],popt[14],popt[15]]
mpl.figure()
mpl.plot(xdata,Gauss4(xdata,*popt),'r-',label="fit")
mpl.plot(xdata,ydata, 'b', label="data")
mpl.xticks(np.arange(0.5,3,0.5))
mpl.axis((0.5,3,0,2.0))
if sep ==1:
	mpl.plot(xdata,Gauss1(xdata,*p1),'k',label="a(x)")
	mpl.plot(xdata,Gauss1(xdata,*p2),'k',label="b(x)")
	mpl.plot(xdata,Gauss1(xdata,*p3),'k',label="c(x)")
	mpl.plot(xdata,Gauss1(xdata,*p4),'k',label="d(x)")
for x in range(0,len(popt)):
	if x%4==0:
		print "%i\t%f\t%f\t%f\t%f" % (numGauss,popt[x],popt[x+1],popt[x+2],popt[x+3])
value = AICc(Gauss4,popt,histPile)
print "%i\t%f" % (numGauss,value)

mpl.show()
assert 1==0

numGauss = 5
initparam = [0.5,1.8,0.5,0,0.5,1.8,0.5,0,0.5,1.8,0.5,0,0.5,1.8,0.5,0,0.5,1.8,0.5,0]
lowerbounds = [0,0.5,0,-1000,0,0.5,0,-1000,0,0.5,0,-1000,0,0.5,0,-1000,0,0.5,0,-1000]
upperbounds = [2,2.2,0.5,1000,2,2.2,0.5,1000,2,2.2,0.5,1000,2,2.2,0.5,1000,2,2.2,0.5,1000]
popt, pcov = curve_fit(Gauss5,xdata,ydata, p0=initparam, bounds=(lowerbounds,upperbounds))
p1 = [popt[0],popt[1],popt[2],popt[3]]
p2 = [popt[4],popt[5],popt[6],popt[7]]
p3 = [popt[8],popt[9],popt[10],popt[11]]
p4 = [popt[12],popt[13],popt[14],popt[15]]
p5 = [popt[16],popt[17],popt[18],popt[19]]
for x in range(0,len(popt)):
	if x%4==0:
		print "%i\t%f\t%f\t%f\t%f" % (numGauss,popt[x],popt[x+1],popt[x+2],popt[x+3])
value = AICc(Gauss5,popt,histPile)
print "%i\t%f" % (numGauss,value)

numGauss = 6
initparam = [0.5,1.8,0.5,0,0.5,1.8,0.5,0,0.5,1.8,0.5,0,0.5,1.8,0.5,0,0.5,1.8,0.5,0,0.5,1.8,0.5,0]
lowerbounds = [0,0.5,0,-1000,0,0.5,0,-1000,0,0.5,0,-1000,0,0.5,0,-1000,0,0.5,0,-1000,0,0.5,0,-1000]
upperbounds = [2,2.2,0.5,1000,2,2.2,0.5,1000,2,2.2,0.5,1000,2,2.2,0.5,1000,2,2.2,0.5,1000,2,2.2,0.5,1000]
popt, pcov = curve_fit(Gauss6,xdata,ydata, p0=initparam, bounds=(lowerbounds,upperbounds))
p1 = [popt[0],popt[1],popt[2],popt[3]]
p2 = [popt[4],popt[5],popt[6],popt[7]]
p3 = [popt[8],popt[9],popt[10],popt[11]]
p4 = [popt[12],popt[13],popt[14],popt[15]]
p5 = [popt[16],popt[17],popt[18],popt[19]]
p6 = [popt[20],popt[21],popt[22],popt[23]]
for x in range(0,len(popt)):
	if x%4==0:
		print "%i\t%f\t%f\t%f\t%f" % (numGauss,popt[x],popt[x+1],popt[x+2],popt[x+3])
value = AICc(Gauss6,popt,histPile)
print "%i\t%f" % (numGauss,value)

numGauss = 7
initparam = [0.5,1.8,0.5,0,0.5,1.8,0.5,0,0.5,1.8,0.5,0,0.5,1.8,0.5,0,0.5,1.8,0.5,0,0.5,1.8,0.5,0,0.5,1.8,0.5,0]
lowerbounds = [0,0.5,0,-1000,0,0.5,0,-1000,0,0.5,0,-1000,0,0.5,0,-1000,0,0.5,0,-1000,0,0.5,0,-1000,0,0.5,0,-1000]
upperbounds = [2,2.2,0.5,1000,2,2.2,0.5,1000,2,2.2,0.5,1000,2,2.2,0.5,1000,2,2.2,0.5,1000,2,2.2,0.5,1000,2,2.2,0.5,1000]
popt, pcov = curve_fit(Gauss7,xdata,ydata, p0=initparam, bounds=(lowerbounds,upperbounds))
p1 = [popt[0],popt[1],popt[2],popt[3]]
p2 = [popt[4],popt[5],popt[6],popt[7]]
p3 = [popt[8],popt[9],popt[10],popt[11]]
p4 = [popt[12],popt[13],popt[14],popt[15]]
p5 = [popt[16],popt[17],popt[18],popt[19]]
p6 = [popt[20],popt[21],popt[22],popt[23]]
p7 = [popt[24],popt[25],popt[26],popt[27]]
for x in range(0,len(popt)):
	if x%4==0:
		print "%i\t%f\t%f\t%f\t%f" % (numGauss,popt[x],popt[x+1],popt[x+2],popt[x+3])
value = AICc(Gauss7,popt,histPile)
print "%i\t%f" % (numGauss,value)
