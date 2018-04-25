import numpy as np
#import matplotlib.pyplot as mpl
from scipy.optimize import curve_fit
from scipy.optimize import least_squares
from scipy.optimize import fmin
from multiprocessing import Pool

#separated or no?
sep = 1

from os import path
#file_path = path.relpath("../for_meredith/CH5plus/Data/5Deut_20K10K/5DeutDDPile.txt")
#file_path = path.relpath("../for_meredith/CH5plus/Data/CH5HHPile.txt")
with open("./finalhistpiles/3DD.txt",'r') as f:
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

def numGauss(x,theta):
	print len(theta)
	assert len(theta)%4 == 0
	num = len(theta)/4
	sum = 0
	for i in xrange(0,num):
		sum += theta[4*i]*np.exp(-((x-theta[4*i+1])**2)/(2*theta[4*i+2]**2))*(1+erf(theta[4*i+3]*(x-theta[4*i+1])/2))
	return sum

def Gauss1(x,A,c1,sig1,k1):
	return A*np.exp(-((x-c1)**2)/(2*sig1**2))*(1+erf(k1*(x-c1)/2))

def logGauss(x,theta):
	assert len(theta)%3 == 0
	num = len(theta)/3
	sum = 0
	for i in xrange(0,num):
		sum += theta[3*i]*np.exp(-((x-theta[3*i+1])**2)/(2*theta[3*i+2]**2))
	return np.log(sum)

count = 0
def lnL(theta):
	global count
	sum = 0
	for i in xrange(0,len(histPile)):
		sum += np.log(Gauss(histPile[i],theta))
	count += 1
	print count
	return sum

def lnLmulti(theta):
	sum = 0
	p = Pool(10)
	extra = lambda x : np.log(Gauss(x,theta))
	array = p.map(logGauss, histPile)
	sum = np.sum(array)
	return sum

def AICc(Gauss,inputs,histPile):
	lnL = 0
	k = len(inputs)
	for HHdist in histPile:
		lnL += np.log(Gauss(HHdist,*inputs))
	AICvalue = 2*k-2*lnL+2*k*(k+1)/(len(histPile)-k-1)
	return AICvalue

gaussguess3 = [ 1.26170482,  1.92307714,  0.20292108,  0.17721031,  1.04439748, 0.18135109,  0.40331624,  1.59993816,  0.27587171]
gaussguess1 = [1.412796, 1.859040, 0.259415]
#neglnL = lambda theta : -lnLmulti(theta)
#max = fmin(neglnL,gaussguess1)
#print max
#initparam = [0.5,1,0.5]
#lowerbounds = [0,0.5,0]
#upperbounds = [2,2.2,0.5]



popt, pcov = curve_fit(numGauss,xdata,ydata)#, p0=initparam, bounds=(lowerbounds,upperbounds))

# numGauss = 7
# initparam = [0.5,1,0.5,0.5,1.4,0.5,0.5,1.7,0.5,0.5,1.8,0.5,0.5,1.9,0.5,0.5,2,0.5,0.5,2,0.5]
# lowerbounds = [0,0.5,0,0,0.5,0,0,0.5,0,0,0.5,0,0,0.5,0,0,0.5,0,0,0.5,0]
# upperbounds = [2,2.2,0.5,2,2.2,0.5,2,2.2,0.5,2,2.2,0.5,2,2.2,0.5,2,2.2,0.5,2,2.2,0.5]
# popt, pcov = curve_fit(Gauss7,xdata,ydata, p0=initparam, bounds=(lowerbounds,upperbounds))
# p1 = [popt[0],popt[1],popt[2]]
# p2 = [popt[3],popt[4],popt[5]]
# p3 = [popt[6],popt[7],popt[8]]
# p4 = [popt[9],popt[10],popt[11]]
# p5 = [popt[12],popt[13],popt[14]]
# p6 = [popt[15],popt[16],popt[17]]
# p7 = [popt[18],popt[19],popt[20]]
# #for x in range(0,len(popt)):
# #	if x%3==1:
# #		print "%i\t%f\t%f\t%f\t%f" % (numGauss,popt[x-1],popt[x],popt[x+1],popt[x-1]*popt[x+1]*np.sqrt(2*np.pi))
# value = AICc(Gauss7,popt,histPile)
# print "%i\t%f" % (numGauss,value)

for j in range(0,len(popt)):
	if j%4==1:
		print "%i\t%f\t%f\t%f\t%f" % (len(popt)/4,popt[x-1],popt[x],popt[x+1],popt[x-1]*popt[x+1]*np.sqrt(2*np.pi))

mpl.figure()
mpl.plot(xdata,numGauss(xdata,*popt),'r-',label="fit")
mpl.plot(xdata,ydata, 'b', label="data")
mpl.xticks(np.arange(0.5,3,0.5))
#mpl.axis((0.5,3,0,2.0))
mpl.axis((0.5,3,0,2.5))
#if sep ==1:
#	mpl.plot(xdata,Gauss1(xdata,*p1),label="a(x)")
#	mpl.plot(xdata,Gauss1(xdata,*p2),label="b(x)")
mpl.show()
