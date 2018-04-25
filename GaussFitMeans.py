from sys import argv
import numpy as np
import matplotlib.pyplot as mpl
from scipy.optimize import curve_fit
from scipy.optimize import least_squares

#how many Gaussians?
#argv = np.array(argv)
#numGauss = argv.astype(int)
numGauss = [2]
#separated or no?
sep = 1

with open('./GNUData/GNU2DeutDD.txt',"r") as f:
	data = f.readlines()

xdata = np.zeros(len(data)-1)
ydata = np.zeros(len(data)-1)

for i in range(1,len(data)):
	current = data[i].split()
	xdata[i-1] = current[0]
	ydata[i-1] = current[1]

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

if 1 in numGauss:
	initparam = [0.5,1,0.5]

	lowerbounds = [0,0.5,0]
	upperbounds = [2,2.2,0.5]
	popt, pcov = curve_fit(Gauss1,xdata,ydata, p0=initparam, bounds=(lowerbounds,upperbounds))
	mpl.plot(xdata,Gauss1(xdata,*popt),'r-',label="fit")
	mpl.plot(xdata,ydata, 'b', label="data")
	mpl.xticks(np.arange(0.5,3,0.5))
	mpl.axis((0.5,3,0,2.0))
	for x in range(0,len(popt)):
		if x%3==1:
			print "%f\t%f\t%f" % (popt[x-1],popt[x],popt[x+1])
	
if 2 in numGauss:
	initparam = [0.5,1.65,0.5,0.5,1.95,0.5]

	lowerbounds = [0,0.5,0,0,0.5,0]
	upperbounds = [2,2.2,0.5,2,2.2,0.5]
	popt, pcov = curve_fit(Gauss2,xdata,ydata, p0=initparam, bounds=(lowerbounds,upperbounds))
	p1 = [popt[0],popt[1],popt[2]]
	p2 = [popt[3],popt[4],popt[5]]
	mpl.figure()
	mpl.plot(xdata,Gauss2(xdata,*popt),'r-',label="fit")
 	mpl.plot(xdata,ydata, 'b', label="data")
	mpl.xticks(np.arange(0.5,3,0.5))
	#mpl.axis((0.5,3,0,2.0))
        mpl.axis((0.5,3,0,2.5))
        if sep ==1:
		mpl.plot(xdata,Gauss1(xdata,*p1),label="a(x)")
		mpl.plot(xdata,Gauss1(xdata,*p2),label="b(x)")
	for x in range(0,len(popt)):
		if x%3==1:
			print "%f\t%f\t%f" % (popt[x-1],popt[x],popt[x+1])
			#print "2\t%f\tarea\t%f" % (popt[x],popt[x-1]*popt[x+1]*np.sqrt(2*np.pi))
	
if 3 in numGauss:
	initparam = [0.5,1.8,0.5,0.5,1.8,0.5,0.5,1.8,0.5]

	lowerbounds = [0,0.5,0,0,0.5,0,0,0.5,0]
	upperbounds = [2,2.2,0.5,2,2.2,0.5,2,2.2,0.5]
	popt, pcov = curve_fit(Gauss3,xdata,ydata, p0=initparam, bounds=(lowerbounds,upperbounds))
	p1 = [popt[0],popt[1],popt[2]]
	p2 = [popt[3],popt[4],popt[5]]
	p3 = [popt[6],popt[7],popt[8]]
	mpl.figure()
	mpl.plot(xdata,Gauss3(xdata,*popt),'r-',label="fit")
	mpl.plot(xdata,ydata, 'b', label="data")
	mpl.xticks(np.arange(0.5,3,0.5))
	mpl.axis((0.5,3,0,2.0))
	if sep ==1:
		mpl.plot(xdata,Gauss1(xdata,*p1),label="a(x)")
		mpl.plot(xdata,Gauss1(xdata,*p2),label="b(x)")
		mpl.plot(xdata,Gauss1(xdata,*p3),label="c(x)")
	for x in range(0,len(popt)):
		if x%3==1:
			#print "3\t%f\tarea\t%f" % (popt[x],popt[x-1]*popt[x+1]*np.sqrt(2*np.pi))
			print "%f\t%f\t%f" % (popt[x-1],popt[x],popt[x+1])
			
if 4 in numGauss:
	initparam = [0.5,1,0.5,0.5,1,0.5,0.5,2,0.5,0.5,2,0.5]

	lowerbounds = [0,0.5,0,0,0.5,0,0,0.5,0,0,0.5,0]
	upperbounds = [2,2.2,0.5,2,2.2,0.5,2,2.2,0.5,2,2.2,0.5]
	popt, pcov = curve_fit(Gauss4,xdata,ydata, p0=initparam, bounds=(lowerbounds,upperbounds))
	p1 = [popt[0],popt[1],popt[2]]
	p2 = [popt[3],popt[4],popt[5]]
	p3 = [popt[6],popt[7],popt[8]]
	p4 = [popt[9],popt[10],popt[11]]
	mpl.figure()
	mpl.plot(xdata,Gauss4(xdata,*popt),'r-',label="fit")
	mpl.plot(xdata,ydata, 'b', label="data")
	mpl.xticks(np.arange(0.5,3,0.5))
	mpl.axis((0.5,3,0,2.0))
	if sep ==1:
		mpl.plot(xdata,Gauss1(xdata,*p1),label="a(x)")
		mpl.plot(xdata,Gauss1(xdata,*p2),label="b(x)")
		mpl.plot(xdata,Gauss1(xdata,*p3),label="c(x)")
		mpl.plot(xdata,Gauss1(xdata,*p4),label="d(x)")
	for x in range(0,len(popt)):
		if x%3==1:
			print "%f\t%f\t%f" % (popt[x-1],popt[x],popt[x+1])
			#print "4\t%f\tarea\t%f" % (popt[x],popt[x-1]*popt[x+1]*np.sqrt(2*np.pi))
	
if 5 in numGauss:
	initparam = [0.5,1,0.5,0.5,1,0.5,0.5,1.7,0.5,0.5,1.9,0.5,0.5,2,0.5]

	lowerbounds = [0,0.5,0,0,0.5,0,0,0.5,0,0,0.5,0,0,0.5,0]
	upperbounds = [2,2.2,0.5,2,2.2,0.5,2,2.2,0.5,2,2.2,0.5,2,2.2,0.5]
	popt, pcov = curve_fit(Gauss5,xdata,ydata, p0=initparam, bounds=(lowerbounds,upperbounds))
	p1 = [popt[0],popt[1],popt[2]]
	p2 = [popt[3],popt[4],popt[5]]
	p3 = [popt[6],popt[7],popt[8]]
	p4 = [popt[9],popt[10],popt[11]]
	p5 = [popt[12],popt[13],popt[14]]
	mpl.figure()	
	mpl.plot(xdata,Gauss5(xdata,*popt),'r-',label="fit")
	mpl.plot(xdata,ydata, 'b', label="data")
	mpl.xticks(np.arange(0.5,3,0.5))
	mpl.axis((0.5,3,0,2.0))
	if sep ==1:
		mpl.plot(xdata,Gauss1(xdata,*p1),label="a(x)")
		mpl.plot(xdata,Gauss1(xdata,*p2),label="b(x)")
		mpl.plot(xdata,Gauss1(xdata,*p3),label="c(x)")
		mpl.plot(xdata,Gauss1(xdata,*p4),label="d(x)")
		mpl.plot(xdata,Gauss1(xdata,*p5),label="e(x)")
	for x in range(0,len(popt)):
		if x%3==1:
			print "%f\t%f\t%f" % (popt[x-1],popt[x],popt[x+1])
			#print "5\t%f\tarea\t%f" % (popt[x],popt[x-1]*popt[x+1]*np.sqrt(2*np.pi))

if 6 in numGauss:
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
	mpl.figure()	
	mpl.plot(xdata,Gauss6(xdata,*popt),'r-',label="fit")
	mpl.plot(xdata,ydata, 'b', label="data")
	mpl.xticks(np.arange(0.5,3,0.5))
	mpl.axis((0.5,3,0,2.0))
	if sep ==1:
		mpl.plot(xdata,Gauss1(xdata,*p1),label="a(x)")
		mpl.plot(xdata,Gauss1(xdata,*p2),label="b(x)")
		mpl.plot(xdata,Gauss1(xdata,*p3),label="c(x)")
		mpl.plot(xdata,Gauss1(xdata,*p4),label="d(x)")
		mpl.plot(xdata,Gauss1(xdata,*p5),label="e(x)")
		mpl.plot(xdata,Gauss1(xdata,*p6),label="f(x)")
	for x in range(0,len(popt)):
		if x%3==1:
			print "%f\t%f\t%f" % (popt[x-1],popt[x],popt[x+1])
			#print "6\t%f\tarea\t%f" % (popt[x],popt[x-1]*popt[x+1]*np.sqrt(2*np.pi))
	
if 7 in numGauss:
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
	mpl.figure()	
	mpl.plot(xdata,Gauss7(xdata,*popt),'r-',label="fit")
	mpl.plot(xdata,ydata, 'b', label="data")
	mpl.xticks(np.arange(0.5,3,0.5))
	mpl.axis((0.5,3,0,2.0))
	if sep ==1:
		mpl.plot(xdata,Gauss1(xdata,*p1),label="a(x)")
		mpl.plot(xdata,Gauss1(xdata,*p2),label="b(x)")
		mpl.plot(xdata,Gauss1(xdata,*p3),label="c(x)")
		mpl.plot(xdata,Gauss1(xdata,*p4),label="d(x)")
		mpl.plot(xdata,Gauss1(xdata,*p5),label="e(x)")
		mpl.plot(xdata,Gauss1(xdata,*p6),label="f(x)")
		mpl.plot(xdata,Gauss1(xdata,*p7),label="f(x)")
	for x in range(0,len(popt)):
		if x%3==1:
			print "7\t%f\tarea\t%f" % (popt[x],popt[x-1]*popt[x+1]*np.sqrt(2*np.pi))

mpl.show()
