import numpy as np

I = np.array([[1,0,0],[0,1,0],[0,0,1]])

zero = np.array([[0,0,0],[0,0,0],[0,0,0]])

A = np.block([[I,zero,zero,zero],[zero,I,zero,zero],[zero,zero,I,zero],[zero,zero,zero,I]])
B = np.block([[zero,I,zero,zero],[I,zero,zero,zero],[zero,zero,I,zero],[zero,zero,zero,I]])
C = np.block([[zero,I,zero,zero],[I,zero,zero,zero],[zero,zero,zero,I],[zero,zero,I,zero]])
D = np.block([[I,zero,zero,zero],[zero,I,zero,zero],[zero,zero,zero,I],[zero,zero,I,zero]])
E = np.block([[zero,zero,I,zero],[zero,zero,zero,I],[I,zero,zero,zero],[zero,I,zero,zero]])
F = np.block([[zero,zero,zero,I],[zero,zero,I,zero],[I,zero,zero,zero],[zero,I,zero,zero]])
G = np.block([[zero,zero,zero,I],[zero,zero,I,zero],[zero,I,zero,zero],[I,zero,zero,zero]])
H = np.block([[zero,zero,I,zero],[zero,zero,zero,I],[zero,I,zero,zero],[I,zero,zero,zero]])
l = np.array([A,B,C,D,E,F,G,H])
k = ['A','B','C','D','E','F','G','H']

for i in xrange(0,len(l)):
	arr = np.linalg.eigvals(l[i])
	print k[i], arr


