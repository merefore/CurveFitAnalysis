import numpy as np
from scipy.spatial.distance import cdist

#H5+ style structure
A = [-1,0,3]
B = [1,0,3]
C = [0,0,2]
D = [0,0,1]

#CD3H2+ style structure
E = [0,1,2]
F = [np.sqrt(3)/2,-1/2,2]
G = [-np.sqrt(3)/2,-1/2,2]

angles = range(10,370,10)
angles = np.radians(angles)
xy1 = []
xy2 = []

for i in xrange(0,len(angles)):
    x = np.cos(angles[i])
    y = np.sin(angles[i])
    xy1.append([x,y,0])
    xy2.append([-x,-y,0])

fixdists = [A,B,C,D]
fixdists2 = [E,F,G]
ans1 = cdist(fixdists,xy1) #H5
ans2 = cdist(fixdists,xy2) #H5
ans3 = cdist(fixdists2,xy1) #CD3H2
ans4 = cdist(fixdists2,xy2) #CD3H2

for i in xrange(0,len(ans1[0])):
    for j in [ans1,ans2]:
        print "%f\t%f\t%f\t%f" % (j[0][i],j[1][i],j[2][i],j[3][i])

#for i in xrange(0,len(ans3[0])):
#    for j in [ans3,ans4]:
#        print "%f\t%f\t%f" % (j[0][i],j[1][i],j[2][i])
