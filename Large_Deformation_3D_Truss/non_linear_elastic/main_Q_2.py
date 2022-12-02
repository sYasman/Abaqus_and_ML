from math import pi
import numpy as np
import pandas as pd
from Truss3D_functions import get_K_NL
from matplotlib import pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
import time
t = time.time()
# Read data, okurken ilk satir data bilgisi oalcak.
X = pd.read_csv("n_cords.txt").to_numpy()
NC = pd.read_csv("el_nodes.txt").to_numpy()
BC = pd.read_csv("restrained_nodes.txt").to_numpy()

x = np.zeros( (int(X.shape[0]), int(X.shape[1])) )
#material properties
A = 1.0
E = 1.0E+7

# node and element numbers
nnode = X.shape[0]
numel = NC.shape[0]
n_el_col = NC.shape[1]
n_BC = BC.shape[0]

MG = np.zeros( (numel,2) )
MG[:,0]=A
MG[:,1]=E
# update DOF numbering
DOF = np.zeros([nnode, 3])
k = 1

for i in range(0, nnode):
    bcflag = 0
    for j in range(0, BC.shape[0]):
        if i + 1 == BC[j, 0]:
            bcflag = 1
            row_id = j

    if bcflag == 0:
        DOF[i, 0] = k
        DOF[i, 1] = k + 1
        DOF[i, 2] = k + 2
        k += 3

    if bcflag == 1 and BC[row_id, 1] == 0:
        DOF[i, 0] = k
        k += 1

    if bcflag == 1 and BC[row_id, 2] == 0:
        DOF[i, 1] = k
        k += 1

    if bcflag == 1 and BC[row_id, 3] == 0:
        DOF[i, 2] = k
        k += 1
        
nfrdof = int(DOF.max())
r, c = DOF.shape

for i in range(0, r):
    for j in range(0, c):
        if DOF[i, j] == 0:
            DOF[i, j] = k
            k += 1

n_dof = int(DOF.max())    

#load = -9.0

tfin = 1.0
numstep = 100
deltime = tfin/numstep
tol = 1.0E-8
itermax = 200
force_int = np.zeros( (numel, 1) )

#rb = int( BC[:,1:4].sum() ) #number of prescribed DOF
#ra = int( DOF.max()-rb ) #number of active DOF

ra = int( nfrdof ) #number of active DOF
rb = n_dof-ra #number of prescribed DOF

#x = X
U = np.zeros( (int(DOF.max()) ,1) )
Fext = np.zeros( (int(DOF.max()) ,1) )        
U_hist = []   
f_hist_x = []
f_hist_y = []
f_hist_z = []
for i in range(1, numstep+1):
    tt = i*deltime        
    for el in BC:
        node_w_displ = el[0]
        if el[1] == 1:
            relatedDof = int(DOF[int(el[0])-1,0])
            U[relatedDof-1] = el[4]*tt/tfin
            
        if el[2] == 1:
            relatedDof = int(DOF[int(el[0])-1,1])
            U[relatedDof-1] = el[5]*tt/tfin    
            
        if el[3] == 1:
            relatedDof = int(DOF[int(el[0])-1,2])
            U[relatedDof-1] = el[6]*tt/tfin              
            
    for j in range(nnode):
        dof_id_x = int(DOF[j,0]-1)
        dof_id_y = int(DOF[j,1]-1)
        dof_id_z = int(DOF[j,2]-1)
        x[j,1] = X[j,1]+U[dof_id_x]
        x[j,2] = X[j,2]+U[dof_id_y]
        x[j,3] = X[j,3]+U[dof_id_z]
        
        
    K, fint = get_K_NL(MG, numel, DOF, NC, x, X)    

    Ra = fint[ 0:ra, 0]
    nrm = np.linalg.norm(Ra)
    
    iter_num = 1
    while(nrm > tol and iter_num < itermax):
        Kaa = K[0:nfrdof, 0:nfrdof]
        Kab = K[0:nfrdof, nfrdof:int(DOF.max())]
        Kba = K[nfrdof:int(DOF.max()), 0:nfrdof]
        Kbb = K[nfrdof:int(DOF.max()), nfrdof:int(DOF.max())]
        
        delUa = np.zeros( (ra,1) ) 
        delUb = np.zeros( (rb,1) )
        Kaa_inv = np.linalg.pinv(Kaa)
        RHS = Ra-np.transpose(np.matmul(Kab,delUb))
        delUa = np.matmul(-Kaa_inv, np.transpose(RHS))
        delU = np.vstack((delUa, delUb))
        U = U + delU
        
        for j in range(nnode):
            dof_id_x = int(DOF[j,0]-1)
            dof_id_y = int(DOF[j,1]-1)
            dof_id_z = int(DOF[j,2]-1)
            x[j,1] = X[j,1]+U[dof_id_x]
            x[j,2] = X[j,2]+U[dof_id_y]
            x[j,3] = X[j,3]+U[dof_id_z]
            
        x[:,0] = X[:,0]    
        K, fint = get_K_NL(MG, numel, DOF, NC, x, X)
        Ra = fint[ 0:ra, 0] - Fext[ 0:ra, 0]
        
        nrm = np.linalg.norm(Ra)
        print('Increment: {1}. Iteration: {2}\
              Norm: {0:.15f}'.format(nrm, i, iter_num))
        iter_num += 1
    f_hist_x.append(fint[27])  
    f_hist_y.append(fint[28]) 
    f_hist_z.append(fint[29]) 
    U_hist.append(U[20])

plt.plot(U_hist, f_hist_x, 'k--o', Label='Internal force x-dir')
plt.grid()
plt.plot(U_hist, f_hist_y, 'b--x',Label='Internal force y-dir')
plt.grid()
plt.plot(U_hist, f_hist_z, 'r--d',Label='Internal force z-dir')
plt.grid()
plt.xlabel('Displecement in y-dir at node 1 (mm)')
plt.ylabel('Internal force (N)')
plt.title('Internal force vs. applied displacement')
plt.axis('tight')
plt.xlim([0,-6.3])
plt.legend()
plt.legend(fontsize=20)
fig = plt.figure()
ax = plt.axes(projection='3d')
kk =1  
for el in NC:
    bn_id = int( el[1]-1 )
    en_id = int( el[2]-1 )
    XX = [X[bn_id, 1], X[en_id, 1] ]
    YY = [X[bn_id, 2], X[en_id, 2] ]
    ZZ = [X[bn_id, 3], X[en_id, 3] ]
    if kk == 1:
        ax.plot3D(XX, YY, ZZ, 'g--d', Label='Undeformed') 
        kk += 1
    ax.plot3D(XX, YY, ZZ, 'g--d') 

kk =1         
for el in NC:
    bn_id = int( el[1]-1 )
    en_id = int( el[2]-1 )
    xx = [x[bn_id, 1], x[en_id, 1] ]
    yy = [x[bn_id, 2], x[en_id, 2] ]
    zz = [x[bn_id, 3], x[en_id, 3] ]
    if kk == 1:
        ax.plot3D(xx, yy, zz, 'r-o',  Label='Deformed')
        kk += 1
    ax.plot3D(xx, yy, zz, 'r-o')
ax.legend()  
plt.legend(fontsize=20)
ax.set_xlabel('X-label')  
ax.set_ylabel('y-label')  
ax.set_zlabel('z-label')  

elapsed = time.time() - t
plt.title('3D truss example, elapsed time {0:.3f} seconds.'.format(elapsed))
    