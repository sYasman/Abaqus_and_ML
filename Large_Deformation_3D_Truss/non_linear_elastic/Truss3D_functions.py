from math import log,sqrt
import numpy as np

def get_K_NL(MG, numel, DOF, NC, x, X):
    dim = int(DOF.max())
    KK = np.zeros( (dim, dim) )
    fint = np.zeros( (dim,1) )
    
    for i in range(numel):
        bn_id = int( NC[i,1]-1 )
        en_id = int( NC[i,2]-1 )
        dist_vec = x[en_id,1:]-x[bn_id,1:]
        Dist_vec = X[en_id,1:]-X[bn_id,1:]
        l = (np.matmul(dist_vec, np.transpose(dist_vec)))**0.50
        L = (np.matmul(Dist_vec, np.transpose(Dist_vec)))**0.50
        n = 1.0/l * dist_vec
        stretch =  l / L
        V = MG[i,0] * L
        tau_KR = MG[i,1] * log(stretch)
        Tb = tau_KR * V / l
        Ta = -Tb
#        T = [ [Ta], [Tb] ]
#        Trans = [[n, 0, 0, 0],
#                 [0, 0, 0, n]]
        T_part1 = np.transpose(Ta*n)
        T_part2 = np.transpose(Tb*n)
        T_global = np.concatenate((T_part1, T_part2))
        
        k = V/(l**2) * (MG[i,1] - 2 * tau_KR)
        n_dyadic_n = [ [n[0]* n[0], n[0]* n[1], n[0]* n[2]],
        [n[1]* n[0], n[1]* n[1], n[1]* n[2]],
        [n[2]* n[0], n[2]* n[1], n[2]* n[2]] ]
        
#        np.multiply(k, n_dyadic_n)
#        np.multiply(Tb/l, np.eye( 3 ))
        Kaa = np.multiply(k, n_dyadic_n) + np.multiply(Tb/l, np.eye( 3 ))
        Kbb = Kaa
        Kba = -Kbb
        Kab = Kba
        
        K1 = np.hstack((Kaa, Kab))
        K2 = np.hstack((Kba, Kbb))
        
        K = np.vstack((K1, K2))
        
        eldof = [ int(DOF[NC[i,1]-1,0]), int(DOF[NC[i,1]-1,1]), int(DOF[NC[i,1]-1,2]), \
                 int(DOF[NC[i,2]-1,0]), int(DOF[NC[i,2]-1,1]), int(DOF[NC[i,2]-1,2]) ]
        
        for ii in range(6):
            for jj in range(6):
                KK[eldof[ii]-1,eldof[jj]-1]+=K[ii,jj]
                
        for ii in range(6):
            fint[eldof[ii]-1,0]+= T_global[ii]
        
    return KK, fint