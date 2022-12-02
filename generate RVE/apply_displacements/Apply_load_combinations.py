import numpy as np
import csv
# 1 2 4 8
# applied strain
# e11 = -3.0e-2# keep this way. Only modify others
# e12 = 1.5e-4
# e22 = 8.00E-02*0
# E11 = [-3.00E-02, -3.00E-02, -3.00E-02, -3.00E-02, -3.00E-02, -3.00E-02, -3.00E-02,\
       # -3.00E-02, -3.00E-02, -3.00E-02, -3.00E-02, -2.40E-02, -2.00E-02, -1.71E-02,\
       # -1.50E-02, -1.33E-02, -1.20E-02, -1.00E-02, -8.57E-03, -7.50E-03, -6.67E-03,\
       # -6.00E-03, -3.00E-03, -1.20E-03, -6.00E-04]

# E12 = [3.00E-05, 1.50E-04, 7.50E-04, 1.50E-03, 3.00E-03, 7.50E-03, 1.50E-02,\
       # 2.25E-02, 3.00E-02, 4.50E-02, 6.00E-02, 6.00E-02, 6.00E-02, 6.00E-02,\
       # 6.00E-02, 6.00E-02, 6.00E-02, 6.00E-02, 6.00E-02, 6.00E-02, 6.00E-02,\
       # 6.00E-02, 6.00E-02, 6.00E-02, 6.00E-02]
    

# E11 = []
# E12 = []
# with open('ND_E22_E12_data.csv', 'r') as csvfile:
with open('ND_E22_E12_data.csv', 'r') as csvfile:    
    csvfile.readline()
    for r in csvfile:
        e12, e11 = r.split(',')
        e12 = e12.replace('\n', '')
        e11 = e11.replace('\n', '')
        print(float(e12), float(e11))
        E11.append(float(e11))
        E12.append(float(e12))
        

# E12 = [3e-3, 6e-3, 1.2e-2, 1.5e-2, 3.0e-2, 4.5e-2, 6e-2, 9e-2, 1.2e-1, 1.8e-1, 2.4e-1]
# E11 = [3e-2, 3e-2, 3e-2, 3e-2, 3e-2, 3e-2, 3e-2, 3e-2, 3e-2, 3e-2, 3e-2,]

assert len(E12) == len(E11), 'E12 and E11 must have same number of elements'
p_file_names = []
k = 0
wait_in_seconds = 6000
max_analysis = 2
n_cpu = 8

names = []
for f_i, (e11, e12) in enumerate(zip(E11, E12)):    
    # rat = e12 / e11
    eps_ = np.array([[e11, e12/2, 0],
    				[e12/2, e22, 0],
    				[0, 0, 0]])
    
    # distance b/w corner nodes
    # required for imposing displacements
    l1 = np.array([0.0489999987185001, 0.0, 0.0], 'd')
    l2 = np.array([0.0, 0.0489999987185001, 0.0], 'd')
    l3 = np.array([0.0, 0.0, 0.00104999996256083], 'd')
    
    # displacement at node c1, c2, c4 and c5
    u1 = [0, 0, 0]
    u2 = np.dot(eps_, l2)
    u4 = np.dot(eps_, l1)
    u5 = np.dot(eps_, l3)
    
    statement=f'''** 
** 
** BOUNDARY CONDITIONS
** 
** e11 = {e11}
** e12 = {e12}
** e22 = {e22}
** rat = {e12/(e11+1e-9)}
** Name: C1 Type: Displacement/Rotation
*Boundary
CORNER_1, 1, 1, {u1[0]:.6f}
CORNER_1, 2, 2, {u1[1]:.6f}
CORNER_1, 3, 3, {u1[2]:.6f}
** Name: C2 Type: Displacement/Rotation
*Boundary
CORNER_2, 1, 1, {u2[0]:.6f}
CORNER_2, 2, 2, {u2[1]:.6f}
CORNER_2, 3, 3, {u2[2]:.6f}
** Name: C4 Type: Displacement/Rotation
*Boundary
CORNER_4, 1, 1, {u4[0]:.6f}
CORNER_4, 2, 2, {u4[1]:.6f}
CORNER_4, 3, 3, {u4[2]:.6f}
** Name: C5 Type: Displacement/Rotation
*Boundary
CORNER_5, 1, 1, {u5[0]:.6f}
CORNER_5, 2, 2, {u5[1]:.6f}
CORNER_5, 3, 3, {u5[2]:.6f}
    '''
         
    # rat = str(rat).replace('.', '_').replace('-', '')
    # f_name = f'c_new_dist_E12_E22_neg_{rat[:6]}.inp'
    f_name = f'p6_more_gap_{str(f_i)}.inp'
    names.append(f_name)
    f_new = open(f_name, 'w')    
    print(f'{f_name} is being processing..')       
    with open('p6_main_NON_LOC.inp', 'r') as f_read:
        for l in f_read:
            if l == '__REPLACE__\n':
                f_new.writelines(statement)
            else:
                f_new.writelines(l)
        
    f_new.close()
    
    t_ = f'call abaqus j={f_name} user=main cpus={n_cpu} gpus=2\n'
    with open('easy_analysis.bat', 'a') as fout:
        fout.writelines(t_)
        k += 1
        p_file_names.append(f_name.split('.')[0])
        if k > max_analysis:
            fout.writelines(f'TIMEOUT /T {wait_in_seconds}\n\n')
            # fout.writelines(f'TIMEOUT /T {wait_in_seconds}\n\n')
            # fout.writelines(f'TIMEOUT /T {wait_in_seconds}\n\n')
            n_text = ' '.join(p_file_names)
            fout.writelines(f'call abaqus python v3_ALL_EXTRACTER.py {n_text}\n\n\n')
            k = 0
            p_file_names = []
            
            for el in names:
                e = el.split('.')[0]
                fout.writelines(f'call del {e}.odb\n')
            names = []
    
    
    
    
    
    