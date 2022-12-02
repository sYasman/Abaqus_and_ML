### Place dummy element on top of Abaqus Elements
### INP file must be created accordingly.
### Be aware of naming conventions
import numpy as np
np.set_printoptions(precision=0)
np.set_printoptions(suppress=True)


# functions
def increase_el_nums(UEL):
    n_UEL = len(UEL)
    max_ = -1
    for i in range(n_UEL):
        UEL[i] = UEL[i].split(',')
        UEL[i] = np.asarray(UEL[i], dtype=np.float64, order='C')
        UEL[i][0] += n_UEL
        for el in UEL[i]:
            if el > max_:
                max_ = int(el)
    return UEL, n_UEL, max_   

# Constants
name_ = 'p6_main'
file_name = name_+'.inp'
file_name_2 = name_+'_NON_LOC.inp'

dummy_el ='type=C3D8T'
dummy_el_replace_1 = '''*USER ELEMENT, NODES=8, TYPE=U1, COORDINATES=3, 
 PROPERTIES=10, VARIABLES=112, UNSYMM
1, 2, 3, 11, 12   
*Element, type=U1, ELSET=USET1\n''' #+3 rows

dummy_element_header = '''*Element, type=C3D8T, ELSET=dummy_set\n'''#+1 row
# material properties
E = 3760.000
nu = 0.390
lc = 0.00150
lc = 0.00300
lc = 0.00075
lc = 0.00020
lc = 1e-3
aa = 0.950
bb = 10.000
# bb = 25.000
nu_p = 0.33

ALL_RIGHT = input('Did you update material properties (y/n)')

if ALL_RIGHT == 'y':
    with open(file_name, 'r') as input_uel:
        lines = input_uel.readlines()
      
    num_lines = len(lines)
    start_found = False
    for i in range(num_lines):
        # replace C3D8T with user element data
        line = lines[i]
        if(line.split(',')[-1].strip() == dummy_el):
            lines[i]=dummy_el_replace_1
            start_found=True
            start_line=i+1
        if(start_found and line.split(',')[0].strip() == '*Nset'):
            end_line=i
            start_found = False
            
    dummy_elements, n_UEL, max_el_num = increase_el_nums(lines[start_line: end_line])
    dummy_elements = np.round(dummy_elements, decimals=0)
    lines.insert(end_line, dummy_element_header)
    j = 0
    ch_to_remove = ['[', ']', ' ']
    str_to_print = ''
    for i in range(n_UEL):
        my_text = np.array_str(dummy_elements[j])
        for ch in ch_to_remove:
            my_text = my_text.replace(ch, "")
        my_text = my_text.split('.')[0:-1]
        for el in my_text:
            str_to_print += str(f'{el}') + ','
        lines.insert(end_line+1+i, str_to_print[0:-1]+'\n')
        str_to_print = ''
        j += 1
    
    soldi_part_1 = f'''*Elset, elset=dummy_set, internal, generate   
    {n_UEL+1}, {2*n_UEL}, 1
    *UEL PROPERTY, ELSET=USET1
    ** N_U: not used
    ** OLD: E, nu, sigma_c_0 (N_U), sigma_t_0 (N_U), sigma_c_inf(N_U), sigma_t_inf(N_U), nu_p, e1p(N_U), e2p(N_U), e3p(N_U)
    ** NEW: E, nu, lc, aa, bb, sigma_t_inf(N_U), nu_p, e1p(N_U), e2p(N_U), e3p(N_U)
    ** keep last three PROPS. They may be required for rotation.
    {E}, {nu}, {lc},  {aa}, {bb},  90., 0.33, 1.0, 0.0, 0.0\n''' 
       
    # modifiy section assignmenet
    num_lines = len(lines)    
    for i in range(num_lines):
        if lines[i].split(' ')[-1].strip() == 'material=MATRIX_MAT':
            start_line = i
    lines.insert(start_line-1, soldi_part_1)
    el_set = lines[start_line+1].split(',')
    new_el_set = lines[start_line+1].split(',')[0] + ', elset=dummy_set, ' + lines[start_line+1].split(',')[-1]
    lines[start_line+1] = new_el_set

    # create new inpu file    
    with open(file_name_2, 'w') as output_uel:
        output_uel.writelines(lines)

    print(f'Original inp file: {file_name}')
    print(f'New inp file: {file_name_2}')
    print(f'Number of elements: {n_UEL}')        
else:
    print('Please update material properties')
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    