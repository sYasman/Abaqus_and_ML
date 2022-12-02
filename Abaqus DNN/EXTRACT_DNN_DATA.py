from odbAccess import *
from abaqusConstants import *

from odbMaterial import *
from odbSection import *

import numpy as np
import sys

### FUNCTIONS

def extract_one_data(frame, field_output, elem_type='C3D8T', step_name='Step-1'):
    E11 = []

    STRAIN_11 = odb.steps[step_name].frames[frame].fieldOutputs[field_output]
    STRAIN_11_field = STRAIN_11.getSubset(position=INTEGRATION_POINT, elementType = elem_type)
    STRAIN_11_values = STRAIN_11_field.values
    for v in STRAIN_11_values:
      E11.append(v.data)
      
    return E11
def get_stress(frame):
    # '''
    # A function to get values below for a specific frame.
    # 1. STRAIN(6) := Total strain tensor
    # 2. PLASTIC_STRAIN(6) := Plastic strain tensor
    # 3. eps_pl_bar_pos := Non-local positive equivalent plastic strain (SCALAR)
    # 4. eps_pl_bar_neg := Non-local negative equivalent plastic strain (SCALAR)
    # 5. D_T := Tension damage (SCALAR)
    # 6. D_C := Compression damage (SCALAR)
    # 7. DSDDE(6,6) := Derivative of stress w.r.t. strain tensor
    # 8. dEps_pl_dEps_bar_pos := Derivative of Equivalent plastic strain w.r.t. Non-local positive equivalent plastic strain
    # 9. dEps_pl_dEps_bar_neg := Derivative of Equivalent plastic strain w.r.t. Non-local negative equivalent plastic strain
    # outputs:
    # DDSDDE[comp_1][comp_2][n_ip]: For example 27th integration points D12 tangent term -> DDSDDE[0][1][26]
    # print DDSDDE for 27th integration point is as follows,
    # ip = 26
    # for i in range(6):
    #   for j in range(6):
    #     print(DDSDDE[i][j][ip])
    # '''
## Extract STRAIN values for each integration point
## UVARM 53-58
    STRAIN_11 = extract_one_data(frame, 'UVARM53') # STRAIN_11 component for each integration point
    STRAIN_22 = extract_one_data(frame, 'UVARM54') # STRAIN_22 component for each integration point
    STRAIN_33 = extract_one_data(frame, 'UVARM55') # STRAIN_33 component for each integration point
    STRAIN_12 = extract_one_data(frame, 'UVARM56') # STRAIN_12 component for each integration point
    STRAIN_13 = extract_one_data(frame, 'UVARM57') # STRAIN_13 component for each integration point
    STRAIN_23 = extract_one_data(frame, 'UVARM58') # STRAIN_23 component for each integration point
    STRAIN = [STRAIN_11, STRAIN_22, STRAIN_33, STRAIN_12, STRAIN_13, STRAIN_23]
#
## Extract PLASTIC STRAIN values for each integration point
## UVARM 59-64
    EPS_PL_STRAIN_11 = extract_one_data(frame, 'UVARM59') # EPS_PL_STRAIN_11 component for each integration point
    EPS_PL_STRAIN_22 = extract_one_data(frame, 'UVARM60') # EPS_PL_STRAIN_22 component for each integration point
    EPS_PL_STRAIN_33 = extract_one_data(frame, 'UVARM61') # EPS_PL_STRAIN_33 component for each integration point
    EPS_PL_STRAIN_12 = extract_one_data(frame, 'UVARM62') # EPS_PL_STRAIN_12 component for each integration point
    EPS_PL_STRAIN_13 = extract_one_data(frame, 'UVARM63') # EPS_PL_STRAIN_13 component for each integration point
    EPS_PL_STRAIN_23 = extract_one_data(frame, 'UVARM64') # EPS_PL_STRAIN_23 component for each integration point
    EPS_PL_STRAIN = [EPS_PL_STRAIN_11, EPS_PL_STRAIN_22, EPS_PL_STRAIN_33, EPS_PL_STRAIN_12, EPS_PL_STRAIN_13, EPS_PL_STRAIN_23]
#
## Extract equivalent plastic strain for tension and compression
    EPS_EQ_POS = extract_one_data(frame, 'UVARM3') 
    EPS_EQ_NEG = extract_one_data(frame, 'UVARM4')
## Extract tension and compression damage    
    DM_TENS = extract_one_data(frame, 'UVARM5') 
    DM_COMP = extract_one_data(frame, 'UVARM6')
## Extract DDSDDE values for each integration point
    DDSDDE_11 = extract_one_data(frame, 'UVARM15') 
    DDSDDE_12 = extract_one_data(frame, 'UVARM16') 
    DDSDDE_13 = extract_one_data(frame, 'UVARM17') 
    DDSDDE_14 = extract_one_data(frame, 'UVARM18') 
    DDSDDE_15 = extract_one_data(frame, 'UVARM19') 
    DDSDDE_16 = extract_one_data(frame, 'UVARM20') 
    
    DDSDDE_21 = extract_one_data(frame, 'UVARM21') 
    DDSDDE_22 = extract_one_data(frame, 'UVARM22') 
    DDSDDE_23 = extract_one_data(frame, 'UVARM23') 
    DDSDDE_24 = extract_one_data(frame, 'UVARM24') 
    DDSDDE_25 = extract_one_data(frame, 'UVARM25') 
    DDSDDE_26 = extract_one_data(frame, 'UVARM26') 

    DDSDDE_31 = extract_one_data(frame, 'UVARM27') 
    DDSDDE_32 = extract_one_data(frame, 'UVARM28') 
    DDSDDE_33 = extract_one_data(frame, 'UVARM29') 
    DDSDDE_34 = extract_one_data(frame, 'UVARM30') 
    DDSDDE_35 = extract_one_data(frame, 'UVARM31') 
    DDSDDE_36 = extract_one_data(frame, 'UVARM32') 
    
    DDSDDE_41 = extract_one_data(frame, 'UVARM33') 
    DDSDDE_42 = extract_one_data(frame, 'UVARM34') 
    DDSDDE_43 = extract_one_data(frame, 'UVARM35') 
    DDSDDE_44 = extract_one_data(frame, 'UVARM36') 
    DDSDDE_45 = extract_one_data(frame, 'UVARM37') 
    DDSDDE_46 = extract_one_data(frame, 'UVARM38')     
    
    DDSDDE_51 = extract_one_data(frame, 'UVARM39') 
    DDSDDE_52 = extract_one_data(frame, 'UVARM40') 
    DDSDDE_53 = extract_one_data(frame, 'UVARM41') 
    DDSDDE_54 = extract_one_data(frame, 'UVARM42') 
    DDSDDE_55 = extract_one_data(frame, 'UVARM43') 
    DDSDDE_56 = extract_one_data(frame, 'UVARM44')     
    
    DDSDDE_61 = extract_one_data(frame, 'UVARM45') 
    DDSDDE_62 = extract_one_data(frame, 'UVARM46') 
    DDSDDE_63 = extract_one_data(frame, 'UVARM47') 
    DDSDDE_64 = extract_one_data(frame, 'UVARM48') 
    DDSDDE_65 = extract_one_data(frame, 'UVARM49') 
    DDSDDE_66 = extract_one_data(frame, 'UVARM50')     
    
    DDSDDE = [[DDSDDE_11, DDSDDE_12, DDSDDE_13, DDSDDE_14, DDSDDE_15, DDSDDE_16],[DDSDDE_21, DDSDDE_22, DDSDDE_23, DDSDDE_24, DDSDDE_25, DDSDDE_26],[DDSDDE_31, DDSDDE_32, DDSDDE_33, DDSDDE_34, DDSDDE_35, DDSDDE_36],[DDSDDE_41, DDSDDE_42, DDSDDE_43, DDSDDE_44, DDSDDE_45, DDSDDE_46],[DDSDDE_51, DDSDDE_52, DDSDDE_53, DDSDDE_54, DDSDDE_55, DDSDDE_56],[DDSDDE_61, DDSDDE_62, DDSDDE_63, DDSDDE_64, DDSDDE_65, DDSDDE_66]]
# Extract derivative of equivalent plastic strain w.r.t non-local counter part
    DEPS_PL_DEPS_PL_BAR_POS = extract_one_data(frame, 'UVARM51') 
    DEPS_PL_DEPS_PL_BAR_POS = extract_one_data(frame, 'UVARM52')  
    return STRAIN, EPS_PL_STRAIN, EPS_EQ_POS, EPS_EQ_NEG, DM_TENS, DM_COMP, DDSDDE, DEPS_PL_DEPS_PL_BAR_POS, DEPS_PL_DEPS_PL_BAR_POS
### (1.05000e-003*4.90000e-002)

###
# no more mask
# session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)

NAME = 'SIMPLE_CUBE_NON_LOC'
print(NAME, ' is started.')
NAME_ODB = NAME +'.odb'
odb = openOdb(path=NAME_ODB)



# reach the last step
lastFrame = odb.steps['Step-1'].frames[-1]
E, EP, ep_eq_pos, ep_eq_neg, dm_tens, dm_comp, DDSDDE, deps_pl_deps_pl_bar_pos, deps_pl_deps_pl_bar_neg = get_stress(100)
n_ip = 8000
# n_ip = 3
with open('INPUT-2.csv', 'a') as t1, open('OUTPUT-2.csv', 'a') as t2:
    for f in range(lastFrame.frameId):
        E, EP, ep_eq_pos, ep_eq_neg, dm_tens, dm_comp, DDSDDE, deps_pl_deps_pl_bar_pos, deps_pl_deps_pl_bar_neg = get_stress(f)
        for i in range(n_ip):
            for j in range(6):
                tout = '{}, '.format(E[j][i])
                t1.write(tout) #ith ip, nth component
                
                tout = '{}, '.format(EP[j][i])
                t1.write(tout) #ith ip, nth component    

            tout = '{}, '.format(ep_eq_pos[i])
            t1.write(tout) #ith ip, nth component    

            tout = '{}, '.format(ep_eq_neg[i])
            t1.write(tout) #ith ip, nth component    

            tout = '{}, '.format(dm_tens[i])
            t1.write(tout) #ith ip, nth component    

            tout = '{}, '.format(dm_comp[i])
            t1.write(tout) #ith ip, nth component                
            t1.write('\n')
            for j in range(6):
                for k in range(6):
                    tout = '{}, '.format(DDSDDE[j][k][i])
                    t2.write(tout) #ith ip, nth component
                tout = '{}, '.format(deps_pl_deps_pl_bar_pos[i])
                t2.write(tout) #ith ip, nth component
                
                tout = '{}, '.format(deps_pl_deps_pl_bar_neg[i])
                t2.write(tout) #ith ip, nth component
            t2.write('\n')
                
# for i in range(int(lastFrame.frameId)):
    # get_stress(i)
