from odbAccess import *
from abaqusConstants import *

from odbMaterial import *
from odbSection import *

import numpy as np
import sys

### Extract ip data from odb file and store a CSV file
### Be aware of naming convention (file names etc.)
# no more mask
session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)

n_files = len(sys.argv)

for i in range(1, n_files):
    NAMES.append(sys.argv[i])

for NAME in NAMES:
	print(NAME, ' is started.')
	NAME_ODB = NAME +'.odb'
	odb = openOdb(path=NAME_ODB)

	# reach assembly
	myAssembly = odb.rootAssembly

	# parts in assembly
	MATRIX_PART = odb.rootAssembly.instances.keys()[0]
	FIBER_PART = odb.rootAssembly.instances.keys()[1]

	MATRIX_ASSEMBLY = myAssembly.instances[MATRIX_PART]
	FIBER_ASSEMBLY = myAssembly.instances[FIBER_PART]
	### MATRIX PART calculations
	Matrix_elements = MATRIX_ASSEMBLY.elementSets['DUMMY_SET']
	# Fiber_elements = FIBER_ASSEMBLY.elementSets['FIBER_PROP_PART']

	# reach the last step
	lastFrame = odb.steps['LOADING_STEP'].frames[-1]

	# variables avaliable
	fields = lastFrame.fieldOutputs.keys()

	# Calculate VTOTAL
	# TOTAL volume of RVE
	# simply summation of the volume of all individual integration points
	# Error check
	# volume of the RVE is 4.9e-2 * 4.9e-2 * 1.05e-3 = VTOT
	VTOT=lastFrame.fieldOutputs['IVOL']
	VTOT_val = VTOT.values
	TOTAL_VOLUME = 0.00
	for v in VTOT_val:
		TOTAL_VOLUME += v

	print('Total volume is: ', TOTAL_VOLUME.data)	

	frame=-1

	def get_stress(frame):
		MATRIX_SET = Matrix_elements
		## get UVARM7-UVARM12 := S11,S22,S33,S12,S13,S23 (For UEL)
		## Whole Model
		UVARM1_FIELD = odb.steps['LOADING_STEP'].frames[frame].fieldOutputs['UVARM7']
		UVARM2_FIELD = odb.steps['LOADING_STEP'].frames[frame].fieldOutputs['UVARM8']
		UVARM3_FIELD = odb.steps['LOADING_STEP'].frames[frame].fieldOutputs['UVARM9']
		UVARM4_FIELD = odb.steps['LOADING_STEP'].frames[frame].fieldOutputs['UVARM10']
		UVARM5_FIELD = odb.steps['LOADING_STEP'].frames[frame].fieldOutputs['UVARM11']
		UVARM6_FIELD = odb.steps['LOADING_STEP'].frames[frame].fieldOutputs['UVARM12']

		
		## get Stress values S11, S22, S33, S12, S13, S23 (For fibers)
		## Whole Model
		STRESS_FIELD = odb.steps['LOADING_STEP'].frames[frame].fieldOutputs['S']
		
		## get volume related to Integration Point (Whole model)
		IVOL_FIELD = odb.steps['LOADING_STEP'].frames[frame].fieldOutputs['IVOL']

		### Extract the value of UVARM for each integration point in C3D8T (Dummy elements)
		## Stress from Matrix part
		field = UVARM1_FIELD.getSubset(position=INTEGRATION_POINT, elementType = 'C3D8T')
		UVARM1_VALUES = field.values
		field = UVARM2_FIELD.getSubset(position=INTEGRATION_POINT, elementType = 'C3D8T')
		UVARM2_VALUES = field.values
		field = UVARM3_FIELD.getSubset(position=INTEGRATION_POINT, elementType = 'C3D8T')
		UVARM3_VALUES = field.values
		field = UVARM4_FIELD.getSubset(position=INTEGRATION_POINT, elementType = 'C3D8T')
		UVARM4_VALUES = field.values
		field = UVARM5_FIELD.getSubset(position=INTEGRATION_POINT, elementType = 'C3D8T')
		UVARM5_VALUES = field.values
		field = UVARM6_FIELD.getSubset(position=INTEGRATION_POINT, elementType = 'C3D8T')
		UVARM6_VALUES = field.values

		### Extract Volume data from Matrix Part
		## Integration Point Volume related to matrix part
		field = IVOL_FIELD.getSubset(region=MATRIX_SET,position=INTEGRATION_POINT, elementType = 'C3D8T')
		IVOL_VALUES = field.values

		### FIBER CALCULATIONS
		## Stress values in fiber integration points
		## If element type changes change related parts
		field = STRESS_FIELD.getSubset(position=INTEGRATION_POINT, elementType = 'C3D8R')
		STRESS_VALUES = field.values

		## Extract integration point volume for fiber part
		field = IVOL_FIELD.getSubset(position=INTEGRATION_POINT, elementType = 'C3D8R')
		IVOL_FIBER_VALUES = field.values

		### Create hom stress fields
		S11 = 0.0
		S22 = 0.0
		S33 = 0.0
		S12 = 0.0
		S13 = 0.0
		S23 = 0.0

	# S_hom = 1/V * w * S * detJ 
	# S := Stress value at integration point S11, S22, S33, S12, S13, S23
	# w = 1 for all integrations points. Hence it is neglected in code
	# detJ = volume of each integtraion point. Stored in IVOL. Components of IVOL are related to integration pt.s
		for i, v in enumerate(STRESS_VALUES):
			vol = IVOL_FIBER_VALUES[i].data #volume of each integraion point (equals detJ)
			S11_temp = v.data[0]
			S22_temp = v.data[1]
			S33_temp = v.data[2]
			S12_temp = v.data[3]
			S13_temp = v.data[4]
			S23_temp = v.data[5]
			S11 +=  vol * S11_temp
			S22 +=  vol * S22_temp
			S33 +=  vol * S33_temp
			S12 +=  vol * S12_temp
			S13 +=  vol * S13_temp
			S23 +=  vol * S23_temp	
			

		for i, v in enumerate(UVARM1_VALUES):
			vol = IVOL_VALUES[i].data
			S11 +=  vol * UVARM1_VALUES[i].data
			S22 +=  vol * UVARM2_VALUES[i].data
			S33 +=  vol * UVARM3_VALUES[i].data
			S12 +=  vol * UVARM4_VALUES[i].data
			S13 +=  vol * UVARM5_VALUES[i].data
			S23 +=  vol * UVARM6_VALUES[i].data

		### S_hom = 1/V * S_ALL
		S11 = S11/TOTAL_VOLUME.data
		S22 = S22/TOTAL_VOLUME.data
		S33 = S33/TOTAL_VOLUME.data
		S12 = S12/TOTAL_VOLUME.data
		S13 = S13/TOTAL_VOLUME.data
		S23 = S23/TOTAL_VOLUME.data

		return S11, S22, S33, S12, S13, S23
	### (1.05000e-003*4.90000e-002)

	S = []
	last_frame = 45
	last_frame = int(odb.steps['LOADING_STEP'].frames[-1].description.split()[1][0:3].replace(':',''))
	n_incr=1
	n_read=2
	for i in range(1, int(last_frame/n_read), n_incr):
		S.append(get_stress(i))
		print('step', i, 'is processed')

	max_S11 = -10.0
	max_S22 = -10.0
	max_S33 = -10.0
	max_S12 = -10.0
	max_S13 = -10.0
	max_S23 = -10.0

	for s in S:
		if s[0] > max_S11:
			max_S11 = s[0]
		if s[1] > max_S22:
			max_S22 = s[1]
		if s[2] > max_S33:
			max_S33 = s[2]   
		if s[3] > max_S12:
			max_S12 = s[3]
		if s[4] > max_S13:
			max_S13 = s[4]
		if s[5] > max_S23:
			max_S23 = s[5]           

	print('#########')
	print('This is', NAME_ODB)
	print( max_S11, max_S22, max_S33, max_S12, max_S13, max_S23) 
	odb.close()
	print('#########')
	f_name = NAME+'.csv'
	np.savetxt(f_name, S, delimiter=',')
