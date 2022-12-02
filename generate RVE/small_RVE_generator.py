from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *

import numpy as np
#
# IMPORTANT
# If the geometry s changed then sets must be redefined
#
# some parameters to define the model
INCR_INIT = 0.001
INCR_MAX = 0.001
INCR_MIN = 1e-9
INCR_MAX_NUM = 5000
STAB_MAGN = 0.002
displ_y = 0.002

mesh_size = 0.0010 # Default mesh size
mesh_size = 0.00075 # Default mesh size
# mesh_size = 0.00025
# mesh_size = 0.02 # Default mesh size
# mesh_size = 0.00025 # Default mesh size
# mesh_size = 0.0002 # Default mesh size
# fiber_mesh_size = mesh_size*2.00

UD_FIBERS_reduced = False
mesh_size = 0.00050 # MODEL_II
# mesh_size = 0.000250 # MODEL_III
# mesh_size = 0.0003250 # MODEL_III
# mesh_size = 0.00050 # MODEL_III

large_num = 25
med_num = 20
small_num = 10
small_num_mat = 5
med_num_mat = 25
large_num_mat = 30	

# large_num = int(25*0.80)
# med_num = int(20*0.80)
# small_num = int(10*0.80)
# small_num_mat = int(5*0.80)
# med_num_mat = int(25*0.80)
# large_num_mat = int(30*0.80)	

FIXED_SIZE = mesh_size
# FIXED_SIZE = 1.05000e-003 #half of the depth
# FIXED_SIZE = mesh_size
fiber_mesh_size = FIXED_SIZE*2.00
matrix_depth_el = 2

remove_other_assemblies = True
transverse_multipler = 5.00
long_multiplier = 0.15
fiber_diameter = 5e-3
matrix_part_depth = long_multiplier*fiber_diameter
# matrix_part_depth = fiber_diameter
# matrix rectangle pt1 (bot_left), pt2 (top_right)
pt_1 = (0.0, 0.0)
pt_2 = (50e-3, 50e-3)

# matrix part params
almost_zero = 2e-12
user_mat_outputs = 7
E_fiber = 74e3
pois_fiber = 0.20

session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)

# matrix geom
MODEL = mdb.models['Model-1']
MODEL.ConstrainedSketch(name='__profile__', sheetSize=0.01)
MODEL.sketches['__profile__'].sketchOptions.setValues(decimalPlaces=4)
MODEL.sketches['__profile__'].rectangle(point1=pt_1, point2=pt_2)
MODEL.Part(dimensionality=THREE_D, name='MATRIX', type=DEFORMABLE_BODY)
MODEL.parts['MATRIX'].BaseSolidExtrude(depth=matrix_part_depth, sketch=MODEL.sketches['__profile__'])
del MODEL.sketches['__profile__']

## define UD_FIBERS
MODEL.ConstrainedSketch(name='__profile__', sheetSize=0.01)
MODEL.sketches['__profile__'].sketchOptions.setValues(decimalPlaces=4)

# dummy rectangle    
MODEL.sketches['__profile__'].rectangle(point1=pt_1, point2=pt_2)

coords = [[-0.0018651575676875264, 0.051198402130492676],
[0.0007143364403018205, 0.03581890812250333],
[0.0037660896582334662, 0.02463382157123835],
[0.008216511318242343, 0.012450066577896136],
[0.01013368841544607, 0.05033288948069241],
[0.012784110075454944, 0.0381491344873502],
[0.01570155348424323, 0.026231691078561915],
[0.020018020417221477, 0.013581890812250336],
[0.0220023080337328, 0.05179760319573901],
[0.024784110075454953, 0.0381491344873502],
[0.027633821571238344, 0.02543275632490013],
[0.032016156236129606, 0.012183754993342212],
[0.03393306702174878, 0.049866844207723034],
[0.0365167332445628, 0.03761651131824234],
[0.03956564580559254, 0.024300932090545936],
[0.04401633377718597, 0.012316910785619181]]


fiber_diameter = 4.80e-3*2
MODEL.sketches['__profile__'].CircleByCenterPerimeter(center=(0.0, 0.0), point1=(0.0, fiber_diameter/2.00))
    
for coord in coords:
  MODEL.sketches['__profile__'].copyMove(objectList=(
    MODEL.sketches['__profile__'].geometry[6], ), vector=coord)
# BURAYA linear pattern gelecek
mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((0.0, -0.0048))
mdb.models['Model-1'].sketches['__profile__'].delete(objectList=(
    mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((0.0, 
    -0.0048), ), ))
mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((-0.001865, 
    0.046398))
mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((0.010134, 
    0.045533))
mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((0.022002, 
    0.046998))
mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((0.033933, 
    0.045067))
mdb.models['Model-1'].sketches['__profile__'].linearPattern(angle1=0.0, angle2=
    270.0, geomList=(
    mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((-0.001865, 
    0.046398), ), 
    mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((0.010134, 
    0.045533), ), 
    mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((0.022002, 
    0.046998), ), 
    mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((0.033933, 
    0.045067), )), number1=1, number2=2, spacing1=0.05, spacing2=0.05, 
    vertexList=())
mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((-0.001865, 
    0.046398))
mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((0.000714, 
    0.031019))
mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((0.003766, 
    0.019834))
mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((-0.001865, 
    -0.003602))
mdb.models['Model-1'].sketches['__profile__'].linearPattern(angle1=0.0, angle2=
    270.0, geomList=(
    mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((-0.001865, 
    0.046398), ), 
    mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((0.000714, 
    0.031019), ), 
    mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((0.003766, 
    0.019834), ), 
    mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((-0.001865, 
    -0.003602), )), number1=2, number2=1, spacing1=0.05, spacing2=0.05, 
    vertexList=())
mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((0.0, 0.025))
mdb.models['Model-1'].sketches['__profile__'].delete(objectList=(
    mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((0.0, 0.025), 
    ), ))
mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((0.025, 0.05))
mdb.models['Model-1'].sketches['__profile__'].delete(objectList=(
    mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((0.025, 
    0.05), ), ))
mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((0.05, 0.025))
mdb.models['Model-1'].sketches['__profile__'].delete(objectList=(
    mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((0.05, 
    0.025), ), ))
mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((0.025, 0.0))
mdb.models['Model-1'].sketches['__profile__'].delete(objectList=(
    mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((0.025, 0.0), 
    ), ))
mdb.models['Model-1'].Part(dimensionality=THREE_D, name='FIBER_PART', type=
    DEFORMABLE_BODY)
mdb.models['Model-1'].parts['FIBER_PART'].BaseSolidExtrude(depth=matrix_part_depth, sketch=
    mdb.models['Model-1'].sketches['__profile__'])
del mdb.models['Model-1'].sketches['__profile__']

## Define materials and sections
MODEL.Material(name='FIBER_MAT')
MODEL.materials['FIBER_MAT'].Elastic(table=((E_fiber, pois_fiber), ))
MODEL.Material(name='MATRIX_MAT')
MODEL.materials['MATRIX_MAT'].UserOutputVariables(n=user_mat_outputs)
MODEL.materials['MATRIX_MAT'].Elastic(table=((almost_zero, 0.0), ))
MODEL.materials['MATRIX_MAT'].Conductivity(table=((almost_zero, ), ))
MODEL.HomogeneousSolidSection(material='FIBER_MAT', name='FIBER_SECTION', thickness=None)
MODEL.HomogeneousSolidSection(material='MATRIX_MAT', name='MATRIX_SECTION', thickness=None)

# TODO: REDINE ONCE THE GEOMETRY IS CHANGED
# fiber faces on XY plane
mdb.models['Model-1'].parts['FIBER_PART'].Set(cells=
    mdb.models['Model-1'].parts['FIBER_PART'].cells.findAt(((0.000714, 
    0.040375, 0.0), ), ((0.008217, 0.017006, 0.0), ), ((0.012784, 0.042706, 
    0.0), ), ((0.020018, 0.018138, 0.0), ), ((0.024784, 0.042706, 0.0), ), ((
    0.032016, 0.01674, 0.0), ), ((0.036517, 0.042173, 0.0), ), ((0.044016, 
    0.016873, 0.0), ), ((0.010134, 0.004889, 0.0), ), ((0.033933, 0.004423, 
    0.0), ), ((0.050714, 0.040375, 0.0), ), ((0.048135, 0.005755, 0.0), ), ((
    0.053766, 0.02919, 0.0), ), ((0.048135, 0.055755, 0.0), ), ((0.022002, 
    0.006354, 0.0), ), ((-0.001865, 0.005755, 0.0), ), ((0.039566, 0.028857, 
    0.0), ), ((0.033933, 0.054423, 0.0), ), ((0.027634, 0.029989, 0.0), ), ((
    0.022002, 0.056354, 0.0), ), ((0.015702, 0.030788, 0.0), ), ((0.010134, 
    0.054889, 0.0), ), ((0.003766, 0.02919, 0.0), ), ((-0.001865, 0.055755, 
    0.0), ), ), name='FIBER_PROP_PART')

aa = MODEL.parts['FIBER_PART'].sets['FIBER_PROP_PART'].faces.pointsOn
for i in range(len(aa)):
    MODEL.parts['FIBER_PART'].SectionAssignment(offset=0.0, 
        offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
        cells=MODEL.parts['FIBER_PART'].cells.findAt(aa[i])     
        ), sectionName='FIBER_SECTION', 
        thicknessAssignment=FROM_SECTION)

MODEL.parts['MATRIX'].SectionAssignment(offset=0.0, 
        offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
        cells=MODEL.parts['MATRIX'].cells.findAt(((0,0,0,), ), )), sectionName='MATRIX_SECTION', thicknessAssignment=
        FROM_SECTION)

# Assembly
# 1. Matrix with holes

# shorten the lines
r = MODEL.rootAssembly

r.DatumCsysByDefault(CARTESIAN)
r.Instance(dependent=ON, name='FIBER_PART-1', part=MODEL.parts['FIBER_PART'])
r.Instance(dependent=ON, name='MATRIX-1', part=MODEL.parts['MATRIX'])
r.InstanceFromBooleanCut(cuttingInstances=(
    r.instances['FIBER_PART-1'], ), 
    instanceToBeCut=r.instances['MATRIX-1'], 
    name='MATRIX_WITH_HOLES', originalInstances=SUPPRESS)    

r.deleteFeatures(('MATRIX-1', 'FIBER_PART-1'))
r.Instance(dependent=ON, name='FIBER_PART-1', 
    part=MODEL.parts['FIBER_PART'])
r.Instance(dependent=ON, name='MATRIX-1', 
    part=MODEL.parts['MATRIX'])
r.InstanceFromBooleanCut(cuttingInstances=(
    r.instances['MATRIX-1'], ), 
    instanceToBeCut=
    r.instances['FIBER_PART-1'], name=
    'Fiber_negative', originalInstances=SUPPRESS)
r.deleteFeatures(('MATRIX-1', 'FIBER_PART-1'))
r.Instance(dependent=ON, name='FIBER_PART-1', 
    part=MODEL.parts['FIBER_PART'])
r.InstanceFromBooleanCut(cuttingInstances=(
    r.instances['Fiber_negative-1'], ), 
    instanceToBeCut=
    r.instances['FIBER_PART-1'], name='UD_FIBERS'
    , originalInstances=SUPPRESS)
r.deleteFeatures(('FIBER_PART-1', 
    'Fiber_negative-1'))

if(remove_other_assemblies):    
  del MODEL.parts['FIBER_PART']
  del MODEL.parts['Fiber_negative']
  del MODEL.parts['MATRIX']

# Remove old parts and assemblies 
# Makes easy to modify inp file    

# # Mesh
# r.regenerate()  
# # 1. UD_FIBERS
# #reduced integration in UD_FIBERS
r.regenerate()

# # # TODO: REDINE ONCE THE GEOMETRY IS CHANGED
# # # Fiber face on XY plane
mdb.models['Model-1'].parts['UD_FIBERS'].Set(cells=
    mdb.models['Model-1'].parts['UD_FIBERS'].cells.findAt(((0.000889, 0.040487, 
    0.00075), ), ((0.008217, 0.017006, 0.0), ), ((0.012784, 0.042706, 0.0), ), 
    ((0.020018, 0.018138, 0.0), ), ((0.024784, 0.042706, 0.0), ), ((0.032016, 
    0.01674, 0.0), ), ((0.036517, 0.042173, 0.0), ), ((0.044016, 0.016873, 
    0.0), ), ((0.008534, 0.000111, 0.00075), ), ((0.032468, 0.000599, 0.00075), 
    ), ((0.049464, 0.037219, 0.00075), ), ((0.04959, 0.003859, 0.00075), ), ((
    0.049744, 0.025159, 0.00075), ), ((0.045883, 0.049502, 0.00075), ), ((
    0.020402, 0.000599, 0.00075), ), ((0.000475, 0.003438, 0.00075), ), ((
    0.039566, 0.028857, 0.0), ), ((0.035448, 0.049441, 0.00075), ), ((0.027634, 
    0.029989, 0.0), ), ((0.023202, 0.049541, 0.00075), ), ((0.015702, 0.030788, 
    0.0), ), ((0.011578, 0.049423, 0.00075), ), ((0.005646, 0.02338, 0.00075), 
    ), ((0.001663, 0.049551, 0.00075), ), ), name='Set-1')

aa2 = MODEL.parts['UD_FIBERS'].sets['Set-1'].faces.pointsOn

mdb.models['Model-1'].parts['UD_FIBERS'].SectionAssignment(offset=0.0, 
    offsetField='', offsetType=MIDDLE_SURFACE, region=
    mdb.models['Model-1'].parts['UD_FIBERS'].sets['FIBER_PROP_PART'], 
    sectionName='FIBER_SECTION', thicknessAssignment=FROM_SECTION)
mdb.models['Model-1'].rootAssembly.regenerate()

for i in range(len(aa2)) :
    MODEL.parts['UD_FIBERS'].setElementType(elemTypes=(ElemType(
        elemCode=C3D8, elemLibrary=STANDARD), ElemType(elemCode=C3D6, 
        elemLibrary=STANDARD), ElemType(elemCode=C3D4, elemLibrary=STANDARD)), 
        regions=(MODEL.parts['UD_FIBERS'].cells.findAt(aa2[i]), )
        )
        
    MODEL.parts['UD_FIBERS'].setMeshControls(algorithm=MEDIAL_AXIS, 
        minTransition=OFF, regions=
        MODEL.parts['UD_FIBERS'].cells.findAt(aa2[i]), technique=SWEEP)

MODEL.parts['UD_FIBERS'].seedPart(deviationFactor=0.1, 
    minSizeFactor=0.1, size=fiber_mesh_size)
	
## if mesh sizes are not the same divide fibers	
mdb.models['Model-1'].rootAssembly.regenerate()

mdb.models['Model-1'].parts['UD_FIBERS'].PartitionCellByPlanePointNormal(cells=
    mdb.models['Model-1'].parts['UD_FIBERS'].cells.findAt(((0.008217, 0.017006, 
    0.0), )), normal=mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt((
    0.0, 0.026122, 0.00075), ), point=
    mdb.models['Model-1'].parts['UD_FIBERS'].InterestingPoint(
    mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt((0.013017, 0.01245, 
    0.00075), ), CENTER))
mdb.models['Model-1'].parts['UD_FIBERS'].PartitionCellByPlanePointNormal(cells=
    mdb.models['Model-1'].parts['UD_FIBERS'].cells.findAt(((0.020018, 0.018138, 
    0.0), )), normal=mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt((
    0.0, 0.026122, 0.00075), ), point=
    mdb.models['Model-1'].parts['UD_FIBERS'].InterestingPoint(
    mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt((0.024818, 0.013582, 
    0.00075), ), CENTER))
mdb.models['Model-1'].parts['UD_FIBERS'].PartitionCellByPlanePointNormal(cells=
    mdb.models['Model-1'].parts['UD_FIBERS'].cells.findAt(((0.032016, 0.01674, 
    0.0), )), normal=mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt((
    0.05, 0.023146, 0.00075), ), point=
    mdb.models['Model-1'].parts['UD_FIBERS'].InterestingPoint(
    mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt((0.036816, 0.012184, 
    0.00075), ), CENTER))
mdb.models['Model-1'].parts['UD_FIBERS'].PartitionCellByPlanePointNormal(cells=
    mdb.models['Model-1'].parts['UD_FIBERS'].cells.findAt(((0.044016, 0.016873, 
    0.0), )), normal=mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt((
    0.05, 0.023146, 0.00075), ), point=
    mdb.models['Model-1'].parts['UD_FIBERS'].InterestingPoint(
    mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt((0.048816, 0.012317, 
    0.00075), ), CENTER))
mdb.models['Model-1'].parts['UD_FIBERS'].PartitionCellByPlanePointNormal(cells=
    mdb.models['Model-1'].parts['UD_FIBERS'].cells.findAt(((0.039566, 0.028857, 
    0.0), )), normal=mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt((
    0.05, 0.023146, 0.00075), ), point=
    mdb.models['Model-1'].parts['UD_FIBERS'].InterestingPoint(
    mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt((0.044366, 0.024301, 
    0.00075), ), CENTER))
mdb.models['Model-1'].parts['UD_FIBERS'].PartitionCellByPlanePointNormal(cells=
    mdb.models['Model-1'].parts['UD_FIBERS'].cells.findAt(((0.027634, 0.029989, 
    0.0), )), normal=mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt((
    0.05, 0.023146, 0.00075), ), point=
    mdb.models['Model-1'].parts['UD_FIBERS'].InterestingPoint(
    mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt((0.032434, 0.025433, 
    0.00075), ), CENTER))
mdb.models['Model-1'].parts['UD_FIBERS'].PartitionCellByPlanePointNormal(cells=
    mdb.models['Model-1'].parts['UD_FIBERS'].cells.findAt(((0.015702, 0.030788, 
    0.0), )), normal=mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt((
    0.05, 0.023146, 0.00075), ), point=
    mdb.models['Model-1'].parts['UD_FIBERS'].InterestingPoint(
    mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt((0.020502, 0.026232, 
    0.00075), ), CENTER))
mdb.models['Model-1'].parts['UD_FIBERS'].PartitionCellByPlanePointNormal(cells=
    mdb.models['Model-1'].parts['UD_FIBERS'].cells.findAt(((0.005646, 0.02338, 
    0.00075), )), normal=mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt(
    (0.0, 0.026122, 0.00075), ), point=
    mdb.models['Model-1'].parts['UD_FIBERS'].InterestingPoint(
    mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt((0.000761, 0.028377, 
    0.00075), ), CENTER))
mdb.models['Model-1'].parts['UD_FIBERS'].PartitionCellByPlanePointNormal(cells=
    mdb.models['Model-1'].parts['UD_FIBERS'].cells.findAt(((0.000889, 0.040487, 
    0.00075), )), normal=mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt(
    (0.0, 0.038192, 0.00075), ), point=
    mdb.models['Model-1'].parts['UD_FIBERS'].InterestingPoint(
    mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt((0.000178, 0.040589, 
    0.00075), ), CENTER))
mdb.models['Model-1'].parts['UD_FIBERS'].PartitionCellByPlanePointNormal(cells=
    mdb.models['Model-1'].parts['UD_FIBERS'].cells.findAt(((0.012784, 0.042706, 
    0.0), )), normal=mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt((
    0.0, 0.039379, 0.00075), ), point=
    mdb.models['Model-1'].parts['UD_FIBERS'].InterestingPoint(
    mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt((0.017584, 0.038149, 
    0.00075), ), CENTER))
mdb.models['Model-1'].parts['UD_FIBERS'].PartitionCellByPlanePointNormal(cells=
    mdb.models['Model-1'].parts['UD_FIBERS'].cells.findAt(((0.024784, 0.042706, 
    0.0), )), normal=mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt((
    0.05, 0.033446, 0.00075), ), point=
    mdb.models['Model-1'].parts['UD_FIBERS'].InterestingPoint(
    mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt((0.029584, 0.038149, 
    0.00075), ), CENTER))
mdb.models['Model-1'].parts['UD_FIBERS'].PartitionCellByPlanePointNormal(cells=
    mdb.models['Model-1'].parts['UD_FIBERS'].cells.findAt(((0.036517, 0.042173, 
    0.0), )), normal=mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt((
    0.05, 0.033446, 0.00075), ), point=
    mdb.models['Model-1'].parts['UD_FIBERS'].InterestingPoint(
    mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt((0.041317, 0.037617, 
    0.00075), ), CENTER))
mdb.models['Model-1'].parts['UD_FIBERS'].PartitionCellByPlanePointNormal(cells=
    mdb.models['Model-1'].parts['UD_FIBERS'].cells.findAt(((0.008534, 0.000111, 
    0.00075), )), normal=mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt(
    (0.007739, 0.0, 0.00075), ), point=
    mdb.models['Model-1'].parts['UD_FIBERS'].InterestingPoint(
    mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt((0.005609, 0.001937, 
    0.00075), ), CENTER))
mdb.models['Model-1'].parts['UD_FIBERS'].PartitionCellByPlanePointNormal(cells=
    mdb.models['Model-1'].parts['UD_FIBERS'].cells.findAt(((0.020402, 0.000599, 
    0.00075), )), normal=mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt(
    (0.019777, 0.0, 0.00075), ), point=
    mdb.models['Model-1'].parts['UD_FIBERS'].InterestingPoint(
    mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt((0.017229, 0.0023, 
    0.00075), ), CENTER))
mdb.models['Model-1'].parts['UD_FIBERS'].PartitionCellByPlanePointNormal(cells=
    mdb.models['Model-1'].parts['UD_FIBERS'].cells.findAt(((0.032468, 0.000599, 
    0.00075), )), normal=mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt(
    (0.031534, 0.0, 0.00075), ), point=
    mdb.models['Model-1'].parts['UD_FIBERS'].InterestingPoint(
    mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt((0.031534, 0.0, 
    0.00075), ), MIDDLE))
mdb.models['Model-1'].parts['UD_FIBERS'].PartitionCellByPlanePointNormal(cells=
    mdb.models['Model-1'].parts['UD_FIBERS'].cells.findAt(((0.045495, 0.012929, 
    0.0), ), ((0.042538, 0.011705, 0.0), ), ), normal=
    mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt((0.046416, 0.012317, 
    0.00075), ), point=
    mdb.models['Model-1'].parts['UD_FIBERS'].InterestingPoint(
    mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt((0.039582, 0.014154, 
    0.00075), ), CENTER))
mdb.models['Model-1'].parts['UD_FIBERS'].PartitionCellByPlanePointNormal(cells=
    mdb.models['Model-1'].parts['UD_FIBERS'].cells.findAt(((0.033494, 0.012796, 
    0.0), ), ((0.030538, 0.011571, 0.0), ), ), normal=
    mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt((0.034416, 0.012184, 
    0.00075), ), point=
    mdb.models['Model-1'].parts['UD_FIBERS'].InterestingPoint(
    mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt((0.027582, 0.014021, 
    0.00075), ), CENTER))
mdb.models['Model-1'].parts['UD_FIBERS'].PartitionCellByPlanePointNormal(cells=
    mdb.models['Model-1'].parts['UD_FIBERS'].cells.findAt(((0.021496, 0.014194, 
    0.0), ), ((0.01854, 0.01297, 0.0), ), ), normal=
    mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt((0.017618, 0.013582, 
    0.00075), ), point=
    mdb.models['Model-1'].parts['UD_FIBERS'].InterestingPoint(
    mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt((0.015583, 0.015419, 
    0.00075), ), CENTER))
mdb.models['Model-1'].parts['UD_FIBERS'].PartitionCellByPlanePointNormal(cells=
    mdb.models['Model-1'].parts['UD_FIBERS'].cells.findAt(((0.009695, 0.013062, 
    0.0), ), ((0.006738, 0.011838, 0.0), ), ), normal=
    mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt((0.005817, 0.01245, 
    0.00075), ), point=
    mdb.models['Model-1'].parts['UD_FIBERS'].InterestingPoint(
    mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt((0.003782, 0.014287, 
    0.00075), ), CENTER))
mdb.models['Model-1'].parts['UD_FIBERS'].PartitionCellByPlanePointNormal(cells=
    mdb.models['Model-1'].parts['UD_FIBERS'].cells.findAt(((0.01718, 0.026844, 
    0.0), ), ((0.014223, 0.025619, 0.0), ), ), normal=
    mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt((0.018102, 0.026232, 
    0.00075), ), point=
    mdb.models['Model-1'].parts['UD_FIBERS'].InterestingPoint(
    mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt((0.011267, 0.028069, 
    0.00075), ), CENTER))
mdb.models['Model-1'].parts['UD_FIBERS'].PartitionCellByPlanePointNormal(cells=
    mdb.models['Model-1'].parts['UD_FIBERS'].cells.findAt(((0.029112, 0.026045, 
    0.0), ), ((0.026156, 0.02482, 0.0), ), ), normal=
    mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt((0.030034, 0.025433, 
    0.00075), ), point=
    mdb.models['Model-1'].parts['UD_FIBERS'].InterestingPoint(
    mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt((0.023199, 0.02727, 
    0.00075), ), CENTER))
mdb.models['Model-1'].parts['UD_FIBERS'].PartitionCellByPlanePointNormal(cells=
    mdb.models['Model-1'].parts['UD_FIBERS'].cells.findAt(((0.041044, 0.024913, 
    0.0), ), ((0.038087, 0.023689, 0.0), ), ), normal=
    mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt((0.041966, 0.024301, 
    0.00075), ), point=
    mdb.models['Model-1'].parts['UD_FIBERS'].InterestingPoint(
    mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt((0.035131, 0.026138, 
    0.00075), ), CENTER))
mdb.models['Model-1'].parts['UD_FIBERS'].PartitionCellByPlanePointNormal(cells=
    mdb.models['Model-1'].parts['UD_FIBERS'].cells.findAt(((0.037995, 0.038229, 
    0.0), ), ((0.035039, 0.037004, 0.0), ), ), normal=
    mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt((0.038917, 0.037617, 
    0.00075), ), point=
    mdb.models['Model-1'].parts['UD_FIBERS'].InterestingPoint(
    mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt((0.032082, 0.039453, 
    0.00075), ), CENTER))
mdb.models['Model-1'].parts['UD_FIBERS'].PartitionCellByPlanePointNormal(cells=
    mdb.models['Model-1'].parts['UD_FIBERS'].cells.findAt(((0.026262, 0.038761, 
    0.0), ), ((0.023306, 0.037537, 0.0), ), ), normal=
    mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt((0.027184, 0.038149, 
    0.00075), ), point=
    mdb.models['Model-1'].parts['UD_FIBERS'].InterestingPoint(
    mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt((0.020349, 0.039986, 
    0.00075), ), CENTER))
mdb.models['Model-1'].parts['UD_FIBERS'].PartitionCellByPlanePointNormal(cells=
    mdb.models['Model-1'].parts['UD_FIBERS'].cells.findAt(((0.014262, 0.038761, 
    0.0), ), ((0.011306, 0.037537, 0.0), ), ), normal=
    mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt((0.010384, 0.038149, 
    0.00075), ), point=
    mdb.models['Model-1'].parts['UD_FIBERS'].InterestingPoint(
    mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt((0.008349, 0.039986, 
    0.00075), ), CENTER))
mdb.models['Model-1'].parts['UD_FIBERS'].PartitionCellByPlanePointNormal(cells=
    mdb.models['Model-1'].parts['UD_FIBERS'].cells.findAt(((0.011578, 0.049423, 
    0.00075), )), normal=mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt(
    (0.012528, 0.05, 0.00075), ), point=
    mdb.models['Model-1'].parts['UD_FIBERS'].InterestingPoint(
    mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt((0.013408, 0.046823, 
    0.00075), ), CENTER))
mdb.models['Model-1'].parts['UD_FIBERS'].PartitionCellByPlanePointNormal(cells=
    mdb.models['Model-1'].parts['UD_FIBERS'].cells.findAt(((0.023202, 0.049541, 
    0.00075), )), normal=mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt(
    (0.024228, 0.05, 0.00075), ), point=
    mdb.models['Model-1'].parts['UD_FIBERS'].InterestingPoint(
    mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt((0.024228, 0.05, 
    0.00075), ), MIDDLE))
mdb.models['Model-1'].parts['UD_FIBERS'].PartitionCellByPlanePointNormal(cells=
    mdb.models['Model-1'].parts['UD_FIBERS'].cells.findAt(((0.035448, 0.049441, 
    0.00075), )), normal=mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt(
    (0.036332, 0.05, 0.00075), ), point=
    mdb.models['Model-1'].parts['UD_FIBERS'].InterestingPoint(
    mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt((0.037374, 0.04652, 
    0.00075), ), CENTER))
mdb.models['Model-1'].parts['UD_FIBERS'].PartitionCellByPlanePointNormal(cells=
    mdb.models['Model-1'].parts['UD_FIBERS'].cells.findAt(((0.049464, 0.037219, 
    0.00075), )), normal=mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt(
    (0.05, 0.033446, 0.00075), ), point=
    mdb.models['Model-1'].parts['UD_FIBERS'].InterestingPoint(
    mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt((0.047076, 0.032688, 
    0.00075), ), CENTER))
mdb.models['Model-1'].parts['UD_FIBERS'].PartitionCellByPlanePointNormal(cells=
    mdb.models['Model-1'].parts['UD_FIBERS'].cells.findAt(((0.049744, 0.025159, 
    0.00075), )), normal=mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt(
    (0.05, 0.023146, 0.00075), ), point=
    mdb.models['Model-1'].parts['UD_FIBERS'].InterestingPoint(
    mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt((0.049232, 0.023059, 
    0.00075), ), CENTER))

mdb.models['Model-1'].parts['UD_FIBERS'].setMeshControls(regions=
    mdb.models['Model-1'].parts['UD_FIBERS'].cells.findAt(((0.003554, 0.036431, 
    0.00075), ), ((0.003582, 0.035279, 0.00075), ), ((0.008829, 0.009372, 0.0), 
    ), ((0.007604, 0.015528, 0.0), ), ((0.008829, 0.015528, 0.0), ), ((
    0.007604, 0.009372, 0.0), ), ((0.013396, 0.035071, 0.0), ), ((0.012172, 
    0.041227, 0.0), ), ((0.013396, 0.041227, 0.0), ), ((0.012172, 0.035071, 
    0.0), ), ((0.02063, 0.010504, 0.0), ), ((0.019406, 0.01666, 0.0), ), ((
    0.02063, 0.01666, 0.0), ), ((0.019406, 0.010504, 0.0), ), ((0.025396, 
    0.035071, 0.0), ), ((0.024172, 0.041227, 0.0), ), ((0.025396, 0.041227, 
    0.0), ), ((0.024172, 0.035071, 0.0), ), ((0.032628, 0.009106, 0.0), ), ((
    0.031404, 0.015262, 0.0), ), ((0.032628, 0.015262, 0.0), ), ((0.031404, 
    0.009106, 0.0), ), ((0.037129, 0.034538, 0.0), ), ((0.035904, 0.040695, 
    0.0), ), ((0.037129, 0.040695, 0.0), ), ((0.035904, 0.034538, 0.0), ), ((
    0.044629, 0.009239, 0.0), ), ((0.043404, 0.015395, 0.0), ), ((0.044629, 
    0.015395, 0.0), ), ((0.043404, 0.009239, 0.0), ), ((0.009618, 0.003337, 
    0.00075), ), ((0.010649, 0.003337, 0.00075), ), ((0.034535, 0.002994, 
    0.00075), ), ((0.033331, 0.002994, 0.00075), ), ((0.047376, 0.036376, 
    0.00075), ), ((0.047376, 0.035262, 0.00075), ), ((0.04959, 0.003859, 
    0.00075), ), ((0.04935, 0.024988, 0.00075), ), ((0.04935, 0.02428, 
    0.00075), ), ((0.045883, 0.049502, 0.00075), ), ((0.021393, 0.004278, 
    0.00075), ), ((0.022612, 0.004278, 0.00075), ), ((0.000475, 0.003438, 
    0.00075), ), ((0.040178, 0.021223, 0.0), ), ((0.038953, 0.027379, 0.0), ), 
    ((0.040178, 0.027379, 0.0), ), ((0.038953, 0.021223, 0.0), ), ((0.03343, 
    0.046792, 0.00075), ), ((0.034436, 0.046792, 0.00075), ), ((0.028246, 
    0.022355, 0.0), ), ((0.027022, 0.028511, 0.0), ), ((0.028246, 0.028511, 
    0.0), ), ((0.027022, 0.022355, 0.0), ), ((0.021534, 0.048068, 0.00075), ), 
    ((0.02247, 0.048068, 0.00075), ), ((0.016314, 0.023153, 0.0), ), ((
    0.015089, 0.02931, 0.0), ), ((0.016314, 0.02931, 0.0), ), ((0.015089, 
    0.023153, 0.0), ), ((0.009547, 0.047133, 0.00075), ), ((0.01072, 0.047133, 
    0.00075), ), ((0.005589, 0.025246, 0.00075), ), ((0.005612, 0.02408, 
    0.00075), ), ((0.001663, 0.049551, 0.00075), ), ), technique=SWEEP)
	
mdb.models['Model-1'].parts['UD_FIBERS'].seedEdgeByNumber(constraint=FIXED, 
    edges=mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt(((0.0, 
    0.035819, 0.000187), ), ((0.005514, 0.035819, 0.000563), ), ((0.0, 
    0.040565, 0.000187), ), ((0.0, 0.031072, 0.000563), ), ((0.008217, 0.01245, 
    0.000187), ), ((0.008217, 0.01725, 0.000188), ), ((0.008217, 0.00765, 
    0.000563), ), ((0.003417, 0.01245, 0.000188), ), ((0.013017, 0.01245, 
    0.000563), ), ((0.012784, 0.038149, 0.000187), ), ((0.012784, 0.042949, 
    0.000188), ), ((0.012784, 0.033349, 0.000563), ), ((0.007984, 0.038149, 
    0.000188), ), ((0.017584, 0.038149, 0.000563), ), ((0.020018, 0.013582, 
    0.000187), ), ((0.020018, 0.018382, 0.000188), ), ((0.020018, 0.008782, 
    0.000563), ), ((0.015218, 0.013582, 0.000188), ), ((0.024818, 0.013582, 
    0.000563), ), ((0.024784, 0.038149, 0.000187), ), ((0.024784, 0.042949, 
    0.000563), ), ((0.024784, 0.033349, 0.000188), ), ((0.019984, 0.038149, 
    0.000563), ), ((0.029584, 0.038149, 0.000188), ), ((0.032016, 0.012184, 
    0.000187), ), ((0.032016, 0.016984, 0.000563), ), ((0.032016, 0.007384, 
    0.000188), ), ((0.027216, 0.012184, 0.000563), ), ((0.036816, 0.012184, 
    0.000188), ), ((0.036517, 0.037617, 0.000187), ), ((0.036517, 0.042417, 
    0.000563), ), ((0.036517, 0.032817, 0.000188), ), ((0.031717, 0.037617, 
    0.000563), ), ((0.041317, 0.037617, 0.000188), ), ((0.044016, 0.012317, 
    0.000187), ), ((0.044016, 0.017117, 0.000563), ), ((0.044016, 0.007517, 
    0.000188), ), ((0.039216, 0.012317, 0.000563), ), ((0.048816, 0.012317, 
    0.000188), ), ((0.010134, 0.0, 0.000187), ), ((0.010134, 0.005133, 
    0.000563), ), ((0.005345, 0.0, 0.000187), ), ((0.014922, 0.0, 0.000563), ), 
    ((0.033933, 0.0, 0.000187), ), ((0.033933, 0.004667, 0.000563), ), ((
    0.029135, 0.0, 0.000187), ), ((0.038731, 0.0, 0.000563), ), ((0.05, 
    0.035819, 0.000187), ), ((0.045914, 0.035819, 0.000563), ), ((0.05, 
    0.031072, 0.000187), ), ((0.05, 0.040565, 0.000563), ), ((0.05, 0.005621, 
    0.000563), ), ((0.05, 0.0, 0.000188), ), ((0.043487, 0.0, 0.000187), ), ((
    0.05, 0.024634, 0.000187), ), ((0.048966, 0.024634, 0.000563), ), ((0.05, 
    0.021658, 0.000187), ), ((0.05, 0.02761, 0.000563), ), ((0.043487, 0.05, 
    0.000563), ), ((0.05, 0.05, 0.000188), ), ((0.05, 0.046776, 0.000187), ), (
    (0.022002, 0.0, 0.000187), ), ((0.022002, 0.006598, 0.000563), ), ((
    0.017552, 0.0, 0.000187), ), ((0.026453, 0.0, 0.000563), ), ((0.0, 
    0.005621, 0.000187), ), ((0.0, 0.0, 0.000188), ), ((0.002783, 0.0, 
    0.000563), ), ((0.039566, 0.024301, 0.000187), ), ((0.039566, 0.029101, 
    0.000563), ), ((0.039566, 0.019501, 0.000188), ), ((0.034766, 0.024301, 
    0.000563), ), ((0.044366, 0.024301, 0.000188), ), ((0.033933, 0.05, 
    0.000187), ), ((0.033933, 0.045067, 0.000563), ), ((0.038731, 0.05, 
    0.000187), ), ((0.029135, 0.05, 0.000563), ), ((0.027634, 0.025433, 
    0.000187), ), ((0.027634, 0.030233, 0.000563), ), ((0.027634, 0.020633, 
    0.000188), ), ((0.022834, 0.025433, 0.000563), ), ((0.032434, 0.025433, 
    0.000188), ), ((0.022002, 0.05, 0.000187), ), ((0.022002, 0.046998, 
    0.000563), ), ((0.026453, 0.05, 0.000187), ), ((0.017552, 0.05, 0.000563), 
    ), ((0.015702, 0.026232, 0.000187), ), ((0.015702, 0.031032, 0.000563), ), 
    ((0.015702, 0.021432, 0.000188), ), ((0.010902, 0.026232, 0.000563), ), ((
    0.020502, 0.026232, 0.000188), ), ((0.010134, 0.05, 0.000187), ), ((
    0.010134, 0.045533, 0.000563), ), ((0.014922, 0.05, 0.000187), ), ((
    0.005345, 0.05, 0.000563), ), ((0.0, 0.024634, 0.000187), ), ((0.008566, 
    0.024634, 0.000563), ), ((0.0, 0.02761, 0.000187), ), ((0.0, 0.021658, 
    0.000563), ), ((0.002783, 0.05, 0.000187), ), ((0.0, 0.05, 0.000188), ), ((
    0.0, 0.046776, 0.000563), ), ), number=matrix_depth_el)

mdb.models['Model-1'].parts['UD_FIBERS'].seedEdgeBySize(constraint=FIXED, 
    deviationFactor=0.1, edges=
    mdb.models['Model-1'].parts['UD_FIBERS'].edges.findAt(((0.0, 0.032259, 
    0.0), ), ((0.0, 0.039379, 0.00075), ), ((0.0, 0.037006, 0.0), ), ((0.0, 
    0.034632, 0.00075), ), ((0.013725, 0.0, 0.0), ), ((0.006542, 0.0, 0.00075), 
    ), ((0.008937, 0.0, 0.0), ), ((0.011331, 0.0, 0.00075), ), ((0.035133, 0.0, 
    0.00075), ), ((0.032734, 0.0, 0.0), ), ((0.037532, 0.0, 0.0), ), ((
    0.030334, 0.0, 0.00075), ), ((0.05, 0.037006, 0.00075), ), ((0.05, 
    0.034632, 0.0), ), ((0.05, 0.039379, 0.0), ), ((0.05, 0.032259, 0.00075), 
    ), ((0.05, 0.001405, 0.00075), ), ((0.05, 0.004216, 0.0), ), ((0.048372, 
    0.0, 0.0), ), ((0.045115, 0.0, 0.00075), ), ((0.05, 0.025378, 0.00075), ), 
    ((0.05, 0.02389, 0.0), ), ((0.05, 0.026866, 0.0), ), ((0.05, 0.022402, 
    0.00075), ), ((0.048372, 0.05, 0.00075), ), ((0.045115, 0.05, 0.0), ), ((
    0.05, 0.049194, 0.0), ), ((0.05, 0.047582, 0.00075), ), ((0.02534, 0.0, 
    0.0), ), ((0.018664, 0.0, 0.00075), ), ((0.02089, 0.0, 0.0), ), ((0.023115, 
    0.0, 0.00075), ), ((0.0, 0.001405, 0.0), ), ((0.0, 0.004216, 0.00075), ), (
    (0.000696, 0.0, 0.00075), ), ((0.002087, 0.0, 0.0), ), ((0.032734, 0.05, 
    0.00075), ), ((0.035133, 0.05, 0.0), ), ((0.030334, 0.05, 0.0), ), ((
    0.037532, 0.05, 0.00075), ), ((0.02089, 0.05, 0.00075), ), ((0.023115, 
    0.05, 0.0), ), ((0.018664, 0.05, 0.0), ), ((0.02534, 0.05, 0.00075), ), ((
    0.008937, 0.05, 0.00075), ), ((0.011331, 0.05, 0.0), ), ((0.006542, 0.05, 
    0.0), ), ((0.013725, 0.05, 0.00075), ), ((0.0, 0.022402, 0.0), ), ((0.0, 
    0.026866, 0.00075), ), ((0.0, 0.025378, 0.0), ), ((0.0, 0.02389, 0.00075), 
    ), ((0.000696, 0.05, 0.0), ), ((0.002087, 0.05, 0.00075), ), ((0.0, 
    0.049194, 0.00075), ), ((0.0, 0.047582, 0.0), ), ), minSizeFactor=0.1, 
    size=fiber_mesh_size)
	
MODEL.parts['UD_FIBERS'].generateMesh()

r.regenerate()
MODEL.parts['MATRIX_WITH_HOLES'].seedPart(deviationFactor=0.1, 
    minSizeFactor=0.1, size=mesh_size)

# mdb.models['Model-1'].parts['MATRIX_WITH_HOLES'].generateMesh()
	
# TODO: for new models, This set must be generated    
r.regenerate()
mdb.models['Model-1'].parts['MATRIX_WITH_HOLES'].Set(faces=
    mdb.models['Model-1'].parts['MATRIX_WITH_HOLES'].faces.findAt(((0.044898, 
    0.031857, 0.00075), )), name='Set-1')

aa3 = MODEL.parts['MATRIX_WITH_HOLES'].sets['Set-1'].faces[0].pointOn     

MODEL.parts['MATRIX_WITH_HOLES'].setElementType(elemTypes=(
    ElemType(elemCode=C3D8T, elemLibrary=STANDARD, secondOrderAccuracy=OFF, 
    distortionControl=DEFAULT), ElemType(elemCode=C3D6T, elemLibrary=STANDARD), 
    ElemType(elemCode=C3D4T, elemLibrary=STANDARD)), regions=(
    MODEL.parts['MATRIX_WITH_HOLES'].cells.findAt(aa3, ), ))
	
mdb.models['Model-1'].parts['MATRIX_WITH_HOLES'].seedEdgeByNumber(constraint=
    FIXED, edges=mdb.models['Model-1'].parts['MATRIX_WITH_HOLES'].edges.findAt(
    ((0.0, 0.040565, 0.000187), ), ((0.0, 0.031072, 0.000563), ), ((0.014922, 
    0.0, 0.000563), ), ((0.005345, 0.0, 0.000187), ), ((0.029135, 0.0, 
    0.000187), ), ((0.038731, 0.0, 0.000563), ), ((0.05, 0.040565, 0.000563), 
    ), ((0.05, 0.031072, 0.000187), ), ((0.05, 0.005621, 0.000563), ), ((
    0.043487, 0.0, 0.000187), ), ((0.05, 0.02761, 0.000563), ), ((0.05, 
    0.021658, 0.000187), ), ((0.043487, 0.05, 0.000563), ), ((0.05, 0.046776, 
    0.000187), ), ((0.026453, 0.0, 0.000563), ), ((0.017552, 0.0, 0.000187), ), 
    ((0.0, 0.005621, 0.000187), ), ((0.002783, 0.0, 0.000563), ), ((0.029135, 
    0.05, 0.000563), ), ((0.038731, 0.05, 0.000187), ), ((0.026453, 0.05, 
    0.000187), ), ((0.017552, 0.05, 0.000563), ), ((0.014922, 0.05, 0.000187), 
    ), ((0.005345, 0.05, 0.000563), ), ((0.0, 0.02761, 0.000187), ), ((0.0, 
    0.021658, 0.000563), ), ((0.0, 0.046776, 0.000563), ), ((0.002783, 0.05, 
    0.000187), ), ), number=matrix_depth_el)
	
mdb.models['Model-1'].parts['MATRIX_WITH_HOLES'].seedEdgeBySize(constraint=
    FIXED, deviationFactor=0.1, edges=
    mdb.models['Model-1'].parts['MATRIX_WITH_HOLES'].edges.findAt(((0.05, 
    0.017649, 0.0), ), ((0.05, 0.017649, 0.00075), ), ((0.05, 0.030207, 
    0.00075), ), ((0.05, 0.030207, 0.0), ), ((0.05, 0.045223, 0.0), ), ((0.05, 
    0.045223, 0.00075), ), ((0.042298, 0.0, 0.0), ), ((0.042298, 0.0, 0.00075), 
    ), ((0.028464, 0.0, 0.00075), ), ((0.028464, 0.0, 0.0), ), ((0.004705, 0.0, 
    0.00075), ), ((0.004705, 0.0, 0.0), ), ((0.0, 0.00963, 0.00075), ), ((0.0, 
    0.00963, 0.0), ), ((0.0, 0.028475, 0.0), ), ((0.0, 0.028475, 0.00075), ), (
    (0.003423, 0.05, 0.0), ), ((0.003423, 0.05, 0.00075), ), ((0.01558, 0.05, 
    0.0), ), ((0.01558, 0.05, 0.00075), ), ((0.027123, 0.05, 0.0), ), ((
    0.027123, 0.05, 0.00075), ), ((0.03992, 0.05, 0.0), ), ((0.03992, 0.05, 
    0.00075), ), ((0.016894, 0.0, 0.00075), ), ((0.0, 0.042118, 0.00075), ), ((
    0.0, 0.042118, 0.0), ), ((0.016894, 0.0, 0.0), ), ), size=mesh_size)
	
mdb.models['Model-1'].parts['MATRIX_WITH_HOLES'].generateMesh()
r.regenerate()   

# Step
MODEL.CoupledTempDisplacementStep(adaptiveDampingRatio=None, 
    amplitude=RAMP, cetol=None, continueDampingFactors=False, creepIntegration=
    None, deltmx=None, initialInc=INCR_INIT, maxInc=INCR_MAX, maxNumInc=INCR_MAX_NUM, minInc=
    INCR_MIN, name='LOADING_STEP', previous='Initial', response=STEADY_STATE, 
    stabilizationMagnitude=STAB_MAGN, stabilizationMethod=DAMPING_FACTOR)

r.regenerate()    

## Interactions & Contact
mdb.models['Model-1'].ContactProperty('IntProp-1')
mdb.models['Model-1'].interactionProperties['IntProp-1'].TangentialBehavior(
    dependencies=0, directionality=ISOTROPIC, elasticSlipStiffness=None, 
    formulation=PENALTY, fraction=0.005, maximumElasticSlip=FRACTION, 
    pressureDependency=OFF, shearStressLimit=None, slipRateDependency=OFF, 
    table=((0.3, ), ), temperatureDependency=OFF)
mdb.models['Model-1'].interactionProperties['IntProp-1'].NormalBehavior(
    allowSeparation=ON, clearanceAtZeroContactPressure=0.0, 
    constraintEnforcementMethod=PENALTY, contactStiffness=DEFAULT, 
    contactStiffnessScaleFactor=1.0, pressureOverclosure=HARD, 
    stiffnessBehavior=LINEAR)
mdb.models['Model-1'].interactionProperties['IntProp-1'].Damage(evolTable=((
    0.2, 0.6, 0.6), ), evolutionType=ENERGY, exponent=1.0, initTable=((80.0, 
    40.0, 40.0), ), mixedModeType=BK, softening=EXPONENTIAL, useEvolution=ON, 
    useMixedMode=ON, useStabilization=ON, viscosityCoef=0.001)

	

## define contact b/w matrix and fibers
mdb.models['Model-1'].rootAssembly.Surface(name='MASTER_SURFACE', side1Faces=
    mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].faces.findAt(((
    0.005483, 0.035272, 0.0005), ), ((0.005474, 0.036442, 0.00025), ), ((
    0.007594, 0.007691, 0.0005), ), ((0.007594, 0.017209, 0.00025), ), ((
    0.008839, 0.017209, 0.0005), ), ((0.008839, 0.007691, 0.00025), ), ((
    0.012161, 0.03339, 0.0005), ), ((0.012161, 0.042909, 0.00025), ), ((
    0.013407, 0.042909, 0.0005), ), ((0.013407, 0.03339, 0.00025), ), ((
    0.019395, 0.008822, 0.0005), ), ((0.019395, 0.018341, 0.00025), ), ((
    0.020641, 0.018341, 0.0005), ), ((0.020641, 0.008822, 0.00025), ), ((
    0.024161, 0.03339, 0.0005), ), ((0.024161, 0.042909, 0.00025), ), ((
    0.025407, 0.03339, 0.00025), ), ((0.025407, 0.042909, 0.0005), ), ((
    0.031393, 0.007424, 0.0005), ), ((0.031393, 0.016943, 0.00025), ), ((
    0.032639, 0.007424, 0.00025), ), ((0.032639, 0.016943, 0.0005), ), ((
    0.035894, 0.032857, 0.0005), ), ((0.035894, 0.042376, 0.00025), ), ((
    0.03714, 0.032857, 0.00025), ), ((0.03714, 0.042376, 0.0005), ), ((
    0.043393, 0.007558, 0.0005), ), ((0.043393, 0.017076, 0.00025), ), ((
    0.044639, 0.007558, 0.00025), ), ((0.044639, 0.017076, 0.0005), ), ((
    0.010655, 0.005104, 0.0005), ), ((0.009612, 0.005104, 0.00025), ), ((
    0.033321, 0.004628, 0.00025), ), ((0.034545, 0.004628, 0.0005), ), ((
    0.045948, 0.035254, 0.00025), ), ((0.045948, 0.036383, 0.0005), ), ((
    0.049802, 0.005699, 0.00025), ), ((0.048979, 0.024278, 0.0005), ), ((
    0.048979, 0.024989, 0.00025), ), ((0.043653, 0.049479, 0.00025), ), ((
    0.022622, 0.006557, 0.0005), ), ((0.021382, 0.006557, 0.00025), ), ((
    0.000508, 0.005371, 0.0005), ), ((0.038943, 0.019542, 0.0005), ), ((
    0.038943, 0.02906, 0.00025), ), ((0.040189, 0.019542, 0.00025), ), ((
    0.040189, 0.02906, 0.0005), ), ((0.034442, 0.045094, 0.00025), ), ((
    0.033424, 0.045094, 0.0005), ), ((0.027011, 0.020673, 0.0005), ), ((
    0.027011, 0.030192, 0.00025), ), ((0.028257, 0.020673, 0.00025), ), ((
    0.028257, 0.030192, 0.0005), ), ((0.022475, 0.047021, 0.00025), ), ((
    0.02153, 0.047021, 0.0005), ), ((0.015079, 0.021472, 0.0005), ), ((
    0.015079, 0.030991, 0.00025), ), ((0.016324, 0.021472, 0.00025), ), ((
    0.016324, 0.030991, 0.0005), ), ((0.01073, 0.04557, 0.00025), ), ((
    0.009538, 0.04557, 0.0005), ), ((0.008533, 0.024073, 0.0005), ), ((
    0.008525, 0.025257, 0.00025), ), ((0.002637, 0.049534, 0.0005), ), ))
mdb.models['Model-1'].rootAssembly.Surface(name='SLAVE_SURFACE', side1Faces=
    mdb.models['Model-1'].rootAssembly.instances['MATRIX_WITH_HOLES-1'].faces.findAt(
    ((0.000606, 0.040618, 0.00025), ), ((0.007594, 0.017209, 0.0005), ), ((
    0.012161, 0.042909, 0.0005), ), ((0.019395, 0.018341, 0.0005), ), ((
    0.024161, 0.042909, 0.0005), ), ((0.031393, 0.016943, 0.0005), ), ((
    0.035894, 0.042376, 0.0005), ), ((0.043393, 0.017076, 0.0005), ), ((
    0.014927, 0.00058, 0.0005), ), ((0.029191, 0.000611, 0.00025), ), ((
    0.049447, 0.040448, 0.0005), ), ((0.049444, 0.005816, 0.0005), ), ((
    0.049694, 0.027174, 0.0005), ), ((0.043653, 0.049479, 0.0005), ), ((
    0.026648, 0.00059, 0.0005), ), ((0.000508, 0.005371, 0.00025), ), ((
    0.038943, 0.02906, 0.0005), ), ((0.029153, 0.049435, 0.00025), ), ((
    0.027011, 0.030192, 0.0005), ), ((0.026223, 0.049511, 0.00025), ), ((
    0.015079, 0.030991, 0.0005), ), ((0.014844, 0.049408, 0.00025), ), ((
    0.000404, 0.02806, 0.00025), ), ((0.00044, 0.046988, 0.0005), ), ))
	
mdb.models['Model-1'].SurfaceToSurfaceContactStd(adjustMethod=TOLERANCE, 
    adjustTolerance=1e-08, clearanceRegion=None, createStepName='LOADING_STEP', 
    datumAxis=None, initialClearance=1e-08, interactionProperty='IntProp-1', 
    master=mdb.models['Model-1'].rootAssembly.surfaces['MASTER_SURFACE'], name=
    'Int-1', slave=mdb.models['Model-1'].rootAssembly.surfaces['SLAVE_SURFACE']
    , sliding=SMALL, supplementaryContact=ALWAYS, surfaceSmoothing=AUTOMATIC, 
    thickness=ON, tied=OFF)
   

## Define X_POS_MATRIX, X_POS_FIBERS etc.
mdb.models['Model-1'].rootAssembly.Set(faces=
    mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].faces.findAt(((
    0.05, 0.037401, 0.00025), ), ((0.05, 0.034237, 0.0005), ), ((0.05, 
    0.003747, 0.0005), ), ((0.05, 0.025626, 0.00025), ), ((0.05, 0.023642, 
    0.0005), ), ((0.05, 0.04785, 0.00025), ), ), name='X_POS_FIBERS')
mdb.models['Model-1'].rootAssembly.Set(faces=
    mdb.models['Model-1'].rootAssembly.instances['MATRIX_WITH_HOLES-1'].faces.findAt(
    ((0.05, 0.010967, 0.00025), ), ((0.05, 0.029918, 0.0005), ), ((0.05, 
    0.042636, 0.00025), ), ), name='X_POS_MATRIX')
mdb.models['Model-1'].rootAssembly.Set(faces=
    mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].faces.findAt(((
    0.0, 0.037401, 0.0005), ), ((0.0, 0.034237, 0.00025), ), ((0.0, 0.003747, 
    0.00025), ), ((0.0, 0.025626, 0.0005), ), ((0.0, 0.023642, 0.00025), ), ((
    0.0, 0.04785, 0.0005), ), ), name='X_NEG_FIBERS')
mdb.models['Model-1'].rootAssembly.Set(faces=
    mdb.models['Model-1'].rootAssembly.instances['MATRIX_WITH_HOLES-1'].faces.findAt(
    ((0.0, 0.010967, 0.0005), ), ((0.0, 0.029918, 0.00025), ), ((0.0, 0.042636, 
    0.0005), ), ), name='X_NEG_MATRIX')
mdb.models['Model-1'].rootAssembly.Set(faces=
    mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].faces.findAt(((
    0.045658, 0.05, 0.0005), ), ((0.032334, 0.05, 0.00025), ), ((0.035532, 
    0.05, 0.0005), ), ((0.020519, 0.05, 0.00025), ), ((0.023486, 0.05, 0.0005), 
    ), ((0.008538, 0.05, 0.00025), ), ((0.01173, 0.05, 0.0005), ), ((0.001855, 
    0.05, 0.00025), ), ), name='Y_POS_FIBERS')
mdb.models['Model-1'].rootAssembly.Set(faces=
    mdb.models['Model-1'].rootAssembly.instances['MATRIX_WITH_HOLES-1'].faces.findAt(
    ((0.004491, 0.05, 0.00025), ), ((0.016675, 0.05, 0.00025), ), ((0.028241, 
    0.05, 0.00025), ), ((0.041902, 0.05, 0.00025), ), ), name='Y_POS_MATRIX')
mdb.models['Model-1'].rootAssembly.Set(faces=
    mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].faces.findAt(((
    0.008538, 0.0, 0.0005), ), ((0.01173, 0.0, 0.00025), ), ((0.035532, 0.0, 
    0.00025), ), ((0.032334, 0.0, 0.0005), ), ((0.045658, 0.0, 0.00025), ), ((
    0.020519, 0.0, 0.0005), ), ((0.023486, 0.0, 0.00025), ), ((0.001855, 0.0, 
    0.0005), ), ), name='Y_NEG_FIBERS')
mdb.models['Model-1'].rootAssembly.Set(faces=
    mdb.models['Model-1'].rootAssembly.instances['MATRIX_WITH_HOLES-1'].faces.findAt(
    ((0.040316, 0.0, 0.00025), ), ((0.028241, 0.0, 0.0005), ), ((0.004491, 0.0, 
    0.0005), ), ((0.015799, 0.0, 0.00025), ), ), name='Y_NEG_MATRIX')
mdb.models['Model-1'].rootAssembly.Set(faces=
    mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].faces.findAt(((
    0.003554, 0.036431, 0.00075), ), ((0.003582, 0.035279, 0.00075), ), ((
    0.008829, 0.015528, 0.00075), ), ((0.007604, 0.009372, 0.00075), ), ((
    0.008829, 0.009372, 0.00075), ), ((0.007604, 0.015528, 0.00075), ), ((
    0.013396, 0.041227, 0.00075), ), ((0.012172, 0.035071, 0.00075), ), ((
    0.013396, 0.035071, 0.00075), ), ((0.012172, 0.041227, 0.00075), ), ((
    0.02063, 0.01666, 0.00075), ), ((0.019406, 0.010504, 0.00075), ), ((
    0.02063, 0.010504, 0.00075), ), ((0.019406, 0.01666, 0.00075), ), ((
    0.025396, 0.041227, 0.00075), ), ((0.024172, 0.035071, 0.00075), ), ((
    0.025396, 0.035071, 0.00075), ), ((0.024172, 0.041227, 0.00075), ), ((
    0.032628, 0.015262, 0.00075), ), ((0.031404, 0.009106, 0.00075), ), ((
    0.032628, 0.009106, 0.00075), ), ((0.031404, 0.015262, 0.00075), ), ((
    0.037129, 0.040695, 0.00075), ), ((0.035904, 0.034538, 0.00075), ), ((
    0.037129, 0.034538, 0.00075), ), ((0.035904, 0.040695, 0.00075), ), ((
    0.044629, 0.015395, 0.00075), ), ((0.043404, 0.009239, 0.00075), ), ((
    0.044629, 0.009239, 0.00075), ), ((0.043404, 0.015395, 0.00075), ), ((
    0.009618, 0.003337, 0.00075), ), ((0.010649, 0.003337, 0.00075), ), ((
    0.034535, 0.002994, 0.00075), ), ((0.033331, 0.002994, 0.00075), ), ((
    0.047376, 0.036376, 0.00075), ), ((0.047376, 0.035262, 0.00075), ), ((
    0.04959, 0.003859, 0.00075), ), ((0.04935, 0.024988, 0.00075), ), ((
    0.04935, 0.02428, 0.00075), ), ((0.045883, 0.049502, 0.00075), ), ((
    0.021393, 0.004278, 0.00075), ), ((0.022612, 0.004278, 0.00075), ), ((
    0.000475, 0.003438, 0.00075), ), ((0.040178, 0.027379, 0.00075), ), ((
    0.038953, 0.021223, 0.00075), ), ((0.040178, 0.021223, 0.00075), ), ((
    0.038953, 0.027379, 0.00075), ), ((0.03343, 0.046792, 0.00075), ), ((
    0.034436, 0.046792, 0.00075), ), ((0.028246, 0.028511, 0.00075), ), ((
    0.027022, 0.022355, 0.00075), ), ((0.028246, 0.022355, 0.00075), ), ((
    0.027022, 0.028511, 0.00075), ), ((0.021534, 0.048068, 0.00075), ), ((
    0.02247, 0.048068, 0.00075), ), ((0.016314, 0.02931, 0.00075), ), ((
    0.015089, 0.023153, 0.00075), ), ((0.016314, 0.023153, 0.00075), ), ((
    0.015089, 0.02931, 0.00075), ), ((0.009547, 0.047133, 0.00075), ), ((
    0.01072, 0.047133, 0.00075), ), ((0.005589, 0.025246, 0.00075), ), ((
    0.005612, 0.02408, 0.00075), ), ((0.001663, 0.049551, 0.00075), ), ), name=
    'Z_POS_FIBERS')
mdb.models['Model-1'].rootAssembly.Set(faces=
    mdb.models['Model-1'].rootAssembly.instances['MATRIX_WITH_HOLES-1'].faces.findAt(
    ((0.044898, 0.031857, 0.00075), )), name='Z_POS_MATRIX')
mdb.models['Model-1'].rootAssembly.Set(faces=
    mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].faces.findAt(((
    0.003582, 0.035279, 0.0), ), ((0.003554, 0.036431, 0.0), ), ((0.008829, 
    0.009372, 0.0), ), ((0.007604, 0.015528, 0.0), ), ((0.008829, 0.015528, 
    0.0), ), ((0.007604, 0.009372, 0.0), ), ((0.013396, 0.035071, 0.0), ), ((
    0.012172, 0.041227, 0.0), ), ((0.013396, 0.041227, 0.0), ), ((0.012172, 
    0.035071, 0.0), ), ((0.02063, 0.010504, 0.0), ), ((0.019406, 0.01666, 0.0), 
    ), ((0.02063, 0.01666, 0.0), ), ((0.019406, 0.010504, 0.0), ), ((0.025396, 
    0.035071, 0.0), ), ((0.024172, 0.041227, 0.0), ), ((0.025396, 0.041227, 
    0.0), ), ((0.024172, 0.035071, 0.0), ), ((0.032628, 0.009106, 0.0), ), ((
    0.031404, 0.015262, 0.0), ), ((0.032628, 0.015262, 0.0), ), ((0.031404, 
    0.009106, 0.0), ), ((0.037129, 0.034538, 0.0), ), ((0.035904, 0.040695, 
    0.0), ), ((0.037129, 0.040695, 0.0), ), ((0.035904, 0.034538, 0.0), ), ((
    0.044629, 0.009239, 0.0), ), ((0.043404, 0.015395, 0.0), ), ((0.044629, 
    0.015395, 0.0), ), ((0.043404, 0.009239, 0.0), ), ((0.010649, 0.003337, 
    0.0), ), ((0.009618, 0.003337, 0.0), ), ((0.033331, 0.002994, 0.0), ), ((
    0.034535, 0.002994, 0.0), ), ((0.047376, 0.035262, 0.0), ), ((0.047376, 
    0.036376, 0.0), ), ((0.049799, 0.003817, 0.0), ), ((0.04935, 0.02428, 0.0), 
    ), ((0.04935, 0.024988, 0.0), ), ((0.045883, 0.049502, 0.0), ), ((0.022612, 
    0.004278, 0.0), ), ((0.021393, 0.004278, 0.0), ), ((0.000475, 0.003438, 
    0.0), ), ((0.040178, 0.021223, 0.0), ), ((0.038953, 0.027379, 0.0), ), ((
    0.040178, 0.027379, 0.0), ), ((0.038953, 0.021223, 0.0), ), ((0.034436, 
    0.046792, 0.0), ), ((0.03343, 0.046792, 0.0), ), ((0.028246, 0.022355, 
    0.0), ), ((0.027022, 0.028511, 0.0), ), ((0.028246, 0.028511, 0.0), ), ((
    0.027022, 0.022355, 0.0), ), ((0.02247, 0.048068, 0.0), ), ((0.021534, 
    0.048068, 0.0), ), ((0.016314, 0.023153, 0.0), ), ((0.015089, 0.02931, 
    0.0), ), ((0.016314, 0.02931, 0.0), ), ((0.015089, 0.023153, 0.0), ), ((
    0.01072, 0.047133, 0.0), ), ((0.009547, 0.047133, 0.0), ), ((0.005612, 
    0.02408, 0.0), ), ((0.005589, 0.025246, 0.0), ), ((0.001663, 0.049551, 
    0.0), ), ), name='Z_NEG_FIBERS')
mdb.models['Model-1'].rootAssembly.Set(faces=
    mdb.models['Model-1'].rootAssembly.instances['MATRIX_WITH_HOLES-1'].faces.findAt(
    ((0.000458, 0.029159, 0.0), )), name='Z_NEG_MATRIX')
	
	
## Construct Sets for P.B.C
# |y
# 2 3
# 1 4-->x
CORNERS = ['CORNER_1', 'CORNER_2', 'CORNER_3', 'CORNER_4', 'CORNER_5', 'CORNER_6',\
'CORNER_7', 'CORNER_8']
r.Set(name=CORNERS[0], vertices=
    r.instances['UD_FIBERS-1'].vertices.findAt(
    ((0.0, 0.0, 0.0), )))
	
r.Set(name=CORNERS[1], vertices=
    r.instances['UD_FIBERS-1'].vertices.findAt(
    ((0, pt_2[1], 0.0), )))	
	
r.Set(name=CORNERS[2], vertices=
    r.instances['UD_FIBERS-1'].vertices.findAt(
    ((pt_2[0], pt_2[1], 0), )))		
	
r.Set(name=CORNERS[3], vertices=
    r.instances['UD_FIBERS-1'].vertices.findAt(
    ((pt_2[0], 0, 0), )))		

# |y
# 6 7
# 5 8-->x
r.Set(name=CORNERS[4], vertices=
    r.instances['UD_FIBERS-1'].vertices.findAt(
    ((0.0, 0.0, matrix_part_depth), )))
	
r.Set(name=CORNERS[5], vertices=
    r.instances['UD_FIBERS-1'].vertices.findAt(
    ((0, pt_2[1], matrix_part_depth), )))	
	
r.Set(name=CORNERS[6], vertices=
    r.instances['UD_FIBERS-1'].vertices.findAt(
    ((pt_2[0], pt_2[1], matrix_part_depth), )))		
	
r.Set(name=CORNERS[7], vertices=
    r.instances['UD_FIBERS-1'].vertices.findAt(
    ((pt_2[0], 0, matrix_part_depth), )))	

# # define edge sets
EDGES = ['EDGE_1', 'EDGE_2', 'EDGE_3', 'EDGE_4', 'EDGE_5', 'EDGE_6',\
'EDGE_7', 'EDGE_8', 'EDGE_9', 'EDGE_10', 'EDGE_11', 'EDGE_12']

mdb.models['Model-1'].rootAssembly.Set(edges=
    mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].edges.findAt(((
    0.0, 0.039379, 0.00075), ), ((0.0, 0.034632, 0.00075), ), ((0.0, 0.004216, 
    0.00075), ), ((0.0, 0.026866, 0.00075), ), ((0.0, 0.02389, 0.00075), ), ((
    0.0, 0.049194, 0.00075), ), ), name='edge_1_fibers')
mdb.models['Model-1'].rootAssembly.Set(edges=
    mdb.models['Model-1'].rootAssembly.instances['MATRIX_WITH_HOLES-1'].edges.findAt(
    ((0.0, 0.00963, 0.00075), ), ((0.0, 0.028475, 0.00075), ), ((0.0, 0.042118, 
    0.00075), ), ), name='edge_1_matrix')
mdb.models['Model-1'].rootAssembly.Set(edges=
    mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].edges.findAt(((
    0.05, 0.037006, 0.00075), ), ((0.05, 0.032259, 0.00075), ), ((0.05, 
    0.001405, 0.00075), ), ((0.05, 0.025378, 0.00075), ), ((0.05, 0.022402, 
    0.00075), ), ((0.05, 0.047582, 0.00075), ), ), name='edge_2_fibers')
mdb.models['Model-1'].rootAssembly.Set(edges=
    mdb.models['Model-1'].rootAssembly.instances['MATRIX_WITH_HOLES-1'].edges.findAt(
    ((0.05, 0.017649, 0.00075), ), ((0.05, 0.030207, 0.00075), ), ((0.05, 
    0.045223, 0.00075), ), ), name='edge_2_matrix')
mdb.models['Model-1'].rootAssembly.Set(edges=
    mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].edges.findAt(((
    0.0, 0.032259, 0.0), ), ((0.0, 0.037006, 0.0), ), ((0.0, 0.001405, 0.0), ), 
    ((0.0, 0.022402, 0.0), ), ((0.0, 0.025378, 0.0), ), ((0.0, 0.047582, 0.0), 
    ), ), name='edge_3_fibers')
mdb.models['Model-1'].rootAssembly.Set(edges=
    mdb.models['Model-1'].rootAssembly.instances['MATRIX_WITH_HOLES-1'].edges.findAt(
    ((0.0, 0.00963, 0.0), ), ((0.0, 0.028475, 0.0), ), ((0.0, 0.042118, 0.0), 
    ), ), name='edge_3_matrix')
mdb.models['Model-1'].rootAssembly.Set(edges=
    mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].edges.findAt(((
    0.05, 0.034632, 0.0), ), ((0.05, 0.039379, 0.0), ), ((0.05, 0.004216, 0.0), 
    ), ((0.05, 0.02389, 0.0), ), ((0.05, 0.026866, 0.0), ), ((0.05, 0.049194, 
    0.0), ), ), name='edge_4_fibers')
mdb.models['Model-1'].rootAssembly.Set(edges=
    mdb.models['Model-1'].rootAssembly.instances['MATRIX_WITH_HOLES-1'].edges.findAt(
    ((0.05, 0.017649, 0.0), ), ((0.05, 0.030207, 0.0), ), ((0.05, 0.045223, 
    0.0), ), ), name='edge_4_matrix')
mdb.models['Model-1'].rootAssembly.Set(edges=
    mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].edges.findAt(((
    0.0, 0.0, 0.000188), )), name='edge_5_fibers')
mdb.models['Model-1'].rootAssembly.Set(edges=
    mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].edges.findAt(((
    0.05, 0.0, 0.000188), )), name='edge_6_fibers')
mdb.models['Model-1'].rootAssembly.Set(edges=
    mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].edges.findAt(((
    0.013725, 0.0, 0.0), ), ((0.008937, 0.0, 0.0), ), ((0.032734, 0.0, 0.0), ), 
    ((0.037532, 0.0, 0.0), ), ((0.048372, 0.0, 0.0), ), ((0.02534, 0.0, 0.0), 
    ), ((0.02089, 0.0, 0.0), ), ((0.002087, 0.0, 0.0), ), ), name=
    'edge_7_fibers')
mdb.models['Model-1'].rootAssembly.Set(edges=
    mdb.models['Model-1'].rootAssembly.instances['MATRIX_WITH_HOLES-1'].edges.findAt(
    ((0.042298, 0.0, 0.0), ), ((0.028464, 0.0, 0.0), ), ((0.004705, 0.0, 0.0), 
    ), ((0.016894, 0.0, 0.0), ), ), name='edge_7_matrix')
mdb.models['Model-1'].rootAssembly.Set(edges=
    mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].edges.findAt(((
    0.006542, 0.0, 0.00075), ), ((0.011331, 0.0, 0.00075), ), ((0.035133, 0.0, 
    0.00075), ), ((0.030334, 0.0, 0.00075), ), ((0.045115, 0.0, 0.00075), ), ((
    0.018664, 0.0, 0.00075), ), ((0.023115, 0.0, 0.00075), ), ((0.000696, 0.0, 
    0.00075), ), ), name='edge_8_fibers')
mdb.models['Model-1'].rootAssembly.Set(edges=
    mdb.models['Model-1'].rootAssembly.instances['MATRIX_WITH_HOLES-1'].edges.findAt(
    ((0.042298, 0.0, 0.00075), ), ((0.028464, 0.0, 0.00075), ), ((0.004705, 
    0.0, 0.00075), ), ((0.016894, 0.0, 0.00075), ), ), name='edge_8_matrix')
mdb.models['Model-1'].rootAssembly.Set(edges=
    mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].edges.findAt(((
    0.0, 0.05, 0.000188), )), name='edge_9_fibers')
mdb.models['Model-1'].rootAssembly.Set(edges=
    mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].edges.findAt(((
    0.05, 0.05, 0.000188), )), name='edge_10_fibers')
mdb.models['Model-1'].rootAssembly.Set(edges=
    mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].edges.findAt(((
    0.045115, 0.05, 0.0), ), ((0.035133, 0.05, 0.0), ), ((0.030334, 0.05, 0.0), 
    ), ((0.023115, 0.05, 0.0), ), ((0.018664, 0.05, 0.0), ), ((0.011331, 0.05, 
    0.0), ), ((0.006542, 0.05, 0.0), ), ((0.000696, 0.05, 0.0), ), ), name=
    'edge_11_fibers')
mdb.models['Model-1'].rootAssembly.Set(edges=
    mdb.models['Model-1'].rootAssembly.instances['MATRIX_WITH_HOLES-1'].edges.findAt(
    ((0.003423, 0.05, 0.0), ), ((0.01558, 0.05, 0.0), ), ((0.027123, 0.05, 
    0.0), ), ((0.03992, 0.05, 0.0), ), ), name='edge_11_matrix')
mdb.models['Model-1'].rootAssembly.Set(edges=
    mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].edges.findAt(((
    0.048372, 0.05, 0.00075), ), ((0.032734, 0.05, 0.00075), ), ((0.037532, 
    0.05, 0.00075), ), ((0.02089, 0.05, 0.00075), ), ((0.02534, 0.05, 0.00075), 
    ), ((0.008937, 0.05, 0.00075), ), ((0.013725, 0.05, 0.00075), ), ((
    0.002087, 0.05, 0.00075), ), ), name='edge_12_fibers')
mdb.models['Model-1'].rootAssembly.Set(edges=
    mdb.models['Model-1'].rootAssembly.instances['MATRIX_WITH_HOLES-1'].edges.findAt(
    ((0.003423, 0.05, 0.00075), ), ((0.01558, 0.05, 0.00075), ), ((0.027123, 
    0.05, 0.00075), ), ((0.03992, 0.05, 0.00075), ), ), name='edge_12_matrix')
	
# extract nodes
r = r

### Dobule check abaove parts. 
### Investigate if all those sets were defined correctly.
def get_nodes(node_arr):
	el_nodes = np.zeros((len(node_arr), 3))	
	k = 0	
	for i in node_arr:		
		el_nodes[k,:] = [i.coordinates[0], i.coordinates[1], i.coordinates[2]]
		k = k +1
	return el_nodes
	
# corner nodes
c1_nodes = get_nodes(r.sets['CORNER_1'].nodes)
c2_nodes = get_nodes(r.sets['CORNER_2'].nodes)
c3_nodes = get_nodes(r.sets['CORNER_3'].nodes)
c4_nodes = get_nodes(r.sets['CORNER_4'].nodes)
c5_nodes = get_nodes(r.sets['CORNER_5'].nodes)
c6_nodes = get_nodes(r.sets['CORNER_6'].nodes)
c7_nodes = get_nodes(r.sets['CORNER_7'].nodes)
c8_nodes = get_nodes(r.sets['CORNER_8'].nodes)
corners = [c1_nodes, c2_nodes, c3_nodes, c4_nodes, c5_nodes, c6_nodes, c7_nodes, c8_nodes]
# edge nodes on matrix		
e1_matrix_nodes = get_nodes(r.sets['edge_1_matrix'].nodes)	
e2_matrix_nodes = get_nodes(r.sets['edge_2_matrix'].nodes)
e3_matrix_nodes = get_nodes(r.sets['edge_3_matrix'].nodes)
e4_matrix_nodes = get_nodes(r.sets['edge_4_matrix'].nodes)
e5_matrix_nodes = get_nodes(r.sets['edge_5_fibers'].nodes)	
e6_matrix_nodes = get_nodes(r.sets['edge_6_fibers'].nodes)
e7_matrix_nodes = get_nodes(r.sets['edge_7_matrix'].nodes)
e8_matrix_nodes = get_nodes(r.sets['edge_8_matrix'].nodes)
e9_matrix_nodes = get_nodes(r.sets['edge_9_fibers'].nodes)	
e10_matrix_nodes = get_nodes(r.sets['edge_10_fibers'].nodes)
e11_matrix_nodes = get_nodes(r.sets['edge_11_matrix'].nodes)
e12_matrix_nodes = get_nodes(r.sets['edge_12_matrix'].nodes)
# edge nodes on fibers
e1_fiber_nodes = get_nodes(r.sets['edge_1_fibers'].nodes)	
e2_fiber_nodes = get_nodes(r.sets['edge_2_fibers'].nodes)
e3_fiber_nodes = get_nodes(r.sets['edge_3_fibers'].nodes)
e4_fiber_nodes = get_nodes(r.sets['edge_4_fibers'].nodes)
e5_fiber_nodes = get_nodes(r.sets['edge_5_fibers'].nodes)	
e6_fiber_nodes = get_nodes(r.sets['edge_6_fibers'].nodes)
e7_fiber_nodes = get_nodes(r.sets['edge_7_fibers'].nodes)
e8_fiber_nodes = get_nodes(r.sets['edge_8_fibers'].nodes)
e9_fiber_nodes = get_nodes(r.sets['edge_9_fibers'].nodes)	
e10_fiber_nodes = get_nodes(r.sets['edge_10_fibers'].nodes)
e11_fiber_nodes = get_nodes(r.sets['edge_11_fibers'].nodes)
e12_fiber_nodes = get_nodes(r.sets['edge_12_fibers'].nodes)
	
# extract nodes on faces (X_POS, X_NEG, ..., Z_NEG)
XPOS_FIBER_nodes = get_nodes(r.sets['X_POS_FIBERS'].nodes)
XPOS_MATRIX_nodes = get_nodes(r.sets['X_POS_MATRIX'].nodes)

XNEG_FIBER_nodes = get_nodes(r.sets['X_NEG_FIBERS'].nodes)
XNEG_MATRIX_nodes = get_nodes(r.sets['X_NEG_MATRIX'].nodes)

YPOS_FIBER_nodes = get_nodes(r.sets['Y_POS_FIBERS'].nodes)
YPOS_MATRIX_nodes = get_nodes(r.sets['Y_POS_MATRIX'].nodes)

YNEG_FIBER_nodes = get_nodes(r.sets['Y_NEG_FIBERS'].nodes)
YNEG_MATRIX_nodes = get_nodes(r.sets['Y_NEG_MATRIX'].nodes)

ZPOS_FIBER_nodes = get_nodes(r.sets['Z_POS_FIBERS'].nodes)
ZPOS_MATRIX_nodes = get_nodes(r.sets['Z_POS_MATRIX'].nodes)

ZNEG_FIBER_nodes = get_nodes(r.sets['Z_NEG_FIBERS'].nodes)
ZNEG_MATRIX_nodes = get_nodes(r.sets['Z_NEG_MATRIX'].nodes)

def is_equal(A, B, name):
    print('\nThis is ', name)
    if len(A) == len(B):
	    print('Lengths are equal.')
    else:
	    print('check the length of POS AND NEG nodes')

is_equal(XPOS_FIBER_nodes, XNEG_FIBER_nodes, 'X_FIBERS')
is_equal(XPOS_MATRIX_nodes, XNEG_MATRIX_nodes, 'X_MATRIX')
is_equal(YPOS_FIBER_nodes, YNEG_FIBER_nodes, 'Y_FIBERS')
is_equal(YPOS_MATRIX_nodes, YNEG_MATRIX_nodes, 'Y_MATRIX')
is_equal(ZPOS_FIBER_nodes, ZNEG_FIBER_nodes, 'Z_FIBERS')
is_equal(ZPOS_MATRIX_nodes, ZNEG_MATRIX_nodes, 'Z_MATRIX')

# applied strain
eps_ = np.array([[1e-1, 0, 0],
				[0, -1e-1, 0],
				[0, 0, 0]])

# distance b/w corner nodes
# required for imposing displacements
l1 = c4_nodes[0] - c1_nodes[0]
l2 = c2_nodes[0] - c1_nodes[0]
l3 = c5_nodes[0] - c1_nodes[0]

# displacement at node c1, c2, c4 and c5
u1 = [0, 0, 0]
u2 = np.dot(eps_, l2)
u4 = np.dot(eps_, l1)
u5 = np.dot(eps_, l3)

# apply BC at corner nodes
mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName=
    'LOADING_STEP', distributionType=UNIFORM, fieldName='', fixed=OFF, 
    localCsys=None, name='C1', region=
    mdb.models['Model-1'].rootAssembly.sets['CORNER_1'], u1=0.0, u2=0.0, u3=0.0
    , ur1=UNSET, ur2=UNSET, ur3=UNSET)
	
mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName=
    'LOADING_STEP', distributionType=UNIFORM, fieldName='', fixed=OFF, 
    localCsys=None, name='C2', region=
    mdb.models['Model-1'].rootAssembly.sets['CORNER_2'], u1=u2[0], u2=u2[1], u3=u2[2]
    , ur1=UNSET, ur2=UNSET, ur3=UNSET)	
	
mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName=
    'LOADING_STEP', distributionType=UNIFORM, fieldName='', fixed=OFF, 
    localCsys=None, name='C4', region=
    mdb.models['Model-1'].rootAssembly.sets['CORNER_4'], u1=u4[0], u2=u4[1], u3=u4[2]
    , ur1=UNSET, ur2=UNSET, ur3=UNSET)		
	
mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName=
    'LOADING_STEP', distributionType=UNIFORM, fieldName='', fixed=OFF, 
    localCsys=None, name='C5', region=
    mdb.models['Model-1'].rootAssembly.sets['CORNER_5'], u1=u5[0], u2=u5[1], u3=u5[2]
    , ur1=UNSET, ur2=UNSET, ur3=UNSET)		
	
# assign equation constraints
for i in range(1,4):
#1. u8 - u5 - u4 = 0
	mdb.models['Model-1'].Equation(name="Constraint-1-%s" % (i), terms=((1.0, 'CORNER_8', 
		i), (-1.0, 'CORNER_5', i), (-1.0, 'CORNER_4', i)))	
#2. u7 - u5 - u4 - u2  = 0		
	mdb.models['Model-1'].Equation(name="Constraint-2-%s" % (i), terms=((1.0, 'CORNER_7', 
		i), (-1.0, 'CORNER_5', i), (-1.0, 'CORNER_4', i), (-1.0, 'CORNER_2', i)))			
#3. u6 - u2 - u5 = 0
	mdb.models['Model-1'].Equation(name="Constraint-3-%s" % (i), terms=((1.0, 'CORNER_6', 
		i), (-1.0, 'CORNER_5', i), (-1.0, 'CORNER_2', i)))
#4. u3 - u2 - u4 = 0
	mdb.models['Model-1'].Equation(name="Constraint-4-%s" % (i), terms=((1.0, 'CORNER_3', 
		i), (-1.0, 'CORNER_2', i), (-1.0, 'CORNER_4', i)))			

ind_ = 1
def sort_arr(arr_, ind_):
	'''
	sort an array arr_ regarding column ind_
	'''	
	return arr_[arr_[:,ind_].argsort()]

def remove_corners(edge_, coors):
	dist = 1e-6
	edge_2 = []
	for i,el in enumerate(edge_):
		if np.linalg.norm(el-coors[0]) < dist or\
		np.linalg.norm(el-coors[1]) < dist or\
		np.linalg.norm(el-coors[2]) < dist or\
		np.linalg.norm(el-coors[3]) < dist or\
		np.linalg.norm(el-coors[4]) < dist or\
		np.linalg.norm(el-coors[5]) < dist or\
		np.linalg.norm(el-coors[6]) < dist or\
		np.linalg.norm(el-coors[7]) < dist:
			pass
		else:
			edge_2.append(el)
			for i in range(8):
				print('diff:', np.linalg.norm(el-coors[1]) )
			
	return edge_2
r.regenerate()  

# PBC for edge block
dist = 1e-5
def get_pairs(e1, e2):
	n_e1 = len(e1)
	n_e2 = len(e2)	
	if n_e1 != n_e2:
		print('# of elements are not the same')
		print('Reducing element size may solve the problem')
		return -1
		
	val = [-999 for j in range(n_e1)]
	for j in range(n_e1):
		min_ = 1e6
		for i in range(n_e2):
			dist = np.linalg.norm(e1[j] - e2[i])
			if np.abs(dist) <= np.abs(min_):
				min_ = dist
				val[j] = i
	return val
	
# remove corner nodes. Manuel version	
def remove_corners_II(c1, e_):
	for i, e in enumerate(e_):
		if e[0] == c1[0][0] and e[1] == c1[0][1] and e[2] == c1[0][2]:
			e_ = np.delete(e_, i, 0)
	
	return e_
	
for c in corners:	
	e1_fiber_nodes = remove_corners_II(c, e1_fiber_nodes)	
	e1_matrix_nodes = remove_corners_II(c, e1_matrix_nodes)	
	
	e2_fiber_nodes = remove_corners_II(c, e2_fiber_nodes)	
	e2_matrix_nodes = remove_corners_II(c, e2_matrix_nodes)	

	e3_fiber_nodes = remove_corners_II(c, e3_fiber_nodes)	
	e3_matrix_nodes = remove_corners_II(c, e3_matrix_nodes)	

	e4_fiber_nodes = remove_corners_II(c, e4_fiber_nodes)	
	e4_matrix_nodes = remove_corners_II(c, e4_matrix_nodes)	

	e5_fiber_nodes = remove_corners_II(c, e5_fiber_nodes)	
	# e5_matrix_nodes = remove_corners_II(c, e5_matrix_nodes)	

	e6_fiber_nodes = remove_corners_II(c, e6_fiber_nodes)	
	# e6_matrix_nodes = remove_corners_II(c, e6_matrix_nodes)	

	e7_fiber_nodes = remove_corners_II(c, e7_fiber_nodes)	
	e7_matrix_nodes = remove_corners_II(c, e7_matrix_nodes)	

	e8_fiber_nodes = remove_corners_II(c, e8_fiber_nodes)	
	e8_matrix_nodes = remove_corners_II(c, e8_matrix_nodes)	

	e9_fiber_nodes = remove_corners_II(c, e9_fiber_nodes)	
	# e9_matrix_nodes = remove_corners_II(c, e9_matrix_nodes)	

	e10_fiber_nodes = remove_corners_II(c, e10_fiber_nodes)	
	# e10_matrix_nodes = remove_corners_II(c, e10_matrix_nodes)	

	e11_fiber_nodes = remove_corners_II(c, e11_fiber_nodes)	
	e11_matrix_nodes = remove_corners_II(c, e11_matrix_nodes)	

	e12_fiber_nodes = remove_corners_II(c, e12_fiber_nodes)	
	e12_matrix_nodes = remove_corners_II(c, e12_matrix_nodes)		
	
for c in corners:		
	XPOS_FIBER_nodes = remove_corners_II(c, XPOS_FIBER_nodes)	
	XPOS_MATRIX_nodes = remove_corners_II(c, XPOS_MATRIX_nodes)
	XNEG_FIBER_nodes = remove_corners_II(c, XNEG_FIBER_nodes)	
	XNEG_MATRIX_nodes = remove_corners_II(c, XNEG_MATRIX_nodes)	
	
	YPOS_FIBER_nodes = remove_corners_II(c, YPOS_FIBER_nodes)	
	YPOS_MATRIX_nodes = remove_corners_II(c, YPOS_MATRIX_nodes)
	YNEG_FIBER_nodes = remove_corners_II(c, YNEG_FIBER_nodes)	
	YNEG_MATRIX_nodes = remove_corners_II(c, YNEG_MATRIX_nodes)	

	ZPOS_FIBER_nodes = remove_corners_II(c, ZPOS_FIBER_nodes)	
	ZPOS_MATRIX_nodes = remove_corners_II(c, ZPOS_MATRIX_nodes)
	ZNEG_FIBER_nodes = remove_corners_II(c, ZNEG_FIBER_nodes)	
	ZNEG_MATRIX_nodes = remove_corners_II(c, ZNEG_MATRIX_nodes)		
	
	

# val_test = get_pairs(XPOS_nodes, XNEG_nodes)	

# new = np.copy(XPOS_nodes)		
# new2 = np.copy(XPOS_nodes)	
# for i, el in enumerate(XPOS_nodes):
	# for edge in e2_matrix_nodes:
		# if np.linalg.norm(el - edge) < 1e-05:
			# # print('Same nodes')
			# new = np.delete(new, i, 0)

def remove_edges(face, edge, TOL, BIG_NUMBER):
	for i, f in enumerate(face):
		for e in edge:
			if np.linalg.norm(f - e) < TOL:
				face[i][0] = BIG_NUMBER
				face[i][1] = BIG_NUMBER
				face[i][2] = BIG_NUMBER
	return face

TOL = 1e-04
BIG_NUMBER = 1e15
###
### arranging XPOS face	
###
### edges = 2, 4, 6, 10
## FIBERS
XPOS_FIBER_nodes = remove_edges(XPOS_FIBER_nodes, e2_matrix_nodes, TOL, BIG_NUMBER)
XPOS_FIBER_nodes = remove_edges(XPOS_FIBER_nodes, e4_matrix_nodes, TOL, BIG_NUMBER)
# XPOS_nodes = remove_edges(XPOS_nodes, e6_matrix_nodes, TOL, BIG_NUMBER)	
# XPOS_nodes = remove_edges(XPOS_nodes, e10_matrix_nodes, TOL, BIG_NUMBER)

XPOS_FIBER_nodes = remove_edges(XPOS_FIBER_nodes, e2_fiber_nodes, TOL, BIG_NUMBER)
XPOS_FIBER_nodes = remove_edges(XPOS_FIBER_nodes, e4_fiber_nodes, TOL, BIG_NUMBER)
XPOS_FIBER_nodes = remove_edges(XPOS_FIBER_nodes, e6_fiber_nodes, TOL, BIG_NUMBER)
XPOS_FIBER_nodes = remove_edges(XPOS_FIBER_nodes, e10_fiber_nodes, TOL, BIG_NUMBER)

XPOS_FIBER = []
for i in range(len(XPOS_FIBER_nodes)):
	if XPOS_FIBER_nodes[i][0] < BIG_NUMBER:
		XPOS_FIBER.append(XPOS_FIBER_nodes[i])
	
## MATRIX
XPOS_MATRIX_nodes = remove_edges(XPOS_MATRIX_nodes, e2_matrix_nodes, TOL, BIG_NUMBER)
XPOS_MATRIX_nodes = remove_edges(XPOS_MATRIX_nodes, e4_matrix_nodes, TOL, BIG_NUMBER)
# XPOS_nodes = remove_edges(XPOS_nodes, e6_matrix_nodes, TOL, BIG_NUMBER)	
# XPOS_nodes = remove_edges(XPOS_nodes, e10_matrix_nodes, TOL, BIG_NUMBER)

XPOS_MATRIX_nodes = remove_edges(XPOS_MATRIX_nodes, e2_fiber_nodes, TOL, BIG_NUMBER)
XPOS_MATRIX_nodes = remove_edges(XPOS_MATRIX_nodes, e4_fiber_nodes, TOL, BIG_NUMBER)
XPOS_MATRIX_nodes = remove_edges(XPOS_MATRIX_nodes, e6_fiber_nodes, TOL, BIG_NUMBER)
XPOS_MATRIX_nodes = remove_edges(XPOS_MATRIX_nodes, e10_fiber_nodes, TOL, BIG_NUMBER)

XPOS_MATRIX = []
for i in range(len(XPOS_MATRIX_nodes)):
	if XPOS_MATRIX_nodes[i][0] < BIG_NUMBER:
		XPOS_MATRIX.append(XPOS_MATRIX_nodes[i])
print('End of XPOS')
###		
### arranging XNEG face	
###
### edges = 1, 3, 5*, 9*
## FIBERS
XNEG_FIBER_nodes = remove_edges(XNEG_FIBER_nodes, e1_matrix_nodes, TOL, BIG_NUMBER)
XNEG_FIBER_nodes = remove_edges(XNEG_FIBER_nodes, e3_matrix_nodes, TOL, BIG_NUMBER)
# XPOS_nodes = remove_edges(XPOS_nodes, e6_matrix_nodes, TOL, BIG_NUMBER)	
# XPOS_nodes = remove_edges(XPOS_nodes, e10_matrix_nodes, TOL, BIG_NUMBER)

XNEG_FIBER_nodes = remove_edges(XNEG_FIBER_nodes, e1_fiber_nodes, TOL, BIG_NUMBER)
XNEG_FIBER_nodes = remove_edges(XNEG_FIBER_nodes, e3_fiber_nodes, TOL, BIG_NUMBER)
XNEG_FIBER_nodes = remove_edges(XNEG_FIBER_nodes, e5_fiber_nodes, TOL, BIG_NUMBER)
XNEG_FIBER_nodes = remove_edges(XNEG_FIBER_nodes, e9_fiber_nodes, TOL, BIG_NUMBER)

XNEG_FIBER = []
for i in range(len(XNEG_FIBER_nodes)):
	if XNEG_FIBER_nodes[i][0] < BIG_NUMBER:
		XNEG_FIBER.append(XNEG_FIBER_nodes[i])

## MATRIX
XNEG_MATRIX_nodes = remove_edges(XNEG_MATRIX_nodes, e1_matrix_nodes, TOL, BIG_NUMBER)
XNEG_MATRIX_nodes = remove_edges(XNEG_MATRIX_nodes, e3_matrix_nodes, TOL, BIG_NUMBER)
# XPOS_nodes = remove_edges(XPOS_nodes, e6_matrix_nodes, TOL, BIG_NUMBER)	
# XPOS_nodes = remove_edges(XPOS_nodes, e10_matrix_nodes, TOL, BIG_NUMBER)

XNEG_MATRIX_nodes = remove_edges(XNEG_MATRIX_nodes, e1_fiber_nodes, TOL, BIG_NUMBER)
XNEG_MATRIX_nodes = remove_edges(XNEG_MATRIX_nodes, e3_fiber_nodes, TOL, BIG_NUMBER)
XNEG_MATRIX_nodes = remove_edges(XNEG_MATRIX_nodes, e5_fiber_nodes, TOL, BIG_NUMBER)
XNEG_MATRIX_nodes = remove_edges(XNEG_MATRIX_nodes, e9_fiber_nodes, TOL, BIG_NUMBER)

XNEG_MATRIX = []
for i in range(len(XNEG_MATRIX_nodes)):
	if XNEG_MATRIX_nodes[i][0] < BIG_NUMBER:
		XNEG_MATRIX.append(XNEG_MATRIX_nodes[i])
print('End of XNEG')		
###
### arranging YPOS face	
###
### edges = 11, 12, 9*, 10*
## FIBERS
YPOS_FIBER_nodes = remove_edges(YPOS_FIBER_nodes, e11_matrix_nodes, TOL, BIG_NUMBER)
YPOS_FIBER_nodes = remove_edges(YPOS_FIBER_nodes, e12_matrix_nodes, TOL, BIG_NUMBER)
# XPOS_nodes = remove_edges(XPOS_nodes, e6_matrix_nodes, TOL, BIG_NUMBER)	
# XPOS_nodes = remove_edges(XPOS_nodes, e10_matrix_nodes, TOL, BIG_NUMBER)

YPOS_FIBER_nodes = remove_edges(YPOS_FIBER_nodes, e11_fiber_nodes, TOL, BIG_NUMBER)
YPOS_FIBER_nodes = remove_edges(YPOS_FIBER_nodes, e12_fiber_nodes, TOL, BIG_NUMBER)
YPOS_FIBER_nodes = remove_edges(YPOS_FIBER_nodes, e9_fiber_nodes, TOL, BIG_NUMBER)
YPOS_FIBER_nodes = remove_edges(YPOS_FIBER_nodes, e10_fiber_nodes, TOL, BIG_NUMBER)

YPOS_FIBER = []
for i in range(len(YPOS_FIBER_nodes)):
	if YPOS_FIBER_nodes[i][0] < BIG_NUMBER:
		YPOS_FIBER.append(YPOS_FIBER_nodes[i])
		
## MATRIX
YPOS_MATRIX_nodes = remove_edges(YPOS_MATRIX_nodes, e11_matrix_nodes, TOL, BIG_NUMBER)
YPOS_MATRIX_nodes = remove_edges(YPOS_MATRIX_nodes, e12_matrix_nodes, TOL, BIG_NUMBER)
# XPOS_nodes = remove_edges(XPOS_nodes, e6_matrix_nodes, TOL, BIG_NUMBER)	
# XPOS_nodes = remove_edges(XPOS_nodes, e10_matrix_nodes, TOL, BIG_NUMBER)

YPOS_MATRIX_nodes = remove_edges(YPOS_MATRIX_nodes, e11_fiber_nodes, TOL, BIG_NUMBER)
YPOS_MATRIX_nodes = remove_edges(YPOS_MATRIX_nodes, e12_fiber_nodes, TOL, BIG_NUMBER)
YPOS_MATRIX_nodes = remove_edges(YPOS_MATRIX_nodes, e9_fiber_nodes, TOL, BIG_NUMBER)
YPOS_MATRIX_nodes = remove_edges(YPOS_MATRIX_nodes, e10_fiber_nodes, TOL, BIG_NUMBER)

YPOS_MATRIX = []
for i in range(len(YPOS_MATRIX_nodes)):
	if YPOS_MATRIX_nodes[i][0] < BIG_NUMBER:
		YPOS_MATRIX.append(YPOS_MATRIX_nodes[i])
print('End of YPOS')		
###		
### arranging YNEG face	
###
### edges = 7, 8, 5*, 6*
## FIBERS
YNEG_FIBER_nodes = remove_edges(YNEG_FIBER_nodes, e7_matrix_nodes, TOL, BIG_NUMBER)
YNEG_FIBER_nodes = remove_edges(YNEG_FIBER_nodes, e8_matrix_nodes, TOL, BIG_NUMBER)
# XPOS_nodes = remove_edges(XPOS_nodes, e6_matrix_nodes, TOL, BIG_NUMBER)	
# XPOS_nodes = remove_edges(XPOS_nodes, e10_matrix_nodes, TOL, BIG_NUMBER)

YNEG_FIBER_nodes = remove_edges(YNEG_FIBER_nodes, e7_fiber_nodes, TOL, BIG_NUMBER)
YNEG_FIBER_nodes = remove_edges(YNEG_FIBER_nodes, e8_fiber_nodes, TOL, BIG_NUMBER)
YNEG_FIBER_nodes = remove_edges(YNEG_FIBER_nodes, e5_fiber_nodes, TOL, BIG_NUMBER)
YNEG_FIBER_nodes = remove_edges(YNEG_FIBER_nodes, e6_fiber_nodes, TOL, BIG_NUMBER)

YNEG_FIBER = []
for i in range(len(YNEG_FIBER_nodes)):
	if YNEG_FIBER_nodes[i][0] < BIG_NUMBER:
		YNEG_FIBER.append(YNEG_FIBER_nodes[i])
		
## MATRIX
YNEG_MATRIX_nodes = remove_edges(YNEG_MATRIX_nodes, e7_matrix_nodes, TOL, BIG_NUMBER)
YNEG_MATRIX_nodes = remove_edges(YNEG_MATRIX_nodes, e8_matrix_nodes, TOL, BIG_NUMBER)
# XPOS_nodes = remove_edges(XPOS_nodes, e6_matrix_nodes, TOL, BIG_NUMBER)	
# XPOS_nodes = remove_edges(XPOS_nodes, e10_matrix_nodes, TOL, BIG_NUMBER)

YNEG_MATRIX_nodes = remove_edges(YNEG_MATRIX_nodes, e7_fiber_nodes, TOL, BIG_NUMBER)
YNEG_MATRIX_nodes = remove_edges(YNEG_MATRIX_nodes, e8_fiber_nodes, TOL, BIG_NUMBER)
YNEG_MATRIX_nodes = remove_edges(YNEG_MATRIX_nodes, e5_fiber_nodes, TOL, BIG_NUMBER)
YNEG_MATRIX_nodes = remove_edges(YNEG_MATRIX_nodes, e6_fiber_nodes, TOL, BIG_NUMBER)

YNEG_MATRIX = []
for i in range(len(YNEG_MATRIX_nodes)):
	if YNEG_MATRIX_nodes[i][0] < BIG_NUMBER:
		YNEG_MATRIX.append(YNEG_MATRIX_nodes[i])	
print('End of YNEG')
###
### arranging ZPOS face	
###
### edges = 1, 2, 8, 12
## FIBERS
ZPOS_FIBER_nodes = remove_edges(ZPOS_FIBER_nodes, e1_matrix_nodes, TOL, BIG_NUMBER)
ZPOS_FIBER_nodes = remove_edges(ZPOS_FIBER_nodes, e2_matrix_nodes, TOL, BIG_NUMBER)
ZPOS_FIBER_nodes = remove_edges(ZPOS_FIBER_nodes, e8_matrix_nodes, TOL, BIG_NUMBER)
ZPOS_FIBER_nodes = remove_edges(ZPOS_FIBER_nodes, e12_matrix_nodes, TOL, BIG_NUMBER)

ZPOS_FIBER_nodes = remove_edges(ZPOS_FIBER_nodes, e1_fiber_nodes, TOL, BIG_NUMBER)
ZPOS_FIBER_nodes = remove_edges(ZPOS_FIBER_nodes, e2_fiber_nodes, TOL, BIG_NUMBER)
ZPOS_FIBER_nodes = remove_edges(ZPOS_FIBER_nodes, e8_fiber_nodes, TOL, BIG_NUMBER)
ZPOS_FIBER_nodes = remove_edges(ZPOS_FIBER_nodes, e12_fiber_nodes, TOL, BIG_NUMBER)

ZPOS_FIBER = []
for i in range(len(ZPOS_FIBER_nodes)):
	if ZPOS_FIBER_nodes[i][0] < BIG_NUMBER:
		ZPOS_FIBER.append(ZPOS_FIBER_nodes[i])
		
## MATRIX
ZPOS_MATRIX_nodes = remove_edges(ZPOS_MATRIX_nodes, e1_matrix_nodes, TOL, BIG_NUMBER)
ZPOS_MATRIX_nodes = remove_edges(ZPOS_MATRIX_nodes, e2_matrix_nodes, TOL, BIG_NUMBER)
ZPOS_MATRIX_nodes = remove_edges(ZPOS_MATRIX_nodes, e8_matrix_nodes, TOL, BIG_NUMBER)
ZPOS_MATRIX_nodes = remove_edges(ZPOS_MATRIX_nodes, e12_matrix_nodes, TOL, BIG_NUMBER)

ZPOS_MATRIX_nodes = remove_edges(ZPOS_MATRIX_nodes, e1_fiber_nodes, TOL, BIG_NUMBER)
ZPOS_MATRIX_nodes = remove_edges(ZPOS_MATRIX_nodes, e2_fiber_nodes, TOL, BIG_NUMBER)
ZPOS_MATRIX_nodes = remove_edges(ZPOS_MATRIX_nodes, e8_fiber_nodes, TOL, BIG_NUMBER)
ZPOS_MATRIX_nodes = remove_edges(ZPOS_MATRIX_nodes, e12_fiber_nodes, TOL, BIG_NUMBER)

ZPOS_MATRIX = []
for i in range(len(ZPOS_MATRIX_nodes)):
	if ZPOS_MATRIX_nodes[i][0] < BIG_NUMBER:
		ZPOS_MATRIX.append(ZPOS_MATRIX_nodes[i])
print('End of ZPOS')
###
### arranging ZNEG face	
###
### edges = 3, 4, 11, 7
## FIBERS
ZNEG_FIBER_nodes = remove_edges(ZNEG_FIBER_nodes, e3_matrix_nodes, TOL, BIG_NUMBER)
ZNEG_FIBER_nodes = remove_edges(ZNEG_FIBER_nodes, e4_matrix_nodes, TOL, BIG_NUMBER)
ZNEG_FIBER_nodes = remove_edges(ZNEG_FIBER_nodes, e11_matrix_nodes, TOL, BIG_NUMBER)
ZNEG_FIBER_nodes = remove_edges(ZNEG_FIBER_nodes, e7_matrix_nodes, TOL, BIG_NUMBER)

ZNEG_FIBER_nodes = remove_edges(ZNEG_FIBER_nodes, e3_fiber_nodes, TOL, BIG_NUMBER)
ZNEG_FIBER_nodes = remove_edges(ZNEG_FIBER_nodes, e4_fiber_nodes, TOL, BIG_NUMBER)
ZNEG_FIBER_nodes = remove_edges(ZNEG_FIBER_nodes, e11_fiber_nodes, TOL, BIG_NUMBER)
ZNEG_FIBER_nodes = remove_edges(ZNEG_FIBER_nodes, e7_fiber_nodes, TOL, BIG_NUMBER)

ZNEG_FIBER = []
for i in range(len(ZNEG_FIBER_nodes)):
	if ZNEG_FIBER_nodes[i][0] < BIG_NUMBER:
		ZNEG_FIBER.append(ZNEG_FIBER_nodes[i])
		
## MATRIX
ZNEG_MATRIX_nodes = remove_edges(ZNEG_MATRIX_nodes, e3_matrix_nodes, TOL, BIG_NUMBER)
ZNEG_MATRIX_nodes = remove_edges(ZNEG_MATRIX_nodes, e4_matrix_nodes, TOL, BIG_NUMBER)
ZNEG_MATRIX_nodes = remove_edges(ZNEG_MATRIX_nodes, e11_matrix_nodes, TOL, BIG_NUMBER)
ZNEG_MATRIX_nodes = remove_edges(ZNEG_MATRIX_nodes, e7_matrix_nodes, TOL, BIG_NUMBER)

ZNEG_MATRIX_nodes = remove_edges(ZNEG_MATRIX_nodes, e3_fiber_nodes, TOL, BIG_NUMBER)
ZNEG_MATRIX_nodes = remove_edges(ZNEG_MATRIX_nodes, e4_fiber_nodes, TOL, BIG_NUMBER)
ZNEG_MATRIX_nodes = remove_edges(ZNEG_MATRIX_nodes, e11_fiber_nodes, TOL, BIG_NUMBER)
ZNEG_MATRIX_nodes = remove_edges(ZNEG_MATRIX_nodes, e7_fiber_nodes, TOL, BIG_NUMBER)

ZNEG_MATRIX = []
for i in range(len(ZNEG_MATRIX_nodes)):
	if ZNEG_MATRIX_nodes[i][0] < BIG_NUMBER:
		ZNEG_MATRIX.append(ZNEG_MATRIX_nodes[i])
print('End of ZNEG')
		
### Relation b/w faces
### 1. XPOS - XNEG - CORNER4
val_XPOS_XNEG = [-999 for j in range(len(XPOS_FIBER))]
val_XPOS_XNEG = get_pairs(XPOS_FIBER, XNEG_FIBER)	
for i in range(len(XPOS_FIBER)):
	# chose elements by sphere
	dd = mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].nodes.getByBoundingSphere(XNEG_FIBER[i], dist)
	dd2 = mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].nodes.getByBoundingSphere(XPOS_FIBER[val_XPOS_XNEG[i]], dist)	
	mdb.models['Model-1'].rootAssembly.Set(name="f_XNEG_%s" % (i), nodes=dd)
	mdb.models['Model-1'].rootAssembly.Set(name="f_XPOS_%s" % (i), nodes=dd2)
	mdb.models['Model-1'].Equation(name="f_EQ_x_XPOS_XNEG_%s" % (i), terms=((1.0, "f_XPOS_%s" % (i), 1), (-1.0, "f_XNEG_%s" % (i), 1), (-1.0, 'CORNER_4', 1) ))
	mdb.models['Model-1'].Equation(name="f_EQ_y_XPOS_XNEG_%s" % (i), terms=((1.0, "f_XPOS_%s" % (i), 2), (-1.0, "f_XNEG_%s" % (i), 2), (-1.0, 'CORNER_4', 2) ))
	mdb.models['Model-1'].Equation(name="f_EQ_z_XPOS_XNEG_%s" % (i), terms=((1.0, "f_XPOS_%s" % (i), 3), (-1.0, "f_XNEG_%s" % (i), 3), (-1.0, 'CORNER_4', 3) ))
	
val_XPOS_XNEG = [-999 for j in range(len(XPOS_MATRIX))]
val_XPOS_XNEG = get_pairs(XPOS_MATRIX, XNEG_MATRIX)	
for i in range(len(XPOS_MATRIX)):
	# chose elements by sphere
	dd = mdb.models['Model-1'].rootAssembly.instances['MATRIX_WITH_HOLES-1'].nodes.getByBoundingSphere(XNEG_MATRIX[i], dist)
	dd2 = mdb.models['Model-1'].rootAssembly.instances['MATRIX_WITH_HOLES-1'].nodes.getByBoundingSphere(XPOS_MATRIX[val_XPOS_XNEG[i]], dist)	
	mdb.models['Model-1'].rootAssembly.Set(name="m_XNEG_%s" % (i), nodes=dd)
	mdb.models['Model-1'].rootAssembly.Set(name="m_XPOS_%s" % (i), nodes=dd2)
	mdb.models['Model-1'].Equation(name="m_EQ_x_XPOS_XNEG_%s" % (i), terms=((1.0, "m_XPOS_%s" % (i), 1), (-1.0, "m_XNEG_%s" % (i), 1), (-1.0, 'CORNER_4', 1) ))
	mdb.models['Model-1'].Equation(name="m_EQ_y_XPOS_XNEG_%s" % (i), terms=((1.0, "m_XPOS_%s" % (i), 2), (-1.0, "m_XNEG_%s" % (i), 2), (-1.0, 'CORNER_4', 2) ))
	mdb.models['Model-1'].Equation(name="m_EQ_z_XPOS_XNEG_%s" % (i), terms=((1.0, "m_XPOS_%s" % (i), 3), (-1.0, "m_XNEG_%s" % (i), 3), (-1.0, 'CORNER_4', 3) ))	
print('END OF 1. XPOS - XNEG - CORNER4')

### 2. YPOS - YNEG - CORNER2
val_YPOS_YNEG = [-999 for j in range(len(YPOS_FIBER))]
val_YPOS_YNEG = get_pairs(YPOS_FIBER, YNEG_FIBER)	
for i in range(len(YPOS_FIBER)):
	# chose elements by sphere
	dd = mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].nodes.getByBoundingSphere(YNEG_FIBER[val_YPOS_YNEG[i]], dist)
	dd2 = mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].nodes.getByBoundingSphere(YPOS_FIBER[i], dist)	
	mdb.models['Model-1'].rootAssembly.Set(name="f_YNEG_%s" % (i), nodes=dd)
	mdb.models['Model-1'].rootAssembly.Set(name="f_YPOS_%s" % (i), nodes=dd2)
	mdb.models['Model-1'].Equation(name="f_EQ_x_YPOS_YNEG_%s" % (i), terms=((1.0, "f_YPOS_%s" % (i), 1), (-1.0, "f_YNEG_%s" % (i), 1), (-1.0, 'CORNER_2', 1) ))
	mdb.models['Model-1'].Equation(name="f_EQ_y_YPOS_YNEG_%s" % (i), terms=((1.0, "f_YPOS_%s" % (i), 2), (-1.0, "f_YNEG_%s" % (i), 2), (-1.0, 'CORNER_2', 2) ))
	mdb.models['Model-1'].Equation(name="f_EQ_z_YPOS_YNEG_%s" % (i), terms=((1.0, "f_YPOS_%s" % (i), 3), (-1.0, "f_YNEG_%s" % (i), 3), (-1.0, 'CORNER_2', 3) ))

val_YPOS_YNEG = [-999 for j in range(len(YPOS_MATRIX))]
val_YPOS_YNEG = get_pairs(YPOS_MATRIX, YNEG_MATRIX)	
for i in range(len(YPOS_MATRIX)):
	# chose elements by sphere
	dd = mdb.models['Model-1'].rootAssembly.instances['MATRIX_WITH_HOLES-1'].nodes.getByBoundingSphere(YNEG_MATRIX[val_YPOS_YNEG[i]], dist)
	dd2 = mdb.models['Model-1'].rootAssembly.instances['MATRIX_WITH_HOLES-1'].nodes.getByBoundingSphere(YPOS_MATRIX[i], dist)	
	mdb.models['Model-1'].rootAssembly.Set(name="m_YNEG_%s" % (i), nodes=dd)
	mdb.models['Model-1'].rootAssembly.Set(name="m_YPOS_%s" % (i), nodes=dd2)
	mdb.models['Model-1'].Equation(name="m_EQ_x_YPOS_YNEG_%s" % (i), terms=((1.0, "m_YPOS_%s" % (i), 1), (-1.0, "m_YNEG_%s" % (i), 1), (-1.0, 'CORNER_2', 1) ))
	mdb.models['Model-1'].Equation(name="m_EQ_y_YPOS_YNEG_%s" % (i), terms=((1.0, "m_YPOS_%s" % (i), 2), (-1.0, "m_YNEG_%s" % (i), 2), (-1.0, 'CORNER_2', 2) ))
	mdb.models['Model-1'].Equation(name="m_EQ_z_YPOS_YNEG_%s" % (i), terms=((1.0, "m_YPOS_%s" % (i), 3), (-1.0, "m_YNEG_%s" % (i), 3), (-1.0, 'CORNER_2', 3) ))
print('END OF 2. YPOS - YNEG - CORNER2')

### 3. ZPOS - ZNEG - CORNER5
val_ZPOS_ZNEG = [-999 for j in range(len(ZPOS_FIBER))]
val_ZPOS_ZNEG = get_pairs(ZPOS_FIBER, ZNEG_FIBER)	
for i in range(len(ZPOS_FIBER)):
	# chose elements by sphere
	dd = mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].nodes.getByBoundingSphere(ZNEG_FIBER[i], dist)
	dd2 = mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].nodes.getByBoundingSphere(ZPOS_FIBER[val_ZPOS_ZNEG[i]], dist)	
	mdb.models['Model-1'].rootAssembly.Set(name="f_ZNEG_%s" % (i), nodes=dd)
	mdb.models['Model-1'].rootAssembly.Set(name="f_ZPOS_%s" % (i), nodes=dd2)
	mdb.models['Model-1'].Equation(name="f_EQ_x_ZPOS_ZNEG_%s" % (i), terms=((1.0, "f_ZPOS_%s" % (i), 1), (-1.0, "f_ZNEG_%s" % (i), 1), (-1.0, 'CORNER_5', 1) ))
	mdb.models['Model-1'].Equation(name="f_EQ_y_ZPOS_ZNEG_%s" % (i), terms=((1.0, "f_ZPOS_%s" % (i), 2), (-1.0, "f_ZNEG_%s" % (i), 2), (-1.0, 'CORNER_5', 2) ))
	mdb.models['Model-1'].Equation(name="f_EQ_z_ZPOS_ZNEG_%s" % (i), terms=((1.0, "f_ZPOS_%s" % (i), 3), (-1.0, "f_ZNEG_%s" % (i), 3), (-1.0, 'CORNER_5', 3) ))
	
val_ZPOS_ZNEG = [-999 for j in range(len(ZPOS_MATRIX))]
val_ZPOS_ZNEG = get_pairs(ZPOS_MATRIX, ZNEG_MATRIX)	
for i in range(len(ZPOS_MATRIX)):
	# chose elements by sphere
	dd = mdb.models['Model-1'].rootAssembly.instances['MATRIX_WITH_HOLES-1'].nodes.getByBoundingSphere(ZNEG_MATRIX[i], dist)
	dd2 = mdb.models['Model-1'].rootAssembly.instances['MATRIX_WITH_HOLES-1'].nodes.getByBoundingSphere(ZPOS_MATRIX[val_ZPOS_ZNEG[i]], dist)	
	mdb.models['Model-1'].rootAssembly.Set(name="m_ZNEG_%s" % (i), nodes=dd)
	mdb.models['Model-1'].rootAssembly.Set(name="m_ZPOS_%s" % (i), nodes=dd2)
	mdb.models['Model-1'].Equation(name="m_EQ_x_ZPOS_ZNEG_%s" % (i), terms=((1.0, "m_ZPOS_%s" % (i), 1), (-1.0, "m_ZNEG_%s" % (i), 1), (-1.0, 'CORNER_5', 1) ))
	mdb.models['Model-1'].Equation(name="m_EQ_y_ZPOS_ZNEG_%s" % (i), terms=((1.0, "m_ZPOS_%s" % (i), 2), (-1.0, "m_ZNEG_%s" % (i), 2), (-1.0, 'CORNER_5', 2) ))
	mdb.models['Model-1'].Equation(name="m_EQ_z_ZPOS_ZNEG_%s" % (i), terms=((1.0, "m_ZPOS_%s" % (i), 3), (-1.0, "m_ZNEG_%s" % (i), 3), (-1.0, 'CORNER_5', 3) ))
print('END OF 3. ZPOS - ZNEG - CORNER5')
	
## eq 1:= e1 - e3 - u5 = 0
val_e3_e1 = [-999 for j in range(len(e3_matrix_nodes))]
val_e3_e1 = get_pairs(e3_matrix_nodes, e1_matrix_nodes)		
for i in range(len(e3_matrix_nodes)):
	# chose elements by sphere
	dd = mdb.models['Model-1'].rootAssembly.instances['MATRIX_WITH_HOLES-1'].nodes.getByBoundingSphere(e1_matrix_nodes[i], dist)
	dd2 = mdb.models['Model-1'].rootAssembly.instances['MATRIX_WITH_HOLES-1'].nodes.getByBoundingSphere(e3_matrix_nodes[val_e3_e1[i]], dist)
	mdb.models['Model-1'].rootAssembly.Set(name="m_eq1_e1_%s_1" % (i), nodes=dd)
	mdb.models['Model-1'].rootAssembly.Set(name="m_eq1_e3_%s_2" % (i), nodes=dd2)
	mdb.models['Model-1'].Equation(name="m_EQ1_x_e1_e3_%s" % (i), terms=((1.0, "m_eq1_e1_%s_1" % (i), 1), (-1.0, "m_eq1_e3_%s_2" % (i), 1), (-1.0, 'CORNER_5', 1))) #
	mdb.models['Model-1'].Equation(name="m_EQ1_y_e1_e3_%s" % (i), terms=((1.0, "m_eq1_e1_%s_1" % (i), 2), (-1.0, "m_eq1_e3_%s_2" % (i), 2), (-1.0, 'CORNER_5', 2))) #	
	mdb.models['Model-1'].Equation(name="m_EQ1_z_e1_e3_%s" % (i), terms=((1.0, "m_eq1_e1_%s_1" % (i), 3), (-1.0, "m_eq1_e3_%s_2" % (i), 3), (-1.0, 'CORNER_5', 3)))		

val_e3_e1 = [-999 for j in range(len(e3_fiber_nodes))]
val_e3_e1 = get_pairs(e3_fiber_nodes, e1_fiber_nodes)	
		
for i in range(len(e3_fiber_nodes)):	
	dd3 = mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].nodes.getByBoundingSphere(e1_fiber_nodes[val_e3_e1[i]], dist)
	dd4 = mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].nodes.getByBoundingSphere(e3_fiber_nodes[i], dist)
	mdb.models['Model-1'].rootAssembly.Set(name="f_eq1_e1_%s_1" % (i), nodes=dd3)
	mdb.models['Model-1'].rootAssembly.Set(name="f_eq1_e3_%s_2" % (i), nodes=dd4)
	mdb.models['Model-1'].Equation(name="f_EQ1_x_e1_e3_%s" % (i), terms=((1.0, "f_eq1_e1_%s_1" % (i), 1), (-1.0, "f_eq1_e3_%s_2" % (i), 1), (-1.0, 'CORNER_5', 1)))	#
	mdb.models['Model-1'].Equation(name="f_EQ1_y_e1_e3_%s" % (i), terms=((1.0, "f_eq1_e1_%s_1" % (i), 2), (-1.0, "f_eq1_e3_%s_2" % (i), 2), (-1.0, 'CORNER_5', 2)))	#	
	mdb.models['Model-1'].Equation(name="f_EQ1_z_e1_e3_%s" % (i), terms=((1.0, "f_eq1_e1_%s_1" % (i), 3), (-1.0, "f_eq1_e3_%s_2" % (i), 3), (-1.0, 'CORNER_5', 3)))			
# end of block	
print('END OF 1:= e1 - e3 - u5 = 0')

## eq 2:= e4 - e3 - u4 = 0
val_e4_e3 = [-999 for j in range(len(e4_matrix_nodes))]
val_e4_e3 = get_pairs(e4_matrix_nodes, e3_matrix_nodes)		
for i in range(len(e4_matrix_nodes)):
	# chose elements by sphere
	dd = mdb.models['Model-1'].rootAssembly.instances['MATRIX_WITH_HOLES-1'].nodes.getByBoundingSphere(e3_matrix_nodes[i], dist)
	dd2 = mdb.models['Model-1'].rootAssembly.instances['MATRIX_WITH_HOLES-1'].nodes.getByBoundingSphere(e4_matrix_nodes[val_e4_e3[i]], dist)
	mdb.models['Model-1'].rootAssembly.Set(name="m_eq2_e3_%s_1" % (i), nodes=dd)
	mdb.models['Model-1'].rootAssembly.Set(name="m_eq2_e4_%s_2" % (i), nodes=dd2)
	mdb.models['Model-1'].Equation(name="m_EQ2_x_e1_e4_%s" % (i), terms=((1.0, "m_eq2_e4_%s_2" % (i), 1), (-1.0, "m_eq2_e3_%s_1" % (i), 1), (-1.0, 'CORNER_4', 1))) #
	mdb.models['Model-1'].Equation(name="m_EQ2_y_e1_e4_%s" % (i), terms=((1.0, "m_eq2_e4_%s_2" % (i), 2), (-1.0, "m_eq2_e3_%s_1" % (i), 2), (-1.0, 'CORNER_4', 2)))	
	mdb.models['Model-1'].Equation(name="m_EQ2_z_e1_e4_%s" % (i), terms=((1.0, "m_eq2_e4_%s_2" % (i), 3), (-1.0, "m_eq2_e3_%s_1" % (i), 3), (-1.0, 'CORNER_4', 3)))	#

val_e4_e3 = [-999 for j in range(len(e4_fiber_nodes))]
val_e4_e3 = get_pairs(e4_fiber_nodes, e3_fiber_nodes)		
for i in range(len(e4_fiber_nodes)):
	# chose elements by sphere
	dd = mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].nodes.getByBoundingSphere(e3_fiber_nodes[i], dist)
	dd2 = mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].nodes.getByBoundingSphere(e4_fiber_nodes[val_e4_e3[i]], dist)
	mdb.models['Model-1'].rootAssembly.Set(name="f_eq2_e3_%s_1" % (i), nodes=dd)
	mdb.models['Model-1'].rootAssembly.Set(name="f_eq2_e4_%s_2" % (i), nodes=dd2)
	mdb.models['Model-1'].Equation(name="f_EQ2_x_e1_e4_%s" % (i), terms=((1.0, "f_eq2_e4_%s_2" % (i), 1), (-1.0, "f_eq2_e3_%s_1" % (i), 1), (-1.0, 'CORNER_4', 1))) #
	mdb.models['Model-1'].Equation(name="f_EQ2_y_e1_e4_%s" % (i), terms=((1.0, "f_eq2_e4_%s_2" % (i), 2), (-1.0, "f_eq2_e3_%s_1" % (i), 2), (-1.0, 'CORNER_4', 2)))	
	mdb.models['Model-1'].Equation(name="f_EQ2_z_e1_e4_%s" % (i), terms=((1.0, "f_eq2_e4_%s_2" % (i), 3), (-1.0, "f_eq2_e3_%s_1" % (i), 3), (-1.0, 'CORNER_4', 3)))	#	
## end of eq. 2
print('END OF 2:= e4 - e3 - u4 = 0')
	
## eq 3:= e9 - e5 - u2 = 0
val_e9_e5 = [-999 for j in range(len(e9_fiber_nodes))]
val_e9_e5 = get_pairs(e9_fiber_nodes, e5_fiber_nodes)		
for i in range(len(e9_fiber_nodes)):
	# chose elements by sphere
	dd = mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].nodes.getByBoundingSphere(e5_fiber_nodes[i], dist)
	dd2 = mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].nodes.getByBoundingSphere(e9_fiber_nodes[val_e9_e5[i]], dist)
	mdb.models['Model-1'].rootAssembly.Set(name="f_eq3_e5_%s_1" % (i), nodes=dd)
	mdb.models['Model-1'].rootAssembly.Set(name="f_eq3_e9_%s_2" % (i), nodes=dd2)
	mdb.models['Model-1'].Equation(name="f_EQ3_x_e5_e9_%s" % (i), terms=((1.0, "f_eq3_e9_%s_2" % (i), 1), (-1.0, "f_eq3_e5_%s_1" % (i), 1), (-1.0, 'CORNER_2', 1))) #
	mdb.models['Model-1'].Equation(name="f_EQ3_y_e5_e9_%s" % (i), terms=((1.0, "f_eq3_e9_%s_2" % (i), 2), (-1.0, "f_eq3_e5_%s_1" % (i), 2), (-1.0, 'CORNER_2', 2)))	
	mdb.models['Model-1'].Equation(name="f_EQ3_z_e5_e9_%s" % (i), terms=((1.0, "f_eq3_e9_%s_2" % (i), 3), (-1.0, "f_eq3_e5_%s_1" % (i), 3), (-1.0, 'CORNER_2', 3)))	#
## end of eq. 3
print('END OF 3:= e9 - e5 - u2 = 0')

## eq4 := e11 - e7 - u2 = 0
val_e11_e7 = [-999 for j in range(len(e11_matrix_nodes))]
val_e11_e7 = get_pairs(e11_matrix_nodes, e7_matrix_nodes)	
for i in range(len(e11_matrix_nodes)):
	# chose elements by sphere
	dd = mdb.models['Model-1'].rootAssembly.instances['MATRIX_WITH_HOLES-1'].nodes.getByBoundingSphere(e7_matrix_nodes[val_e11_e7[i]], dist)
	dd2 = mdb.models['Model-1'].rootAssembly.instances['MATRIX_WITH_HOLES-1'].nodes.getByBoundingSphere(e11_matrix_nodes[i], dist)
	mdb.models['Model-1'].rootAssembly.Set(name="m_eq4_e7_%s_1" % (i), nodes=dd)
	mdb.models['Model-1'].rootAssembly.Set(name="m_eq4_e11_%s_2" % (i), nodes=dd2)
	mdb.models['Model-1'].Equation(name="m_EQ4_x_e7_e11_%s" % (i), terms=((1.0, "m_eq4_e11_%s_2" % (i), 1), (-1.0, "m_eq4_e7_%s_1" % (i), 1), (-1.0, 'CORNER_2', 1))) #
	mdb.models['Model-1'].Equation(name="m_EQ4_y_e7_e11_%s" % (i), terms=((1.0, "m_eq4_e11_%s_2" % (i), 2), (-1.0, "m_eq4_e7_%s_1" % (i), 2), (-1.0, 'CORNER_2', 2)))	
	mdb.models['Model-1'].Equation(name="m_EQ4_z_e7_e11_%s" % (i), terms=((1.0, "m_eq4_e11_%s_2" % (i), 3), (-1.0, "m_eq4_e7_%s_1" % (i), 3), (-1.0, 'CORNER_2', 3))) #
	
val_e11_e7 = [-999 for j in range(len(e11_fiber_nodes))]
val_e11_e7 = get_pairs(e11_fiber_nodes, e7_fiber_nodes)		
for i in range(len(e11_fiber_nodes)):
	# chose elements by sphere
	dd = mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].nodes.getByBoundingSphere(e7_fiber_nodes[val_e11_e7[i]], dist)
	dd2 = mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].nodes.getByBoundingSphere(e11_fiber_nodes[i], dist)
	mdb.models['Model-1'].rootAssembly.Set(name="f_eq4_e7_%s_1" % (i), nodes=dd)
	mdb.models['Model-1'].rootAssembly.Set(name="f_eq4_e11_%s_2" % (i), nodes=dd2)
	mdb.models['Model-1'].Equation(name="f_EQ4_x_e7_e11_%s" % (i), terms=((1.0, "f_eq4_e11_%s_2" % (i), 1), (-1.0, "f_eq4_e7_%s_1" % (i), 1), (-1.0, 'CORNER_2', 1))) #
	mdb.models['Model-1'].Equation(name="f_EQ4_y_e7_e11_%s" % (i), terms=((1.0, "f_eq4_e11_%s_2" % (i), 2), (-1.0, "f_eq4_e7_%s_1" % (i), 2), (-1.0, 'CORNER_2', 2)))	
	mdb.models['Model-1'].Equation(name="f_EQ4_z_e7_e11_%s" % (i), terms=((1.0, "f_eq4_e11_%s_2" % (i), 3), (-1.0, "f_eq4_e7_%s_1" % (i), 3), (-1.0, 'CORNER_2', 3)))		#
## end of eq. 4
print('END OF eq4 := e11 - e7 - u2 = 0')

### eq5 := e8 - e7 - u5 = 0
val_e8_e7 = [-999 for j in range(len(e8_matrix_nodes))]
val_e8_e7 = get_pairs(e8_matrix_nodes, e7_matrix_nodes)		
for i in range(len(e8_matrix_nodes)):
	# chose elements by sphere
	dd = mdb.models['Model-1'].rootAssembly.instances['MATRIX_WITH_HOLES-1'].nodes.getByBoundingSphere(e7_matrix_nodes[i], dist)
	dd2 = mdb.models['Model-1'].rootAssembly.instances['MATRIX_WITH_HOLES-1'].nodes.getByBoundingSphere(e8_matrix_nodes[val_e8_e7[i]], dist)
	mdb.models['Model-1'].rootAssembly.Set(name="m_eq5_e7_%s_1" % (i), nodes=dd)
	mdb.models['Model-1'].rootAssembly.Set(name="m_eq5_e8_%s_2" % (i), nodes=dd2)
	mdb.models['Model-1'].Equation(name="m_EQ5_x_e7_e8_%s" % (i), terms=((1.0, "m_eq5_e8_%s_2" % (i), 1), (-1.0, "m_eq5_e7_%s_1" % (i), 1), (-1.0, 'CORNER_5', 1))) #
	mdb.models['Model-1'].Equation(name="m_EQ5_y_e7_e8_%s" % (i), terms=((1.0, "m_eq5_e8_%s_2" % (i), 2), (-1.0, "m_eq5_e7_%s_1" % (i), 2), (-1.0, 'CORNER_5', 2)))	#
	mdb.models['Model-1'].Equation(name="m_EQ5_z_e7_e8_%s" % (i), terms=((1.0, "m_eq5_e8_%s_2" % (i), 3), (-1.0, "m_eq5_e7_%s_1" % (i), 3), (-1.0, 'CORNER_5', 3)))	
	
val_e8_e7 = [-999 for j in range(len(e8_fiber_nodes))]
val_e8_e7 = get_pairs(e8_fiber_nodes, e7_fiber_nodes)		
for i in range(len(e8_fiber_nodes)):
	# chose elements by sphere
	dd = mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].nodes.getByBoundingSphere(e7_fiber_nodes[i], dist)
	dd2 = mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].nodes.getByBoundingSphere(e8_fiber_nodes[val_e8_e7[i]], dist)
	mdb.models['Model-1'].rootAssembly.Set(name="f_eq5_e7_%s_1" % (i), nodes=dd)
	mdb.models['Model-1'].rootAssembly.Set(name="f_eq5_e8_%s_2" % (i), nodes=dd2)
	mdb.models['Model-1'].Equation(name="f_EQ5_x_e7_e8_%s" % (i), terms=((1.0, "f_eq5_e8_%s_2" % (i), 1), (-1.0, "f_eq5_e7_%s_1" % (i), 1), (-1.0, 'CORNER_5', 1))) #
	mdb.models['Model-1'].Equation(name="f_EQ5_y_e7_e8_%s" % (i), terms=((1.0, "f_eq5_e8_%s_2" % (i), 2), (-1.0, "f_eq5_e7_%s_1" % (i), 2), (-1.0, 'CORNER_5', 2))) #	
	mdb.models['Model-1'].Equation(name="f_EQ5_z_e7_e8_%s" % (i), terms=((1.0, "f_eq5_e8_%s_2" % (i), 3), (-1.0, "f_eq5_e7_%s_1" % (i), 3), (-1.0, 'CORNER_5', 3)))		
### end of eq. 5
print('END OF eq5 := e8 - e7 - u5 = 0')

### eq 6:= e6 - e5 - u4 = 0
val_e6_e5 = [-999 for j in range(len(e6_fiber_nodes))]
val_e6_e5 = get_pairs(e6_fiber_nodes, e5_fiber_nodes)		
for i in range(len(e6_fiber_nodes)):
	# chose elements by sphere
	dd = mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].nodes.getByBoundingSphere(e5_fiber_nodes[i], dist)
	dd2 = mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].nodes.getByBoundingSphere(e6_fiber_nodes[val_e9_e5[i]], dist)
	mdb.models['Model-1'].rootAssembly.Set(name="f_eq6_e5_%s_1" % (i), nodes=dd)
	mdb.models['Model-1'].rootAssembly.Set(name="f_eq6_e6_%s_2" % (i), nodes=dd2)
	mdb.models['Model-1'].Equation(name="f_EQ6_x_e5_e6_%s" % (i), terms=((1.0, "f_eq6_e6_%s_2" % (i), 1), (-1.0, "f_eq6_e5_%s_1" % (i), 1), (-1.0, 'CORNER_4', 1)))
	mdb.models['Model-1'].Equation(name="f_EQ6_y_e5_e6_%s" % (i), terms=((1.0, "f_eq6_e6_%s_2" % (i), 2), (-1.0, "f_eq6_e5_%s_1" % (i), 2), (-1.0, 'CORNER_4', 2)))	#
	mdb.models['Model-1'].Equation(name="f_EQ6_z_e5_e6_%s" % (i), terms=((1.0, "f_eq6_e6_%s_2" % (i), 3), (-1.0, "f_eq6_e5_%s_1" % (i), 3), (-1.0, 'CORNER_4', 3)))	#	
### end of eq. 6
print('END OF 6:= e6 - e5 - u4 = 0')

### eq7 := e2 - e3 - u4 - u5 = 0
val_e3_e2 = [-999 for j in range(len(e3_matrix_nodes))]
val_e3_e2 = get_pairs(e3_matrix_nodes, e2_matrix_nodes)	
for i in range(len(e3_matrix_nodes)):
	# chose elements by sphere
	dd = mdb.models['Model-1'].rootAssembly.instances['MATRIX_WITH_HOLES-1'].nodes.getByBoundingSphere(e2_matrix_nodes[i], dist)
	dd2 = mdb.models['Model-1'].rootAssembly.instances['MATRIX_WITH_HOLES-1'].nodes.getByBoundingSphere(e3_matrix_nodes[val_e3_e2[i]], dist)	
	mdb.models['Model-1'].rootAssembly.Set(name="m_eq7_e2_%s_1" % (i), nodes=dd)
	mdb.models['Model-1'].rootAssembly.Set(name="m_eq7_e3_%s_2" % (i), nodes=dd2)
	mdb.models['Model-1'].Equation(name="m_EQ7_x_e2_e3_%s" % (i), terms=((1.0, "m_eq7_e2_%s_1" % (i), 1), (-1.0, "m_eq7_e3_%s_2" % (i), 1), (-1.0, 'CORNER_4', 1), (-1.0, 'CORNER_5', 1) ))
	mdb.models['Model-1'].Equation(name="m_EQ7_y_e2_e3_%s" % (i), terms=((1.0, "m_eq7_e2_%s_1" % (i), 2), (-1.0, "m_eq7_e3_%s_2" % (i), 2), (-1.0, 'CORNER_4', 2), (-1.0, 'CORNER_5', 2) ))
	mdb.models['Model-1'].Equation(name="m_EQ7_z_e2_e3_%s" % (i), terms=((1.0, "m_eq7_e2_%s_1" % (i), 3), (-1.0, "m_eq7_e3_%s_2" % (i), 3), (-1.0, 'CORNER_4', 3), (-1.0, 'CORNER_5', 3) )) #
	
val_e3_e2 = [-999 for j in range(len(e3_fiber_nodes))]
val_e3_e2 = get_pairs(e3_fiber_nodes, e2_fiber_nodes)	
for i in range(len(e3_fiber_nodes)):
	# chose elements by sphere
	dd = mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].nodes.getByBoundingSphere(e2_fiber_nodes[i], dist)
	dd2 = mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].nodes.getByBoundingSphere(e3_fiber_nodes[val_e3_e2[i]], dist)	
	mdb.models['Model-1'].rootAssembly.Set(name="f_eq7_e2_%s_1" % (i), nodes=dd)
	mdb.models['Model-1'].rootAssembly.Set(name="f_eq7_e3_%s_2" % (i), nodes=dd2)
	mdb.models['Model-1'].Equation(name="f_EQ7_x_e2_e3_%s" % (i), terms=((1.0, "f_eq7_e2_%s_1" % (i), 1), (-1.0, "f_eq7_e3_%s_2" % (i), 1), (-1.0, 'CORNER_4', 1), (-1.0, 'CORNER_5', 1) ))
	mdb.models['Model-1'].Equation(name="f_EQ7_y_e2_e3_%s" % (i), terms=((1.0, "f_eq7_e2_%s_1" % (i), 2), (-1.0, "f_eq7_e3_%s_2" % (i), 2), (-1.0, 'CORNER_4', 2), (-1.0, 'CORNER_5', 2) ))
	mdb.models['Model-1'].Equation(name="f_EQ7_z_e2_e3_%s" % (i), terms=((1.0, "f_eq7_e2_%s_1" % (i), 3), (-1.0, "f_eq7_e3_%s_2" % (i), 3), (-1.0, 'CORNER_4', 3), (-1.0, 'CORNER_5', 3) ))	#
### end of eq. 7
print('END OF eq7 := e2 - e3 - u4 - u5 = 0')

### eq8 := e12 - e7 - u2 - u5 = 0
val_e12_e7 = [-999 for j in range(len(e12_matrix_nodes))]
val_e12_e7 = get_pairs(e12_matrix_nodes, e7_matrix_nodes)	
for i in range(len(e12_matrix_nodes)):
	# chose elements by sphere
	dd = mdb.models['Model-1'].rootAssembly.instances['MATRIX_WITH_HOLES-1'].nodes.getByBoundingSphere(e7_matrix_nodes[val_e12_e7[i]], dist)
	dd2 = mdb.models['Model-1'].rootAssembly.instances['MATRIX_WITH_HOLES-1'].nodes.getByBoundingSphere(e12_matrix_nodes[i], dist)	
	mdb.models['Model-1'].rootAssembly.Set(name="m_eq8_e7_%s_2" % (i), nodes=dd)
	mdb.models['Model-1'].rootAssembly.Set(name="m_eq8_e12_%s_1" % (i), nodes=dd2)
	mdb.models['Model-1'].Equation(name="m_EQ8_x_e2_e3_%s" % (i), terms=((1.0, "m_eq8_e12_%s_1" % (i), 1), (-1.0, "m_eq8_e7_%s_2" % (i), 1), (-1.0, 'CORNER_2', 1), (-1.0, 'CORNER_5', 1) )) #
	mdb.models['Model-1'].Equation(name="m_EQ8_y_e2_e3_%s" % (i), terms=((1.0, "m_eq8_e12_%s_1" % (i), 2), (-1.0, "m_eq8_e7_%s_2" % (i), 2), (-1.0, 'CORNER_2', 2), (-1.0, 'CORNER_5', 2) ))
	mdb.models['Model-1'].Equation(name="m_EQ8_z_e2_e3_%s" % (i), terms=((1.0, "m_eq8_e12_%s_1" % (i), 3), (-1.0, "m_eq8_e7_%s_2" % (i), 3), (-1.0, 'CORNER_2', 3), (-1.0, 'CORNER_5', 3) ))
	
val_e12_e7 = [-999 for j in range(len(e12_fiber_nodes))]
val_e12_e7 = get_pairs(e12_fiber_nodes, e7_fiber_nodes)	
for i in range(len(e12_fiber_nodes)):
	# chose elements by sphere
	dd = mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].nodes.getByBoundingSphere(e7_fiber_nodes[val_e12_e7[i]], dist)
	dd2 = mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].nodes.getByBoundingSphere(e12_fiber_nodes[i], dist)	
	mdb.models['Model-1'].rootAssembly.Set(name="f_eq8_e7_%s_2" % (i), nodes=dd)
	mdb.models['Model-1'].rootAssembly.Set(name="f_eq8_e12_%s_1" % (i), nodes=dd2)
	mdb.models['Model-1'].Equation(name="f_EQ8_x_e2_e3_%s" % (i), terms=((1.0, "f_eq8_e12_%s_1" % (i), 1), (-1.0, "f_eq8_e7_%s_2" % (i), 1), (-1.0, 'CORNER_2', 1), (-1.0, 'CORNER_5', 1) )) #
	mdb.models['Model-1'].Equation(name="f_EQ8_y_e2_e3_%s" % (i), terms=((1.0, "f_eq8_e12_%s_1" % (i), 2), (-1.0, "f_eq8_e7_%s_2" % (i), 2), (-1.0, 'CORNER_2', 2), (-1.0, 'CORNER_5', 2) ))
	mdb.models['Model-1'].Equation(name="f_EQ8_z_e2_e3_%s" % (i), terms=((1.0, "f_eq8_e12_%s_1" % (i), 3), (-1.0, "f_eq8_e7_%s_2" % (i), 3), (-1.0, 'CORNER_2', 3), (-1.0, 'CORNER_5', 3) ))
### end of eq. 8
print('END OF eq8 := e12 - e7 - u2 - u5 = 0')

### eq9 := e10 - e5 - u2 - u4 = 0
val_e10_e5 = [-999 for j in range(len(e10_fiber_nodes))]
val_e10_e5 = get_pairs(e10_fiber_nodes, e5_fiber_nodes)	
for i in range(len(e10_fiber_nodes)):
	# chose elements by sphere
	dd = mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].nodes.getByBoundingSphere(e5_fiber_nodes[i], dist)
	dd2 = mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].nodes.getByBoundingSphere(e10_fiber_nodes[val_e10_e5[i]], dist)	
	mdb.models['Model-1'].rootAssembly.Set(name="f_eq9_e5_%s_2" % (i), nodes=dd)
	mdb.models['Model-1'].rootAssembly.Set(name="f_eq9_e10_%s_1" % (i), nodes=dd2)
	mdb.models['Model-1'].Equation(name="f_EQ9_x_e2_e3_%s" % (i), terms=((1.0, "f_eq9_e10_%s_1" % (i), 1), (-1.0, "f_eq9_e5_%s_2" % (i), 1), (-1.0, 'CORNER_2', 1), (-1.0, 'CORNER_4', 1) ))
	mdb.models['Model-1'].Equation(name="f_EQ9_y_e2_e3_%s" % (i), terms=((1.0, "f_eq9_e10_%s_1" % (i), 2), (-1.0, "f_eq9_e5_%s_2" % (i), 2), (-1.0, 'CORNER_2', 2), (-1.0, 'CORNER_4', 2) ))
	mdb.models['Model-1'].Equation(name="f_EQ9_z_e2_e3_%s" % (i), terms=((1.0, "f_eq9_e10_%s_1" % (i), 3), (-1.0, "f_eq9_e5_%s_2" % (i), 3), (-1.0, 'CORNER_2', 3), (-1.0, 'CORNER_4', 3) )) #
### end of eq. 9	
print('END OF eq9 := e10 - e5 - u2 - u4 = 0')	
