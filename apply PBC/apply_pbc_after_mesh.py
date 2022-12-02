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

# Apply PBC after mesh generation.
# Number of elements at edges must be the same.
#  Use fixed element size for this purpose.
session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)

# ### new ....	
# # Define boundaries and loading
# # TODO: define sets for boundaries

r = mdb.models['Model-1'].rootAssembly
r.regenerate()
mdb.models['Model-1'].rootAssembly.Set(faces=
    mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].faces.findAt(((
    0.025, 0.024027, 0.0005), ), ((0.025, 0.00562, 0.0005), ), ((0.025, 
    0.001265, 0.0005), ), ((0.025, 0.011659, 0.0005), ), ), name=
    'X_POS_FIBERS')
mdb.models['Model-1'].rootAssembly.Set(faces=
    mdb.models['Model-1'].rootAssembly.instances['MATRIX_WITH_HOLES-1'].faces.findAt(
    ((0.025, 0.018921, 0.0005), ), ((0.025, 0.008812, 0.0005), ), ((0.025, 
    0.002679, 0.0005), ), ), name='X_POS_MATRIX')
mdb.models['Model-1'].rootAssembly.Set(faces=
    mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].faces.findAt(((
    0.001039, 0.025, 0.0005), ), ((0.009504, 0.025, 0.0005), ), ((0.023816, 
    0.025, 0.0005), ), ((0.020078, 0.025, 0.0005), ), ((0.004987, 0.025, 
    0.0005), ), ), name='Y_POS_FIBERS')
mdb.models['Model-1'].rootAssembly.Set(faces=
    mdb.models['Model-1'].rootAssembly.instances['MATRIX_WITH_HOLES-1'].faces.findAt(
    ((0.014729, 0.025, 0.0005), ), ((0.023046, 0.025, 0.0005), ), ), name=
    'Y_POS_MATRIX')
mdb.models['Model-1'].rootAssembly.Set(faces=
    mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].faces.findAt(((
    0.007036, 0.006739, 0.00075), ), ((0.005496, 0.007376, 0.00075), ), ((
    0.011617, 0.009357, 0.00075), ), ((0.010077, 0.009994, 0.00075), ), ((
    0.012518, 0.00423, 0.00075), ), ((0.010979, 0.004868, 0.00075), ), ((
    0.016236, 0.022143, 0.00075), ), ((0.014696, 0.02278, 0.00075), ), ((
    0.018301, 0.016217, 0.00075), ), ((0.016761, 0.016855, 0.00075), ), ((
    0.021681, 0.020798, 0.00075), ), ((0.020141, 0.021436, 0.00075), ), ((
    0.000304, 0.023027, 0.00075), ), ((0.011148, 0.024947, 0.00075), ), ((
    0.023798, 0.02483, 0.00075), ), ((0.024815, 0.005398, 0.00075), ), ((
    0.02158, 0.000133, 0.00075), ), ((0.002986, 0.011482, 0.00075), ), ((
    0.002986, 0.010891, 0.00075), ), ((0.00027, 0.001294, 0.00075), ), ((
    0.009606, 0.000266, 0.00075), ), ((0.002681, 0.00529, 0.00075), ), ((
    0.002681, 0.004674, 0.00075), ), ((0.006374, 0.000161, 0.00075), ), ((
    0.02393, 0.00028, 0.00075), ), ((0.024891, 0.011434, 0.00075), ), ((
    0.021332, 0.02478, 0.00075), ), ((0.006046, 0.024812, 0.00075), ), ((
    0.022582, 0.013173, 0.00075), ), ((0.021042, 0.013811, 0.00075), ), ((
    0.019503, 0.005963, 0.00075), ), ((0.017963, 0.006601, 0.00075), ), ((
    0.017888, 0.011014, 0.00075), ), ((0.016348, 0.011651, 0.00075), ), ((
    0.01357, 0.014123, 0.00075), ), ((0.01203, 0.014761, 0.00075), ), ((
    0.011767, 0.019014, 0.00075), ), ((0.010228, 0.019652, 0.00075), ), ((
    0.007787, 0.01314, 0.00075), ), ((0.006247, 0.013778, 0.00075), ), ((
    0.005046, 0.017581, 0.00075), ), ((0.003506, 0.018218, 0.00075), ), ), 
    name='Z_POS_FIBERS')
mdb.models['Model-1'].rootAssembly.Set(faces=
    mdb.models['Model-1'].rootAssembly.instances['MATRIX_WITH_HOLES-1'].faces.findAt(
    ((0.024478, 0.013515, 0.00075), )), name='Z_POS_MATRIX')
mdb.models['Model-1'].rootAssembly.Set(faces=
    mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].faces.findAt(((
    0.0, 0.023054, 0.0005), ), ((0.0, 0.011659, 0.0005), ), ((0.0, 0.010715, 
    0.00025), ), ((0.0, 0.000632, 0.0005), ), ((0.0, 0.00562, 0.0005), ), ((
    0.0, 0.004344, 0.00025), ), ), name='X_NEG_FIBERS')
mdb.models['Model-1'].rootAssembly.Set(faces=
    mdb.models['Model-1'].rootAssembly.instances['MATRIX_WITH_HOLES-1'].faces.findAt(
    ((0.0, 0.002288, 0.0005), ), ((0.0, 0.007853, 0.0005), ), ((0.0, 0.015762, 
    0.0005), ), ), name='X_NEG_MATRIX')
mdb.models['Model-1'].rootAssembly.Set(faces=
    mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].faces.findAt(((
    0.021518, 0.0, 0.0005), ), ((0.002079, 0.0, 0.0005), ), ((0.011139, 0.0, 
    0.0005), ), ((0.006269, 0.0, 0.0005), ), ((0.024408, 0.0, 0.0005), ), ), 
    name='Y_NEG_FIBERS')
mdb.models['Model-1'].rootAssembly.Set(faces=
    mdb.models['Model-1'].rootAssembly.instances['MATRIX_WITH_HOLES-1'].faces.findAt(
    ((0.023135, 0.0, 0.0005), ), ((0.016684, 0.0, 0.0005), ), ((0.007763, 0.0, 
    0.0005), ), ((0.003509, 0.0, 0.0005), ), ), name='Y_NEG_MATRIX')
mdb.models['Model-1'].rootAssembly.Set(faces=
    mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].faces.findAt(((
    0.007036, 0.007376, 0.0), ), ((0.005496, 0.006739, 0.0), ), ((0.011617, 
    0.009994, 0.0), ), ((0.010077, 0.009357, 0.0), ), ((0.012518, 0.004868, 
    0.0), ), ((0.010979, 0.00423, 0.0), ), ((0.016236, 0.02278, 0.0), ), ((
    0.014696, 0.022143, 0.0), ), ((0.018301, 0.016855, 0.0), ), ((0.016761, 
    0.016217, 0.0), ), ((0.021681, 0.021436, 0.0), ), ((0.020141, 0.020798, 
    0.0), ), ((0.002112, 0.024716, 0.0), ), ((0.010531, 0.022138, 0.0), ), ((
    0.024721, 0.023192, 0.0), ), ((0.024815, 0.005398, 0.0), ), ((0.021681, 
    0.000407, 0.0), ), ((0.002986, 0.010891, 0.0), ), ((0.002986, 0.011482, 
    0.0), ), ((0.001981, 0.000253, 0.0), ), ((0.011037, 0.000266, 0.0), ), ((
    0.002681, 0.004674, 0.0), ), ((0.002681, 0.00529, 0.0), ), ((0.006479, 
    0.001153, 0.0), ), ((0.024729, 0.001132, 0.0), ), ((0.024891, 0.011434, 
    0.0), ), ((0.020265, 0.02478, 0.0), ), ((0.005209, 0.024812, 0.0), ), ((
    0.022582, 0.013811, 0.0), ), ((0.021042, 0.013173, 0.0), ), ((0.019503, 
    0.006601, 0.0), ), ((0.017963, 0.005963, 0.0), ), ((0.017888, 0.011651, 
    0.0), ), ((0.016348, 0.011014, 0.0), ), ((0.01357, 0.014761, 0.0), ), ((
    0.01203, 0.014123, 0.0), ), ((0.011767, 0.019652, 0.0), ), ((0.010228, 
    0.019014, 0.0), ), ((0.007787, 0.013778, 0.0), ), ((0.006247, 0.01314, 
    0.0), ), ((0.005046, 0.018218, 0.0), ), ((0.003506, 0.017581, 0.0), ), ), 
    name='Z_NEG_FIBERS')
mdb.models['Model-1'].rootAssembly.Set(faces=
    mdb.models['Model-1'].rootAssembly.instances['MATRIX_WITH_HOLES-1'].faces.findAt(
    ((0.024478, 0.012877, 0.0), )), name='Z_NEG_MATRIX')
# Save by suley on 2022_09_10-15.54.26; build 6.14-1 2014_06_05-01.11.02 134264

# ## Construct Sets for P.B.C
# # |y
# # 2 3
# # 1 4-->x
CORNERS = ['CORNER_1', 'CORNER_2', 'CORNER_3', 'CORNER_4', 'CORNER_5', 'CORNER_6',\
'CORNER_7', 'CORNER_8']
r.Set(name='CORNER_1', vertices=
    r.instances['UD_FIBERS-1'].vertices.findAt(
    ((0.0, 0.0, 0.0), )))
	
mdb.models['Model-1'].rootAssembly.Set(name='CORNER_2', vertices=
    mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].vertices.findAt(
    ((0.0, 0.025, 0.0), )))
mdb.models['Model-1'].rootAssembly.Set(name='CORNER_3', vertices=
    mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].vertices.findAt(
    ((0.025, 0.025, 0.0), )))
mdb.models['Model-1'].rootAssembly.Set(name='CORNER_4', vertices=
    mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].vertices.findAt(
    ((0.025, 0.0, 0.0), )))
mdb.models['Model-1'].rootAssembly.Set(name='CORNER_5', vertices=
    mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].vertices.findAt(
    ((0.0, 0.0, 0.00075), )))
mdb.models['Model-1'].rootAssembly.Set(name='CORNER_6', vertices=
    mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].vertices.findAt(
    ((0.0, 0.025, 0.00075), )))
mdb.models['Model-1'].rootAssembly.Set(name='CORNER_7', vertices=
    mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].vertices.findAt(
    ((0.025, 0.025, 0.00075), )))
mdb.models['Model-1'].rootAssembly.Set(name='CORNER_8', vertices=
    mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].vertices.findAt(
    ((0.025, 0.0, 0.00075), )))

# # # define edge sets
# EDGES = ['EDGE_1', 'EDGE_2', 'EDGE_3', 'EDGE_4', 'EDGE_5', 'EDGE_6',\
# 'EDGE_7', 'EDGE_8', 'EDGE_9', 'EDGE_10', 'EDGE_11', 'EDGE_12']
mdb.models['Model-1'].rootAssembly.Set(edges=
    mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].edges.findAt(((
    0.0, 0.02281, 0.00075), ), ((0.0, 0.011541, 0.00075), ), ((0.0, 0.010125, 
    0.00075), ), ((0.0, 0.000474, 0.00075), ), ((0.0, 0.00546, 0.00075), ), ((
    0.0, 0.003547, 0.00075), ), ), name='edge_1_fibers')
mdb.models['Model-1'].rootAssembly.Set(edges=
    mdb.models['Model-1'].rootAssembly.instances['MATRIX_WITH_HOLES-1'].edges.findAt(
    ((0.0, 0.00219, 0.00075), ), ((0.0, 0.007614, 0.00075), ), ((0.0, 0.014972, 
    0.00075), ), ), name='edge_1_matrix')
mdb.models['Model-1'].rootAssembly.Set(edges=
    mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].edges.findAt(((
    0.025, 0.02427, 0.00075), ), ((0.025, 0.005938, 0.00075), ), ((0.025, 
    0.001423, 0.00075), ), ((0.025, 0.011895, 0.00075), ), ), name=
    'edge_2_fibers')
mdb.models['Model-1'].rootAssembly.Set(edges=
    mdb.models['Model-1'].rootAssembly.instances['MATRIX_WITH_HOLES-1'].edges.findAt(
    ((0.025, 0.019711, 0.00075), ), ((0.025, 0.009052, 0.00075), ), ((0.025, 
    0.002776, 0.00075), ), ), name='edge_2_matrix')
mdb.models['Model-1'].rootAssembly.Set(edges=
    mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].edges.findAt(((
    0.0, 0.02281, 0.0), ), ((0.0, 0.010125, 0.0), ), ((0.0, 0.011541, 0.0), ), 
    ((0.0, 0.000474, 0.0), ), ((0.0, 0.003547, 0.0), ), ((0.0, 0.00546, 0.0), 
    ), ), name='edge_3_fibers')
mdb.models['Model-1'].rootAssembly.Set(edges=
    mdb.models['Model-1'].rootAssembly.instances['MATRIX_WITH_HOLES-1'].edges.findAt(
    ((0.0, 0.00219, 0.0), ), ((0.0, 0.007614, 0.0), ), ((0.0, 0.014972, 0.0), 
    ), ), name='edge_3_matrix')
mdb.models['Model-1'].rootAssembly.Set(edges=
    mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].edges.findAt(((
    0.025, 0.02427, 0.0), ), ((0.025, 0.005938, 0.0), ), ((0.025, 0.001423, 
    0.0), ), ((0.025, 0.011895, 0.0), ), ), name='edge_4_fibers')
mdb.models['Model-1'].rootAssembly.Set(edges=
    mdb.models['Model-1'].rootAssembly.instances['MATRIX_WITH_HOLES-1'].edges.findAt(
    ((0.025, 0.019711, 0.0), ), ((0.025, 0.009052, 0.0), ), ((0.025, 0.002776, 
    0.0), ), ), name='edge_4_matrix')
mdb.models['Model-1'].rootAssembly.Set(edges=
    mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].edges.findAt(((
    0.0, 0.0, 0.000188), )), name='edge_5_fibers')
mdb.models['Model-1'].rootAssembly.Set(edges=
    mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].edges.findAt(((
    0.025, 0.0, 0.000188), )), name='edge_6_fibers')
mdb.models['Model-1'].rootAssembly.Set(edges=
    mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].edges.findAt(((
    0.021878, 0.0, 0.0), ), ((0.002339, 0.0, 0.0), ), ((0.011548, 0.0, 0.0), ), 
    ((0.006589, 0.0, 0.0), ), ((0.024556, 0.0, 0.0), ), ), name=
    'edge_7_fibers')
mdb.models['Model-1'].rootAssembly.Set(edges=
    mdb.models['Model-1'].rootAssembly.instances['MATRIX_WITH_HOLES-1'].edges.findAt(
    ((0.023157, 0.0, 0.0), ), ((0.017173, 0.0, 0.0), ), ((0.007789, 0.0, 0.0), 
    ), ((0.003558, 0.0, 0.0), ), ), name='edge_7_matrix')
mdb.models['Model-1'].rootAssembly.Set(edges=
    mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].edges.findAt(((
    0.021878, 0.0, 0.00075), ), ((0.002339, 0.0, 0.00075), ), ((0.011548, 0.0, 
    0.00075), ), ((0.006589, 0.0, 0.00075), ), ((0.024556, 0.0, 0.00075), ), ), 
    name='edge_8_fibers')
mdb.models['Model-1'].rootAssembly.Set(edges=
    mdb.models['Model-1'].rootAssembly.instances['MATRIX_WITH_HOLES-1'].edges.findAt(
    ((0.023157, 0.0, 0.00075), ), ((0.017173, 0.0, 0.00075), ), ((0.007789, 
    0.0, 0.00075), ), ((0.003558, 0.0, 0.00075), ), ), name='edge_8_matrix')
mdb.models['Model-1'].rootAssembly.Set(edges=
    mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].edges.findAt(((
    0.0, 0.025, 0.000188), )), name='edge_9_fibers')
mdb.models['Model-1'].rootAssembly.Set(edges=
    mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].edges.findAt(((
    0.025, 0.025, 0.000188), )), name='edge_10_fibers')
mdb.models['Model-1'].rootAssembly.Set(edges=
    mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].edges.findAt(((
    0.00078, 0.025, 0.0), ), ((0.009095, 0.025, 0.0), ), ((0.023668, 0.025, 
    0.0), ), ((0.019719, 0.025, 0.0), ), ((0.004666, 0.025, 0.0), ), ), name=
    'edge_11_fibers')
mdb.models['Model-1'].rootAssembly.Set(edges=
    mdb.models['Model-1'].rootAssembly.instances['MATRIX_WITH_HOLES-1'].edges.findAt(
    ((0.003265, 0.025, 0.0), ), ((0.00763, 0.025, 0.0), ), ((0.01424, 0.025, 
    0.0), ), ((0.023024, 0.025, 0.0), ), ), name='edge_11_matrix')
mdb.models['Model-1'].rootAssembly.Set(edges=
    mdb.models['Model-1'].rootAssembly.instances['UD_FIBERS-1'].edges.findAt(((
    0.00078, 0.025, 0.00075), ), ((0.009095, 0.025, 0.00075), ), ((0.023668, 
    0.025, 0.00075), ), ((0.019719, 0.025, 0.00075), ), ((0.004666, 0.025, 
    0.00075), ), ), name='edge_12_fibers')
mdb.models['Model-1'].rootAssembly.Set(edges=
    mdb.models['Model-1'].rootAssembly.instances['MATRIX_WITH_HOLES-1'].edges.findAt(
    ((0.003265, 0.025, 0.00075), ), ((0.00763, 0.025, 0.00075), ), ((0.01424, 
    0.025, 0.00075), ), ((0.023024, 0.025, 0.00075), ), ), name=
    'edge_12_matrix')

## extract nodes

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
# e5_matrix_nodes = get_nodes(r.sets['edge_5_fibers'].nodes)	
# e6_matrix_nodes = get_nodes(r.sets['edge_6_fibers'].nodes)
e7_matrix_nodes = get_nodes(r.sets['edge_7_matrix'].nodes)
e8_matrix_nodes = get_nodes(r.sets['edge_8_matrix'].nodes)
# e9_matrix_nodes = get_nodes(r.sets['edge_9_fibers'].nodes)	
# e10_matrix_nodes = get_nodes(r.sets['edge_10_fibers'].nodes)
e11_matrix_nodes = get_nodes(r.sets['edge_11_matrix'].nodes)
e12_matrix_nodes = get_nodes(r.sets['edge_12_matrix'].nodes)
#edge nodes on fibers
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
    'Step-1', distributionType=UNIFORM, fieldName='', fixed=OFF, 
    localCsys=None, name='C1', region=
    mdb.models['Model-1'].rootAssembly.sets['CORNER_1'], u1=0.0, u2=0.0, u3=0.0
    , ur1=UNSET, ur2=UNSET, ur3=UNSET)
	
mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName=
    'Step-1', distributionType=UNIFORM, fieldName='', fixed=OFF, 
    localCsys=None, name='C2', region=
    mdb.models['Model-1'].rootAssembly.sets['CORNER_2'], u1=u2[0], u2=u2[1], u3=u2[2]
    , ur1=UNSET, ur2=UNSET, ur3=UNSET)	
	
mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName=
    'Step-1', distributionType=UNIFORM, fieldName='', fixed=OFF, 
    localCsys=None, name='C4', region=
    mdb.models['Model-1'].rootAssembly.sets['CORNER_4'], u1=u4[0], u2=u4[1], u3=u4[2]
    , ur1=UNSET, ur2=UNSET, ur3=UNSET)		
	
mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName=
    'Step-1', distributionType=UNIFORM, fieldName='', fixed=OFF, 
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

TOL = 1e-05
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
	
print('nodes are arranged sucsessfully.')	
### Relation b/w faces
print('Relations b/w faces.')

print('1: XPOS - XNEG - CORNER4 is started')
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

print('1: XPOS - XNEG - CORNER4 is completed')

### 2. YPOS - YNEG - CORNER2
print('2:YPOS - YNEG - CORNER2 is started')
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

print('2:YPOS - YNEG - CORNER2 is completed')

### 3. ZPOS - ZNEG - CORNER5
print('3. ZPOS - ZNEG - CORNER5 is started.')
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

print('3.1')
val_ZPOS_ZNEG = [-999 for j in range(len(ZPOS_MATRIX))]
val_ZPOS_ZNEG = get_pairs(ZPOS_MATRIX, ZNEG_MATRIX)	
for i in range(len(ZNEG_MATRIX)):
	# chose elements by sphere
	dd = mdb.models['Model-1'].rootAssembly.instances['MATRIX_WITH_HOLES-1'].nodes.getByBoundingSphere(ZNEG_MATRIX[i], dist)
	dd2 = mdb.models['Model-1'].rootAssembly.instances['MATRIX_WITH_HOLES-1'].nodes.getByBoundingSphere(ZPOS_MATRIX[val_ZPOS_ZNEG[i]], dist)	
	mdb.models['Model-1'].rootAssembly.Set(name="m_ZNEG_%s" % (i), nodes=dd)
	mdb.models['Model-1'].rootAssembly.Set(name="m_ZPOS_%s" % (i), nodes=dd2)
	mdb.models['Model-1'].Equation(name="m_EQ_x_ZPOS_ZNEG_%s" % (i), terms=((1.0, "m_ZPOS_%s" % (i), 1), (-1.0, "m_ZNEG_%s" % (i), 1), (-1.0, 'CORNER_5', 1) ))
	mdb.models['Model-1'].Equation(name="m_EQ_y_ZPOS_ZNEG_%s" % (i), terms=((1.0, "m_ZPOS_%s" % (i), 2), (-1.0, "m_ZNEG_%s" % (i), 2), (-1.0, 'CORNER_5', 2) ))
	mdb.models['Model-1'].Equation(name="m_EQ_z_ZPOS_ZNEG_%s" % (i), terms=((1.0, "m_ZPOS_%s" % (i), 3), (-1.0, "m_ZNEG_%s" % (i), 3), (-1.0, 'CORNER_5', 3) ))

print('3. ZPOS - ZNEG - CORNER5 is completed')
	
## eq 1:= e1 - e3 - u5 = 0
print('eq 1:= e1 - e3 - u5 = 0 is started')
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
print('eq 1:= e1 - e3 - u5 = 0 is completed')

## eq 2:= e4 - e3 - u4 = 0
print('eq 2:= e4 - e3 - u4 = 0 is started')
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
print('eq 2:= e4 - e3 - u4 = 0 is completed')
	
## eq 3:= e9 - e5 - u2 = 0
print('eq 3:= e9 - e5 - u2 = 0 is started')
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
print('eq 3:= e9 - e5 - u2 = 0 is completed')

## eq4 := e11 - e7 - u2 = 0
print('eq4 := e11 - e7 - u2 = 0 is started')
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
print('eq4 := e11 - e7 - u2 = 0 is completed')

### eq5 := e8 - e7 - u5 = 0
print('eq5 := e8 - e7 - u5 = 0 is started')
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
print('eq5 := e8 - e7 - u5 = 0 is completed')

### eq 6:= e6 - e5 - u4 = 0
print('eq 6:= e6 - e5 - u4 = 0 is started')
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
print('eq 6:= e6 - e5 - u4 = 0 is completed')

### eq7 := e2 - e3 - u4 - u5 = 0
print('eq7 := e2 - e3 - u4 - u5 = 0 is started')
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
print('eq7 := e2 - e3 - u4 - u5 = 0 is completed')

### eq8 := e12 - e7 - u2 - u5 = 0
print('eq8 := e12 - e7 - u2 - u5 = 0 is started')
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
print('eq8 := e12 - e7 - u2 - u5 = 0 is completed')

### eq9 := e10 - e5 - u2 - u4 = 0
print('eq9 := e10 - e5 - u2 - u4 = 0 is started')
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
print('eq9 := e10 - e5 - u2 - u4 = 0 is completed')
### end of eq. 9
		

r.regenerate()    

# Create job
name_ = 'THICK_1'
mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
    explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
    memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF, 
    multiprocessingMode=DEFAULT, name=name_, 
    nodalOutputPrecision=SINGLE, numCpus=1, numGPUs=0, queue=None, 
    resultsFormat=ODB, scratch='', type=ANALYSIS, userSubroutine='', waitHours=
    0, waitMinutes=0)    

# # set field-output
r.regenerate()
MODEL.fieldOutputRequests['F-Output-1'].setValues(variables=(
    'S', 'E', 'U', 'RF', 'CSTRESS', 'CDISP', 'CFORCE', 'CSTATUS', 'CSDMG', 
    'CSMAXSCRT', 'CSMAXUCRT', 'CSQUADSCRT', 'CSQUADUCRT', 'NT', 'UVARM'))    

