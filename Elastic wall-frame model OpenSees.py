# -*- coding: utf-8 -*-
"""
Created on Wed May 22 11:45:11 2024

@author: Liam Pledger : liam.pledger@pg.canterbury.ac.nz

The following python script provides the user the mode-shape and fundamental 
period of the structure. The structure dimensions such as column and beam sizes,
wall length, number of storeys, can all be varied based on user inputs to 
obtain an initial period and mode-shape. 

Plots of the nodes and elements of the structure are created during the 
running of the script to help the user visualise the output.

"""

import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt
import opsvis as opsv
import vfo.vfo as vfo

ops.wipe()

ops.model('basic', '-ndm', 2, '-ndf', 3)

###################################################################################################
###################################################################################################
"""PARAMETERS TO BE VARIED BY THE USER"""

# Define section properties and elastic beam column elements
E = 25 * 10**6        # kN / m^2
NStories = 10         # No. of stories
HStory1 = 4.0         # 1st floor height (m)
HStoryTyp = 3.6       # Typical storey height (m)
HBuilding = HStory1 + (NStories - 1) * HStoryTyp  # height of building

NBays = 4            # No. of bays
 
WBay = 6.0            # unit m     Bay width

dc = 0.6              # unit m     Column depth

wb = 0.3              # unit m     Beam width
db = 0.4              # unit m     Beam depth

Lw = 6.0              # unit m     Wall length
tw = 0.25             # unit m     Wall thickness

seismic_weight = 8    # kPa        Seismic weight of the structure expressing in kPa

###################################################################################################
###################################################################################################


mass_frame = seismic_weight * WBay**2 * (NBays+ 1) / 9.81          # in Tonnes

Ic = dc**4 / 12        # m^4   Square column second moment of area
Ac = dc**2             # m^2   Column area

Ib = wb * db**3 / 12   # m^4   Beam second moment of area
Ab = wb*db             # m^2   Beam area

Iw = Lw**3 * tw / 12   # m^4   Wall second moment of area
Aw = Lw*tw          # m^2   Wall area

###################################################################################################
"Define nodes and cordinates"
###################################################################################################

build_height_array = np.zeros(NStories + 1)
for i in range(NStories + 1):
    build_height_array[i] = HStory1 + HStoryTyp * (i - 1)
    build_height_array[0] = 0
    build_height_array[1] = HStory1

"Saving plots based on type of analyses used"

ops.geomTransf('Linear', 1)    # Uses the natural deformation scheme, accounting for PDelta and rigid body motion - typically the most accurate

# Element linking matrix
# Main nodes -  format ij = i:xdirec j: ydirec
# j:1-5 floors, i: 1,2 main beam column nodes, + 2 centre nodes for the braces

"Main nodes for the structure"

""" Defining the main nodes and coordinates for the model at the ends of each beam and column for the 2D frame
    
    Nodes are defined as follows: numbering = ij where i is the column in the x direction starting at 1,2,3...n
    and j is the row / storey number starting at the ground which is the 1st floor  01, 02, 03...n.
    
    eg. the node at the bottom left is '101' with the coordinates [0.0, 0.0] """

main_coords_x = []
main_coords_y = []
for j in range(NStories + 1):
    for i in range(NBays + 1):
        main_coords_x.append(WBay * i)
        main_coords_y.append(build_height_array[j])
main_coords_y.sort()
main_coords = np.array([main_coords_x, main_coords_y])   # represents cor1

main_nodes = []
for i in range(1, NStories + 2):
    for j in range(1, NBays + 2):
        if i < 10:
            no = str(j) + str(0) + str(i)
        else:
            no = str(j) + str(i)
        main_nodes.append(int(no))

# Define nodes
ntag = main_nodes
cor = main_coords

# Nodes
for i in range(0, len(cor[0, :])):
    ops.node(int(ntag[i]), cor[0, i], cor[1, i])

###################################################################################################
"Define boundary conditions"
###################################################################################################
#fixity
"fix the columns at ground floor - assuming pinned connections for the system (including the leaning column)"

for i in range(NBays + 1):
    ops.fix(main_nodes[i], 1, 1, 1) # fixed

main_nodes_vert = main_nodes[::]
main_nodes_vert.sort()

opsv.plot_model(nodes_only=True, axis_off=1)
plt.title('Plot of Nodes - including zero-length beam and column hinge nodes')

###################################################################################################
" Define mass distribution of nodes "
###################################################################################################

#Mass distribution
mass_nodes = []
for i in range(len(main_nodes_vert)):
    mass_nodes.append(mass_frame / (NBays + 1))
    if main_nodes_vert[i] % 100 == 1:
        mass_nodes[i] = 0

for i in range(0, len(main_nodes)):
    "mass command is used to set the mass at each node"
    ops.mass(int(main_nodes_vert[i]), mass_nodes[i], mass_nodes[i], 0)

###################################################################################################
" Define beam, column, and wall elements between nodes"
###################################################################################################

"creating the element tags for the columns and walls"
# eleID convention:  "1xy" where 1 = col,  x = Pier #, y = Story #
col_eleTag = main_nodes[: -(NBays+1)]
for i in range(len(col_eleTag)):
    no = str(1) + str(col_eleTag[i])
    col_eleTag[i] = int(no)

"The wall is defined as the RHS column line"
n_col_ele = (NBays+1)*NStories
for i in range(0, n_col_ele):
    if (i) % (NBays + 1) == 0:
        ops.element('elasticBeamColumn', int(col_eleTag[i]), int(main_nodes[i]), int(main_nodes[i + NBays+1]), float(Aw), E, float(Iw), int(1))
    else:
        ops.element('elasticBeamColumn', int(col_eleTag[i]), int(main_nodes[i]), int(main_nodes[i + NBays+1]), float(Ac), E, float(Ic), int(1))
     
opsv.plot_model(node_labels=0, axis_off=1)
plt.title('Plot of Elements')

"creating the element tags for the beams"
# eleID convention:  "2xy" where 2 = beam, x = Bay #, y = Floor #
beam_eleTag = []
for j in range(1, NBays + 1):
    for i in range(1, NStories + 1):
        no = str(2) + str(j) + str(i)
        beam_eleTag.append(int(no))
        
mn = main_nodes[NBays + 1 :]  # removing the base nodes because they don't have beams between them
mn.sort()

"Defining the beam elements for the structure"
n_beam_ele = NBays * NStories
for i in range(0, n_beam_ele):
    ops.element('elasticBeamColumn', int(beam_eleTag[i]), int(mn[i]), int(mn[i + NStories]), float(Ab), E, float(Ib), int(1))
       
###################################################################################################
"Plotting the model"
###################################################################################################

opsv.plot_model(node_labels=0, axis_off=1)
plt.title('Plot of Elements')

###################################################################################################
"Applying vertical load to frame to determine the building period"
###################################################################################################
# Gravity load
vfo.createODB("Nonlin_RCWall", "Gravity", Nmodes=3)
ops.timeSeries('Linear', 1)  # applies the load in a linear manner (not all at once)
ops.pattern("Plain", 1, 1) # create a plain load pattern - similar to ELF

def create_cntrlnodes(n):
    return [i for i in range(101, 101 + n)]
cntrlnodes = create_cntrlnodes(NStories + 1)

"Create the nodal load - command: load nodeID xForce yForce"
for i in range(0, len(main_nodes)):
    "mass command is used to set the mass at each node"
    ops.load(int(main_nodes_vert[i]), 0, -mass_nodes[i] * 9.81, 0)
    
# create DOF number
ops.numberer("RCM")
# create SOE
ops.system('BandGeneral')
# create constraint handler
ops.constraints("Transformation")
# create number of steps
nsteps=1
# create integrator
ops.integrator('LoadControl', 1/nsteps)
# create algorithm
ops.test('RelativeEnergyIncr', 1e-1, 200, 0) # convergence scheme applied to the model
ops.algorithm("Newton")                      # solution scheme - Newton-Raphson method of solving non-linear equations
# create analysis object
ops.analysis("Static")                       # analysis type - ie, static, transient etc..
# perform the analysis
ops.recorder('Node', '-file', "results/modal/eigen.out",'-closeOnWrite','-dof',1,2,3,'eigen')

ops.analyze(nsteps)

vfo.createODB("Nonlin_RCWall", "Gravity", Nmodes=3)

# printModel()
# opsv.plot_model(fig_wi_he=(20., 14.))
ops.eigen('-genBandArpack', 1)
ops.record
ops.loadConst('-time', 0.0)
ops.wipeAnalysis()

a=ops.eigen(10)
w1=sqrt(a[0])                       # angular frequency of first mode
w2=sqrt(a[1])                       # angular frequency of second mode
w3=sqrt(a[2])

zeta=0.02                                  #Assumed damping ratio (currently 2%)
a0    =zeta*2.0*w1*w3/(w1 + w3);  	        # mass damping coefficient based on first and third modes
a1    =zeta*2.0/(w1 + w3);		        # stiffness damping coefficient based on first and third modes
ops.rayleigh(a0, 0, 0, a1)

ms_1 = np.zeros([len(cntrlnodes)])
for i in range(0, len(cntrlnodes)):
    ms_1[i] = ops.nodeEigenvector(cntrlnodes[i], 1, 1)   # obtains the mode shape of the first mode
ms_1 = ms_1/(ms_1[-1])                               # normalised the mode shape (value at the roof equal to 1)
ms_1 = np.round(ms_1, 4)

print('T1 (s) = ' +  str(np.round(2*np.pi/w1, 3)))
print('Fundamental mode-shape = ' + str(np.round(ms_1, 2)))

