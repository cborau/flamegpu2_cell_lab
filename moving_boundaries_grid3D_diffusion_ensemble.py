# +====================================================================+
# | Author: Carlos Borau, University of Zaragoza [cborau@unizar.es]    |
# | Version: 1.00                                                      |
# | Last update: 2023-07-17                                            |
# +====================================================================+
# +--------------------------------------------------------------------+
# | This code aims to simulate the cellular microenvironment incluing: |
# |                                                                    |
# |    ----------------------------------------------------------------|
# |    | Boundary conditions                                           |
# |    ----------------------------------------------------------------|
# |        - mechanical:                                               |
# |              . restricted movement                                 |
# |              . free movement                                       |
# |              . imposed displacements (linear, oscilatory)          |
# |              . sticky                                              |
# |              . parallel sliding                                    |
# |                                                                    |
# |        - chemical:                                                 |
# |              . initial concentration                               |
# |              . fixed concentration                                 |
# |                                                                    |
# |    ----------------------------------------------------------------|
# |    | Agents                                                        |
# |    ----------------------------------------------------------------|
# |        - Corner agents: visual purpose (they show domain limits)   |
# |                                                                    |
# |        - ECM agents: agents representing the extracelluar matrix   |
# |                                                                    |
# |              . mechanical relationship: dumped system              |
# |                                                                    |
# |                       +---------|D|--------+                       |
# |                 a1 ---|                    |--- a2                 |
# |                       +-/\/K1\/\--/\/K2\/\-+                       |
# |                                                                    |
# |              . interchange of substances:                          |
# |                                                                    |
# |                 a1 <--|diffusion|--> a2                            |
# |                                                                    |
# |              . orientation of fibers that:                         |
# |                                                                    |
# |                 affects the elastic constant (K)                   |
# |                 re-aligns over time towards higher traction dir    |
# |                                                                    |
# |                        /   ↑                               |       |
# |                 a1 ---/--- | traction  --> time -->  a1 ---|---    |
# |                      /     ↓                               |       |
# |                                                                    |
# |        - Vascularization agents: agents representing the vessels   |
# |                                                                    |
# |              . transport of substances:                            |
# |                                                                    |
# |                 v1 ---|diffusion|--> a1                            |
# |                                                                    |
# +--------------------------------------------------------------------+
# +====================================================================+


from pyflamegpu import *
import sys, random, math
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pathlib
from dataclasses import make_dataclass
import time

sns.set()
start_time = time.time()

# +====================================================================+
# | GLOBAL PARAMETERS                                                  |
# +====================================================================+
# Set whether to run single model or ensemble, agent population size, and simulation steps 
ENSEMBLE = False
ENSEMBLE_RUNS = 0
VISUALISATION = False  # Change to false if pyflamegpu has not been built with visualisation support
DEBUG_PRINTING = False
PAUSE_EVERY_STEP = False  # If True, the visualization stops every step until P is pressed
SAVE_PICKLE = True  # If True, dumps agent and boudary force data into a pickle file for post-processing
SHOW_PLOTS = False  # Show plots at the end of the simulation
SAVE_DATA_TO_FILE = True  # If true, agent data is exported to .vtk file every SAVE_EVERY_N_STEPS steps
SAVE_EVERY_N_STEPS = 20  # Affects both the .vtk files and the Dataframes storing boundary data

CURR_PATH = pathlib.Path().absolute()
RES_PATH = CURR_PATH / 'result_files'
RES_PATH.mkdir(parents=True, exist_ok=True)
EPSILON = 0.0000000001

print("Executing in ", CURR_PATH)
# Minimum number of agents per direction (x,y,z). 
# If domain is not cubical, N is asigned to the shorter dimension and more agents are added to the longer ones
# +--------------------------------------------------------------------+
N = 16


# Time simulation parameters
# +--------------------------------------------------------------------+
TIME_STEP = 0.025  # time. WARNING: diffusionn and cell migration events might need different scales
STEPS = 5000

# Boundary interactions and mechanical parameters
# +--------------------------------------------------------------------+
ECM_K_ELAST = 2.0  # [N/units/kg]
ECM_D_DUMPING = 0.4  # [N*s/units/kg]
ECM_MASS = 1.0  # [dimensionless to make K and D mass dependent]
ECM_GEL_CONCENTRATION = 1.0  # [dimensionless, 1.0 represents 2.5mg/ml]
#BOUNDARY_COORDS = [0.5, -0.5, 0.5, -0.5, 0.5, -0.5]  # +X,-X,+Y,-Y,+Z,-Z
BOUNDARY_COORDS = [1000.0, -1000.0, 650.0, -650.0, 150.0, -150.0] # microdevice dimensions in um
BOUNDARY_COORDS = [coord / 1000.0 for coord in BOUNDARY_COORDS] # in mm
BOUNDARY_DISP_RATES = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]  # perpendicular to each surface (+X,-X,+Y,-Y,+Z,-Z) [units/time]
#BOUNDARY_DISP_RATES_PARALLEL = [0.0, 0.0, 0.0, 0.0, 0.0025, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] # parallel to each surface (+X_y,+X_z,-X_y,-X_z,+Y_x,+Y_z,-Y_x,-Y_z,+Z_x,+Z_y,-Z_x,-Z_y)[units/time]
BOUNDARY_DISP_RATES_PARALLEL = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]  # parallel to each surface (+X_y,+X_z,-X_y,-X_z,+Y_x,+Y_z,-Y_x,-Y_z,+Z_x,+Z_y,-Z_x,-Z_y)[units/time]

POISSON_DIRS = [0, 1]  # 0: xdir, 1:ydir, 2:zdir. poisson_ratio ~= -incL(dir1)/incL(dir2) dir2 is the direction in which the load is applied
ALLOW_BOUNDARY_ELASTIC_MOVEMENT = [0, 0, 0, 0, 0, 0]  # [bool]
RELATIVE_BOUNDARY_STIFFNESS = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
BOUNDARY_STIFFNESS_VALUE = 10.0  # N/units
BOUNDARY_DUMPING_VALUE = 5.0
BOUNDARY_STIFFNESS = [BOUNDARY_STIFFNESS_VALUE * x for x in RELATIVE_BOUNDARY_STIFFNESS]
BOUNDARY_DUMPING = [BOUNDARY_DUMPING_VALUE * x for x in RELATIVE_BOUNDARY_STIFFNESS]
#CLAMP_AGENT_TOUCHING_BOUNDARY = [0, 0, 1, 1, 0, 0]  # +X,-X,+Y,-Y,+Z,-Z [bool] - shear assay
CLAMP_AGENT_TOUCHING_BOUNDARY = [1, 1, 1, 1, 1, 1]  # +X,-X,+Y,-Y,+Z,-Z [bool]
ALLOW_AGENT_SLIDING = [0, 0, 0, 0, 0, 0]  # +X,-X,+Y,-Y,+Z,-Z [bool]

# Adjusting number of agents if domain is not cubical
# +--------------------------------------------------------------------+
# Calculate the differences between opposite pairs along each axis
diff_x = abs(BOUNDARY_COORDS[0] - BOUNDARY_COORDS[1])
diff_y = abs(BOUNDARY_COORDS[2] - BOUNDARY_COORDS[3])
diff_z = abs(BOUNDARY_COORDS[4] - BOUNDARY_COORDS[5])

# Check if the differences are equal
if diff_x == diff_y == diff_z:
    ECM_AGENTS_PER_DIR = [N, N, N]
    print("The domain is cubical.")    
else:
    print("The domain is not cubical.")
    min_length = min(diff_x, diff_y, diff_z)
    dist_agents = min_length / (N - 1)
    ECM_AGENTS_PER_DIR = [int(diff_x / dist_agents) + 1, int(diff_y / dist_agents) + 1, int(diff_z / dist_agents) + 1]
    # Redefine BOUNDARY_COORDS due to rounding values
    diff_x = dist_agents * (ECM_AGENTS_PER_DIR[0] - 1)
    diff_y = dist_agents * (ECM_AGENTS_PER_DIR[1] - 1)
    diff_z = dist_agents * (ECM_AGENTS_PER_DIR[2] - 1)
    BOUNDARY_COORDS = [round(diff_x / 2, 2), -round(diff_x / 2, 2), round(diff_y / 2, 2), -round(diff_y / 2, 2), round(diff_z / 2, 2), -round(diff_z / 2, 2)] 
    
print('DOMAIN SIZE: {0},{1},{2}'.format(
    abs(BOUNDARY_COORDS[0] - BOUNDARY_COORDS[1]),
    abs(BOUNDARY_COORDS[3] - BOUNDARY_COORDS[2]),
    abs(BOUNDARY_COORDS[5] - BOUNDARY_COORDS[4])
))
print('ECM_AGENTS_PER_DIR: {0},{1},{2}'.format(ECM_AGENTS_PER_DIR[0], ECM_AGENTS_PER_DIR[1], ECM_AGENTS_PER_DIR[2]))
ECM_POPULATION_SIZE = ECM_AGENTS_PER_DIR[0] * ECM_AGENTS_PER_DIR[1] * ECM_AGENTS_PER_DIR[2]
print('Total number of agents: {0}'.format(ECM_POPULATION_SIZE))
ECM_ECM_EQUILIBRIUM_DISTANCE = (BOUNDARY_COORDS[0] - BOUNDARY_COORDS[1]) / (ECM_AGENTS_PER_DIR[0] - 1) 
print("ECM_ECM_EQUILIBRIUM_DISTANCE [units]: ", ECM_ECM_EQUILIBRIUM_DISTANCE)
ECM_BOUNDARY_INTERACTION_RADIUS = 0.05
ECM_BOUNDARY_EQUILIBRIUM_DISTANCE = 0.0

INCLUDE_FIBER_ALIGNMENT = True
ECM_ORIENTATION_RATE = 0.1 * 0.1 / (
            ECM_ECM_EQUILIBRIUM_DISTANCE * ECM_K_ELAST)  # [1/time] This is adjusted to aprox 1.0/max_force so that the reorientation is not too slow with small forces
print("ECM_ORIENTATION_RATE [1/s]: ", ECM_ORIENTATION_RATE)
MAX_SEARCH_RADIUS_VASCULARIZATION = ECM_ECM_EQUILIBRIUM_DISTANCE  # this strongly affects the number of bins and therefore the memory allocated for simulations (more bins -> more memory -> faster (in theory))
MAX_SEARCH_RADIUS_CELLS = 2 * ECM_ECM_EQUILIBRIUM_DISTANCE # this radius is used to find other agents
MAX_SEARCH_RADIUS_CELL_CELL_INTERACTION = ECM_ECM_EQUILIBRIUM_DISTANCE # this radius is used to check if cells interact with each other
print("MAX_SEARCH_RADIUS for VASCULARIZATION [units]: ", MAX_SEARCH_RADIUS_VASCULARIZATION)
print("MAX_SEARCH_RADIUS for CELLS [units]: ", MAX_SEARCH_RADIUS_CELLS)
OSCILLATORY_SHEAR_ASSAY = False  # if true, BOUNDARY_DISP_RATES_PARALLEL options are overrun but used to make the boundaries oscillate in their corresponding planes following a sin() function
OSCILLATORY_AMPLITUDE = 0.25  # range [0-1]
OSCILLATORY_FREQ = 0.1  # strain oscillation frequency [time^-1]
OSCILLATORY_W = 2 * math.pi * OSCILLATORY_FREQ * TIME_STEP

# Fitting parameters for the fiber strain-stiffening phenomena
# Ref: https://bio.physik.fau.de/publications/Steinwachs%20Nat%20Meth%202016.pdf
# +--------------------------------------------------------------------+
BUCKLING_COEFF_D0 = 0.1
STRAIN_STIFFENING_COEFF_DS = 0.25
CRITICAL_STRAIN = 0.1

# Parallel disp rate values are overrun in oscillatory assays
# +--------------------------------------------------------------------+
if OSCILLATORY_SHEAR_ASSAY:
    for d in range(12):
        if abs(BOUNDARY_DISP_RATES_PARALLEL[d]) > 0.0:
            BOUNDARY_DISP_RATES_PARALLEL[d] = OSCILLATORY_AMPLITUDE * math.cos(
                OSCILLATORY_W * 0.0) * OSCILLATORY_W / TIME_STEP  # cos(w*t)*w is used because the slope of the sin(w*t) function is needed. Expressed in units/sec

# Diffusion related paramenters
# +--------------------------------------------------------------------+
INCLUDE_DIFFUSION = False
N_SPECIES = 2  # number of diffusing species.WARNING: make sure that the value coincides with the one declared in ecm_output_grid_location_data.cpp, ecm_boundary_concentration_conditions.cpp, ecm_ecm_interaction_grid3D.cpp
DIFFUSION_COEFF_MULTI = [0.02, 0.02]  # diffusion coefficient in [units^2/s] per specie
BOUNDARY_CONC_INIT_MULTI = [[-1.0, -1.0, -1.0, -1.0, -1.0, -1.0],
                            # initial concentration at each surface (+X,-X,+Y,-Y,+Z,-Z) [units^2/s]. -1.0 means no condition assigned. All agents are assigned 0 by default.
                            [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0]]  # add as many lines as different species

BOUNDARY_CONC_FIXED_MULTI = [[-1.0, -1.0, -1.0, -1.0, -1.0, -1.0],
                             # concentration boundary conditions at each surface. WARNING: -1.0 means initial condition prevails. Don't use 0.0 as initial condition if that value is not fixed. Use -1.0 instead
                             [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0]]  # add as many lines as different species

INIT_ECM_CONCENTRATION_VALS = [0.0, 0.0]  # initial concentration of each species on the ECM agents
INCLUDE_VASCULARIZATION = False  # if True, vascularization is taken into account
INIT_VASCULARIZATION_CONCENTRATION_VALS = [-1.0,
                                           -1.0]  # initial concentration of each species on the VASCULARIZATION agents
N_VASCULARIZATION_POINTS = 0  # declared here. Points are loaded from file
VASCULARIZATION_POINTS_COORDS = None  # declared here. Coords loaded from file

# Cell agent related paramenters
# +--------------------------------------------------------------------+
INCLUDE_CELLS = True
INCLUDE_CELL_ORIENTATION = True
INCLUDE_CELL_CELL_INTERACTION = False
INCLUDE_CELL_CYCLE = False
PERIODIC_BOUNDARIES_FOR_CELLS = False
CELL_ORIENTATION_RATE = 1.0 * ECM_ORIENTATION_RATE # [1/time] TODO: check whether cell reorient themselves faster than ECM
N_CELLS = 1000
CELL_K_ELAST = 2.0  # [N/units/kg]
CELL_D_DUMPING = 0.4  # [N*time/units/kg]
CELL_RADIUS = ECM_ECM_EQUILIBRIUM_DISTANCE / 2 # [units]
CELL_SPEED_REF = ECM_ECM_EQUILIBRIUM_DISTANCE / TIME_STEP / 10.0 # [units/time]
CYCLE_PHASE_G1_DURATION = 10.0 #[h]
CYCLE_PHASE_S_DURATION = 8.0
CYCLE_PHASE_G2_DURATION = 4.0
CYCLE_PHASE_M_DURATION = 2.0
CYCLE_PHASE_G1_START = 0.0 #[h]
CYCLE_PHASE_S_START = CYCLE_PHASE_G1_DURATION
CYCLE_PHASE_G2_START = CYCLE_PHASE_G1_DURATION + CYCLE_PHASE_S_DURATION
CYCLE_PHASE_M_START = CYCLE_PHASE_G1_DURATION + CYCLE_PHASE_S_DURATION + CYCLE_PHASE_G2_DURATION
CELL_CYCLE_DURATION = CYCLE_PHASE_G1_DURATION + CYCLE_PHASE_S_DURATION + CYCLE_PHASE_G2_DURATION + CYCLE_PHASE_M_DURATION # typically 24h [h]

# Other simulation parameters: TODO: INCLUDE PARALLEL DISP RATES
# +--------------------------------------------------------------------+
MAX_EXPECTED_BOUNDARY_POS = max(BOUNDARY_DISP_RATES) * STEPS * TIME_STEP + max(diff_x, diff_y, diff_z) / 2
MIN_EXPECTED_BOUNDARY_POS = min(BOUNDARY_DISP_RATES) * STEPS * TIME_STEP - max(diff_x, diff_y, diff_z) / 2
print("Max expected boundary position: ", MAX_EXPECTED_BOUNDARY_POS)
print("Min expected boundary position: ", MIN_EXPECTED_BOUNDARY_POS)

# Dataframe initialization data storage
# +--------------------------------------------------------------------+
BPOS = make_dataclass("BPOS", [("xpos", float), ("xneg", float), ("ypos", float), ("yneg", float), ("zpos", float),
                               ("zneg", float)])
# Use a dataframe to store boundary positions over time
BPOS_OVER_TIME = pd.DataFrame([BPOS(BOUNDARY_COORDS[0], BOUNDARY_COORDS[1], BOUNDARY_COORDS[2], BOUNDARY_COORDS[3],
                                    BOUNDARY_COORDS[4], BOUNDARY_COORDS[5])])
OSOT = make_dataclass("OSOT", [("strain", float)])
OSCILLATORY_STRAIN_OVER_TIME = pd.DataFrame([OSOT(0)])

# Checking for incompatible conditions
# +--------------------------------------------------------------------+
critical_error = False
msg_poisson = "WARNING: poisson ratio directions are not well defined or might not have sense due to boundary conditions \n"
if (BOUNDARY_DISP_RATES[0] != 0.0 or BOUNDARY_DISP_RATES[1] != 0.0) and POISSON_DIRS[1] != 0:
    print(msg_poisson)
if (BOUNDARY_DISP_RATES[2] != 0.0 or BOUNDARY_DISP_RATES[3] != 0.0) and POISSON_DIRS[1] != 1:
    print(msg_poisson)
if (BOUNDARY_DISP_RATES[4] != 0.0 or BOUNDARY_DISP_RATES[5] != 0.0) and POISSON_DIRS[1] != 2:
    print(msg_poisson)

msg_incompatible_conditions = "ERROR: CLAMP_AGENT_TOUCHING_BOUNDARY condition is incompatible with ALLOW_BOUNDARY_ELASTIC_MOVEMENT in position [{}]"
for i in range(6):
    if CLAMP_AGENT_TOUCHING_BOUNDARY[i] > 0 and ALLOW_BOUNDARY_ELASTIC_MOVEMENT[i] > 0:
        print(msg_incompatible_conditions.format(i))
        critical_error = True

if INCLUDE_DIFFUSION:
    if (len(DIFFUSION_COEFF_MULTI) != N_SPECIES) or (len(BOUNDARY_CONC_INIT_MULTI) != N_SPECIES) or (
            len(BOUNDARY_CONC_FIXED_MULTI) != N_SPECIES):
        print('ERROR: you must define a diffusion coefficient and the boundary conditions for each species simulated')
        critical_error = True
    # Check diffusion values for numerical stability
    dx = 1.0 / (ECM_AGENTS_PER_DIR[0] - 1)
    for i in range(N_SPECIES):
        Fi_x = 3 * (DIFFUSION_COEFF_MULTI[i] * TIME_STEP / (dx * dx))  # this value should be < 0.5
        print('Fi_x value: {0} for species {1}'.format(Fi_x, i + 1))
        if Fi_x > 0.5:
            print(
                'ERROR: diffusion problem is ill conditioned (Fi_x should be < 0.5), check parameters and consider decreasing time step')
            critical_error = True
    dy = 1.0 / (ECM_AGENTS_PER_DIR[1] - 1)
    for i in range(N_SPECIES):
        Fi_y = 3 * (DIFFUSION_COEFF_MULTI[i] * TIME_STEP / (dy * dy))  # this value should be < 0.5
        print('Fi_y value: {0} for species {1}'.format(Fi_y, i + 1))
        if Fi_y > 0.5:
            print(
                'ERROR: diffusion problem is ill conditioned (Fi_y should be < 0.5), check parameters and consider decreasing time step')
            critical_error = True
    dz = 1.0 / (ECM_AGENTS_PER_DIR[2] - 1)
    for i in range(N_SPECIES):
        Fi_z = 3 * (DIFFUSION_COEFF_MULTI[i] * TIME_STEP / (dz * dz))  # this value should be < 0.5
        print('Fi_z value: {0} for species {1}'.format(Fi_z, i + 1))
        if Fi_z > 0.5:
            print(
                'ERROR: diffusion problem is ill conditioned (Fi_z should be < 0.5), check parameters and consider decreasing time step')
            critical_error = True
elif INCLUDE_VASCULARIZATION:
    print('ERROR: Diffusion is deactivated. Vascularization cannot be included')
    critical_error = True

if INCLUDE_CELLS:
    if MAX_SEARCH_RADIUS_CELL_CELL_INTERACTION < (2 * CELL_RADIUS):
        print('MAX_SEARCH_RADIUS_CELL_CELL_INTERACTION: {0} must be higher than 2 * CELL_RADIUS: 2 * {1}'.format(MAX_SEARCH_RADIUS_CELL_CELL_INTERACTION, CELL_RADIUS))
        critical_error = True

if critical_error:
    quit()
else:
    if SAVE_PICKLE:
        from helper_module import ModelParameterConfig

        # Store all parameters in a class to be stored in the pickle file
        model_config = ModelParameterConfig(SAVE_EVERY_N_STEPS, ECM_AGENTS_PER_DIR, TIME_STEP, STEPS,
                                            ECM_K_ELAST, ECM_D_DUMPING, ECM_MASS, ECM_GEL_CONCENTRATION,
                                            BOUNDARY_COORDS, BOUNDARY_DISP_RATES, BOUNDARY_DISP_RATES_PARALLEL,
                                            POISSON_DIRS,
                                            ALLOW_BOUNDARY_ELASTIC_MOVEMENT,
                                            BOUNDARY_STIFFNESS, BOUNDARY_DUMPING,
                                            CLAMP_AGENT_TOUCHING_BOUNDARY, ALLOW_AGENT_SLIDING,
                                            ECM_ECM_EQUILIBRIUM_DISTANCE, ECM_BOUNDARY_INTERACTION_RADIUS,
                                            ECM_BOUNDARY_EQUILIBRIUM_DISTANCE,
                                            INCLUDE_FIBER_ALIGNMENT, ECM_ORIENTATION_RATE,
                                            BUCKLING_COEFF_D0, STRAIN_STIFFENING_COEFF_DS, CRITICAL_STRAIN,
                                            OSCILLATORY_SHEAR_ASSAY, OSCILLATORY_AMPLITUDE, OSCILLATORY_FREQ,
                                            OSCILLATORY_W,
                                            INCLUDE_DIFFUSION, N_SPECIES, DIFFUSION_COEFF_MULTI,
                                            BOUNDARY_CONC_INIT_MULTI, BOUNDARY_CONC_FIXED_MULTI,
                                            INIT_ECM_CONCENTRATION_VALS,
                                            INCLUDE_VASCULARIZATION, INIT_VASCULARIZATION_CONCENTRATION_VALS)

# +====================================================================+
# | FLAMEGPU2 IMPLEMENTATION                                           |
# +====================================================================+

"""
  FLAME GPU 2 implementation of mechanical assays via moving boundaries and extracellular matrix (ecm) agents by using spatial3D messaging.
"""

# Files containing agent functions for agents, which outputs publicly visible properties to a message list

"""
  VASCULARIZATION  
"""
vascularization_output_location_data_file = "vascularization_output_location_data.cpp"
vascularization_move_file = "vascularization_move.cpp"

"""
  BCORNER  
"""
bcorner_output_location_data_file = "bcorner_output_location_data.cpp"
bcorner_move_file = "bcorner_move.cpp"

"""
  CELL  
"""
cell_output_location_data_file = "cell_output_location_data.cpp"
cell_move_file = "cell_move.cpp"
cell_cycle_file = "cell_cycle.cpp"

"""
  ECM
"""
ecm_output_grid_location_data_file = "ecm_output_grid_location_data.cpp"
ecm_boundary_interaction_file = "ecm_boundary_interaction.cpp"
ecm_ecm_interaction_file = "ecm_ecm_interaction_grid3D.cpp"
ecm_vascularization_interaction_file = "ecm_vascularization_interaction.cpp"
ecm_cell_interaction_file = "ecm_cell_interaction.cpp"
cell_ecm_interaction_file = "cell_ecm_interaction.cpp"
cell_cell_interaction_file = "cell_cell_interaction.cpp"
ecm_boundary_concentration_conditions_file = "ecm_boundary_concentration_conditions.cpp"

"""
  ecm_move agent function for ECM agents
"""
ecm_move_file = "ecm_move.cpp"

model = pyflamegpu.ModelDescription("ECM_Moving_Boundaries")

"""
  GLOBALS
"""
env = model.Environment()
# Population size to generate, if no agents are loaded from disk
env.newPropertyUInt("ECM_POPULATION_TO_GENERATE", ECM_POPULATION_SIZE)
env.newPropertyUInt("CURRENT_ID", 0)
env.newPropertyArrayUInt("ECM_AGENTS_PER_DIR", ECM_AGENTS_PER_DIR)

# Number of steps to simulate
env.newPropertyUInt("STEPS", STEPS)
# Time increment 
env.newPropertyFloat("DELTA_TIME", TIME_STEP)
# Diffusion coefficient
env.newPropertyUInt("INCLUDE_DIFFUSION", INCLUDE_DIFFUSION)
env.newPropertyArrayFloat("DIFFUSION_COEFF_MULTI", DIFFUSION_COEFF_MULTI)
# Number of diffusing species
env.newPropertyUInt("N_SPECIES", N_SPECIES)
env.newPropertyUInt("N_VASCULARIZATION_POINTS", N_VASCULARIZATION_POINTS)

# ------------------------------------------------------
# ECM BEHAVIOUR 
# ------------------------------------------------------
# Equilibrium radius at which elastic force is 0. 
# If ECM_ECM_INTERACTION_RADIUS > ECM_ECM_EQUILIBRIUM_DISTANCE: both repulsion/atraction can occur
# If ECM_ECM_INTERACTION_RADIUS <= ECM_ECM_EQUILIBRIUM_DISTANCE: only repulsion can occur
env.newPropertyFloat("ECM_ECM_EQUILIBRIUM_DISTANCE", ECM_ECM_EQUILIBRIUM_DISTANCE)
# Mechanical parameters
env.newPropertyFloat("ECM_K_ELAST", ECM_K_ELAST)  # initial K_ELAST for agents
env.newPropertyFloat("ECM_D_DUMPING", ECM_D_DUMPING)
env.newPropertyFloat("ECM_MASS", ECM_MASS)
env.newPropertyUInt("INCLUDE_FIBER_ALIGNMENT", INCLUDE_FIBER_ALIGNMENT)
env.newPropertyFloat("ECM_ORIENTATION_RATE", ECM_ORIENTATION_RATE)
env.newPropertyFloat("ECM_GEL_CONCENTRATION", ECM_GEL_CONCENTRATION)
env.newPropertyFloat("BUCKLING_COEFF_D0", BUCKLING_COEFF_D0)
env.newPropertyFloat("STRAIN_STIFFENING_COEFF_DS", STRAIN_STIFFENING_COEFF_DS)
env.newPropertyFloat("CRITICAL_STRAIN", CRITICAL_STRAIN)

# ------------------------------------------------------
# BOUNDARY BEHAVIOUR 
# ------------------------------------------------------
# Boundaries position
bcs = [BOUNDARY_COORDS[0], BOUNDARY_COORDS[1], BOUNDARY_COORDS[2], BOUNDARY_COORDS[3], BOUNDARY_COORDS[4],
       BOUNDARY_COORDS[5]]  # +X,-X,+Y,-Y,+Z,-Z
env.newPropertyArrayFloat("COORDS_BOUNDARIES", bcs)
env.newPropertyArrayFloat("INIT_COORDS_BOUNDARIES",
                          bcs)  # this is used to compute elastic forces with respect to initial position

# Boundaries displacement rate (units/time). 
# e.g. DISP_BOUNDARY_X_POS = 0.1 means that this boundary moves 0.1 units per time towards +X
env.newPropertyArrayFloat("DISP_RATES_BOUNDARIES", BOUNDARY_DISP_RATES)
env.newPropertyArrayFloat("DISP_RATES_BOUNDARIES_PARALLEL", BOUNDARY_DISP_RATES_PARALLEL)

# Boundary-Agent behaviour
env.newPropertyArrayUInt("CLAMP_AGENT_TOUCHING_BOUNDARY", CLAMP_AGENT_TOUCHING_BOUNDARY)
env.newPropertyArrayUInt("ALLOW_BOUNDARY_ELASTIC_MOVEMENT", ALLOW_BOUNDARY_ELASTIC_MOVEMENT)
env.newPropertyArrayFloat("BOUNDARY_STIFFNESS", BOUNDARY_STIFFNESS)
env.newPropertyArrayFloat("BOUNDARY_DUMPING", BOUNDARY_DUMPING)
env.newPropertyArrayUInt("ALLOW_AGENT_SLIDING", ALLOW_AGENT_SLIDING)
env.newPropertyFloat("ECM_BOUNDARY_INTERACTION_RADIUS", ECM_BOUNDARY_INTERACTION_RADIUS)
env.newPropertyFloat("ECM_BOUNDARY_EQUILIBRIUM_DISTANCE", ECM_BOUNDARY_EQUILIBRIUM_DISTANCE)

# Boundary diffusion behavior for multiple species. WARNING: as they are Macro properties, need to be initialized in a host function
env.newMacroPropertyFloat("BOUNDARY_CONC_INIT_MULTI", N_SPECIES,
                          6)  # a 2D matrix with the 6 boundary conditions (columns) for each species (rows)
env.newMacroPropertyFloat("BOUNDARY_CONC_FIXED_MULTI", N_SPECIES,
                          6)  # a 2D matrix with the 6 boundary conditions (columns) for each species (rows)

# Cell properties
env.newPropertyUInt("INCLUDE_CELL_ORIENTATION", INCLUDE_CELL_ORIENTATION)
env.newPropertyUInt("INCLUDE_CELL_CELL_INTERACTION", INCLUDE_CELL_CELL_INTERACTION)
env.newPropertyUInt("PERIODIC_BOUNDARIES_FOR_CELLS", PERIODIC_BOUNDARIES_FOR_CELLS)
env.newPropertyUInt("N_CELLS", N_CELLS)
env.newPropertyFloat("CELL_K_ELAST", CELL_K_ELAST)
env.newPropertyFloat("CELL_D_DUMPING", CELL_D_DUMPING)
env.newPropertyFloat("CELL_RADIUS", CELL_RADIUS)
env.newPropertyFloat("CELL_SPEED_REF", CELL_SPEED_REF)
env.newPropertyFloat("CELL_ORIENTATION_RATE", CELL_ORIENTATION_RATE)
env.newPropertyFloat("MAX_SEARCH_RADIUS_CELLS", MAX_SEARCH_RADIUS_CELLS)
env.newPropertyFloat("MAX_SEARCH_RADIUS_CELL_CELL_INTERACTION", MAX_SEARCH_RADIUS_CELL_CELL_INTERACTION)
env.newPropertyFloat("CELL_CYCLE_DURATION", CELL_CYCLE_DURATION)
env.newPropertyFloat("CYCLE_PHASE_G1_DURATION", CYCLE_PHASE_G1_DURATION)
env.newPropertyFloat("CYCLE_PHASE_S_DURATION", CYCLE_PHASE_S_DURATION)
env.newPropertyFloat("CYCLE_PHASE_G2_DURATION", CYCLE_PHASE_G2_DURATION)
env.newPropertyFloat("CYCLE_PHASE_M_DURATION", CYCLE_PHASE_M_DURATION)
env.newPropertyFloat("CYCLE_PHASE_G1_START", CYCLE_PHASE_G1_START)
env.newPropertyFloat("CYCLE_PHASE_S_START", CYCLE_PHASE_S_START)
env.newPropertyFloat("CYCLE_PHASE_G2_START", CYCLE_PHASE_G2_START)
env.newPropertyFloat("CYCLE_PHASE_M_START", CYCLE_PHASE_M_START)

# Other globals
env.newPropertyFloat("PI", 3.1415)
env.newPropertyUInt("DEBUG_PRINTING", DEBUG_PRINTING)
env.newPropertyFloat("EPSILON", EPSILON)

"""
  LOCATION MESSAGES
"""
bcorner_location_message = model.newMessageSpatial3D("bcorner_location_message")
# Set the range and bounds.
bcorner_location_message.setRadius(MAX_EXPECTED_BOUNDARY_POS - MIN_EXPECTED_BOUNDARY_POS)  # corners are not actually interacting with anything
bcorner_location_message.setMin(MIN_EXPECTED_BOUNDARY_POS, MIN_EXPECTED_BOUNDARY_POS, MIN_EXPECTED_BOUNDARY_POS)
bcorner_location_message.setMax(MAX_EXPECTED_BOUNDARY_POS, MAX_EXPECTED_BOUNDARY_POS, MAX_EXPECTED_BOUNDARY_POS)
# A message to hold the location of an agent. WARNING: spatial3D messages already define x,y,z variables internally.
bcorner_location_message.newVariableInt("id")

vascularization_location_message = model.newMessageSpatial3D("vascularization_location_message")
# Set the range and bounds.
vascularization_location_message.setRadius(
    MAX_SEARCH_RADIUS_VASCULARIZATION)  # TODO: PROPERLY DEFINE THE VALUE OF THE RADIUS
vascularization_location_message.setMin(MIN_EXPECTED_BOUNDARY_POS, MIN_EXPECTED_BOUNDARY_POS,
                                        MIN_EXPECTED_BOUNDARY_POS)
vascularization_location_message.setMax(MAX_EXPECTED_BOUNDARY_POS, MAX_EXPECTED_BOUNDARY_POS,
                                        MAX_EXPECTED_BOUNDARY_POS)
vascularization_location_message.newVariableInt("id")
vascularization_location_message.newVariableFloat("vx")
vascularization_location_message.newVariableFloat("vy")
vascularization_location_message.newVariableFloat("vz")
vascularization_location_message.newVariableArrayFloat("concentration_multi", N_SPECIES)

cell_location_message = model.newMessageSpatial3D("cell_location_message")
# Set the range and bounds.
cell_location_message.setRadius(MAX_SEARCH_RADIUS_CELLS)  # TODO: PROPERLY DEFINE THE VALUE OF THE RADIUS
cell_location_message.setMin(MIN_EXPECTED_BOUNDARY_POS, MIN_EXPECTED_BOUNDARY_POS, MIN_EXPECTED_BOUNDARY_POS)
cell_location_message.setMax(MAX_EXPECTED_BOUNDARY_POS, MAX_EXPECTED_BOUNDARY_POS, MAX_EXPECTED_BOUNDARY_POS)
cell_location_message.newVariableInt("id")
cell_location_message.newVariableFloat("fmag")
cell_location_message.newVariableFloat("k_elast")
cell_location_message.newVariableFloat("d_dumping")
cell_location_message.newVariableFloat("vx")
cell_location_message.newVariableFloat("vy")
cell_location_message.newVariableFloat("vz")
cell_location_message.newVariableFloat("orx")
cell_location_message.newVariableFloat("ory")
cell_location_message.newVariableFloat("orz")
cell_location_message.newVariableFloat("alignment")
cell_location_message.newVariableFloat("radius")
cell_location_message.newVariableFloat("cycle_phase")
cell_location_message.newVariableFloat("clock")
cell_location_message.newVariableInt("completed_cycles")

ecm_grid_location_message = model.newMessageArray3D("ecm_grid_location_message")
ecm_grid_location_message.setDimensions(ECM_AGENTS_PER_DIR[0], ECM_AGENTS_PER_DIR[1], ECM_AGENTS_PER_DIR[2])
ecm_grid_location_message.newVariableInt("id")
ecm_grid_location_message.newVariableFloat("x")
ecm_grid_location_message.newVariableFloat("y")
ecm_grid_location_message.newVariableFloat("z")
ecm_grid_location_message.newVariableFloat("vx")
ecm_grid_location_message.newVariableFloat("vy")
ecm_grid_location_message.newVariableFloat("vz")
ecm_grid_location_message.newVariableFloat("orx")
ecm_grid_location_message.newVariableFloat("ory")
ecm_grid_location_message.newVariableFloat("orz")
ecm_grid_location_message.newVariableFloat("alignment")
ecm_grid_location_message.newVariableFloat("gel_conc")
ecm_grid_location_message.newVariableFloat("k_elast")
ecm_grid_location_message.newVariableUInt8("grid_i")
ecm_grid_location_message.newVariableUInt8("grid_j")
ecm_grid_location_message.newVariableUInt8("grid_k")
ecm_grid_location_message.newVariableArrayFloat("concentration_multi", N_SPECIES)

"""
  AGENTS
"""

"""
  Vascularization agent
"""
if INCLUDE_VASCULARIZATION:
    vascularization_agent = model.newAgent("VASCULARIZATION")
    vascularization_agent.newVariableInt("id")
    vascularization_agent.newVariableFloat("x")
    vascularization_agent.newVariableFloat("y")
    vascularization_agent.newVariableFloat("z")
    vascularization_agent.newVariableFloat("vx")
    vascularization_agent.newVariableFloat("vy")
    vascularization_agent.newVariableFloat("vz")
    vascularization_agent.newVariableArrayFloat("concentration_multi", N_SPECIES)
    vascularization_agent.newRTCFunctionFile("vascularization_output_location_data",
                                             vascularization_output_location_data_file).setMessageOutput(
        "vascularization_location_message")
    vascularization_agent.newRTCFunctionFile("vascularization_move", vascularization_move_file).setMessageInput(
        "ecm_grid_location_message")

"""
  Boundary corner agent
"""
bcorner_agent = model.newAgent("BCORNER")
bcorner_agent.newVariableInt("id")
bcorner_agent.newVariableFloat("x")
bcorner_agent.newVariableFloat("y")
bcorner_agent.newVariableFloat("z")

bcorner_agent.newRTCFunctionFile("bcorner_output_location_data", bcorner_output_location_data_file).setMessageOutput(
    "bcorner_location_message")
bcorner_agent.newRTCFunctionFile("bcorner_move", bcorner_move_file)

"""
  Cell agent
"""
if INCLUDE_CELLS:
    cell_agent = model.newAgent("CELL")
    cell_agent.newVariableInt("id")
    cell_agent.newVariableFloat("x")
    cell_agent.newVariableFloat("y")
    cell_agent.newVariableFloat("z")
    cell_agent.newVariableFloat("vx", 0.0)
    cell_agent.newVariableFloat("vy", 0.0)
    cell_agent.newVariableFloat("vz", 0.0)
    cell_agent.newVariableFloat("fx", 0.0)
    cell_agent.newVariableFloat("fy", 0.0)
    cell_agent.newVariableFloat("fz", 0.0)
    cell_agent.newVariableFloat("f_extension", 0.0)
    cell_agent.newVariableFloat("f_compression", 0.0)
    cell_agent.newVariableFloat("elastic_energy", 0.0)
    cell_agent.newVariableFloat("fmag", 0.0)
    cell_agent.newVariableFloat("k_elast")
    cell_agent.newVariableFloat("d_dumping")
    cell_agent.newVariableFloat("orx")
    cell_agent.newVariableFloat("ory")
    cell_agent.newVariableFloat("orz")
    cell_agent.newVariableFloat("alignment", 0.0)
    cell_agent.newVariableFloat("radius", CELL_RADIUS)
    cell_agent.newVariableInt("cycle_phase", 1) # [1:G1] [2:S] [3:G2] [4:M]
    cell_agent.newVariableFloat("clock", 0.0) # internal clock of the cell to switch phases
    cell_agent.newVariableInt("completed_cycles", 0)
    
    cell_agent.newRTCFunctionFile("cell_output_location_data", cell_output_location_data_file).setMessageOutput(
        "cell_location_message")
    cell_agent.newRTCFunctionFile("cell_move", cell_move_file)
    if INCLUDE_CELL_ORIENTATION:
        cell_agent.newRTCFunctionFile("cell_ecm_interaction", cell_ecm_interaction_file).setMessageInput("ecm_grid_location_message")
    if INCLUDE_CELL_CELL_INTERACTION:
        cell_agent.newRTCFunctionFile("cell_cell_interaction", cell_cell_interaction_file).setMessageInput(
            "cell_location_message")
    if INCLUDE_CELL_CYCLE:
        t = cell_agent.newRTCFunctionFile("cell_cycle", cell_cycle_file)
        t.setAgentOutput(cell_agent)
    # t.setAllowAgentDeath(true) 

"""
  ECM agent
"""
ecm_agent = model.newAgent("ECM")
ecm_agent.newVariableInt("id")
ecm_agent.newVariableFloat("x")
ecm_agent.newVariableFloat("y")
ecm_agent.newVariableFloat("z")
ecm_agent.newVariableFloat("vx")
ecm_agent.newVariableFloat("vy")
ecm_agent.newVariableFloat("vz")
ecm_agent.newVariableFloat("fx")
ecm_agent.newVariableFloat("fy")
ecm_agent.newVariableFloat("fz")
ecm_agent.newVariableFloat("orx")  # orientation of the fibers represented by the agent
ecm_agent.newVariableFloat("ory")
ecm_agent.newVariableFloat("orz")
ecm_agent.newVariableFloat("alignment")
ecm_agent.newVariableFloat("gel_conc")
ecm_agent.newVariableFloat("k_elast")
ecm_agent.newVariableFloat("d_dumping")
ecm_agent.newVariableFloat("mass")
ecm_agent.newVariableFloat(
    "boundary_fx")  # boundary_f[A]: normal force coming from boundary [A] when elastic boundaries option is selected.
ecm_agent.newVariableFloat("boundary_fy")
ecm_agent.newVariableFloat("boundary_fz")
ecm_agent.newVariableFloat(
    "f_bx_pos")  # f_b[A]_[B]: normal force transmitted to the boundary [A]_[B] when agent is clamped
ecm_agent.newVariableFloat("f_bx_neg")
ecm_agent.newVariableFloat("f_by_pos")
ecm_agent.newVariableFloat("f_by_neg")
ecm_agent.newVariableFloat("f_bz_pos")
ecm_agent.newVariableFloat("f_bz_neg")
ecm_agent.newVariableFloat(
    "f_bx_pos_y")  # f_b[A]_[B]_[C]: shear force transmitted to the boundary [A]_[B] in the direction [C] when agent is clamped
ecm_agent.newVariableFloat("f_bx_pos_z")
ecm_agent.newVariableFloat("f_bx_neg_y")
ecm_agent.newVariableFloat("f_bx_neg_z")
ecm_agent.newVariableFloat("f_by_pos_x")
ecm_agent.newVariableFloat("f_by_pos_z")
ecm_agent.newVariableFloat("f_by_neg_x")
ecm_agent.newVariableFloat("f_by_neg_z")
ecm_agent.newVariableFloat("f_bz_pos_x")
ecm_agent.newVariableFloat("f_bz_pos_y")
ecm_agent.newVariableFloat("f_bz_neg_x")
ecm_agent.newVariableFloat("f_bz_neg_y")
ecm_agent.newVariableFloat("f_extension")
ecm_agent.newVariableFloat("f_compression")
ecm_agent.newVariableFloat("elastic_energy")
ecm_agent.newVariableArrayFloat("concentration_multi", N_SPECIES)
ecm_agent.newVariableUInt8("clamped_bx_pos")
ecm_agent.newVariableUInt8("clamped_bx_neg")
ecm_agent.newVariableUInt8("clamped_by_pos")
ecm_agent.newVariableUInt8("clamped_by_neg")
ecm_agent.newVariableUInt8("clamped_bz_pos")
ecm_agent.newVariableUInt8("clamped_bz_neg")
ecm_agent.newVariableUInt8("grid_i")
ecm_agent.newVariableUInt8("grid_j")
ecm_agent.newVariableUInt8("grid_k")

ecm_agent.newRTCFunctionFile("ecm_output_grid_location_data", ecm_output_grid_location_data_file).setMessageOutput(
    "ecm_grid_location_message")
ecm_agent.newRTCFunctionFile("ecm_boundary_interaction", ecm_boundary_interaction_file)
ecm_agent.newRTCFunctionFile("ecm_ecm_interaction", ecm_ecm_interaction_file).setMessageInput(
    "ecm_grid_location_message")
if INCLUDE_VASCULARIZATION:
    ecm_agent.newRTCFunctionFile("ecm_vascularization_interaction",
                                 ecm_vascularization_interaction_file).setMessageInput(
        "vascularization_location_message")
if INCLUDE_DIFFUSION:
    ecm_agent.newRTCFunctionFile("ecm_boundary_concentration_conditions", ecm_boundary_concentration_conditions_file)
ecm_agent.newRTCFunctionFile("ecm_move", ecm_move_file)
if INCLUDE_CELLS:
    ecm_agent.newRTCFunctionFile("ecm_cell_interaction", ecm_cell_interaction_file).setMessageInput(
        "cell_location_message")

""" 
  Helper functions 
"""


def getRandomCoords3D(n, minx, maxx, miny, maxy, minz, maxz):
    """
    Generates an array (nx3 matrix) of random numbers with specific ranges for each column.

    Args:
        n (int): Number of rows in the array.
        minx, maxx (float): Range for the values in the first column [minx, maxx].
        miny, maxy (float): Range for the values in the second column [miny, maxy].
        minz, maxz (float): Range for the values in the third column [minz, maxz].

    Returns:
        numpy.ndarray: Array of random numbers with shape (n, 3).
    """
    np.random.seed()
    random_array = np.random.uniform(low=[minx, miny, minz], high=[maxx, maxy, maxz], size=(n, 3))
    return random_array
    

def randomVector3D():
    """
    Generates a random 3D unit vector (direction) with a uniform spherical distribution

    Returns
    -------
    (x,y,z) : tuple
        Coordinates of the vector.
    """
    np.random.seed()
    phi = np.random.uniform(0.0, np.pi * 2.0)
    costheta = np.random.uniform(-1.0, 1.0)
    theta = np.arccos(costheta)
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)
    return (x, y, z)


def getRandomVectors3D(n_vectors: int):
    """
    Generates an array of random 3D unit vectors (directions) with a uniform spherical distribution

    Parameters
    ----------
    n_vectors : int
        Number of vectors to be generated
    Returns
    -------
    v_array : Numpy array
        Coordinates of the vectors. Shape: [n_vectors, 3].
    """
    v_array = np.zeros((n_vectors, 3))
    for i in range(n_vectors):
        vi = randomVector3D()
        v_array[i, :] = np.array(vi, dtype='float')

    return v_array


def getFixedVectors3D(n_vectors: int, v_dir: np.array):
    """
    Generates an array of 3D unit vectors (directions) in the specified direction

    Parameters
    ----------
    n_vectors : int
        Number of vectors to be generated
    v_dir : Numpy array
        Direction of the vectors
    Returns
    -------
    v_array : Numpy array
        Coordinates of the vectors. Shape: [n_vectors, 3].
    """
    v_array = np.tile(v_dir, (n_vectors, 1))

    return v_array
    
    
def getRandomCoordsAroundPoint(n, px, py, pz, radius):
    """
    Generates N random 3D coordinates within a sphere of a specific radius around a central point.

    Parameters
    ----------
    n : int
        The number of random coordinates to generate.
    px : float
        The x-coordinate of the central point.
    py : float
        The y-coordinate of the central point.
    pz : float
        The z-coordinate of the central point.
    radius : float
        The radius of the sphere.

    Returns
    -------
    coords
        A numpy array of randomly generated 3D coordinates with shape (n, 3).
    """
    central_point = np.array([px, py, pz])
    rand_dirs = getRandomVectors3D(n)
    coords = np.zeros((n, 3))
    np.random.seed()
    for i in range(n):
        radius_i = np.random.uniform(0.0, 1.0) * radius        
        coords[i, :] = central_point + np.array(rand_dirs[i, :] * radius_i, dtype='float')
    

    return coords


"""
  Population initialisation functions
"""


# This class is used to ensure that corner agents are assigned the first 8 ids
class initAgentPopulations(pyflamegpu.HostFunction):
    def run(self, FLAMEGPU):
        global INIT_ECM_CONCENTRATION_VALS, N_SPECIES, INCLUDE_DIFFUSION, INCLUDE_VASCULARIZATION
        global N_VASCULARIZATION_POINTS, VASCULARIZATION_POINTS_COORDS, INCLUDE_FIBER_ALIGNMENT, INCLUDE_CELLS, N_CELLS
        # BOUNDARY CORNERS
        current_id = FLAMEGPU.environment.getPropertyUInt("CURRENT_ID")
        coord_boundary = FLAMEGPU.environment.getPropertyArrayFloat("COORDS_BOUNDARIES")
        coord_boundary_x_pos = coord_boundary[0]
        coord_boundary_x_neg = coord_boundary[1]
        coord_boundary_y_pos = coord_boundary[2]
        coord_boundary_y_neg = coord_boundary[3]
        coord_boundary_z_pos = coord_boundary[4]
        coord_boundary_z_neg = coord_boundary[5]
        print("CORNERS:")
        print("current_id:", current_id)

        for i in range(1, 9):
            instance = FLAMEGPU.agent("BCORNER").newAgent()
            instance.setVariableInt("id", current_id + i)
            if i == 1:
                # +x,+y,+z
                instance.setVariableFloat("x", coord_boundary_x_pos)
                instance.setVariableFloat("y", coord_boundary_y_pos)
                instance.setVariableFloat("z", coord_boundary_z_pos)
            elif i == 2:
                # -x,+y,+z
                instance.setVariableFloat("x", coord_boundary_x_neg)
                instance.setVariableFloat("y", coord_boundary_y_pos)
                instance.setVariableFloat("z", coord_boundary_z_pos)
            elif i == 3:
                # -x,-y,+z
                instance.setVariableFloat("x", coord_boundary_x_neg)
                instance.setVariableFloat("y", coord_boundary_y_neg)
                instance.setVariableFloat("z", coord_boundary_z_pos)
            elif i == 4:
                # +x,-y,+z
                instance.setVariableFloat("x", coord_boundary_x_pos)
                instance.setVariableFloat("y", coord_boundary_y_neg)
                instance.setVariableFloat("z", coord_boundary_z_pos)
            elif i == 5:
                # +x,+y,-z
                instance.setVariableFloat("x", coord_boundary_x_pos)
                instance.setVariableFloat("y", coord_boundary_y_pos)
                instance.setVariableFloat("z", coord_boundary_z_neg)
            elif i == 6:
                # -x,+y,-z
                instance.setVariableFloat("x", coord_boundary_x_neg)
                instance.setVariableFloat("y", coord_boundary_y_pos)
                instance.setVariableFloat("z", coord_boundary_z_neg)
            elif i == 7:
                # -x,-y,-z
                instance.setVariableFloat("x", coord_boundary_x_neg)
                instance.setVariableFloat("y", coord_boundary_y_neg)
                instance.setVariableFloat("z", coord_boundary_z_neg)
            elif i == 8:
                # +x,-y,-z
                instance.setVariableFloat("x", coord_boundary_x_pos)
                instance.setVariableFloat("y", coord_boundary_y_neg)
                instance.setVariableFloat("z", coord_boundary_z_neg)
            else:
                sys.exit("Bad initialization of boundary corners!")

        FLAMEGPU.environment.setPropertyUInt("CURRENT_ID", 8)

        # ECM
        populationSize = FLAMEGPU.environment.getPropertyUInt("ECM_POPULATION_TO_GENERATE")
        min_pos = -1.0
        max_pos = 1.0
        min_speed = 0.0
        max_speed = 0.0
        k_elast = FLAMEGPU.environment.getPropertyFloat("ECM_K_ELAST")
        d_dumping = FLAMEGPU.environment.getPropertyFloat("ECM_D_DUMPING")
        mass = FLAMEGPU.environment.getPropertyFloat("ECM_MASS")
        current_id = FLAMEGPU.environment.getPropertyUInt("CURRENT_ID")
        gel_conc = FLAMEGPU.environment.getPropertyFloat("ECM_GEL_CONCENTRATION")
        current_id += 1
        print("ECM:")
        print("current_id:", current_id)
        agents_per_dir = FLAMEGPU.environment.getPropertyArrayUInt("ECM_AGENTS_PER_DIR")
        print("agents per dir", agents_per_dir)
        offset = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]  # +X,-X,+Y,-Y,+Z,-Z
        coords_x = np.linspace(coord_boundary[1] + offset[1], coord_boundary[0] - offset[0], agents_per_dir[0])
        coords_y = np.linspace(coord_boundary[3] + offset[3], coord_boundary[2] - offset[2], agents_per_dir[1])
        coords_z = np.linspace(coord_boundary[5] + offset[5], coord_boundary[4] - offset[4], agents_per_dir[2])
        if INCLUDE_FIBER_ALIGNMENT:
            orientations = getRandomVectors3D(populationSize)
        else:
            fixed_dir = np.zeros((1, 3))
            fixed_dir[0, 1] = 1.0
            orientations = getFixedVectors3D(populationSize, fixed_dir)
        count = -1
        i = -1
        j = -1
        k = -1

        for x in coords_x:
            i += 1
            j = -1
            for y in coords_y:
                j += 1
                k = -1
                for z in coords_z:
                    k += 1
                    count += 1
                    instance = FLAMEGPU.agent("ECM").newAgent()
                    instance.setVariableInt("id", current_id + count)
                    instance.setVariableFloat("x", x)
                    instance.setVariableFloat("y", y)
                    instance.setVariableFloat("z", z)
                    instance.setVariableFloat("vx", 0.0)
                    instance.setVariableFloat("vy", 0.0)
                    instance.setVariableFloat("vz", 0.0)
                    instance.setVariableFloat("fx", 0.0)
                    instance.setVariableFloat("fy", 0.0)
                    instance.setVariableFloat("fz", 0.0)
                    instance.setVariableFloat("orx", orientations[count, 0])
                    instance.setVariableFloat("ory", orientations[count, 1])
                    instance.setVariableFloat("orz", orientations[count, 2])
                    instance.setVariableFloat("alignment", 0.0)
                    instance.setVariableFloat("gel_conc", gel_conc)
                    instance.setVariableFloat("k_elast", k_elast)
                    instance.setVariableFloat("d_dumping", d_dumping)
                    instance.setVariableFloat("mass", mass)
                    instance.setVariableFloat("boundary_fx", 0.0)
                    instance.setVariableFloat("boundary_fy", 0.0)
                    instance.setVariableFloat("boundary_fz", 0.0)
                    instance.setVariableFloat("f_bx_pos", 0.0)
                    instance.setVariableFloat("f_bx_neg", 0.0)
                    instance.setVariableFloat("f_by_pos", 0.0)
                    instance.setVariableFloat("f_by_neg", 0.0)
                    instance.setVariableFloat("f_bz_pos", 0.0)
                    instance.setVariableFloat("f_bz_neg", 0.0)
                    instance.setVariableFloat("f_bx_pos_y", 0.0)
                    instance.setVariableFloat("f_bx_pos_z", 0.0)
                    instance.setVariableFloat("f_bx_neg_y", 0.0)
                    instance.setVariableFloat("f_bx_neg_z", 0.0)
                    instance.setVariableFloat("f_by_pos_x", 0.0)
                    instance.setVariableFloat("f_by_pos_z", 0.0)
                    instance.setVariableFloat("f_by_neg_x", 0.0)
                    instance.setVariableFloat("f_by_neg_z", 0.0)
                    instance.setVariableFloat("f_bz_pos_x", 0.0)
                    instance.setVariableFloat("f_bz_pos_y", 0.0)
                    instance.setVariableFloat("f_bz_neg_x", 0.0)
                    instance.setVariableFloat("f_bz_neg_y", 0.0)
                    instance.setVariableFloat("f_extension", 0.0)
                    instance.setVariableFloat("f_compression", 0.0)
                    instance.setVariableFloat("elastic_energy", 0.0)
                    instance.setVariableArrayFloat("concentration_multi", INIT_ECM_CONCENTRATION_VALS)
                    instance.setVariableUInt8("clamped_bx_pos", 0)
                    instance.setVariableUInt8("clamped_bx_neg", 0)
                    instance.setVariableUInt8("clamped_by_pos", 0)
                    instance.setVariableUInt8("clamped_by_neg", 0)
                    instance.setVariableUInt8("clamped_bz_pos", 0)
                    instance.setVariableUInt8("clamped_bz_neg", 0)
                    instance.setVariableUInt8("grid_i", i)
                    instance.setVariableUInt8("grid_j", j)
                    instance.setVariableUInt8("grid_k", k)

        FLAMEGPU.environment.setPropertyUInt("CURRENT_ID", current_id + count)
        # VASCULARIZATION
        if INCLUDE_VASCULARIZATION:
            current_id = FLAMEGPU.environment.getPropertyUInt("CURRENT_ID")
            current_id += 1
            count = -1
            VASCULARIZATION_POINTS_COORDS = np.genfromtxt('vascularization_points.txt', delimiter=' ')
            # print(VASCULARIZATION_POINTS_COORDS)
            N_VASCULARIZATION_POINTS = len(VASCULARIZATION_POINTS_COORDS)
            for i in range(N_VASCULARIZATION_POINTS):
                count += 1
                instance = FLAMEGPU.agent("VASCULARIZATION").newAgent()
                instance.setVariableInt("id", current_id + count)
                instance.setVariableFloat("x", VASCULARIZATION_POINTS_COORDS[i][0])
                instance.setVariableFloat("y", VASCULARIZATION_POINTS_COORDS[i][1])
                instance.setVariableFloat("z", VASCULARIZATION_POINTS_COORDS[i][2])
                # instance.setVariableFloat("x", -0.25)
                # instance.setVariableFloat("y", VASCULARIZATION_POINTS_COORDS[i][1])
                # instance.setVariableFloat("z", 0.0)
                instance.setVariableFloat("vx", 0.0)
                instance.setVariableFloat("vy", 0.0)
                instance.setVariableFloat("vz", 0.0)
                instance.setVariableArrayFloat("concentration_multi", INIT_VASCULARIZATION_CONCENTRATION_VALS)

            FLAMEGPU.environment.setPropertyUInt("N_VASCULARIZATION_POINTS", N_VASCULARIZATION_POINTS)
            FLAMEGPU.environment.setPropertyUInt("CURRENT_ID", current_id + count)

        # CELLS
        if INCLUDE_CELLS:
            current_id = FLAMEGPU.environment.getPropertyUInt("CURRENT_ID")
            current_id += 1
            count = -1
            cell_orientations = getRandomVectors3D(N_CELLS)
            # cell_orientations = np.array([[0.0, 1.0, 1.0], [1.0, 1.0, 1.0]], dtype='float')
            cell_pos = getRandomCoords3D(N_CELLS,
                                         coord_boundary[0], coord_boundary[1],
                                         coord_boundary[2], coord_boundary[3],
                                         coord_boundary[4], coord_boundary[5])
            #cell_pos = getRandomCoordsAroundPoint(N_CELLS, 0.0, 0.0, 0.0, 0.25)
            k_elast = FLAMEGPU.environment.getPropertyFloat("CELL_K_ELAST")
            d_dumping = FLAMEGPU.environment.getPropertyFloat("CELL_D_DUMPING")
            radius = FLAMEGPU.environment.getPropertyFloat("CELL_RADIUS")
            for i in range(N_CELLS):
                count += 1
                instance = FLAMEGPU.agent("CELL").newAgent()
                instance.setVariableInt("id", current_id + count)
                instance.setVariableFloat("x", cell_pos[i, 0])
                instance.setVariableFloat("y", cell_pos[i, 1])
                instance.setVariableFloat("z", cell_pos[i, 2])
                instance.setVariableFloat("fmag", 0.0)
                instance.setVariableFloat("fx", 0.0)
                instance.setVariableFloat("fy", 0.0)
                instance.setVariableFloat("fz", 0.0)
                instance.setVariableFloat("f_extension", 0.0)
                instance.setVariableFloat("f_compression", 0.0)
                instance.setVariableFloat("elastic_energy", 0.0)
                instance.setVariableFloat("k_elast", k_elast)
                instance.setVariableFloat("d_dumping", d_dumping)
                instance.setVariableFloat("orx", cell_orientations[count, 0])
                instance.setVariableFloat("ory", cell_orientations[count, 1])
                instance.setVariableFloat("orz", cell_orientations[count, 2])
                instance.setVariableFloat("alignment", 0.0)
                instance.setVariableFloat("vx", 0.0)
                instance.setVariableFloat("vy", 0.0)
                instance.setVariableFloat("vz", 0.0)
                instance.setVariableFloat("radius", radius)
                cycle_phase = random.randint(1, 4) # [1:G1] [2:S] [3:G2] [4:M]
                instance.setVariableInt("cycle_phase", cycle_phase)
                cycle_clock = 0.0
                if cycle_phase == 1:
                    cycle_clock = FLAMEGPU.environment.getPropertyFloat("CYCLE_PHASE_G1_START") 
                    + np.random.uniform(0.0, 1.0) * FLAMEGPU.environment.getPropertyFloat("CYCLE_PHASE_G1_DURATION")                
                elif cycle_phase == 2:
                    cycle_clock = FLAMEGPU.environment.getPropertyFloat("CYCLE_PHASE_S_START")
                    + np.random.uniform(0.0, 1.0) * FLAMEGPU.environment.getPropertyFloat("CYCLE_PHASE_S_DURATION")                    
                elif cycle_phase == 3:
                    cycle_clock = FLAMEGPU.environment.getPropertyFloat("CYCLE_PHASE_G2_START")
                    + np.random.uniform(0.0, 1.0) * FLAMEGPU.environment.getPropertyFloat("CYCLE_PHASE_G2_DURATION")                    
                elif cycle_phase == 4:
                    cycle_clock = FLAMEGPU.environment.getPropertyFloat("CYCLE_PHASE_M_START")
                    + np.random.uniform(0.0, 1.0) * FLAMEGPU.environment.getPropertyFloat("CYCLE_PHASE_M_DURATION")                    
                instance.setVariableFloat("clock", cycle_clock)
                instance.setVariableInt("completed_cycles",0)

            FLAMEGPU.environment.setPropertyUInt("CURRENT_ID", current_id + count)

        return


def resetMacroProperties(self, FLAMEGPU):
    global BOUNDARY_CONC_INIT_MULTI, BOUNDARY_CONC_FIXED_MULTI
    bcim = FLAMEGPU.environment.getMacroPropertyFloat("BOUNDARY_CONC_INIT_MULTI")
    bcfm = FLAMEGPU.environment.getMacroPropertyFloat("BOUNDARY_CONC_FIXED_MULTI")
    for i in range(len(BOUNDARY_CONC_INIT_MULTI)):
        for j in range(len(BOUNDARY_CONC_INIT_MULTI[i])):
            bcim[i][j] = BOUNDARY_CONC_INIT_MULTI[i][j]
    for i in range(len(BOUNDARY_CONC_FIXED_MULTI)):
        for j in range(len(BOUNDARY_CONC_FIXED_MULTI[i])):
            bcfm[i][j] = BOUNDARY_CONC_FIXED_MULTI[i][j]
    print("Reseting MacroProperties")
    print(BOUNDARY_CONC_INIT_MULTI)
    print(BOUNDARY_CONC_FIXED_MULTI)
    return


# This class is used to reset the MacroProperties to the values stored in the global variables
class initMacroProperties(pyflamegpu.HostFunction):
    def run(self, FLAMEGPU):
        resetMacroProperties(self, FLAMEGPU)
        return


# Add function callback to INIT functions for population generation
initialAgentPopulation = initAgentPopulations()
initialMacroProperties = initMacroProperties()
model.addInitFunction(initialAgentPopulation)
model.addInitFunction(initialMacroProperties)

"""
  STEP FUNCTIONS
"""


class MoveBoundaries(pyflamegpu.HostFunction):
    """
     pyflamegpu requires step functions to be a class which extends the StepFunction base class.
     This class must extend the handle function
     """

    # Define Python class 'constructor'
    def __init__(self):
        super().__init__()
        self.apply_parallel_disp = list()
        for d in range(12):
            if abs(BOUNDARY_DISP_RATES_PARALLEL[d]) > 0.0:
                self.apply_parallel_disp.append(True)
            else:
                self.apply_parallel_disp.append(False)

    # Override C++ method: virtual void run(FLAMEGPU_HOST_API*)
    def run(self, FLAMEGPU):
        stepCounter = FLAMEGPU.getStepCounter() + 1
        global BOUNDARY_DISP_RATES, ALLOW_BOUNDARY_ELASTIC_MOVEMENT, BOUNDARY_STIFFNESS, BOUNDARY_DUMPING, BPOS_OVER_TIME
        global CLAMP_AGENT_TOUCHING_BOUNDARY, OSCILLATORY_SHEAR_ASSAY, OSCILLATORY_AMPLITUDE, OSCILLATORY_W, OSCILLATORY_STRAIN_OVER_TIME
        global DEBUG_PRINTING, PAUSE_EVERY_STEP, TIME_STEP

        boundaries_moved = False
        if PAUSE_EVERY_STEP:
            input()  # pause everystep
    
        coord_boundary = list(FLAMEGPU.environment.getPropertyArrayFloat("COORDS_BOUNDARIES"))
        if OSCILLATORY_SHEAR_ASSAY:
            if stepCounter % SAVE_EVERY_N_STEPS == 0 or stepCounter == 1:
                new_val = pd.DataFrame([OSOT(OSCILLATORY_AMPLITUDE * math.sin(OSCILLATORY_W * stepCounter))])
                # OSCILLATORY_STRAIN_OVER_TIME = OSCILLATORY_STRAIN_OVER_TIME.append(new_val, ignore_index=True) #TODO: FIX?
                OSCILLATORY_STRAIN_OVER_TIME = pd.concat([OSCILLATORY_STRAIN_OVER_TIME, new_val], ignore_index=True)
            for d in range(12):
                if self.apply_parallel_disp[d]:
                    BOUNDARY_DISP_RATES_PARALLEL[d] = OSCILLATORY_AMPLITUDE * math.cos(
                        OSCILLATORY_W * stepCounter) * OSCILLATORY_W / TIME_STEP  # cos(w*t)*t is used because the slope of the sin(w*t) function is needed

            FLAMEGPU.environment.setPropertyArrayFloat("DISP_RATES_BOUNDARIES_PARALLEL", BOUNDARY_DISP_RATES_PARALLEL)

        if any(catb < 1 for catb in CLAMP_AGENT_TOUCHING_BOUNDARY) or any(
                abem > 0 for abem in ALLOW_BOUNDARY_ELASTIC_MOVEMENT):
            boundaries_moved = True
            agent = FLAMEGPU.agent("ECM")
            minmax_positions = list()
            minmax_positions.append(agent.maxFloat("x"))
            minmax_positions.append(agent.minFloat("x"))
            minmax_positions.append(agent.maxFloat("y"))
            minmax_positions.append(agent.minFloat("y"))
            minmax_positions.append(agent.maxFloat("z"))
            minmax_positions.append(agent.minFloat("z"))
            boundary_equil_distances = list()
            boundary_equil_distances.append(ECM_BOUNDARY_EQUILIBRIUM_DISTANCE)
            boundary_equil_distances.append(-ECM_BOUNDARY_EQUILIBRIUM_DISTANCE)
            boundary_equil_distances.append(ECM_BOUNDARY_EQUILIBRIUM_DISTANCE)
            boundary_equil_distances.append(-ECM_BOUNDARY_EQUILIBRIUM_DISTANCE)
            boundary_equil_distances.append(ECM_BOUNDARY_EQUILIBRIUM_DISTANCE)
            boundary_equil_distances.append(-ECM_BOUNDARY_EQUILIBRIUM_DISTANCE)
            for i in range(6):
                if CLAMP_AGENT_TOUCHING_BOUNDARY[i] < 1:
                    if ALLOW_BOUNDARY_ELASTIC_MOVEMENT[i] > 0:
                        coord_boundary[i] = minmax_positions[i] + boundary_equil_distances[i]
                    else:
                        coord_boundary[i] = minmax_positions[i]

            bcs = [coord_boundary[0], coord_boundary[1], coord_boundary[2], coord_boundary[3], coord_boundary[4],
                   coord_boundary[5]]  # +X,-X,+Y,-Y,+Z,-Z
            FLAMEGPU.environment.setPropertyArrayFloat("COORDS_BOUNDARIES", bcs)

            if stepCounter % SAVE_EVERY_N_STEPS == 0 or stepCounter == 1:
                print("====== MOVING FREE BOUNDARIES  ======")
                print("New boundary positions [+X,-X,+Y,-Y,+Z,-Z]: ", coord_boundary)
                print("=====================================")

        if any(dr > 0.0 or dr < 0.0 for dr in BOUNDARY_DISP_RATES):
            boundaries_moved = True
            for i in range(6):
                coord_boundary[i] += (BOUNDARY_DISP_RATES[i] * TIME_STEP)

            bcs = [coord_boundary[0], coord_boundary[1], coord_boundary[2], coord_boundary[3], coord_boundary[4],
                   coord_boundary[5]]  # +X,-X,+Y,-Y,+Z,-Z
            FLAMEGPU.environment.setPropertyArrayFloat("COORDS_BOUNDARIES", bcs)
            if stepCounter % SAVE_EVERY_N_STEPS == 0 or stepCounter == 1:
                print("====== MOVING BOUNDARIES DUE TO CONDITIONS ======")
                print("New boundary positions [+X,-X,+Y,-Y,+Z,-Z]: ", coord_boundary)
                print("=================================================")

        # if any(abem > 0 for abem in ALLOW_BOUNDARY_ELASTIC_MOVEMENT):
        #   boundaries_moved = True
        #   print ("====== MOVING BOUNDARIES DUE TO FORCES ======")
        #   agent = FLAMEGPU.agent("ECM")
        #   sum_bx_pos = agent.sumFloat("f_bx_pos")
        #   sum_bx_neg = agent.sumFloat("f_bx_neg")
        #   sum_by_pos = agent.sumFloat("f_by_pos")
        #   sum_by_neg = agent.sumFloat("f_by_neg")
        #   sum_bz_pos = agent.sumFloat("f_bz_pos")
        #   sum_bz_neg = agent.sumFloat("f_bz_neg")
        #   print ("Total forces [+X,-X,+Y,-Y,+Z,-Z]: ", sum_bx_pos, sum_bx_neg, sum_by_pos, sum_by_neg, sum_bz_pos, sum_bz_neg)
        #   boundary_forces = [sum_bx_pos, sum_bx_neg, sum_by_pos, sum_by_neg, sum_bz_pos, sum_bz_neg]
        #   for i in range(6):
        #       if BOUNDARY_DISP_RATES[i] < EPSILON and BOUNDARY_DISP_RATES[i] > -EPSILON and ALLOW_BOUNDARY_ELASTIC_MOVEMENT[i]:
        #           #u = boundary_forces[i] / BOUNDARY_STIFFNESS[i]
        #           u = (boundary_forces[i] * TIME_STEP)/ (BOUNDARY_STIFFNESS[i] * TIME_STEP + BOUNDARY_DUMPING[i])
        #           print ("Displacement for boundary {} = {}".format(i,u))
        #           coord_boundary[i] += u

        #   bcs = [coord_boundary[0], coord_boundary[1], coord_boundary[2], coord_boundary[3], coord_boundary[4], coord_boundary[5]]  #+X,-X,+Y,-Y,+Z,-Z
        #   FLAMEGPU.environment.setPropertyArrayFloat("COORDS_BOUNDARIES", bcs)
        #   print ("New boundary positions [+X,-X,+Y,-Y,+Z,-Z]: ", coord_boundary)
        #   print ("=================================================")

        if boundaries_moved:
            if stepCounter % SAVE_EVERY_N_STEPS == 0 or stepCounter == 1:
                new_pos = pd.DataFrame([BPOS(coord_boundary[0], coord_boundary[1], coord_boundary[2],
                                             coord_boundary[3], coord_boundary[4], coord_boundary[5])])
                # BPOS_OVER_TIME = BPOS_OVER_TIME.append(new_pos, ignore_index=True)
                BPOS_OVER_TIME = pd.concat([BPOS_OVER_TIME, new_pos], ignore_index=True)

        # print ("End of step: ", stepCounter)


class SaveDataToFile(pyflamegpu.HostFunction):
    def __init__(self):
        global ECM_AGENTS_PER_DIR, N_VASCULARIZATION_POINTS
        super().__init__()
        self.header = list()
        self.header.append("# vtk DataFile Version 3.0")
        self.header.append("ECM data")
        self.header.append("ASCII")
        self.header.append("DATASET POLYDATA")
        self.header.append("POINTS {} float".format(8 + ECM_AGENTS_PER_DIR[0] * ECM_AGENTS_PER_DIR[1] * ECM_AGENTS_PER_DIR[2]))  # number of ECM agents + 8 corners
        # self.header.append("POINTS {} float".format(8))
        self.domaindata = list()
        self.domaindata.append("POLYGONS 6 30")
        cube_conn = [[4, 0, 3, 7, 4], [4, 1, 2, 6, 5], [4, 1, 0, 4, 5], [4, 2, 3, 7, 6], [4, 0, 1, 2, 3],
                     [4, 4, 5, 6, 7]]
        for i in range(len(cube_conn)):
            for j in range(len(cube_conn[i])):
                if j > 0:
                    cube_conn[i][j] = cube_conn[i][j] + ECM_AGENTS_PER_DIR[0] * ECM_AGENTS_PER_DIR[1] * ECM_AGENTS_PER_DIR[2]
            self.domaindata.append(' '.join(str(x) for x in cube_conn[i]))

        # self.domaindata.append("4 0 3 7 4")
        # self.domaindata.append("4 1 2 6 5")
        # self.domaindata.append("4 1 0 4 5")
        # self.domaindata.append("4 2 3 7 6")
        # self.domaindata.append("4 0 1 2 3")
        # self.domaindata.append("4 4 5 6 7")
        self.domaindata.append("CELL_DATA 6")
        self.domaindata.append("SCALARS boundary_index int 1")
        self.domaindata.append("LOOKUP_TABLE default")
        self.domaindata.append("0")
        self.domaindata.append("1")
        self.domaindata.append("2")
        self.domaindata.append("3")
        self.domaindata.append("4")
        self.domaindata.append("5")
        self.domaindata.append("NORMALS boundary_normals float")
        self.domaindata.append("1 0 0")
        self.domaindata.append("-1 0 0")
        self.domaindata.append("0 1 0")
        self.domaindata.append("0 -1 0")
        self.domaindata.append("0 0 1")
        self.domaindata.append("0 0 -1")
        # VASCULARIZATION
        self.vascularizationdata = list()  # a different file is created to show the position of the vascularization points
        self.vascularizationdata.append("# vtk DataFile Version 3.0")
        self.vascularizationdata.append("Vascularization points")
        self.vascularizationdata.append("ASCII")
        self.vascularizationdata.append("DATASET UNSTRUCTURED_GRID")
        # CELLS
        self.celldata = list()  # a different file is created to show cell agent data
        self.celldata.append("# vtk DataFile Version 3.0")
        self.celldata.append("Cell agents")
        self.celldata.append("ASCII")
        self.celldata.append("DATASET UNSTRUCTURED_GRID")

    def run(self, FLAMEGPU):
        global SAVE_DATA_TO_FILE, SAVE_EVERY_N_STEPS, N_SPECIES
        global RES_PATH, ENSEMBLE
        global fileCounter, INCLUDE_VASCULARIZATION
        global INCLUDE_CELLS
        BUCKLING_COEFF_D0 = FLAMEGPU.environment.getPropertyFloat("BUCKLING_COEFF_D0")
        STRAIN_STIFFENING_COEFF_DS = FLAMEGPU.environment.getPropertyFloat("STRAIN_STIFFENING_COEFF_DS")
        CRITICAL_STRAIN = FLAMEGPU.environment.getPropertyFloat("CRITICAL_STRAIN")
        stepCounter = FLAMEGPU.getStepCounter() + 1
        
        coord_boundary = list(FLAMEGPU.environment.getPropertyArrayFloat("COORDS_BOUNDARIES"))

        if SAVE_DATA_TO_FILE:
            if stepCounter % SAVE_EVERY_N_STEPS == 0 or stepCounter == 1:

                if INCLUDE_VASCULARIZATION:
                    vasc_coords = list()
                    file_name = 'vascularization_points_t{:04d}.vtk'.format(stepCounter)
                    file_path = RES_PATH / file_name
                    vasc_agent = FLAMEGPU.agent("VASCULARIZATION")
                    av = vasc_agent.getPopulationData()  # this returns a DeviceAgentVector
                    for ai in av:
                        coords_ai = (ai.getVariableFloat("x"), ai.getVariableFloat("y"), ai.getVariableFloat("z"))
                        vasc_coords.append(coords_ai)
                    with open(str(file_path), 'w') as file:
                        for line in self.vascularizationdata:
                            file.write(line + '\n')
                        file.write("POINTS {} float \n".format(FLAMEGPU.environment.getPropertyUInt(
                            "N_VASCULARIZATION_POINTS")))  # number of vascularization agents
                        for coords_ai in vasc_coords:
                            file.write("{} {} {} \n".format(coords_ai[0], coords_ai[1], coords_ai[2]))

                if INCLUDE_CELLS:
                    cell_coords = list()
                    cell_velocity = list()
                    cell_orientation = list()
                    cell_alignment = list()
                    cell_radius = list()
                    cell_clock = list()
                    cell_cycle_phase = list()
                    file_name = 'cells_t{:04d}.vtk'.format(stepCounter)
                    file_path = RES_PATH / file_name
                    cell_agent = FLAMEGPU.agent("CELL")
                    cell_agent.sortInt("id", pyflamegpu.HostAgentAPI.Asc); # this is critical to ensure cell ids are kept in order for visualization
                    av = cell_agent.getPopulationData()  # this returns a DeviceAgentVector
                    for ai in av:
                        coords_ai = (ai.getVariableFloat("x"), ai.getVariableFloat("y"), ai.getVariableFloat("z"))
                        velocity_ai = (ai.getVariableFloat("vx"), ai.getVariableFloat("vy"), ai.getVariableFloat("vz"))
                        orientation_ai = (
                        ai.getVariableFloat("orx"), ai.getVariableFloat("ory"), ai.getVariableFloat("orz"))
                        alignment_ai = ai.getVariableFloat("alignment")
                        radius_ai = ai.getVariableFloat("radius")
                        clock_ai = ai.getVariableFloat("clock")
                        cycle_phase_ai = ai.getVariableInt("cycle_phase")
                        cell_coords.append(coords_ai)
                        cell_velocity.append(velocity_ai)
                        cell_orientation.append(orientation_ai)
                        cell_alignment.append(alignment_ai)
                        cell_radius.append(radius_ai)
                        cell_clock.append(clock_ai)
                        cell_cycle_phase.append(cycle_phase_ai)
                    with open(str(file_path), 'w') as file:
                        for line in self.celldata:
                            file.write(line + '\n')
                        file.write("POINTS {} float \n".format(
                            FLAMEGPU.environment.getPropertyUInt("N_CELLS")))  # number of cell agents
                        for coords_ai in cell_coords:
                            file.write("{} {} {} \n".format(coords_ai[0], coords_ai[1], coords_ai[2]))
                        file.write("POINT_DATA {} \n".format(
                            FLAMEGPU.environment.getPropertyUInt("N_CELLS")))  
                        file.write("SCALARS alignment float 1" + '\n')
                        file.write("LOOKUP_TABLE default" + '\n')
                        for a_ai in cell_alignment:
                            file.write("{:.4f} \n".format(a_ai))                            
                        file.write("SCALARS radius float 1" + '\n')
                        file.write("LOOKUP_TABLE default" + '\n')
                        for r_ai in cell_radius:
                            file.write("{:.4f} \n".format(r_ai))                            
                        file.write("SCALARS clock float 1" + '\n')
                        file.write("LOOKUP_TABLE default" + '\n')
                        for c_ai in cell_clock:
                            file.write("{:.4f} \n".format(c_ai))                        
                        file.write("SCALARS cycle_phase int 1" + '\n')
                        file.write("LOOKUP_TABLE default" + '\n')
                        for ccp_ai in cell_cycle_phase:
                            file.write("{} \n".format(ccp_ai))                       
                        file.write("VECTORS velocity float" + '\n')
                        for v_ai in cell_velocity:
                            file.write("{:.4f} {:.4f} {:.4f} \n".format(v_ai[0], v_ai[1], v_ai[2]))
                        file.write("VECTORS orientation float" + '\n')
                        for o_ai in cell_orientation:
                            file.write("{:.4f} {:.4f} {:.4f} \n".format(o_ai[0], o_ai[1], o_ai[2]))

                file_name = 'ecm_data_t{:04d}.vtk'.format(stepCounter)
                if ENSEMBLE:
                    dir_name = f"BUCKLING_COEFF_D0_{BUCKLING_COEFF_D0:.3f}_STRAIN_STIFFENING_COEFF_DS_{STRAIN_STIFFENING_COEFF_DS:.3f}_CRITICAL_STRAIN_{CRITICAL_STRAIN:.3f}"
                    # Combine the base directory with the current directory name
                    file_path = RES_PATH / dir_name / file_name
                else:
                    file_path = RES_PATH / file_name

                agent = FLAMEGPU.agent("ECM")
                # reaction forces, thus, opposite to agent-applied forces
                sum_bx_pos = -agent.sumFloat("f_bx_pos")
                sum_bx_neg = -agent.sumFloat("f_bx_neg")
                sum_by_pos = -agent.sumFloat("f_by_pos")
                sum_by_neg = -agent.sumFloat("f_by_neg")
                sum_bz_pos = -agent.sumFloat("f_bz_pos")
                sum_bz_neg = -agent.sumFloat("f_bz_neg")
                sum_bx_pos_y = -agent.sumFloat("f_bx_pos_y")
                sum_bx_pos_z = -agent.sumFloat("f_bx_pos_z")
                sum_bx_neg_y = -agent.sumFloat("f_bx_neg_y")
                sum_bx_neg_z = -agent.sumFloat("f_bx_neg_z")
                sum_by_pos_x = -agent.sumFloat("f_by_pos_x")
                sum_by_pos_z = -agent.sumFloat("f_by_pos_z")
                sum_by_neg_x = -agent.sumFloat("f_by_neg_x")
                sum_by_neg_z = -agent.sumFloat("f_by_neg_z")
                sum_bz_pos_x = -agent.sumFloat("f_bz_pos_x")
                sum_bz_pos_y = -agent.sumFloat("f_bz_pos_y")
                sum_bz_neg_x = -agent.sumFloat("f_bz_neg_x")
                sum_bz_neg_y = -agent.sumFloat("f_bz_neg_y")

                coords = list()
                velocity = list()
                orientation = list()
                alignment = list()
                gel_conc = list()
                force = list()
                elastic_energy = list()
                concentration_multi = list()  # this is a list of tuples. Each tuple has N_SPECIES elements
                av = agent.getPopulationData()  # this returns a DeviceAgentVector
                for ai in av:
                    coords_ai = (ai.getVariableFloat("x"), ai.getVariableFloat("y"), ai.getVariableFloat("z"))
                    velocity_ai = (ai.getVariableFloat("vx"), ai.getVariableFloat("vy"), ai.getVariableFloat("vz"))
                    force_ai = (ai.getVariableFloat("fx"), ai.getVariableFloat("fy"), ai.getVariableFloat("fz"))
                    orientation_ai = (
                    ai.getVariableFloat("orx"), ai.getVariableFloat("ory"), ai.getVariableFloat("orz"))
                    alignment.append(ai.getVariableFloat("alignment"))
                    gel_conc.append(ai.getVariableFloat("gel_conc"))
                    coords.append(coords_ai)
                    velocity.append(velocity_ai)
                    force.append(force_ai)
                    orientation.append(orientation_ai)
                    elastic_energy.append(ai.getVariableFloat("elastic_energy"))
                    concentration_multi.append(ai.getVariableArrayFloat("concentration_multi"))
                print("====== SAVING DATA FROM Step {:03d} TO FILE ======".format(stepCounter))
                with open(str(file_path), 'w') as file:
                    for line in self.header:
                        file.write(line + '\n')
                    for coords_ai in coords:
                        file.write("{} {} {} \n".format(coords_ai[0], coords_ai[1], coords_ai[2]))
                    # Write boundary positions at the end so that corner points don't cover the points underneath
                    file.write("{} {} {} \n".format(coord_boundary[0], coord_boundary[2], coord_boundary[4]))
                    file.write("{} {} {} \n".format(coord_boundary[1], coord_boundary[2], coord_boundary[4]))
                    file.write("{} {} {} \n".format(coord_boundary[1], coord_boundary[3], coord_boundary[4]))
                    file.write("{} {} {} \n".format(coord_boundary[0], coord_boundary[3], coord_boundary[4]))
                    file.write("{} {} {} \n".format(coord_boundary[0], coord_boundary[2], coord_boundary[5]))
                    file.write("{} {} {} \n".format(coord_boundary[1], coord_boundary[2], coord_boundary[5]))
                    file.write("{} {} {} \n".format(coord_boundary[1], coord_boundary[3], coord_boundary[5]))
                    file.write("{} {} {} \n".format(coord_boundary[0], coord_boundary[3], coord_boundary[5]))
                    for line in self.domaindata:
                        file.write(line + '\n')
                    file.write("SCALARS boundary_normal_forces float 1" + '\n')
                    file.write("LOOKUP_TABLE default" + '\n')
                    file.write(str(sum_bx_pos) + '\n')
                    file.write(str(sum_bx_neg) + '\n')
                    file.write(str(sum_by_pos) + '\n')
                    file.write(str(sum_by_neg) + '\n')
                    file.write(str(sum_bz_pos) + '\n')
                    file.write(str(sum_bz_neg) + '\n')
                    file.write("SCALARS boundary_normal_force_scaling float 1" + '\n')
                    file.write("LOOKUP_TABLE default" + '\n')
                    file.write(str(abs(sum_bx_pos)) + '\n')
                    file.write(str(abs(sum_bx_neg)) + '\n')
                    file.write(str(abs(sum_by_pos)) + '\n')
                    file.write(str(abs(sum_by_neg)) + '\n')
                    file.write(str(abs(sum_bz_pos)) + '\n')
                    file.write(str(abs(sum_bz_neg)) + '\n')
                    file.write("VECTORS boundary_normal_force_dir float" + '\n')
                    file.write("1 0 0 \n" if sum_bx_pos > 0 else "-1 0 0 \n")
                    file.write("1 0 0 \n" if sum_bx_neg > 0 else "-1 0 0 \n")
                    file.write("0 1 0 \n" if sum_by_pos > 0 else "0 -1 0 \n")
                    file.write("0 1 0 \n" if sum_by_neg > 0 else "0 -1 0 \n")
                    file.write("0 0 1 \n" if sum_bz_pos > 0 else "0 0 -1 \n")
                    file.write("0 0 1 \n" if sum_bz_neg > 0 else "0 0 -1 \n")
                    # must be divided in blocks of 6 (one value per face of the cube)
                    file.write("SCALARS boundary_shear_forces_pos float 1" + '\n')
                    file.write("LOOKUP_TABLE default" + '\n')
                    file.write(str(sum_bx_pos_y) + '\n')
                    file.write(str(sum_bx_pos_z) + '\n')
                    file.write(str(sum_by_pos_x) + '\n')
                    file.write(str(sum_by_pos_z) + '\n')
                    file.write(str(sum_bz_pos_x) + '\n')
                    file.write(str(sum_bz_pos_y) + '\n')
                    file.write("SCALARS boundary_shear_forces_neg float 1" + '\n')
                    file.write("LOOKUP_TABLE default" + '\n')
                    file.write(str(sum_bx_neg_y) + '\n')
                    file.write(str(sum_bx_neg_z) + '\n')
                    file.write(str(sum_by_neg_x) + '\n')
                    file.write(str(sum_by_neg_z) + '\n')
                    file.write(str(sum_bz_neg_x) + '\n')
                    file.write(str(sum_bz_neg_y) + '\n')
                    file.write("SCALARS boundary_shear_force_scaling_pos float 1" + '\n')
                    file.write("LOOKUP_TABLE default" + '\n')
                    file.write(str(abs(sum_bx_pos_y)) + '\n')
                    file.write(str(abs(sum_bx_pos_z)) + '\n')
                    file.write(str(abs(sum_by_pos_x)) + '\n')
                    file.write(str(abs(sum_by_pos_z)) + '\n')
                    file.write(str(abs(sum_bz_pos_x)) + '\n')
                    file.write(str(abs(sum_bz_pos_y)) + '\n')
                    file.write("SCALARS boundary_shear_force_scaling_neg float 1" + '\n')
                    file.write("LOOKUP_TABLE default" + '\n')
                    file.write(str(abs(sum_bx_neg_y)) + '\n')
                    file.write(str(abs(sum_bx_neg_z)) + '\n')
                    file.write(str(abs(sum_by_neg_x)) + '\n')
                    file.write(str(abs(sum_by_neg_z)) + '\n')
                    file.write(str(abs(sum_bz_neg_x)) + '\n')
                    file.write(str(abs(sum_bz_neg_y)) + '\n')

                    file.write("VECTORS boundary_shear_force_dir_pos float" + '\n')
                    file.write("0 1 0 \n" if sum_bx_pos_y > 0 else "0 -1 0 \n")
                    file.write("0 0 1 \n" if sum_bx_pos_z > 0 else "0 0 -1 \n")
                    file.write("1 0 0 \n" if sum_by_pos_x > 0 else "-1 0 0 \n")
                    file.write("0 0 1 \n" if sum_by_pos_z > 0 else "0 0 -1 \n")
                    file.write("1 0 0 \n" if sum_bz_pos_x > 0 else "-1 0 0 \n")
                    file.write("0 1 0 \n" if sum_bz_pos_y > 0 else "0 -1 0 \n")
                    file.write("VECTORS boundary_shear_force_dir_neg float" + '\n')
                    file.write("0 1 0 \n" if sum_bx_neg_y > 0 else "0 -1 0 \n")
                    file.write("0 0 1 \n" if sum_bx_neg_z > 0 else "0 0 -1 \n")
                    file.write("1 0 0 \n" if sum_by_neg_x > 0 else "-1 0 0 \n")
                    file.write("0 0 1 \n" if sum_by_neg_z > 0 else "0 0 -1 \n")
                    file.write("1 0 0 \n" if sum_bz_neg_x > 0 else "-1 0 0 \n")
                    file.write("0 1 0 \n" if sum_bz_neg_y > 0 else "0 -1 0 \n")

                    file.write("POINT_DATA {} \n".format(8 + ECM_AGENTS_PER_DIR[0] * ECM_AGENTS_PER_DIR[1] * ECM_AGENTS_PER_DIR[2]))  # 8 corners + number of ECM agents

                    file.write(
                        "SCALARS is_corner int 1" + '\n')  # create this variable to remove them from representations
                    file.write("LOOKUP_TABLE default" + '\n')

                    for ee_ai in elastic_energy:
                        file.write("{0} \n".format(0))
                    for i in range(8):
                        file.write("1 \n")  # boundary corners

                    file.write("SCALARS elastic_energy float 1" + '\n')
                    file.write("LOOKUP_TABLE default" + '\n')
                    for ee_ai in elastic_energy:
                        file.write("{:.4f} \n".format(ee_ai))
                    for i in range(8):
                        file.write("0.0 \n")  # boundary corners

                    file.write("SCALARS alignment float 1" + '\n')
                    file.write("LOOKUP_TABLE default" + '\n')
                    for a_ai in alignment:
                        file.write("{:.4f} \n".format(a_ai))
                    for i in range(8):
                        file.write("0.0 \n")  # boundary corners

                    file.write("SCALARS gel_conc float 1" + '\n')
                    file.write("LOOKUP_TABLE default" + '\n')
                    for gc_ai in gel_conc:
                        file.write("{:.4f} \n".format(gc_ai))
                    for i in range(8):
                        file.write("0.0 \n")  # boundary corners

                    for s in range(N_SPECIES):
                        file.write("SCALARS concentration_species_{0} float 1 \n".format(s))
                        file.write("LOOKUP_TABLE default" + '\n')

                        for c_ai in concentration_multi:
                            file.write("{:.4f} \n".format(c_ai[s]))
                        for i in range(8):
                            file.write("0.0 \n")  # boundary corners

                    file.write("VECTORS velocity float" + '\n')
                    for v_ai in velocity:
                        file.write("{:.4f} {:.4f} {:.4f} \n".format(v_ai[0], v_ai[1], v_ai[2]))
                    for i in range(8):
                        file.write("0.0 0.0 0.0 \n")  # boundary corners

                    file.write("VECTORS force float" + '\n')
                    for f_ai in force:
                        file.write("{:.4f} {:.4f} {:.4f} \n".format(f_ai[0], f_ai[1], f_ai[2]))
                    for i in range(8):
                        file.write("0.0 0.0 0.0 \n")  # boundary corners

                    file.write("VECTORS orientation float" + '\n')
                    for o_ai in orientation:
                        file.write("{:.4f} {:.4f} {:.4f} \n".format(o_ai[0], o_ai[1], o_ai[2]))
                    for i in range(8):
                        file.write("0.0 0.0 0.0 \n")  # boundary corners

                print("... succesful save ")
                print("=================================")


class UpdateBoundaryConcentrationMulti(pyflamegpu.HostFunction):
    def __init__(self):
        super().__init__()

    def run(self, FLAMEGPU):
        global BOUNDARY_CONC_INIT_MULTI, BOUNDARY_CONC_FIXED_MULTI
        stepCounter = FLAMEGPU.getStepCounter() + 1
        if stepCounter == 2:  # after first step BOUNDARY_CONC_INIT_MULTI is removed (set to -1.0) and BOUNDARY_CONC_FIXED_MULTI prevails
            print("====== CONCENTRATION MULTI BOUNDARY CONDITIONS SET  ======")
            print("Initial concentration boundary conditions [+X,-X,+Y,-Y,+Z,-Z]: ", BOUNDARY_CONC_INIT_MULTI)
            print("Fixed concentration boundary conditions [+X,-X,+Y,-Y,+Z,-Z]: ", BOUNDARY_CONC_FIXED_MULTI)
            for i in range(len(BOUNDARY_CONC_INIT_MULTI)):
                for j in range(len(BOUNDARY_CONC_INIT_MULTI[i])):
                    BOUNDARY_CONC_INIT_MULTI[i][j] = -1.0
            resetMacroProperties(self, FLAMEGPU)
            
class UpdateAgentCount(pyflamegpu.HostFunction): # if cells proliferate, N_CELLS must be updated
    def __init__(self):
        super().__init__()

    def run(self, FLAMEGPU):
        FLAMEGPU.environment.setPropertyUInt("N_CELLS", FLAMEGPU.agent("CELL").count())
         
        


if INCLUDE_DIFFUSION:
    ubcm = UpdateBoundaryConcentrationMulti()
    model.addStepFunction(ubcm)

if INCLUDE_CELLS:
    uac = UpdateAgentCount()
    model.addStepFunction(uac)

sdf = SaveDataToFile()
model.addStepFunction(sdf)

mb = MoveBoundaries()
model.addStepFunction(mb)

"""
  END OF STEP FUNCTIONS
"""

"""
  Control flow
"""
# First set of layers: location of agents
layer_count = 1
model.newLayer("L" + str(layer_count)).addAgentFunction("ECM", "ecm_output_grid_location_data")
model.Layer("L" + str(layer_count)).addAgentFunction("BCORNER", "bcorner_output_location_data")
if INCLUDE_VASCULARIZATION:
    model.Layer("L" + str(layer_count)).addAgentFunction("VASCULARIZATION", "vascularization_output_location_data")
if INCLUDE_CELLS:
    model.Layer("L" + str(layer_count)).addAgentFunction("CELL", "cell_output_location_data")

layer_count += 1
# Second set of layers: interactions
model.newLayer("L" + str(layer_count)).addAgentFunction("ECM", "ecm_boundary_interaction")
layer_count += 1
if INCLUDE_DIFFUSION:
    model.newLayer("L" + str(layer_count)).addAgentFunction("ECM", "ecm_boundary_concentration_conditions")
    layer_count += 1
if INCLUDE_VASCULARIZATION:
    model.newLayer("L" + str(layer_count)).addAgentFunction("ECM", "ecm_vascularization_interaction")
    layer_count += 1

model.newLayer("L" + str(layer_count)).addAgentFunction("ECM", "ecm_ecm_interaction")
layer_count += 1

if INCLUDE_CELLS:
    model.newLayer("L" + str(layer_count)).addAgentFunction("ECM", "ecm_cell_interaction")
    if INCLUDE_CELL_ORIENTATION:
        model.Layer("L" + str(layer_count)).addAgentFunction("CELL", "cell_ecm_interaction")
    layer_count += 1
    if INCLUDE_CELL_CELL_INTERACTION:
        model.newLayer("L" + str(layer_count)).addAgentFunction("CELL", "cell_cell_interaction")
        layer_count += 1
    if INCLUDE_CELL_CYCLE:
        model.newLayer("L" + str(layer_count)).addAgentFunction("CELL", "cell_cycle")
        layer_count += 1

# Third set of layers: diffusion from boundaries
if INCLUDE_DIFFUSION:
    model.newLayer("L" + str(layer_count)).addAgentFunction("ECM",
                                                            "ecm_boundary_concentration_conditions")  # called twice to ensure concentration at boundaries is properly shown visually
    layer_count += 1

# Fourth set of layers: agent movement
model.newLayer("L" + str(layer_count)).addAgentFunction("ECM", "ecm_move")
model.Layer("L" + str(layer_count)).addAgentFunction("BCORNER", "bcorner_move")
layer_count += 1
if INCLUDE_VASCULARIZATION:
    model.newLayer("L" + str(layer_count)).addAgentFunction("VASCULARIZATION", "vascularization_move")
    layer_count += 1
if INCLUDE_CELLS:
    model.newLayer("L" + str(layer_count)).addAgentFunction("CELL", "cell_move")
    layer_count += 1

# Create and configure logging details 
logging_config = pyflamegpu.LoggingConfig(model)
logging_config.logEnvironment("CURRENT_ID")
logging_config.logEnvironment("COORDS_BOUNDARIES")
ecm_agent_log = logging_config.agent("ECM")
ecm_agent_log.logCount()
ecm_agent_log.logSumFloat("f_bx_pos")
ecm_agent_log.logSumFloat("f_bx_neg")
ecm_agent_log.logSumFloat("f_by_pos")
ecm_agent_log.logSumFloat("f_by_neg")
ecm_agent_log.logSumFloat("f_bz_pos")
ecm_agent_log.logSumFloat("f_bz_neg")

ecm_agent_log.logSumFloat("f_bx_pos_y")
ecm_agent_log.logSumFloat("f_bx_pos_z")
ecm_agent_log.logSumFloat("f_bx_neg_y")
ecm_agent_log.logSumFloat("f_bx_neg_z")
ecm_agent_log.logSumFloat("f_by_pos_x")
ecm_agent_log.logSumFloat("f_by_pos_z")
ecm_agent_log.logSumFloat("f_by_neg_x")
ecm_agent_log.logSumFloat("f_by_neg_z")
ecm_agent_log.logSumFloat("f_bz_pos_x")
ecm_agent_log.logSumFloat("f_bz_pos_y")
ecm_agent_log.logSumFloat("f_bz_neg_x")
ecm_agent_log.logSumFloat("f_bz_neg_y")

ecm_agent_log.logMeanFloat("f_bx_pos")
ecm_agent_log.logMeanFloat("f_bx_neg")
ecm_agent_log.logMeanFloat("f_by_pos")
ecm_agent_log.logMeanFloat("f_by_neg")
ecm_agent_log.logMeanFloat("f_bz_pos")
ecm_agent_log.logMeanFloat("f_bz_neg")
ecm_agent_log.logStandardDevFloat("f_bx_pos")
ecm_agent_log.logStandardDevFloat("f_bx_neg")
ecm_agent_log.logStandardDevFloat("f_by_pos")
ecm_agent_log.logStandardDevFloat("f_by_neg")
ecm_agent_log.logStandardDevFloat("f_bz_pos")
ecm_agent_log.logStandardDevFloat("f_bz_neg")

step_log = pyflamegpu.StepLoggingConfig(logging_config)
step_log.setFrequency(1)

"""
  Create Model Runner
"""
if ENSEMBLE:
    """
  Create Run Plan Vector
  """
    run_plan_vector = pyflamegpu.RunPlanVector(model, ENSEMBLE_RUNS)
    run_plan_vector.setSteps(env.getPropertyUInt("STEPS"))
    print(env.getPropertyUInt("STEPS"))
    simulation_seed = random.randint(0, 99999)
    run_plan_vector.setRandomSimulationSeed(simulation_seed, 1000)
    simulation = pyflamegpu.CUDAEnsemble(model)

    """
    Create Control Run Plan
  """
    # Create a control run plan, this will define the common properties across all plans
    # https://docs.flamegpu.com/guide/running-multiple-simulations/index.html#creating-a-runplanvector
    run_control = pyflamegpu.RunPlan(model)

    # Ensure that repeated runs use the same Random values within the RunPlans
    # run_control.setRandomPropertySeed(34523) # This method only exists at the vector level, and you're not using setPropertyRandom() so it woud have no effect.
    run_control.setRandomSimulationSeed(
        34523)  # This is possibly what you want, means all simulations will use the same random seed.
    # All runs have the same steps
    run_control.setSteps(STEPS)
    run_control.setPropertyUInt("STEPS", STEPS)

    # Create the first dimension of the parameter sweep
    runs_final = pyflamegpu.RunPlanVector(model, 0)
    for tBUCKLING_COEFF_D0 in np.linspace(0.1, 0.5, 3):
        for tSTRAIN_STIFFENING_COEFF_DS in np.linspace(0.1, 0.75, 5):
            for tCRITICAL_STRAIN in np.linspace(0.0, 0.5, 5):
                run_control.setPropertyFloat("BUCKLING_COEFF_D0", tBUCKLING_COEFF_D0)
                run_control.setPropertyFloat("STRAIN_STIFFENING_COEFF_DS", tSTRAIN_STIFFENING_COEFF_DS)
                run_control.setPropertyFloat("CRITICAL_STRAIN", tCRITICAL_STRAIN)
                runs_final += run_control
                # Create directory names using the parameter values
                dir_name = f"BUCKLING_COEFF_D0_{tBUCKLING_COEFF_D0:.3f}_STRAIN_STIFFENING_COEFF_DS_{tSTRAIN_STIFFENING_COEFF_DS:.3f}_CRITICAL_STRAIN_{tCRITICAL_STRAIN:.3f}"
                # Combine the base directory with the current directory name
                full_path = RES_PATH / dir_name
                # Create the directory if it doesn't exist
                full_path.mkdir(parents=True, exist_ok=True)

    # Create a CUDAEnsemble to execute the RunPlanVector
    ensemble = pyflamegpu.CUDAEnsemble(model)

    # Override config defaults
    ensemble.Config().out_directory = RES_PATH.as_posix()
    ensemble.Config().out_format = "json"
    ensemble.Config().concurrent_runs = 1  # This is concurrent runs per device, higher values may improve performance for "small" models
    ensemble.Config().timing = False
    ensemble.Config().error_level = pyflamegpu.CUDAEnsembleConfig.Fast  # Kills the ensemble as soon as the first error is detected
    # ensemble.Config().devices = pyflamegpu.IntSet([0]) # By default ensemble will use all available CUDA devices

    # Pass any logging configs to the CUDAEnsemble
    # https://docs.flamegpu.com/guide/running-multiple-simulations/index.html#creating-a-logging-configuration
    ensemble.setStepLog(step_log)
    ensemble.setExitLog(logging_config)

else:
    simulation = pyflamegpu.CUDASimulation(model)
    simulation.SimulationConfig().steps = STEPS
    simulation.setStepLog(step_log)
    simulation.setExitLog(logging_config)

"""
  Create Visualisation
"""
if pyflamegpu.VISUALISATION and VISUALISATION and not ENSEMBLE:
    visualisation = simulation.getVisualisation()
    # Configure vis
    envWidth = MAX_EXPECTED_BOUNDARY_POS - MIN_EXPECTED_BOUNDARY_POS
    INIT_CAM = MAX_EXPECTED_BOUNDARY_POS * 4.5
    # Visualisation.setInitialCameraLocation(INIT_CAM * 2, INIT_CAM, INIT_CAM)
    visualisation.setInitialCameraLocation(0.0, 0.0, INIT_CAM)
    visualisation.setCameraSpeed(0.002 * envWidth)
    if DEBUG_PRINTING:
        visualisation.setSimulationSpeed(1)
    visualisation.setBeginPaused(True)
    circ_ecm_agt = visualisation.addAgent("ECM")
    # Position vars are named x, y, z so they are used by default
    circ_ecm_agt.setModel(pyflamegpu.ICOSPHERE)
    # circ_ecm_agt.setModelScale(env.getPropertyFloat("ECM_ECM_INTERACTION_RADIUS")/7.5)
    circ_ecm_agt.setModelScale(0.03 * envWidth)
    circ_ecm_agt.setColor(pyflamegpu.GREEN)
    # circ_ecm_agt.setColor(pyflamegpu.ViridisInterpolation("y", -1.0, 1.0))
    # circ_ecm_agt.setColor(pyflamegpu.HSVInterpolation("y", 0.0, 360.0))
    f_max = ECM_K_ELAST * (ECM_ECM_EQUILIBRIUM_DISTANCE)
    max_energy = 0.5 * (f_max * f_max) / ECM_K_ELAST
    print("max force, max energy: ", f_max, max_energy)
    circ_ecm_agt.setColor(pyflamegpu.HSVInterpolation.GREENRED("elastic_energy", 0.00000001, max_energy * 1.0))
    square_bcorner_agt = visualisation.addAgent("BCORNER")
    square_bcorner_agt.setModel(pyflamegpu.CUBE)
    square_bcorner_agt.setModelScale(0.025 * envWidth)
    square_bcorner_agt.setColor(pyflamegpu.RED)
    if INCLUDE_VASCULARIZATION:
        square_vascularization_agt = visualisation.addAgent("VASCULARIZATION")
        square_vascularization_agt.setModel(pyflamegpu.CUBE)
        square_vascularization_agt.setModelScale(0.018 * envWidth)
        square_vascularization_agt.setColor(pyflamegpu.BLUE)

    if INCLUDE_CELLS:
        circ_cell_agt = visualisation.addAgent("CELL")
        circ_cell_agt.setModel(pyflamegpu.ICOSPHERE)
        circ_cell_agt.setModelScale(0.045 * envWidth)
        circ_cell_agt.setColor(pyflamegpu.Color("#fc03e7"))
        #circ_cell_agt.setColor(pyflamegpu.DiscreteColor("cycle_phase", pyflamegpu.Viridis(10), pyflamegpu.WHITE, 1, 2))
        
    coord_boundary = list(env.getPropertyArrayFloat("COORDS_BOUNDARIES"))
    pen = visualisation.newLineSketch(1, 1, 1, 0.8)
    pen.addVertex(coord_boundary[0], coord_boundary[2], coord_boundary[4])
    pen.addVertex(coord_boundary[0], coord_boundary[2], coord_boundary[5])
    pen.addVertex(coord_boundary[0], coord_boundary[3], coord_boundary[4])
    pen.addVertex(coord_boundary[0], coord_boundary[3], coord_boundary[5])
    pen.addVertex(coord_boundary[1], coord_boundary[2], coord_boundary[4])
    pen.addVertex(coord_boundary[1], coord_boundary[2], coord_boundary[5])
    pen.addVertex(coord_boundary[1], coord_boundary[3], coord_boundary[4])
    pen.addVertex(coord_boundary[1], coord_boundary[3], coord_boundary[5])

    pen.addVertex(coord_boundary[0], coord_boundary[2], coord_boundary[4])
    pen.addVertex(coord_boundary[0], coord_boundary[3], coord_boundary[4])
    pen.addVertex(coord_boundary[0], coord_boundary[2], coord_boundary[5])
    pen.addVertex(coord_boundary[0], coord_boundary[3], coord_boundary[5])
    pen.addVertex(coord_boundary[1], coord_boundary[2], coord_boundary[4])
    pen.addVertex(coord_boundary[1], coord_boundary[3], coord_boundary[4])
    pen.addVertex(coord_boundary[1], coord_boundary[2], coord_boundary[5])
    pen.addVertex(coord_boundary[1], coord_boundary[3], coord_boundary[5])

    pen.addVertex(coord_boundary[0], coord_boundary[2], coord_boundary[4])
    pen.addVertex(coord_boundary[1], coord_boundary[2], coord_boundary[4])
    pen.addVertex(coord_boundary[0], coord_boundary[3], coord_boundary[4])
    pen.addVertex(coord_boundary[1], coord_boundary[3], coord_boundary[4])
    pen.addVertex(coord_boundary[0], coord_boundary[2], coord_boundary[5])
    pen.addVertex(coord_boundary[1], coord_boundary[2], coord_boundary[5])
    pen.addVertex(coord_boundary[0], coord_boundary[3], coord_boundary[5])
    pen.addVertex(coord_boundary[1], coord_boundary[3], coord_boundary[5])

    visualisation.activate()

"""
  Execution
"""
if ENSEMBLE:
    # Execute the ensemble using the specified RunPlans
    errs = ensemble.simulate(runs_final)
else:
    simulation.simulate()

"""
  Export Pop
"""
# simulation.exportData("end.xml")

# Join Visualisation
if pyflamegpu.VISUALISATION and VISUALISATION and not ENSEMBLE:
    visualisation.join()

print("--- EXECUTION TIME: %s seconds ---" % (time.time() - start_time))

incL_dir1 = (BPOS_OVER_TIME.iloc[:, POISSON_DIRS[0] * 2] - BPOS_OVER_TIME.iloc[:, POISSON_DIRS[0] * 2 + 1]) - (
        BPOS_OVER_TIME.iloc[0, POISSON_DIRS[0] * 2] - BPOS_OVER_TIME.iloc[0, POISSON_DIRS[0] * 2 + 1])
incL_dir2 = (BPOS_OVER_TIME.iloc[:, POISSON_DIRS[1] * 2] - BPOS_OVER_TIME.iloc[:, POISSON_DIRS[1] * 2 + 1]) - (
        BPOS_OVER_TIME.iloc[0, POISSON_DIRS[1] * 2] - BPOS_OVER_TIME.iloc[0, POISSON_DIRS[1] * 2 + 1])

print(incL_dir1)
print('/')
print(incL_dir2)

POISSON_RATIO_OVER_TIME = -1 * incL_dir1 / incL_dir2


def manageLogs(steps, is_ensemble, idx):
    global SAVE_EVERY_N_STEPS, SAVE_PICKLE, SHOW_PLOTS, RES_PATH, model_config
    global BPOS_OVER_TIME, BFORCE_OVER_TIME, BFORCE_SHEAR_OVER_TIME, POISSON_RATIO_OVER_TIME, OSCILLATORY_STRAIN_OVER_TIME
    ecm_agent_counts = [None] * len(steps)
    counter = 0
    BFORCE = make_dataclass("BFORCE",
                            [("fxpos", float), ("fxneg", float), ("fypos", float), ("fyneg", float), ("fzpos", float),
                             ("fzneg", float)])
    BFORCE_SHEAR = make_dataclass("BFORCE_SHEAR",
                                  [("fxpos_y", float), ("fxpos_z", float), ("fxneg_y", float), ("fxneg_z", float),
                                   ("fypos_x", float), ("fypos_z", float), ("fyneg_x", float), ("fyneg_z", float),
                                   ("fzpos_x", float), ("fzpos_y", float), ("fzneg_x", float), ("fzneg_y", float)])
    for step in steps:
        stepcount = step.getStepCount()
        if stepcount % SAVE_EVERY_N_STEPS == 0 or stepcount == 1:
            ecm_agents = step.getAgent("ECM")
            ecm_agent_counts[counter] = ecm_agents.getCount()
            f_bx_pos = ecm_agents.getSumFloat("f_bx_pos")
            f_bx_neg = ecm_agents.getSumFloat("f_bx_neg")
            f_by_pos = ecm_agents.getSumFloat("f_by_pos")
            f_by_neg = ecm_agents.getSumFloat("f_by_neg")
            f_bz_pos = ecm_agents.getSumFloat("f_bz_pos")
            f_bz_neg = ecm_agents.getSumFloat("f_bz_neg")
            f_bx_pos_y = ecm_agents.getSumFloat("f_bx_pos_y")
            f_bx_pos_z = ecm_agents.getSumFloat("f_bx_pos_z")
            f_bx_neg_y = ecm_agents.getSumFloat("f_bx_neg_y")
            f_bx_neg_z = ecm_agents.getSumFloat("f_bx_neg_z")
            f_by_pos_x = ecm_agents.getSumFloat("f_by_pos_x")
            f_by_pos_z = ecm_agents.getSumFloat("f_by_pos_z")
            f_by_neg_x = ecm_agents.getSumFloat("f_by_neg_x")
            f_by_neg_z = ecm_agents.getSumFloat("f_by_neg_z")
            f_bz_pos_x = ecm_agents.getSumFloat("f_bz_pos_x")
            f_bz_pos_y = ecm_agents.getSumFloat("f_bz_pos_y")
            f_bz_neg_x = ecm_agents.getSumFloat("f_bz_neg_x")
            f_bz_neg_y = ecm_agents.getSumFloat("f_bz_neg_y")

            step_bforce = pd.DataFrame([BFORCE(f_bx_pos, f_bx_neg, f_by_pos, f_by_neg, f_bz_pos, f_bz_neg)])
            step_bforce_shear = pd.DataFrame([BFORCE_SHEAR(f_bx_pos_y, f_bx_pos_z, f_bx_neg_y, f_bx_neg_z,
                                                           f_by_pos_x, f_by_pos_z, f_by_neg_x, f_by_neg_z,
                                                           f_bz_pos_x, f_bz_pos_y, f_bz_neg_x, f_bz_neg_y)])
            if counter == 0:
                BFORCE_OVER_TIME = pd.DataFrame([BFORCE(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)])
                BFORCE_SHEAR_OVER_TIME = pd.DataFrame(
                    [BFORCE_SHEAR(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)])
            else:
                # BFORCE_OVER_TIME = BFORCE_OVER_TIME.append(step_bforce, ignore_index=True) # deprecated
                BFORCE_OVER_TIME = pd.concat([BFORCE_OVER_TIME, step_bforce], ignore_index=True)
                # BFORCE_SHEAR_OVER_TIME = BFORCE_SHEAR_OVER_TIME.append(step_bforce_shear, ignore_index=True) # deprecated
                BFORCE_SHEAR_OVER_TIME = pd.concat([BFORCE_SHEAR_OVER_TIME, step_bforce_shear], ignore_index=True)
            counter += 1
    if not is_ensemble:
        print()
        print("============================")
        print("BOUNDARY POSITIONS OVER TIME")
        print(BPOS_OVER_TIME)
        print()
        print("============================")
        print("BOUNDARY FORCES OVER TIME")
        print(BFORCE_OVER_TIME)
        print()
        print("============================")
        print("BOUNDARY SHEAR FORCES OVER TIME")
        print(BFORCE_SHEAR_OVER_TIME)
        print()
        print("============================")
        print("POISSON RATIO OVER TIME")
        print(POISSON_RATIO_OVER_TIME)
        print()
        print("============================")
        print("STRAIN OVER TIME")
        print(OSCILLATORY_STRAIN_OVER_TIME)
        print()
    # Saving pickle
    if SAVE_PICKLE:
        import pickle

        file_name = f'output_data_{idx}.pickle'
        file_path = RES_PATH / file_name
        with open(str(file_path), 'wb') as file:
            pickle.dump({'BPOS_OVER_TIME': BPOS_OVER_TIME,
                         'BFORCE_OVER_TIME': BFORCE_OVER_TIME,
                         'BFORCE_SHEAR_OVER_TIME': BFORCE_SHEAR_OVER_TIME,
                         'POISSON_RATIO_OVER_TIME': POISSON_RATIO_OVER_TIME,
                         'OSCILLATORY_STRAIN_OVER_TIME': OSCILLATORY_STRAIN_OVER_TIME,
                         'model_config': model_config},
                        file, protocol=pickle.HIGHEST_PROTOCOL)

            print('Results successfully saved to {0}'.format(file_path))
    # Plotting
    if SHOW_PLOTS and not is_ensemble:
        # fig,ax=plt.subplots(2,3)
        fig = plt.figure()
        gs = fig.add_gridspec(2, 3)
        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[0, 1])
        ax3 = fig.add_subplot(gs[:, 2])
        ax4 = fig.add_subplot(gs[1, 0])
        ax5 = fig.add_subplot(gs[1, 1])
        # BPOS_OVER_TIME.plot()

        BPOS_OVER_TIME.plot(ax=ax1)
        # ax = df['size'].plot(secondary_y=True, color='k', marker='o')
        ax1.set_xlabel('time step')
        ax1.set_ylabel('pos')
        BFORCE_OVER_TIME.plot(ax=ax2)
        ax2.set_ylabel('normal force')
        ax2.set_xlabel('time step')
        BFORCE_SHEAR_OVER_TIME.plot(ax=ax3)
        ax3.set_ylabel('shear force')
        ax3.set_xlabel('time step')
        POISSON_RATIO_OVER_TIME.plot(ax=ax4)
        ax4.set_ylabel('poisson ratio')
        ax4.set_xlabel('time step')
        plt.sca(ax5)
        plt.plot(BPOS_OVER_TIME['ypos'] - 0.5, BFORCE_OVER_TIME['fypos'])
        ax5.set_ylabel('normal force')
        ax5.set_xlabel('disp')

        fig.tight_layout()

        if OSCILLATORY_SHEAR_ASSAY:
            OSCILLATORY_STRAIN_OVER_TIME.plot()
            fig2 = plt.figure()
            colors = np.arange(0, STEPS + 1, 1).tolist()
            plt.scatter(OSCILLATORY_STRAIN_OVER_TIME['strain'].abs(), BFORCE_SHEAR_OVER_TIME['fypos_x'].abs(),
                        marker='o', c=colors, alpha=0.3, cmap='viridis')
            plt.xlabel('strain')
            plt.ylabel('shear force')
            fig3 = plt.figure()
            colors = np.arange(0, STEPS + 1, 1).tolist()
            plt.scatter(OSCILLATORY_STRAIN_OVER_TIME['strain'].abs(), BFORCE_SHEAR_OVER_TIME['fypos_x'], marker='o',
                        c=colors, alpha=0.3, cmap='viridis')
            plt.xlabel('strain')
            plt.ylabel('shear force')

            fig4, ax41 = plt.subplots()
            x = colors
            ax42 = ax41.twinx()
            ax41.plot(x, OSCILLATORY_STRAIN_OVER_TIME['strain'], 'g-')
            ax42.plot(x, BFORCE_SHEAR_OVER_TIME['fypos_x'], 'b-')

            ax41.set_xlabel('steps')
            ax41.set_ylabel('strain', color='g')
            ax42.set_ylabel('shear force', color='b')
            ax42.set_ylim(-35, 35)

        plt.show()


# Deal with logs
if ENSEMBLE:
    logs = simulation.getLogs()
    for i in range(len(logs)):
        steps = logs[i].getStepLog()
        manageLogs(steps, ENSEMBLE, i)
else:
    logs = simulation.getRunLog()
    steps = logs.getStepLog()
    manageLogs(steps, ENSEMBLE, 0)
