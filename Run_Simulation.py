'''
Run_Simulation.py

'''

import sys, pickle
import networkx as nx
from Epidemics_Class_2 import Epidemic_Model
from Epidemics_Class_2 import Save_Info
from Epidemics_Class_2 import Plot_Data

from  Settings import *
import matplotlib.pyplot   as plt
from matplotlib.pyplot import cm
import numpy as np

# ===========================================================================
# 						Initial parameters
# ===========================================================================

# Input parameters
N             = int(sys.argv[1]) #20

# Transmision rate
beta          = float(str(sys.argv[2]))     #0.5

# Run number
run_number    = str(sys.argv[3])  # "1"

# The networkx graph
G             = nx.erdos_renyi_graph(N,0.6) 
# G = nx.watts_strogatz_graph(N,20,0.3)
# G = nx.barabasi_albert_graph(N,10)
# The minimum and maximum value of the links weight
W_min, W_max  = int(sys.argv[4]), int(sys.argv[5])  # 1,10
Model = str(sys.argv[6])

# steps = 50
# ===========================================================================
# 						Run a single simulation
# ===========================================================================

# Epi    = Epidemic_Model(Model, N, G, beta, steps, N_infected = 1, infectious_period=10, inmmunity_period = 30, incubation_period = 6, W_range=(W_min,W_max), C_utilidad =10, x_utilidad = 15, q_vac=0.05, efect_vac=0.25 )
Epi    =Epidemic_Model(Model, N, G, beta, steps, N_infected = 1, infectious_period=6, inmmunity_period = 30, incubation_period = 5, W_range=(W_min,W_max), C_utilidad =1, x_utilidad = 0.07,vaccination_cost=0.20, vaccine_effectiveness = 0.25,utility_matrix='Payoff_Opinion_Preference_Game',strategy='mayority', initial_antivaccine_percentage = 50)
G_evol = Epi.Epidemic_Evolution()

# ===========================================================================
# 						Saving data
# ===========================================================================

outfile = open('Epi_Evolution_Run_' + run_number + '.pkl','wb')
pickle.dump(G_evol,outfile)
outfile.close()

if model_type == 'SEIR':
    dir_path = 'Output_SEIR/beta_' + str(beta) +  '/RUN_' + str(run_number)
else:
    dir_path = 'Output/beta_' + str(beta) +  '/RUN_' + str(run_number)

# Create the directory where the simulation will be stored
Save_Info().MakeDir(dir_path)

# Save the array of graphs with the simulation resutls in the directory 
Save_Info().SaveFiles(dir_path, make_dir=False, file_format='*.pkl')



# for nodo in G_evol[20]:
#     if G_evol[20].nodes[nodo]['state'] == 'V':
#         print(nodo)
    
