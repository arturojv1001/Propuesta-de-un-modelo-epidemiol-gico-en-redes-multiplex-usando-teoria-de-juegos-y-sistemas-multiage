"""
@name  : Run_Ensamble.py
@author: ---
@date:   January, 2022
@update: ----- Sympy
"""

import pickle,csv,sys

from tqdm            import tqdm
from numpy           import linspace
from random          import uniform
from subprocess      import call
from tqdm            import tqdm
from Epidemics_Class_2 import Save_Info
from Settings        import *

# ========================================================================================================================
#                                     Parameters
# =========================================================================================================================

# Run over different values of the beta parameter (transmission rate)
for beta_i in [ round(x * 0.1,1) for x in range(beta_min, beta_max)]:

    # For each beta value, the code performs an ensamble of N_runs_ensamble_max of runs
    for run_i in tqdm(range(1,N_runs_ensamble_max+1)):

        # ========================================================================================================================
        #                                     Run a simulation with parameter beta_i
        # =========================================================================================================================

        # Con la libreria subprocess, compilar el archivo Mapa_Logistico.py con los parámetros de entrada r, xo y run_i
        call(['python','Run_Simulation.py', str(N), str(beta_i), str(run_i), str(W_min), str(W_max), str(model_type) ])
