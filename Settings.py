'''
Settings.py
'''

N             = 50         #20 #100      # Number of nodes
steps         = 100         #15  #150      # Number of steps
W_min, W_max  = 1,10        # The minimum and maximum value of the links weights

# Number of runs
N_runs_ensamble_max = 20

# The range to explore beta parameter
beta_min = 4    # from 0.5
beta_max = 5   # to   1.0
# model_type = 'SEIR_Vaccination_Game'
model_type = 'SEIR'
# model_type = 'SEIR_Risk'