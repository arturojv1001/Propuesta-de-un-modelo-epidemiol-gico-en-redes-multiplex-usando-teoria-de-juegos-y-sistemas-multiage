# -*- coding: utf-8 -*-

'''
----------------------------------------------------------------------------
  Author:      --
  Date:        2022-01-27
  Last:        
  Repo/Bugs:   
  Version:     1.0.1   9
  Codec:       utf_8
  License:     GPLv3
----------------------------------------------------------------------------

	Descripción:

	Código para entrar a cada una de las carpetas generadas con el archivo Run_Simulations.py, leer los archivos hdf5 y pkl, y graficar
	la serie temporal del mapa logístico:  x(t+1) = r*x(t)*(1-x(t))


	Entradas:

			  N_runs_ensamble_max - int - numero total de carpetas ( o número de runs.)
	Salidas:
			  Una figura, en formato jpg, con la gráfica de la serie temporal dentro de la carpeta carpeta /Outputs/Run_<run_number>
			  

	Para compilar este código, es necesario primero compilar el código Run_Simulations, y después, con la instrucción 

	python3 Plot_Results.py 

	se realiza la compilación de este código.

'''

import pickle, os
import networkx            as nx
import pandas              as pd
import matplotlib.pyplot   as plt

from Epidemics_Class_2 import Save_Info
from Epidemics_Class_2 import Plot_Data
from tqdm            import tqdm
from numpy           import linspace, arange, mean, argmax
from Settings        import *



# ===========================================================================================
#           		FUNCTIONS
# ===========================================================================================

def main_menu():

  os.system("clear")

  print('\n \t\t\t\t =================[ MAIN MENU ]=================','\n')

  print('\t\t\t\t\t Options: \n')
  choice = input('''
                  1.- Generate csv file with daily infection for each run and plot.
                  2.- Read csv file with daily infection data, and get mean value of infected and plot.
                  3.- Generate csv file with risk data for each run and plot.
                  4.- Read csv file with the risk data, and get mean value of infected and plot.
                  5.- Generate csv file with daily vaccinated individuals for each run and plot.
                  6.- Read csv file with daily vaccinated data, and get mean value of infected and plot.
                  7.- Generate csv file with daily infection for each run and plot comparing simulations with and without vaccination.
                  8.- Generate csv file with number of each extremists and moderates ( both anti-vaccine & pro-vaccine) in every period of time.
                  9.- Read csv file of number of extemists and plot mean value of extremists and moderates. 
                  10.- Generate csv files with the opinion dynamic of the individuals who begin the simulation as extremists (anti-vaccines)
                  11.- Read csv files and get mean values of opinion dynamics and plot them
                  Enter your choice:

            
Note: Check that the variable 'model_type', in Settings.py, is either set to 'SEIR' or 'SEIR_Risk'  
                ''')
  return choice


def List_of_Nodes_in_State(State,Gc):
	return [ vi for vi,vj in nx.get_node_attributes(Gc, 'state').items()  if vj == State]
def List_of_Nodes_Attribute(State,Gc,type_of_state):
        return [ vi for vi,vj in nx.get_node_attributes(Gc, type_of_state).items()  if vj == State]

# ===========================================================================================
#           		SETTINGS
# ===========================================================================================

# Number of runs
#N_runs_ensamble_max = 1

# Choice and plot option
choice = main_menu()

# ===========================================================================================
#           		READ, PLOT AND SAVE DATA
# ===========================================================================================


if choice == '1':
  if model_type == 'SEIR':
      Out = 'Output_SEIR'
  else:
      Out = 'Output'
  # Run over different values of the beta parameter (transmission rate)
  for beta_i in tqdm([ round(x * 0.1,1) for x in range(beta_min, beta_max)]):

      # For each value of beta_i, a pandas data frame is generated in order to construct a csv file with the daily number of infected nodes for each run.
      Dayly_Infected = pd.DataFrame()

      # For each beta value, the code read each run file
      for run_i in range(1,N_runs_ensamble_max+1):

        # Path directory for beta_i and each run.
        
        dir_path = Out + '/beta_' + str(beta_i) +  '/RUN_' + str(run_i)

        # El nombre del archivo pkl
        pkl_file_name = 'Epi_Evolution_Run_' +  str(run_i) + '.pkl'

        # Abriendo el archivo pkl
        infile         = open(dir_path + '/' + pkl_file_name,'rb')
        G_path         = pickle.load(infile)

        # The number of infected nodes at each time step
        N_infected = [   len( List_of_Nodes_in_State('I',Gi) ) for Gi in G_path ]
        
        # Save data in a data frame 
        Dayly_Infected['Run_'+str(run_i)] = N_infected

        title  = 'Simulación número: ' + str(run_i)
        legend = r'$\beta =$ ' + str(beta_i) 

        plot   = Plot_Data(G_path)
        plot.Infected(N_infected,run_i,title,legend)

        # Guarda los archivos jpg  en la carpeta RUN_<run_i>
        Save_Info().SaveFiles(Out+'/beta_'+str(beta_i), make_dir=False, file_format='*.png')

      # Guardamos la información en un archivo csv.
      Dayly_Infected.to_csv('Daily_Infected_for_beta_' + str(beta_i) + '.csv')
      Save_Info().SaveFiles(Out+'/beta_' + str(beta_i) , make_dir=False, file_format='*.csv')

if choice == '2':

  fontsize  = 15
  linewidth = 2
  labelsize = 15
  alpha     = 0.4

  cont = 1
  if model_type == 'SEIR':
      Out = 'Output_SEIR'
  else:
      Out = 'Output'
  for beta_i in tqdm([ round(x * 0.1,1) for x in range(beta_min, beta_max)] ):

    Daily_I_Data = pd.read_csv( Out+'/beta_' + str(beta_i) + '/Daily_Infected_for_beta_' + str(beta_i) + '.csv' )

    Mean_Daily_I = list( sum(  [ Daily_I_Data['Run_'+str(i)] for i in range(1,N_runs_ensamble_max+1)]  )/N_runs_ensamble_max )

    colors = plt.get_cmap('binary')(linspace(0, 1.0, N_runs_ensamble_max+1))

    plt.clf()
    plt.figure(100+cont,figsize = (8,6))
    for i in range(1,N_runs_ensamble_max+1):
      plt.plot(Daily_I_Data['Run_'+str(i)],'o-',markersize=3,color='lightgray',linewidth=linewidth, alpha = alpha)

    # l, = plt.plot(Mean_Daily_I,'o-',color='red',markersize=3,linewidth=linewidth, alpha = 1., label = 'Average daily number')
    l, = plt.plot(Mean_Daily_I,color='red',linewidth=linewidth, alpha = 1., label = 'Número promedio diario')
    l.set_markerfacecolor("lightcoral")
    l.set_color('red')

    # plt.title(r'Grafo Erd$\"{o}$s-Renyi $\beta = $' + str(beta_i))
    plt.title(r'Grafo Barabasi-Albert $\beta = $' + str(beta_i))

    plt.xlabel("Tiempo (días)",fontsize=fontsize)
    plt.ylabel("Número de casos diarios",fontsize=fontsize)

    plt.tick_params(axis='y',labelsize=labelsize)
    plt.tick_params(axis='x',labelsize=labelsize)

    plt.legend(loc=1, bbox_to_anchor=(1.1,1.),fontsize=10,ncol=1, fancybox=False, shadow=False)
    plt.savefig('Daily_I_beta_' + str(beta_i) + '.png',dpi=320)
    Save_Info().SaveFiles(Out+'/', make_dir=False, file_format='*.png')

    cont += 1

if choice == '3':

  # Run over different values of the beta parameter (transmission rate)
  for beta_i in tqdm([ round(x * 0.1,1) for x in range(beta_min, beta_max)]):

      # For each value of beta_i, a pandas data frame is generated in order to construct a csv file with the daily number of infected nodes for each run.
      Dayly_Risk = pd.DataFrame(columns=['Run','Node','Infected_Neighbors','Exposed','Second_Inf_Neighborhood','Danger','Risk','Risk_exp'])

      # For each beta value, the code read each run file
      for run_i in range(1,N_runs_ensamble_max+1):

        # Path directory for beta_i and each run.
        dir_path = 'Output/beta_' + str(beta_i) +  '/RUN_' + str(run_i)

        # El nombre del archivo pkl
        pkl_file_name = 'Epi_Evolution_Run_' +  str(run_i) + '.pkl'

        # Abriendo el archivo pkl
        infile         = open(dir_path + '/' + pkl_file_name,'rb')
        G_path         = pickle.load(infile)


        Dayly_Risk['Run']                     = [ run_i                                                 for node_i in G_path[0].nodes() for run_i,Gi in enumerate(G_path) ]
        Dayly_Risk['Node']                    = [ node_i                                                for node_i in G_path[0].nodes() for Gi       in G_path ]
        Dayly_Risk['Infected_Neighbors']      = [ round(Gi.nodes[node_i]['Infected_Neighbors'],3)       for node_i in G_path[0].nodes() for Gi       in G_path ]
        Dayly_Risk['Exposed']                 = [ round(Gi.nodes[node_i]['Exposed'],3)                  for node_i in G_path[0].nodes() for Gi       in G_path ]
        Dayly_Risk['Second_Inf_Neighborhood'] = [ round(Gi.nodes[node_i]['Second_Inf_Neighborhood'],3)  for node_i in G_path[0].nodes() for Gi       in G_path ]
        Dayly_Risk['Danger']                  = [ round(Gi.nodes[node_i]['Danger'],3)                   for node_i in G_path[0].nodes() for Gi       in G_path ]
        Dayly_Risk['Risk']                    = [ round(Gi.nodes[node_i]['Risk'],3)                     for node_i in G_path[0].nodes() for Gi       in G_path ]
        Dayly_Risk['Risk_exp']                = [ round(Gi.nodes[node_i]['Risk_exp'],3)                 for node_i in G_path[0].nodes() for Gi       in G_path ]

      # Guardamos la información en un archivo csv.
      Dayly_Risk.to_csv('Risk_for_beta_' + str(beta_i) + '.csv')
      Save_Info().SaveFiles('Output/beta_' + str(beta_i) , make_dir=False, file_format='*.csv')

if choice == '4':
  
  # Plot settings
  fontsize  = 20
  linewidth = 1.
  ms        = 1.2
  labelsize = 10
  alpha     = 0.4

  # Colors for each node
  colors     = plt.get_cmap('jet')(linspace(0, 1.0, N+1))
  colors_exp = plt.get_cmap('binary')(linspace(0, 1.0, N+1))

  cont = 1
  for beta_i in tqdm([ round(x * 0.1,1) for x in range(beta_min, beta_max)] ):

    Risk_Data = pd.read_csv('Output/beta_' + str(beta_i) + '/Risk_for_beta_' + str(beta_i) + '.csv' )

    plt.clf()
    plt.figure(200+cont,figsize = (8,6))
    plt.title(r'Exposed for $\beta = $' + str(beta_i))

    plt.subplot(311)
    for vi in range(N):
      plt.plot(Risk_Data.loc[Risk_Data['Node']==vi]['Exposed'].values,'o-',markersize=ms,color=colors[vi],linewidth=linewidth, alpha = alpha)

    plt.xticks(color='w')
    plt.tick_params(axis='y',labelsize=labelsize)
    plt.ylabel("Exposicion",fontsize=fontsize)

    plt.subplot(312)
    for vi in range(N):
      plt.plot(Risk_Data.loc[Risk_Data['Node']==vi]['Danger'].values,'o-',markersize=ms,color=colors[vi],linewidth=linewidth, alpha = alpha)

    plt.ylabel("Peligro",fontsize=fontsize)
    plt.tick_params(axis='y',labelsize=labelsize)
    plt.xticks(color='w')

    plt.subplot(313)
    for vi in range(N):
      plt.plot(Risk_Data.loc[Risk_Data['Node']==vi]['Risk'].values,'o-',markersize=ms,color=colors[vi],linewidth=linewidth, alpha = alpha)
      plt.plot(Risk_Data.loc[Risk_Data['Node']==vi]['Risk_exp'].values,'o-',markersize=ms,color=colors_exp[vi],linewidth=linewidth, alpha = 0.1)

    plt.ylabel("Riesgo",fontsize=fontsize)
    plt.xlabel("Tiempo (días)",fontsize=fontsize)
    
    plt.tick_params(axis='y',labelsize=labelsize)
    plt.tick_params(axis='x',labelsize=labelsize)

    plt.savefig('Exposed_beta_' + str(beta_i) + '.png',dpi=320)

    plt.clf()
    plt.figure(300+cont,figsize = (8,6))
    plt.title(r'Exposed for $\beta = $' + str(beta_i))

    plt.subplot(211)
    for vi in range(N):
      plt.plot(Risk_Data.loc[Risk_Data['Node']==vi]['Infected_Neighbors'].values,'o-',markersize=ms,color=colors[vi],linewidth=linewidth, alpha = alpha,label = 'Node ' + str(vi+1))

    #plt.legend(loc=1, bbox_to_anchor=(0.95,1.6),fontsize=10,ncol= int(N/2), fancybox=False, shadow=False)

    plt.xticks(color='w')
    plt.tick_params(axis='y',labelsize=labelsize)
    plt.ylabel('$\frac{s(i)}{k(i)}$     ',fontsize=25,rotation=0)
    
    plt.subplot(212)
    for vi in range(N):
      plt.plot(Risk_Data.loc[Risk_Data['Node']==vi]['Second_Inf_Neighborhood'].values,'o-',markersize=ms,color=colors[vi],linewidth=linewidth, alpha = alpha)

    plt.ylabel(r"$H_i$    ",fontsize=35,rotation=0)
    plt.tick_params(axis='y',labelsize=labelsize)
    plt.xticks(color='w')

    plt.xlabel("Tiempo (días)",fontsize=fontsize)
    
    plt.tick_params(axis='y',labelsize=labelsize)
    plt.tick_params(axis='x',labelsize=labelsize)

    plt.savefig('Neighborhood_' + str(beta_i) + '.png',dpi=320)
    
    Save_Info().SaveFiles('Output/', make_dir=False, file_format='*.png')

    
    cont += 1




if choice == '5':

  # Run over different values of the beta parameter (transmission rate)
  for beta_i in tqdm([ round(x * 0.1,1) for x in range(beta_min, beta_max)]):

      # For each value of beta_i, a pandas data frame is generated in order to construct a csv file with the daily number of infected nodes for each run.
      Dayly_Vaccinated = pd.DataFrame()

      # For each beta value, the code read each run file
      for run_i in range(1,N_runs_ensamble_max+1):

        # Path directory for beta_i and each run.
        dir_path = 'Output/beta_' + str(beta_i) +  '/RUN_' + str(run_i) 

        # El nombre del archivo pkl
        pkl_file_name = 'Epi_Evolution_Run_' +  str(run_i) + '.pkl'

        # Abriendo el archivo pkl
        infile         = open(dir_path + '/' + pkl_file_name,'rb')
        G_path         = pickle.load(infile)

        # The number of infected nodes at each time step
        # N_vaccinated = [   len( List_of_Nodes_Any_Attribute(2,Gi,'M_scale') ) for Gi in G_path ]
        N_vaccinated = [   len( List_of_Nodes_Attribute(1,Gi,'Vaccinated_flag') ) for Gi in G_path ]
        # N_vaccinated = [   len( List_of_Nodes_in_State('V',Gi) ) for Gi in G_path ]
        # Save data in a data frame 
        Dayly_Vaccinated['Run_'+str(run_i)] = N_vaccinated

        title  = 'Simulation number: ' + str(run_i)
        legend = r'$\beta =$ ' + str(beta_i) 

        plot   = Plot_Data(G_path)
        plot.Vaccinated(N_vaccinated,run_i,title,legend)

        # Guarda los archivos jpg  en la carpeta RUN_<run_i>
        Save_Info().SaveFiles('Output/beta_'+str(beta_i), make_dir=False, file_format='*.png')

      # Guardamos la información en un archivo csv.
      Dayly_Vaccinated.to_csv('Daily_Vaccinated_for_beta_' + str(beta_i) + '.csv')
      Save_Info().SaveFiles('Output/beta_' + str(beta_i) , make_dir=False, file_format='*.csv')
      
      
    
if choice == '6':

  fontsize  = 15
  linewidth = 2
  labelsize = 15
  alpha     = 0.4

  cont = 1
  for beta_i in tqdm([ round(x * 0.1,1) for x in range(beta_min, beta_max)] ):

    Daily_V_Data = pd.read_csv( 'Output/beta_' + str(beta_i) + '/Daily_Vaccinated_for_beta_' + str(beta_i) + '.csv' )

    Mean_Daily_V = list( sum(  [ Daily_V_Data['Run_'+str(i)] for i in range(1,N_runs_ensamble_max+1)]  )/N_runs_ensamble_max )

    colors = plt.get_cmap('binary')(linspace(0, 1.0, N_runs_ensamble_max+1))

    plt.clf()
    plt.figure(100+cont,figsize = (8,6))
    for i in range(1,N_runs_ensamble_max+1):
      plt.plot(Daily_V_Data['Run_'+str(i)],'o-',markersize=3,color='lightgray',linewidth=linewidth, alpha = alpha)

    l, = plt.plot(Mean_Daily_V,'o-',color='red',markersize=3,linewidth=linewidth, alpha = 1., label = 'Número promedio diario')
    l.set_markerfacecolor("lightcoral")
    l.set_color('red')

    plt.title(r'Serie de tiempo de vacunados para $\beta = $' + str(beta_i))

    plt.xlabel("Tiempo (días)",fontsize=fontsize)
    plt.ylabel("Número de agentes vacunados",fontsize=fontsize)

    plt.tick_params(axis='y',labelsize=labelsize)
    plt.tick_params(axis='x',labelsize=labelsize)

    plt.legend(loc=1, bbox_to_anchor=(1.1,1.),fontsize=10,ncol=1, fancybox=False, shadow=False)
    plt.savefig('Daily_V_beta_' + str(beta_i) + '.png',dpi=320)
    Save_Info().SaveFiles('Output/', make_dir=False, file_format='*.png')

    cont += 1
    
    
    
    
if choice == '7':

  fontsize  = 15
  linewidth = 2
  labelsize = 15
  alpha     = 0.4

  cont = 1
  for beta_i in tqdm([ round(x * 0.1,1) for x in range(beta_min, beta_max)] ):

    Daily_IV_Data = pd.read_csv( 'Output/beta_' + str(beta_i) + '/Daily_Infected_for_beta_' + str(beta_i) + '.csv' )
    Daily_I_Data = pd.read_csv( 'Output_SEIR/beta_' + str(beta_i) + '/Daily_Infected_for_beta_' + str(beta_i) + '.csv' )
    
    Mean_Daily_IV = list( sum(  [ Daily_IV_Data['Run_'+str(i)] for i in range(1,N_runs_ensamble_max+1)]  )/N_runs_ensamble_max )
    Mean_Daily_I = list( sum(  [ Daily_I_Data['Run_'+str(i)] for i in range(1,N_runs_ensamble_max+1)]  )/N_runs_ensamble_max )
    
    colors = plt.get_cmap('binary')(linspace(0, 1.0, N_runs_ensamble_max+1))

    plt.clf()
    plt.figure(100+cont,figsize = (8,6))
    # for i in range(1,N_runs_ensamble_max+1):
      # plt.plot(Daily_I_Data['Run_'+str(i)],'o-',markersize=3,color='lightgray',linewidth=linewidth, alpha = alpha)

    l, = plt.plot(Mean_Daily_I,'o-',color='red',markersize=3,linewidth=linewidth, alpha = 1., label = 'promedio infectados')
    l2, = plt.plot(Mean_Daily_IV,'o-',color='blue',markersize=3,linewidth=linewidth, alpha = 1., label = 'promedio infectados implementando vacunas')
    l.set_markerfacecolor("lightcoral")
    l.set_color('red')
    l2.set_markerfacecolor("lightseagreen")
    l2.set_color('blue')
    

    plt.title(r'Curva promedio de contagios con y sin vacunación')
      
    plt.xlabel("Tiempo (días)",fontsize=fontsize)
    plt.ylabel("Número de casos diarios",fontsize=fontsize)
      
    plt.tick_params(axis='y',labelsize=labelsize)
    plt.tick_params(axis='x',labelsize=labelsize)
      
    plt.legend(loc=1, bbox_to_anchor=(1.1,1.),fontsize=10,ncol=1, fancybox=False, shadow=False)
    plt.savefig('Average_I_beta_' + str(beta_i) + '.png',dpi=320)
    Save_Info().SaveFiles('Output/', make_dir=False, file_format='*.png')
  
    cont += 1
    
    
    
    
if choice == '8':

  # Run over different values of the beta parameter (transmission rate)
  for beta_i in tqdm([ round(x * 0.1,1) for x in range(beta_min, beta_max)]):

      # For each value of beta_i, a pandas data frame is generated in order to construct a csv file with the daily number of infected nodes for each run.
      Dayly_Extremists = pd.DataFrame()
      Dayly_provaccine_moderate = pd.DataFrame()
      Dayly_provaccine_extremists = pd.DataFrame()
      Dayly_antivaccine_moderate = pd.DataFrame()
      
      # For each beta value, the code read each run file
      for run_i in range(1,N_runs_ensamble_max+1):

        # Path directory for beta_i and each run.
        dir_path = 'Output/beta_' + str(beta_i) +  '/RUN_' + str(run_i) 

        # El nombre del archivo pkl
        pkl_file_name = 'Epi_Evolution_Run_' +  str(run_i) + '.pkl'

        # Abriendo el archivo pkl
        infile         = open(dir_path + '/' + pkl_file_name,'rb')
        G_path         = pickle.load(infile)

        # The number of infected nodes at each time step
        N_Extremists = [   len( List_of_Nodes_Attribute(-2,Gi,'M_scale') ) for Gi in G_path ]
        antivaccine_moderate = [   len( List_of_Nodes_Attribute(-1,Gi,'M_scale') ) for Gi in G_path ]
        provaccine_moderate = [   len( List_of_Nodes_Attribute(1,Gi,'M_scale') ) for Gi in G_path ]
        provaccine_extremist = [   len( List_of_Nodes_Attribute(2,Gi,'M_scale') ) for Gi in G_path ]
        
        # Save data in a data frame 
        Dayly_Extremists['Run_'+str(run_i)] = N_Extremists
        Dayly_Extremists.index.name = 'Time'
        
        Dayly_antivaccine_moderate['Run_'+str(run_i)] = antivaccine_moderate
        Dayly_provaccine_moderate['Run_'+str(run_i)]  = provaccine_moderate
        Dayly_provaccine_extremists['Run_'+str(run_i)] = provaccine_extremist
        
        title  = 'Simulation number: ' + str(run_i)
        # legend = r'$\beta =$ ' + str(beta_i) 

        plot   = Plot_Data(G_path)
        plot.Extremists(N_Extremists, antivaccine_moderate, provaccine_moderate, provaccine_extremist, run_i,title)

        # Guarda los archivos jpg  en la carpeta RUN_<run_i>
        Save_Info().SaveFiles('Output/beta_'+str(beta_i), make_dir=False, file_format='*.png')

      # Guardamos la información en un archivo csv.
      Dayly_Extremists.to_csv('Daily_Extremists_for_beta_' + str(beta_i) + '.csv')
      Save_Info().SaveFiles('Output/beta_' + str(beta_i) , make_dir=False, file_format='*.csv')
      
      Dayly_antivaccine_moderate.to_csv('Daily_antivaccine_moderate_for_beta_' + str(beta_i) + '.csv')
      Save_Info().SaveFiles('Output/beta_' + str(beta_i) , make_dir=False, file_format='*.csv')
      
      Dayly_provaccine_moderate.to_csv('Daily_provaccine_moderate_for_beta_' + str(beta_i) + '.csv')
      Save_Info().SaveFiles('Output/beta_' + str(beta_i) , make_dir=False, file_format='*.csv')
      
      Dayly_provaccine_extremists.to_csv('Daily_provaccine_extremists_for_beta_' + str(beta_i) + '.csv')
      Save_Info().SaveFiles('Output/beta_' + str(beta_i) , make_dir=False, file_format='*.csv')
      
    
if choice == '9':

  fontsize  = 15
  linewidth = 2
  labelsize = 15
  alpha     = 0.4

  cont = 1
  for beta_i in tqdm([ round(x * 0.1,1) for x in range(beta_min, beta_max)] ):

    Daily_Ex_Data = pd.read_csv( 'Output/beta_' + str(beta_i) + '/Daily_Extremists_for_beta_' + str(beta_i) + '.csv' )
    Daily_AM_Data = pd.read_csv( 'Output/beta_' + str(beta_i) + '/Daily_antivaccine_moderate_for_beta_' + str(beta_i) + '.csv' )
    Daily_PM_Data = pd.read_csv( 'Output/beta_' + str(beta_i) + '/Daily_provaccine_moderate_for_beta_' + str(beta_i) + '.csv' )
    Daily_PE_Data = pd.read_csv( 'Output/beta_' + str(beta_i) + '/Daily_provaccine_extremists_for_beta_' + str(beta_i) + '.csv' )
    # Provaccine_Opinion_Counter_Data = pd.read_csv( 'Output/beta_' + str(beta_i) + '/Provaccine_Opinion_Counter_for_beta_' + str(beta_i) + '.csv' )

    Mean_Daily_Ex = list( sum(  [ Daily_Ex_Data['Run_'+str(i)] for i in range(1,N_runs_ensamble_max+1)]  )/N_runs_ensamble_max )
    Mean_Daily_AM = list( sum(  [ Daily_AM_Data['Run_'+str(i)] for i in range(1,N_runs_ensamble_max+1)]  )/N_runs_ensamble_max )
    Mean_Daily_PM = list( sum(  [ Daily_PM_Data['Run_'+str(i)] for i in range(1,N_runs_ensamble_max+1)]  )/N_runs_ensamble_max )
    Mean_Daily_PE = list( sum(  [ Daily_PE_Data['Run_'+str(i)] for i in range(1,N_runs_ensamble_max+1)]  )/N_runs_ensamble_max )
    # provaccine_dominance_opinion_counter = sum(  [ Provaccine_Opinion_Counter_Data['Run_'+str(i)] for i in range(1,N_runs_ensamble_max+1)]  )
    provaccine_dominance_opinion_counter = 0
    last_row_PM = Daily_PM_Data.iloc[-1]
    last_row_PE = Daily_PE_Data.iloc[-1]
    for i in range(1,N_runs_ensamble_max+1):
        if last_row_PM['Run_'+str(i)] + last_row_PE['Run_'+str(i)] > int(N/2):
            provaccine_dominance_opinion_counter += 1
    
    colors = plt.get_cmap('binary')(linspace(0, 1.0, N_runs_ensamble_max+1))

    plt.clf()
    plt.figure(100+cont,figsize = (8,6))
    # for i in range(1,N_runs_ensamble_max+1):
    #   plt.plot(Daily_Ex_Data['Run_'+str(i)],'o-',markersize=3,color='lightgray',linewidth=linewidth, alpha = alpha)

    l, = plt.plot(Mean_Daily_Ex,'o-',color='red',markersize=3,linewidth=linewidth, alpha = 1., label = 'extremistas anti-vacunas')
    l2, = plt.plot(Mean_Daily_AM,'o-',color='orange',markersize=3,linewidth=linewidth, alpha = 1., label = 'moderados anti-vacunas')
    l3, = plt.plot(Mean_Daily_PM,'o-',color='blue',markersize=3,linewidth=linewidth, alpha = 1., label = 'moderados pro-vacunas')
    l4, = plt.plot(Mean_Daily_PE,'o-',color='green',markersize=3,linewidth=linewidth, alpha = 1., label = 'extremistas pro-vacunas')
    
    l.set_markerfacecolor("lightcoral")
    l.set_color('red')

    plt.title(r'Dinámica de opiniones de los agentes')

    plt.xlabel("Tiempo (días)",fontsize=fontsize)
    plt.ylabel("Número de agentes",fontsize=fontsize)
    plt.text(40.1,1.0,"Dominancia de opinión pro-vacunas: "+str(100*provaccine_dominance_opinion_counter/N_runs_ensamble_max) +'%',fontsize=8)

    plt.tick_params(axis='y',labelsize=labelsize)
    plt.tick_params(axis='x',labelsize=labelsize)

    plt.legend(loc=1, bbox_to_anchor=(1.1,1.),fontsize=10,ncol=1, fancybox=False, shadow=False)
    plt.savefig('Daily_Ex_beta_' + str(beta_i) + '.png',dpi=320)
    Save_Info().SaveFiles('Output/', make_dir=False, file_format='*.png')

    cont += 1
    
'''    
if choice == '10':

  # Run over different values of the beta parameter (transmission rate)
  for beta_i in tqdm([ round(x * 0.1,1) for x in range(beta_min, beta_max)]):

      # For each value of beta_i, a pandas data frame is generated in order to construct a csv file with the daily number of infected nodes for each run.
      Dayly_Extremists = pd.DataFrame()

      # For each beta value, the code read each run file
      for run_i in range(1,N_runs_ensamble_max+1):

        # Path directory for beta_i and each run.
        dir_path = 'Output/beta_' + str(beta_i) +  '/RUN_' + str(run_i) 

        # El nombre del archivo pkl
        pkl_file_name = 'Epi_Evolution_Run_' +  str(run_i) + '.pkl'

        # Abriendo el archivo pkl
        infile         = open(dir_path + '/' + pkl_file_name,'rb')
        G_path         = pickle.load(infile)

        # The number of infected nodes at each time step
        N_Extremists =     List_of_Nodes_Attribute(1,G_path[-1],'-2_flag') 
        # Dayly_Extremists['Run_'+str(run_i)] = N_Extremists
        Dayly_Extremists['Run'] = run_i
        
        title  = 'Simulation number: ' + str(run_i)
        legend = r'$\beta =$ ' + str(beta_i) 

        plot   = Plot_Data(G_path)
        plot.Extremists_Opinion(N_Extremists,run_i,title,legend,Dayly_Extremists)
        Dayly_Extremists.index.name = 'Time'
        
        # Guarda los archivos jpg  en la carpeta RUN_<run_i>
        Save_Info().SaveFiles('Output/beta_'+str(beta_i), make_dir=False, file_format='*.png')
        
      # Guardamos la información en un archivo csv.
        Dayly_Extremists.to_csv('Daily_Initial_Extremists_for_beta_' + str(run_i) + '.csv')
        Save_Info().SaveFiles('Output/beta_' + str(beta_i) , make_dir=False, file_format='*.csv')
      
        Dayly_Extremists_Total = pd.DataFrame()
      for run_i in range(1,N_runs_ensamble_max+1):
          temp = pd.read_csv( 'Output/beta_' + str(beta_i) + '/Daily_Initial_Extremists_for_beta_' + str(run_i) + '.csv' )    
          Dayly_Extremists_Total = pd.concat([Dayly_Extremists_Total,temp])
            
      Dayly_Extremists_Total.to_csv('Total_Initial_antivaccines.csv',index=False)
      Save_Info().SaveFiles('Output/beta_' + str(beta_i) , make_dir=False, file_format='*.csv')
        
      
        
if choice == '11':

  fontsize  = 15
  linewidth = 2
  labelsize = 15
  alpha     = 0.4

  for beta_i in tqdm([ round(x * 0.1,1) for x in range(beta_min, beta_max)] ):
      Dayly_Extremists_Total = pd.read_csv( 'Output/beta_' + str(beta_i) + '/Total_Initial_antivaccines.csv' )
      Dayly_Extremists_Total['Average'] = Dayly_Extremists_Total.loc[:,[c for c in Dayly_Extremists_Total.columns if (c!= "Time" and c!="Run")]].mean(axis=1)
      # Dayly_Extremists_Total.to_csv('Total_Initial_Extremists.csv',index=False)
      # Save_Info().SaveFiles('Output/beta_' + str(beta_i) , make_dir=False, file_format='*.csv')
      
      Av = Dayly_Extremists_Total.drop('Run',axis=1)
      Av= Dayly_Extremists_Total.groupby("Time").mean()
      Av.to_csv('Av.csv',index=False)
      Save_Info().SaveFiles('Output/beta_' + str(beta_i) , make_dir=False, file_format='*.csv')
      
      plt.plot(Av['Average'])

      plt.title('Opinion dynamics for antivaccines individuals')

      plt.xlabel("Time",fontsize=fontsize)
      plt.ylabel("Initial Extremist's people opinion",fontsize=fontsize)

      plt.tick_params(axis='y',labelsize=labelsize)
      plt.tick_params(axis='x',labelsize=labelsize)

      plt.legend(loc=1, bbox_to_anchor=(1.1,1.),fontsize=10,ncol=1, fancybox=False, shadow=False)
      plt.savefig('Output/Opinion_dynamics_initial_antivaccines.png',dpi=320)
      
'''      
      
      
      
      
      
if choice == '10':

  # Run over different values of the beta parameter (transmission rate)
  for beta_i in tqdm([ round(x * 0.1,1) for x in range(beta_min, beta_max)]):

      # For each value of beta_i, a pandas data frame is generated in order to construct a csv file with the daily number of infected nodes for each run.
      Dayly_Extremists = pd.DataFrame()
      Dayly_Extremists_Total = pd.DataFrame()
      # For each beta value, the code read each run file
      for run_i in range(1,N_runs_ensamble_max+1):

        # Path directory for beta_i and each run.
        dir_path = 'Output/beta_' + str(beta_i) +  '/RUN_' + str(run_i) 

        # El nombre del archivo pkl
        pkl_file_name = 'Epi_Evolution_Run_' +  str(run_i) + '.pkl'

        # Abriendo el archivo pkl
        infile         = open(dir_path + '/' + pkl_file_name,'rb')
        G_path         = pickle.load(infile)

        # The number of infected nodes at each time step
        N_Extremists =     List_of_Nodes_Attribute(1,G_path[-1],'-2_flag') 
        # Dayly_Extremists['Run_'+str(run_i)] = N_Extremists
        Dayly_Extremists['Run'] = str(run_i)
        
        title  = 'Simulation number: ' + str(run_i)
        legend = r'$\beta =$ ' + str(beta_i) 

        plot   = Plot_Data(G_path)
        plot.Extremists_Opinion(N_Extremists,run_i,title,legend,Dayly_Extremists)
        Dayly_Extremists.index.name = 'Time'
        Dayly_Extremists_Total.name = 'Time'
        # Guarda los archivos jpg  en la carpeta RUN_<run_i>
        Save_Info().SaveFiles('Output/beta_'+str(beta_i), make_dir=False, file_format='*.png')
        Dayly_Extremists_Total = pd.concat([Dayly_Extremists_Total,Dayly_Extremists])
      # Guardamos la información en un archivo csv.
        # Dayly_Extremists.to_csv('Daily_Initial_Extremists_for_beta_' + str(run_i) + '.csv')
        # Save_Info().SaveFiles('Output/beta_' + str(beta_i) , make_dir=False, file_format='*.csv')
      
        # Dayly_Extremists_Total = pd.DataFrame()
      # for run_i in range(1,N_runs_ensamble_max+1):
      #     temp = pd.read_csv( 'Output/beta_' + str(beta_i) + '/Daily_Initial_Extremists_for_beta_' + str(run_i) + '.csv' )    
          # Dayly_Extremists_Total = pd.concat([Dayly_Extremists_Total,temp])
            
      Dayly_Extremists_Total.to_csv('Total_Initial_antivaccines.csv')
      Save_Info().SaveFiles('Output/beta_' + str(beta_i) , make_dir=False, file_format='*.csv')
        
      
        
if choice == '11':

  fontsize  = 15
  linewidth = 2
  labelsize = 15
  alpha     = 0.4

  for beta_i in tqdm([ round(x * 0.1,1) for x in range(beta_min, beta_max)] ):
      Dayly_Extremists_Total = pd.read_csv( 'Output/beta_' + str(beta_i) + '/Total_Initial_antivaccines.csv' )
      Dayly_Extremists_Total['Average'] = Dayly_Extremists_Total.loc[:,[c for c in Dayly_Extremists_Total.columns if (c!= "Time" and c!="Run")]].mean(axis=1)
      # Dayly_Extremists_Total.to_csv('Total_Initial_Extremists.csv',index=False)
      # Save_Info().SaveFiles('Output/beta_' + str(beta_i) , make_dir=False, file_format='*.csv')
      
      Av = Dayly_Extremists_Total.drop('Run',axis=1)
      Av= Dayly_Extremists_Total.groupby("Time").mean()
      Av.to_csv('Av.csv',index=False)
      Save_Info().SaveFiles('Output/beta_' + str(beta_i) , make_dir=False, file_format='*.csv')
      
      plt.plot(Av['Average'])

      plt.title('Opinion dynamics for antivaccines individuals')

      plt.xlabel("Time",fontsize=fontsize)
      plt.ylabel("Initial Extremist's people opinion",fontsize=fontsize)

      plt.tick_params(axis='y',labelsize=labelsize)
      plt.tick_params(axis='x',labelsize=labelsize)

      plt.legend(loc=1, bbox_to_anchor=(1.1,1.),fontsize=10,ncol=1, fancybox=False, shadow=False)
      plt.savefig('Output/Opinion_dynamics_initial_antivaccines.png',dpi=320)
  
  cont = 1
  '''
  for beta_i in tqdm([ round(x * 0.1,1) for x in range(beta_min, beta_max)] ):

    Daily_Initial_Ex_Data = pd.read_csv( 'Output/beta_' + str(beta_i) + '/Total_Initial_Extremists.csv' )

    Mean_Initial_Daily_Ex = list( sum(  [ Daily_Initial_Ex_Data['Run_'+str(i)] for i in range(1,N_runs_ensamble_max+1)]  )/N_runs_ensamble_max )

    colors = plt.get_cmap('binary')(linspace(0, 1.0, N_runs_ensamble_max+1))

    plt.clf()
    plt.figure(100+cont,figsize = (8,6))
    for i in range(1,N_runs_ensamble_max+1):
      plt.plot(Daily_Initial_Ex_Data['Run_'+str(i)],'o-',markersize=3,color='lightgray',linewidth=linewidth, alpha = alpha)

    l, = plt.plot(Mean_Initial_Daily_Ex,'o-',color='red',markersize=3,linewidth=linewidth, alpha = 1., label = 'Average daily number')
    l.set_markerfacecolor("lightcoral")
    l.set_color('red')

    plt.title(r'Daily number per run simulation for $\beta = $' + str(beta_i))

    plt.xlabel("Time",fontsize=fontsize)
    plt.ylabel("Initial Extremist's people opinion",fontsize=fontsize)

    plt.tick_params(axis='y',labelsize=labelsize)
    plt.tick_params(axis='x',labelsize=labelsize)

    plt.legend(loc=1, bbox_to_anchor=(1.1,1.),fontsize=10,ncol=1, fancybox=False, shadow=False)
    plt.savefig('Daily_Initial_Ex_beta_' + str(beta_i) + '.png',dpi=320)
    Save_Info().SaveFiles('Output/', make_dir=False, file_format='*.png')

    cont += 1
    '''