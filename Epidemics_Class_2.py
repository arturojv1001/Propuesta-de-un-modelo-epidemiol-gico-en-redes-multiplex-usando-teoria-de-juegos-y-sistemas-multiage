'''
Epidemics_Class.py
'''

# agregar reinfeccion de vacunados (incluir contador de vacunacion)
# juntar grafica de curva de contagios con y sin vacunacion
# variar parametros de costos y efectividad de la vacuna

import shutil,os,glob
import numpy             as np
import networkx          as nx
import matplotlib.pyplot as plt
import random
from matplotlib.pyplot import cm

from random import sample
from tqdm   import tqdm


class Epidemic_Model:

    def __init__(self, Model, N, Go, beta, steps, 
                 N_infected=1, infectious_period=6, inmmunity_period =30, incubation_period = 5, W_range = (1,10),  
                 C_utilidad =1, x_utilidad = 0.07, vaccination_cost=0.4, vaccine_effectiveness = 0.85,utility_matrix='Payoff_Opinion_Preference_Game',strategy='mayority', initial_antivaccine_percentage=50):
        
        self.Model                 = Model                        # (str)   The name of the model. Options: 'SEIR'
        self.N                     = N                            # (int)   Number of nodes
        self.Go                    = Go                           # (networkx graph) Network configuration without initial atributes
        self.beta                  = beta                         # (float) Transmission rate
        self.steps                 = steps                        # (int)   Number of time steps        

        self.N_infected            = N_infected                   # (int) Number of initial infected nodes
        self.infectious_period     = infectious_period            # (int) number of periods of time before before an infected node can be recover
        self.incubation_period     = incubation_period            # (int) Number of periods of time before an infected person is contagious (before it can infect other nodes) syntomatic? - duda, parece que es igual a infectious_period
        self.inmmunity_period      = inmmunity_period             # (int) number of periods of time before before a recovered node can become susceptible again    
        self.W_range               = W_range                      # (tuple) The minimum and maximum weight an edge can have

        self.C_utilidad            = C_utilidad
        self.x_utilidad            = x_utilidad
        self.vaccination_cost      = vaccination_cost                        #
        self.vaccine_effectiveness = vaccine_effectiveness        # (float) vaccine has a certain effectiveness w
        self.utility_matrix        = utility_matrix
        self.strategy              = strategy
        self.initial_antivaccine_percentage = initial_antivaccine_percentage 
        return

    def Initial_Graph_Atributes(self):
        '''
        Function to define the initial atributes for each node of the network as the epedimic state (susceptible or infected), cooperator or detractor, vulnerability; 
        danges, exposition, risk, and also the define the weight of each link as a int number between zero and W_max.

        The list of a node atributes are:
        'state'         - which defines the state of the node in the epidemic (Susceptible,Infected,Recovered,Deceased)
        'days_I'        - amount of days that a node has been infected
        'days_R'        - amount of days that a node has been recovered
        'CoD'           - state that defines wether a node is a detractor or cooperator (for the Prisioner´s Dilemma Game), every node is inicialized as a Coopertator and then the Detractors are selected
        'Vulnerability' - a number that represents how vulnerable an individual is towars infection, the greater the number the more vulnerable a node is
        'Utility'       - the amount of utility (well-being) an individual has, we inicialize every node will null utility
        'Danger'        - the amount (numerical) of danger a node has at a given time in the simulation, also inicialized with 0
        'Exposicion'    - the amount (numerical) of exposure a node has at a given time in the simulation, also inicialized with 0
        'Risk'          - the amount of risk

        Input:
                detractor_prob   - (float) The probability that an individual is a detractor

        Output:
                Go               - (networkx graph) graph with inicialized atributes and edge weights 

        '''
        if self.Model == 'SEIR':
            for node in self.Go.nodes():
                self.Go.nodes[node]['state']         = 'S'
                self.Go.nodes[node]['days_S']        = 0
                self.Go.nodes[node]['days_E']        = 0
                self.Go.nodes[node]['days_I']        = 0
                self.Go.nodes[node]['days_R']        = 0

        if self.Model == 'SEIR_Risk':
            for node in self.Go.nodes():
                self.Go.nodes[node]['state']                   = 'S'
                self.Go.nodes[node]['days_S']                  = 0
                self.Go.nodes[node]['days_E']                  = 0
                self.Go.nodes[node]['days_I']                  = 0
                self.Go.nodes[node]['days_R']                  = 0
                self.Go.nodes[node]['Infected_Neighbors']      = 0
                self.Go.nodes[node]['Exposed']                 = 0
                self.Go.nodes[node]['Second_Inf_Neighborhood'] = 0
                self.Go.nodes[node]['Danger']                  = 0
                self.Go.nodes[node]['Risk']                    = 0
                self.Go.nodes[node]['Risk_exp']                = 0

        if self.Model == 'SEIR_Vaccination_Game':
            for node in self.Go.nodes():
                self.Go.nodes[node]['state']                   = 'S'
                self.Go.nodes[node]['days_S']                  = 0
                self.Go.nodes[node]['days_E']                  = 0
                self.Go.nodes[node]['days_I']                  = 0
                self.Go.nodes[node]['days_R']                  = 0
                self.Go.nodes[node]['Infected_Neighbors']      = 0
                self.Go.nodes[node]['Exposed']                 = 0
                self.Go.nodes[node]['Second_Inf_Neighborhood'] = 0
                self.Go.nodes[node]['Danger']                  = 0
                self.Go.nodes[node]['Risk']                    = 0
                self.Go.nodes[node]['M_scale']                 = 0
                self.Go.nodes[node]['Payoff']                  = 0
                self.Go.nodes[node]['Vaccinated_flag']         = 0
                self.Go.nodes[node]['anti-vaccine']            = 0

            # Randomly select an opinion scale according to the M model.
            for nodo in self.Go.nodes():
                self.Go.nodes[nodo]['M_scale'] = random.choice([1, 2])  
            extremistas = random.sample(self.Go.nodes(), int(self.N*(self.initial_antivaccine_percentage/100)) )
            for nodo in extremistas:
                self.Go.nodes[nodo]['M_scale'] = random.choice([-1, -2]) 
                if self.Go.nodes[nodo]['M_scale'] == -2:
                    self.Go.nodes[nodo]['-2_flag'] = 1

        # Edge weight distribution Section 
        for u,v,w in self.Go.edges(data=True):
            s = np.random.randint(self.W_range[0],self.W_range[1],1)  # np.random.uniform(self.W_range[0],self.W_range[1],1)  #randint(self.W_range[0],self.W_range[1],1)
            w['weight'] = s                                           # we asign the number from the distribution to the weight of the edge    

        # A random sample of size 'N_infected' is selected from the node
        seeds  = sample(list(self.Go.nodes()), self.N_infected)

        # Loop through the intial infected to change their state from 'S' to 'I'
        for vi in seeds: self.Go.nodes[vi]['state'] = 'I'

        # the number of days infected for the first graph are initialized as one (day one for the first infected nodes)
        self.Go = self.add_or_reset_state_days(seeds,self.Go,'I','add')

        return self.Go

    def add_or_reset_state_days(self,Nodes_list, Gi, Epi_State, add_or_reset):
        '''
        This Function adds(+1) or resets(=0) the number of days that a node has been in a certain state (Infected or Recovery)
            Input:
                    Nodes_list    - (list) nodes that will change the number whose number of days in the state State will add or reset
                    Epi_State     - (string) state of the number of days that will be updated, the options are 'days_I' or 'days_R'
                    add_or_reset  - (string) 'add' will add a day to the number of days and 'reset' will make it 0
            Output:
                    G                        - (networkx graph) Graph with the nodes days with a certain state updated 
        '''
        if add_or_reset == 'add':
            for vi in Nodes_list:                       # we make a loop through the list of the nodes that want to be changed
                Gi.nodes[vi]['days_'+Epi_State] +=  1   # the states are changed to the input variable "new_state"
        elif add_or_reset == 'reset':
            for vi in Nodes_list:                       # we make a loop through the list of the nodes that want to be reseted
                Gi.nodes[vi]['days_'+Epi_State] = 0     # the number of days on that certain state resets to 0
        return Gi

    def Count_Days_in_State(self,State,t_period,Gc):
        '''
        This function get the list of nodes that have been t_period in the State = S,I,R

            Input:
                State    - (string) the state of the node to be searching ('I' or 'R')
                t_period - (int) the minimum number of time steps that a node has been in State
                Gc       - (networkx graph) the current graph

            Output
                List of nodes in State = 'I' or 'R' with t_period time steps in such state
        '''
        # Loop over the list of infected or recovered nodes
        return [node_i  for node_i in self.List_of_Nodes_in_State(State,Gc) if Gc.nodes[node_i]['days_' + State] >= t_period ]

    def List_of_Nodes_in_State(self,State,Gc):
        '''
        This function get the list of nodes that currently are in State = I or R

            Input:
                State    - (string) the state of the node to be searching ('I' or 'R')
                t_period - (int) the minimum number of time steps that a node has been in State
                Gc       - (networkx graph) the current graph

            Output
                List of nodes in State = 'I' or 'R' with t_period time steps in such state
        '''
        return [ vi for vi,vj in nx.get_node_attributes(Gc, 'state').items()  if vj == State]
    
    def List_of_Nodes_Attribute(self,State,Gc,type_of_state):
        '''
        This function get the list of nodes that currently are in State = I or R

            Input:
                State    - (string) the state of the node to be searching ('I' or 'R')
                t_period - (int) the minimum number of time steps that a node has been in State
                Gc       - (networkx graph) the current graph

            Output
                List of nodes in State = 'I' or 'R' with t_period time steps in such state
        '''
        return [ vi for vi,vj in nx.get_node_attributes(Gc, type_of_state).items()  if vj == State]

    def Next_Infected_Nodes(self,Gc):
        '''
        Function that simulates infections by looping through every infected node and picking a sample from the list of its susceptible neighbors.
        The selction of the sample. The probability of selecting a susceptible neighbor that will get infected is determined by the weight
        of the edge it has with the infected node.In particular this probability is defined as 
        p_ij = (weight of edge(i,j)) /sum(weights of edges between the infected node and its susceptible neighbors) 
        After selecting the sample, the epidemic state of the nodes in the sample are changed from 'S' (susceptible) to 'E' (exposed).

            Input:
                Gc - (networkx graph) Current graph to be updated
            Output:
                G_next - (networkx graph) The updated graph after the epidmic process.
        '''
        G_next = Gc.copy() # Create a copy of the current graph

        # loop through the list of infected nodes that are currently infected and has been in such state over the infectious_period 
        for infected_node_i in self.List_of_Nodes_in_State('I',G_next): #self.Count_Days_in_State('I', self.infectious_period, Gc):

            # Obtain the list of suceptible neighbors of each infected node by intersecting the list of susceptible nodes and the list of the neighbors of an infected node
            susceptible_neighbors    = list( [  vi for vi in  Gc.neighbors(infected_node_i) if  Gc.nodes[vi]['state'] == 'S' ]  )

            # if the infected node has at least one suscpetible neighbor, we obtain the number of susceptible neighbors that will get infected by the infected node in the current iteration
            if len(susceptible_neighbors) != 0:

                # Create a list with the weights of the edges between infected node and its susceptible neighbors
                W =  [ Gc[infected_node_i][node_i]['weight'][0]  for node_i in susceptible_neighbors ] 

                # a susceptible neighbor is selected to get infected with a probability proportional to its edge weight with the infected node
                new_infected = np.random.choice(susceptible_neighbors, p=W/sum(W))

                # Random process to determine if the infected node will infect the selected node. The process change the state of the node to exposed, during the time interval: incubation_period
                if np.random.uniform(0,1) < self.beta:
                    G_next.nodes[new_infected]['state'] = 'E'
                    self.add_or_reset_state_days([new_infected], G_next, 'E', 'reset') # Reset the clock days_StateB for the nodes in the list_of_nodes 

        return G_next
    
    
    def Next_Infected_Nodes_V(self,Gc):
        '''
        Function that simulates infections by looping through every infected node and picking a sample from the list of its susceptible neighbors.
        The selction of the sample. The probability of selecting a susceptible neighbor that will get infected is determined by the weight
        of the edge it has with the infected node.In particular this probability is defined as 
        p_ij = (weight of edge(i,j)) /sum(weights of edges between the infected node and its susceptible neighbors) 
        After selecting the sample, the epidemic state of the nodes in the sample are changed from 'S' (susceptible) to 'E' (exposed).

            Input:
                Gc - (networkx graph) Current graph to be updated
            Output:
                G_next - (networkx graph) The updated graph after the epidmic process.
        '''
        G_next = Gc.copy() # Create a copy of the current graph

        # loop through the list of infected nodes that are currently infected and has been in such state over the infectious_period 
        for infected_node_i in self.List_of_Nodes_in_State('I',G_next): #self.Count_Days_in_State('I', self.infectious_period, Gc):

            # Obtain the list of suceptible neighbors of each infected node by intersecting the list of susceptible nodes and the list of the neighbors of an infected node
            susceptible_and_vaccinated_neighbors    = list( [  vi for vi in  Gc.neighbors(infected_node_i) if  (Gc.nodes[vi]['state'] == 'S' or Gc.nodes[vi]['state'] == 'V') ]  )

            # if the infected node has at least one suscpetible neighbor, we obtain the number of susceptible neighbors that will get infected by the infected node in the current iteration
            if len(susceptible_and_vaccinated_neighbors) != 0:

                # Create a list with the weights of the edges between infected node and its susceptible neighbors
                W =  [ Gc[infected_node_i][node_i]['weight'][0]  for node_i in susceptible_and_vaccinated_neighbors ] 

                # a susceptible neighbor is selected to get infected with a probability proportional to its edge weight with the infected node
                new_infected = np.random.choice(susceptible_and_vaccinated_neighbors, p=W/sum(W))
                
                if G_next.nodes[new_infected]['state'] == 'V':
                    probability_of_infection = self.beta * (1-self.vaccine_effectiveness)
                else:
                    probability_of_infection = self.beta

                # Random process to determine if the infected node will infect the selected node. The process change the state of the node to exposed, during the time interval: incubation_period
                if np.random.uniform(0,1) < probability_of_infection:
                    G_next.nodes[new_infected]['state'] = 'E'
                    self.add_or_reset_state_days([new_infected], G_next, 'E', 'reset') # Reset the clock days_StateB for the nodes in the list_of_nodes 

        return G_next

    def Vaccinated_become_Infected(self,Gc):
        '''
        Function that simulates infections by looping through every infected node and picking a sample from the list of its Vaccinated neighbors.
        The probability of selecting a vaccinated neighbor that will get infected is determined by the weight
        of the edge it has with the infected node. In particular this probability is defined as 
        p_ij = (weight of edge(i,j)) /sum(weights of edges between the infected node and its susceptible neighbors.
        - The vaccine has a certain effectiveness w = self.vac_effective and as a consequence vaccinated nodes can be infected with probability β(1 − ω) if they are in contact with an infected neighbor. 
        After selecting the sample, the epidemic state of the nodes in the sample are changed from 'S' (susceptible) to 'E' (exposed).

            Input:
                Gc - (networkx graph) Current graph to be updated
            Output:
                G_next - (networkx graph) The updated graph after the epidmic process.
        '''

        G_next = Gc.copy() # Create a copy of the current graph

        # loop through the list of infected nodes that are currently infected and has been in such state over the infectious_period 
        for infected_node_i in self.List_of_Nodes_in_State('I',G_next): 

            # Obtain the list of vaccinated neighbors of each infected node.
            vaccinated_neighbors    = list( [  vi for vi in  Gc.neighbors(infected_node_i) if Gc.nodes[vi]['state'] == 'V' ]  )

            # if the infected node has at least one vaccinated neighbor, we obtain the number of vaccinated neighbors that will get infected in the current iteration
            if len(vaccinated_neighbors) != 0:

                # Create a list with the weights of the edges between infected node and its vaccinated neighbors
                W =  [ Gc[infected_node_i][node_i]['weight'][0]  for node_i in vaccinated_neighbors ] 

                # a vaccinated neighbor is selected to get infected with a probability proportional to its edge weight with the infected node
                new_infected = np.random.choice(vaccinated_neighbors, p=W/sum(W))

                # Random process to determine if the infected node will infect the selected vaccinated node. The process change the state 'V' of the node to infected 'I'.
                if np.random.uniform(0,1) < self.beta*(1-self.vaccine_effectiveness):
                    
                    # The vaccinated agent 'new_infected' change its epidemic state from 'V' to 'E'.
                    G_next.nodes[new_infected]['state'] = 'E'

                    # Reset the clock days_StateB for the nodes in the list_of_nodes 
                    self.add_or_reset_state_days([new_infected], G_next, 'E', 'reset') 

                    # the new infected agent changes his opinion from extremist positive to totally against the vaccine.
                    G_next.nodes[new_infected]['M_scale'] = -2

        return G_next

    def from_StateA_to_StateB(self,StateA, StateB, time_period, Gc):
        '''
        Function to change the state of all the agents whose epidemic state (S,E,V,I or R), has been StateA during the interval time_period in the current graph Gc, and change it to the StateB.

            Input:
                StateA      - (string) the state of the node to be searching ( 'S', 'E', 'I', 'R' or 'V')
                StateB      - (string) the state of the node to be changing  ( 'S', 'E', 'I', 'R' or 'V')
                time_period - (int) the number of iteration that agents have been in state StateA
                Gc          - (networkx graph) Current graph to be updated
        '''

        list_of_nodes = self.Count_Days_in_State(StateA,time_period,Gc)

        if len(list_of_nodes) != 0:
            for node_i in list_of_nodes:
                    
                Gc.nodes[node_i]['state'] = StateB

            # Reset the clock days_StateB for the nodes in the list_of_nodes 
            self.add_or_reset_state_days(list_of_nodes, Gc, StateB, 'reset')
            
        return Gc
    
    def from_StateR_to_StateS(self,StateA,StateB,time_period,Gc):

        list_of_nodes = self.Count_Days_in_State(StateA,time_period,Gc)

        if len(list_of_nodes) != 0:
            for node_i in list_of_nodes:
                if Gc.nodes[node_i]['Vaccinated_flag']>0:
                    # Gc.nodes[node_i]['Vaccinated_flag']+=1
                    Gc.nodes[node_i]['state'] = 'V'
                    continue
                    
                Gc.nodes[node_i]['state'] = StateB

            # Reset the clock days_StateB for the nodes in the list_of_nodes 
            self.add_or_reset_state_days(list_of_nodes, Gc, StateB, 'reset')
            
        return Gc

    def Neighbors_in_State(self,State,vi,Gc):
        '''
        Given the label of the node vi, this function returns the set of neighbors, in the first shell, in State = 'S','E', 'I' or 'R'

            Inputs:
                    State - (string) the state of the node:  'S','E', 'I' or 'R'
                    vi    - (int)    the node lable
                    Gc    - (networkx graph) the current graph
            Output:
                    List of nodes that are neighbors of node i and are in state 'S','E', 'I' or 'R'
        '''
        return list(  set( self.List_of_Nodes_in_State(State,Gc) ) & set( Gc.neighbors(vi) ) )

    def Danger(self,Gc):
        
        G_next = Gc.copy() #Create a copy of the current graph
        
        for node_i in G_next.nodes():
            Hi                                               = sum( [ len(self.Neighbors_in_State('I',vj,Gc)) / Gc.degree(vj)  for vj in self.Neighbors_in_State('S',node_i,Gc)  ]) #self.Second_Infected_Neighborhood(node_i,Gc)
            G_next.nodes[node_i]['Second_Inf_Neighborhood']  = Hi
            G_next.nodes[node_i]['Danger']                   = 1-np.exp( -Hi )
        return G_next
    
    def Exposed(self,Gc,J=1,alpha=1):
        
        G_next = Gc.copy() #Create a copy of the current graph
        
        for node_i in G_next.nodes():    
            Hj                                         = len( self.Neighbors_in_State('I',node_i,Gc)) /  Gc.degree(node_i)
            G_next.nodes[node_i]['Infected_Neighbors'] = Hj
            G_next.nodes[node_i]['Exposed']            = 1-np.exp( - J*pow(Hj,alpha) )
        return G_next

    def Risk(self,Gc,J=1,alpha=1):
        
        '''
            Function to asses the risk of each node in the graph Gc.

            Input:
                Gc - (networkx graph). The current graph to be updated in the atributed Risk.
            Output:
                G_next - (networkx graph). The updated graph with the current value of Risk.
        '''

        G_next = Gc.copy() #Create a copy of the current graph
        
        for node_i in G_next.nodes():
            G_next.nodes[node_i]['Risk']     = self.beta*(1- (1-G_next.nodes[node_i]['Danger'])*(1-G_next.nodes[node_i]['Exposed']))
            G_next.nodes[node_i]['Risk_exp'] = self.beta*np.exp(-J*pow( len( self.Neighbors_in_State('I',node_i,Gc)) /  Gc.degree(node_i),alpha) )
        return G_next

    def Risk_V(self,Gc,J=1,alpha=1):
        '''
            Function to asses the risk of each node in the graph Gc. - For the Vaccination Game

            Input:
                Gc - (networkx graph). The current graph to be updated in the atributed Risk.
            Output:
                G_next - (networkx graph). The updated graph with the current value of Risk.
        '''

        G_next = Gc.copy() #Create a copy of the current graph
        
        for node_i in G_next.nodes():
            G_next.nodes[node_i]['Risk'] = self.beta*(1- (1-G_next.nodes[node_i]['Danger'])*(1-G_next.nodes[node_i]['Exposed']))

        return G_next

    def Update_Epidemic_SEIR(self,G_array):

        '''
        Function to update the epidemic process in a single iteration. The process is the following:

            1.- Takes the last graph of G_array and call the function Next_Infected_nodes to emulate the contagious process.
            2.- Update the state of exposed nodes that satisfy days_E = incubation_period.
            3.- Update the state of infected nodes that satisfy days_I = infectious_period.
            4.- Update the state of recovered nodes that satisfy days_I = inmmunity_period.
            5.- Update the clocks for the states S, E, I and R
            6.- Append the updated graph to G_array

        Input:
            G_array - (array) the set of graphs updated at each iteration
        Output:
            G_array - (array) the input G_array with the new networkx graph, where nodes states are updated. 
        '''

        # Epidemic process (from susceptible to infected)
        G_next = self.Next_Infected_Nodes(G_array[-1])

        # Exposed nodes with t_period = incubation_period, change its state to I
        G_next = self.from_StateA_to_StateB('E','I',self.incubation_period,G_next)

        # Infected nodes with t_period = infectious_period, change its state to R
        G_next = self.from_StateA_to_StateB('I','R',self.infectious_period,G_next)

        # Recovered nodes, with t_period = inmune_period, change its state to S
        # G_next = self.from_StateR_to_StateS('R','S',self.inmmunity_period,G_next)
        G_next = self.from_StateA_to_StateB('R','S',self.inmmunity_period,G_next)

        # we add a day to the atribute 'days_E' to all the nodes currently exposed at the end of the period
        G_next = self.add_or_reset_state_days(self.List_of_Nodes_in_State('S',G_next), G_next, 'S', 'add')

        # we add a day to the atribute 'days_E' to all the nodes currently exposed at the end of the period
        G_next = self.add_or_reset_state_days(self.List_of_Nodes_in_State('E',G_next), G_next, 'E', 'add')

        # we add a day to the atribute 'days_I' to all the nodes currently infected at the end of the period
        G_next = self.add_or_reset_state_days(self.List_of_Nodes_in_State('I',G_next), G_next, 'I', 'add')

        # we add a day to the atribute 'days_R' to all the nodes currently recovered at the end of the period
        G_next = self.add_or_reset_state_days(self.List_of_Nodes_in_State('R',G_next), G_next, 'R', 'add')

        # the Graph is updated by changing the state of the new infected and the recovered nodes
        G_array.append(G_next)

        return G_array

    def Update_Epidemic_SEIR_Risk(self,G_array,J=1,alpha=1):

        '''
        Function to update the epidemic process in a single iteration. The process is the following:

            1.- Takes the last graph of G_array and call the function Next_Infected_nodes to emulate the contagious process.
            2.- Update the state of exposed nodes that satisfy days_E = incubation_period.
            3.- Update the state of infected nodes that satisfy days_I = infectious_period.
            4.- Update the state of recovered nodes that satisfy days_I = inmmunity_period.
            5.- Update the clocks for the states S, E, I and R
            6.- Append the updated graph to G_array

        Input:
            G_array - (array) the set of graphs updated at each iteration
        Output:
            G_array - (array) the input G_array with the new networkx graph, where nodes states are updated. 
        '''

        # Epidemic process (from susceptible to exposed)
        G_next = self.Next_Infected_Nodes(G_array[-1])

        # Exposed nodes with t_period = incubation_period, change its state to I
        G_next = self.from_StateA_to_StateB('E','I',self.incubation_period,G_next)

        # Infected nodes with t_period = infectious_period, change its state to R
        G_next = self.from_StateA_to_StateB('I','R',self.infectious_period,G_next)

        # Recovered nodes, with t_period = inmune_period, change its state to S
        # G_next = self.from_StateR_to_StateS('R','S',self.inmmunity_period,G_next)
        G_next = self.from_StateA_to_StateB('R','S',self.inmmunity_period,G_next)

        # we add a day to the atribute 'days_E' to all the nodes currently exposed at the end of the period
        G_next = self.add_or_reset_state_days(self.List_of_Nodes_in_State('S',G_next), G_next, 'S', 'add')

        # we add a day to the atribute 'days_E' to all the nodes currently exposed at the end of the period
        G_next = self.add_or_reset_state_days(self.List_of_Nodes_in_State('E',G_next), G_next, 'E', 'add')

        # we add a day to the atribute 'days_I' to all the nodes currently infected at the end of the period
        G_next = self.add_or_reset_state_days(self.List_of_Nodes_in_State('I',G_next), G_next, 'I', 'add')

        # we add a day to the atribute 'days_R' to all the nodes currently recovered at the end of the period
        G_next = self.add_or_reset_state_days(self.List_of_Nodes_in_State('R',G_next), G_next, 'R', 'add')

        # For each node, assess the danger.
        G_next = self.Danger(G_next)

        # For each node, assess the exposition.
        G_next = self.Exposed(G_next, J, alpha)

        # For each node, assess the risk.
        G_next = self.Risk(G_next, J, alpha)
        G_array.append(G_next)
        
        
        return G_array

    def Update_Epidemic_SEIR_Vaccination_Game(self,G_array, J=1,alpha=1):

        '''
        Function to update the epidemic process in a single iteration. The process is the following:

            1.- Takes the last graph of G_array and call the function Next_Infected_nodes to emulate the contagious process.
            2.- The vaccinated nodes become infected with probability (1-w)*beta, with w the vaccine effectiveness.
            2.- Update the state of exposed nodes that satisfy days_E = incubation_period.
            3.- Update the state of infected nodes that satisfy days_I = infectious_period.
            4.- Update the state of recovered nodes that satisfy days_I = inmmunity_period.
            5.- Update the clocks for the states S, E, I and R
            6.- Append the updated graph to G_array

        Input:
            G_array - (array) the set of graphs updated at each iteration
        Output:
            G_array - (array) the input G_array with the new networkx graph, where nodes states are updated. 
        '''

        # Epidemic process (from susceptible to exposed)
        G_next = self.Next_Infected_Nodes_V(G_array[-1])

        # Vaccinated nodes becomes Infected with probability (1-w)*beta, with w the vaccine effectiveness
        # G_next = self.Vaccinated_become_Infected(G_next)

        # Exposed nodes with t_period = incubation_period, change its state to I
        G_next = self.from_StateA_to_StateB('E','I',self.incubation_period,G_next)

        # Infected nodes with t_period = infectious_period, change its state to R
        G_next = self.from_StateA_to_StateB('I','R',self.infectious_period,G_next)

        # Recovered nodes, with t_period = inmune_period, change its state to S
        G_next = self.from_StateR_to_StateS('R','S',self.inmmunity_period,G_next)

        # we add a day to the atribute 'days_S' to all the nodes currently exposed at the end of the period
        G_next = self.add_or_reset_state_days(self.List_of_Nodes_in_State('S',G_next), G_next, 'S', 'add')

        # we add a day to the atribute 'days_E' to all the nodes currently exposed at the end of the period
        G_next = self.add_or_reset_state_days(self.List_of_Nodes_in_State('E',G_next), G_next, 'E', 'add')

        # we add a day to the atribute 'days_I' to all the nodes currently infected at the end of the period
        G_next = self.add_or_reset_state_days(self.List_of_Nodes_in_State('I',G_next), G_next, 'I', 'add')

        # we add a day to the atribute 'days_R' to all the nodes currently recovered at the end of the period
        G_next = self.add_or_reset_state_days(self.List_of_Nodes_in_State('R',G_next), G_next, 'R', 'add')

        # For each node, assess the danger.
        G_next = self.Danger(G_next)

        # For each node, assess the exposition.
        G_next = self.Exposed(G_next, J, alpha)

        # For each node, assess the risk.
        G_next = self.Risk_V(G_next, J, alpha)

        if self.strategy == 'mayority':
            # Each node gets its utility by playing with all its neighbors
            G_next = Game.Update_Binary_Opinion_Mayority_Choice(self,G_next) #Update_Payoff(self,G_next)

        if self.strategy == 'replication':
            G_next = Game.Update_Binary_Opinion_Replication_Choice(self,G_next)
        
        # Each node gets its utility by playing with all its neighbors
        #G_next = self.Opinion_Changing_Mayority_Choice(G_next)             # --------- AQUI
        # G_next = self.Opinion_Changing_Deterministic_Replication(G_next)
        # G_next = self.Opinion_Changing_Probabilistic_Replication(G_next)
        
        # Each node who is convinced by the vaccine, decides wether to take the vaccine 
        G_next = Game.Vaccination_Decision(self,G_next)
        # G_next = self.Vaccinated_become_Infected(G_next)
        # the Graph is updated by changing the state of the new infected and the recovered nodes
        G_array.append(G_next)

        return G_array

    def Epidemic_Evolution(self):
    
        """ Function to simulate the course of the epidemic, game and coevolution of the network over a period of time. 
    
            Input: 
                     
            Output:
                Gtrack     - (list) array in which every entry is the updated graph en each step of the epidemic
        """

        self.Go = self.Initial_Graph_Atributes()

        # initialize the array of graphs with the initial graph(introduced in the input)
        Gtrack = [self.Go]   

        if self.Model == 'SEIR':
            for step_i in range(1,self.steps):
                # Simulate the epidemic and update the epidemic state of the nodes in the current iteration
                Gtrack = self.Update_Epidemic_SEIR(Gtrack)

        if self.Model == 'SEIR_Risk':
            for step_i in range(1,self.steps):
            	# Simulate the epidemic and update the epidemic state of the nodes in the current iteration
            	Gtrack = self.Update_Epidemic_SEIR_Risk(Gtrack)

        if self.Model == 'SEIR_Vaccination_Game':
            for step_i in range(1,self.steps):
                Gtrack = self.Update_Epidemic_SEIR_Vaccination_Game(Gtrack)


        return Gtrack
    
class Game(Epidemic_Model):

    def __init__(self):
        return

    def Payoff_Opinion_Preference_Game(self,gamer_i,gamer_j,Gc): 
        '''
        Given the current graph, this function gets the payoff of gamer_i when play with gamer_j based on the opinion preference game, whose utility matrix is:

                                                    Opinion A            Opinion B

        Opinion A (in favor of vaccines)             C + x                 x

        Opion B (against vaccines)                     0                   C

        That is, this game positively rewards players with the same opinion and benefit opinion A.

        Inputs:
            gamer_i - (int) the label of gamer 1
            gamer_j - (int) the label of gamer 2
            Gc      - ((networkx graph) Current graph to be updated.

        Output:
            G_c     - (networkx graph) The updated graph after the opinion formation process.

        '''
        # If both gamers have the same opinion
        if (Gc.nodes[gamer_i]['M_scale']  > 0  and Gc.nodes[gamer_j]['M_scale'] > 0 ):
            Gc.nodes[gamer_i]['Payoff'] += self.C_utilidad + self.x_utilidad

        # If the gamer_i opinion is positive and the opinion of his neighbor is negative  
        elif (Gc.nodes[gamer_i]['M_scale'] > 0 and Gc.nodes[gamer_j]['M_scale'] < 0):
            Gc.nodes[gamer_i]['Payoff'] += self.x_utilidad

        # If the gamer_i opinion is negative and the opinion of his neighbor is negative
        elif (Gc.nodes[gamer_i]['M_scale'] < 0 and Gc.nodes[gamer_j]['M_scale'] < 0):
            Gc.nodes[gamer_i]['Payoff'] +=   self.C_utilidad 

        #  If the gamer_i opinion is negative and the opinion of his neighbor is positive 
        else:
            Gc.nodes[gamer_i]['Payoff'] += 0

        return

    def Update_Payoff(self,Gc):
        '''
        Function to update the payoff of agents according the utility matrix.

        Input:
            Gc             - (networkx graph) Current graph to be updated
            utility_matrix - (string) the name of the utility matrix using in the game

        Output:
            G_next - (networkx graph) The updated graph after the epidmic process.
        '''

        G_next = Gc.copy()

        if self.utility_matrix == 'Payoff_Opinion_Preference_Game':

            # Loop over all the agents (here called gamers)
            for gamer_i in G_next.nodes():

                # Reset the payoff of gamer_i
                G_next.nodes[gamer_i]['Payoff'] = 0

                # Loop over the neighbors of gamer_i and play the game
                for neighbor_j in G_next.neighbors(gamer_i):
                    Game.Payoff_Opinion_Preference_Game(self,gamer_i,neighbor_j,G_next)
                
        return G_next

    def List_of_Nodes_Against_Vaccines(self,Gc):
        '''
        This function get the list of nodes that currently have the opinion scale M != 2 in the graph Gc.

            Input:
                Opinion_Scale   - (int) - the scale of the nodes M = -2,-1,1 or 2.
                Gc              - (networkx graph) the current graph

            Output
                List of nodes with opinion M.
        '''
        return [ vi for vi,vj in nx.get_node_attributes(Gc,'M_scale').items()  if vj != 2 ]

    def Update_Binary_Opinion_Mayority_Choice(self,Gc):
        '''
        Function to update the opinion of the agents according to the mayority choice method. Here, 
        '''

        # Create a copy of the graph Gc to be updated with the next agent opinions.
        G_next   = Gc.copy()

        if self.utility_matrix == 'Payoff_Opinion_Preference_Game':

            for gamer_i in Game.List_of_Nodes_Against_Vaccines(self,Gc):

                # For each gamer in Gc, an auxiliar graph is created from which the M scale of such gamer is changed in oder to obtain its payoff if it is opossite to the current opinion.
                G_opp = Gc.copy()

                # For gamer_i, its opinion is changing in the temporal and auxiliar graph G_opposite_opinion
                G_opp.nodes[gamer_i]['M_scale'] = G_opp.nodes[gamer_i]['M_scale'] * (-1)

                # Reset the payoff of gamer_i in G_next
                G_next.nodes[gamer_i]['Payoff'] = 0

                # Reset the payoff of gamer_i in G_opposite
                G_opp.nodes[gamer_i]['Payoff']  = 0

                # Loop over the neighbors of gamer_i and play the game
                for neighbor_j in G_next.neighbors(gamer_i):
                    # Update gamer_i payoff of and its neighbor in G_next
                    Game.Payoff_Opinion_Preference_Game(self,gamer_i,neighbor_j,G_next)

                    # Update gamer_i payoff and its neighbro in G_opp
                    Game.Payoff_Opinion_Preference_Game(self,gamer_i,neighbor_j,G_opp)

                # Process to update the gamer_i opinion in the scale M = -2,-1, 1, 2.
                # If the payoff of think opposite is bigger than the current opinion scale.  
                if G_opp.nodes[gamer_i]['Payoff'] > G_next.nodes[gamer_i]['Payoff']:
                    
                    # If M = -1 or 1, the M scale change to 1 or -1, respectively
                    if abs(G_next.nodes[gamer_i]['M_scale']) == 1:
                        G_next.nodes[gamer_i]['M_scale'] *= (-1)
                    # If M = -2, it change to M = -1
                    else:
                        G_next.nodes[gamer_i]['M_scale'] = -1
                # If the payoff of perserve the current opinion scale is bigger than the opinion is opposite.
                else:
                    if G_next.nodes[gamer_i]['M_scale'] == 1:
                        G_next.nodes[gamer_i]['M_scale'] = 2
                    else: 
                        G_next.nodes[gamer_i]['M_scale'] = -2

                # Falta borrar G_opp.del()

        return G_next

    def Update_Binary_Opinion_Replication_Choice(self,Gc):

        G_next = Gc.copy()
        for gamer_i in Game.List_of_Nodes_Against_Vaccines(self,Gc):

            # List of neighbors with a current payoff bigger than gamer_i
            Neighbors_with_biggest_payoff = [neighbor_j for neighbor_j in G_next.neighbors(gamer_i) if G_next.nodes[neighbor_j]['Payoff'] > G_next.nodes[gamer_i]['Payoff'] ]

            if len(Neighbors_with_biggest_payoff) != 0:
                
                # Chose a neighbor
                chosen_neighbor = random.choice(Neighbors_with_biggest_payoff)

                # Copy its opinion
                G_next.nodes[gamer_i]['M_scale'] = G_next.nodes[chosen_neighbor]['M_scale']
            else:
                continue
        
        return G_next

    def Vaccination_Decision(self,Gc):

        G_next = Gc.copy()

        candidates_to_vaccinate =  self.List_of_Nodes_in_State('S',G_next) + self.List_of_Nodes_in_State('R',G_next)
        
        for gamer_i in G_next.nodes():
            if G_next.nodes[gamer_i]['M_scale'] == 2 and  gamer_i in candidates_to_vaccinate:
                if self.vaccination_cost < G_next.nodes[gamer_i]['Risk']*(self.vaccine_effectiveness): 
                    G_next.nodes[gamer_i]['state'] = 'V'
                    G_next.nodes[gamer_i]['Vaccinated_flag'] = 1
        return G_next

class Plot_Data(Epidemic_Model):

    def __init__(self,G_path,fontsize = 15, linewidth=2,labelsize=15,alpha=0.6,formato='png'):
        
        self.G_path    = G_path

        self.fontsize  = fontsize
        self.linewidth = linewidth
        self.labelsize = labelsize
        self.alpha     = alpha
        self.formato   = formato

        return

    def Infected(self, N_infected, run_i, title, label,color = 'r', fig_size = (8,6)):
        '''
        Function to plot the number of infected nodes at each time step
        Input:
                run_i  - (int) number of simulation run
                title  - (str) title of the graph
            
        Output:
            prints a graph with the x-axis and y-axis representing time and current infections over time respectively
        '''

        #N_infected = []         # initialize the array with the current infected in every period of time
        #for Gi in self.G_path:
        #    # the number of current infected in each time step
        #    N_infected.append( len(  self.List_of_Nodes_in_State('I',Gi) ) )

        # the number of current infected in each time step
        #N_infected = [   len(  self.List_of_Nodes_in_State('I',Gi) ) for Gi in self.G_path ]
            
        plt.clf()
        plt.figure(1000+run_i,figsize = fig_size)
        plt.plot(N_infected,'o-', color = color, linewidth=self.linewidth, alpha = self.alpha,label=label)
        
        plt.title(title)
        
        plt.xlabel("Tiempo (días)",fontsize=self.fontsize)
        plt.ylabel("Número de casos diarios",fontsize=self.fontsize)

        plt.tick_params(axis='y',labelsize=self.labelsize)
        plt.tick_params(axis='x',labelsize=self.labelsize)

        plt.legend(loc=2, bbox_to_anchor=(0.1,0.98),fontsize=10,ncol=1, fancybox=False, shadow=False)

        plt.savefig('Daily_I_Run' + str(run_i) + '.' + self.formato)
        
        return

    def Vaccinated(self, N_vaccinated, run_i, title, label,color = 'r', fig_size = (8,6)):
        '''
        Function to plot the number of infected nodes at each time step
        Input:
                run_i  - (int) number of simulation run
                title  - (str) title of the graph
            
        Output:
            prints a graph with the x-axis and y-axis representing time and current infections over time respectively
        '''
            
        plt.clf()
        plt.figure(1000+run_i,figsize = fig_size)
        plt.plot(N_vaccinated,'o-', color = color, linewidth=self.linewidth, alpha = self.alpha,label=label)
        
        plt.title(title)
        
        plt.xlabel("Tiempo (días)",fontsize=self.fontsize)
        plt.ylabel("Número de casos diarios",fontsize=self.fontsize)

        plt.tick_params(axis='y',labelsize=self.labelsize)
        plt.tick_params(axis='x',labelsize=self.labelsize)

        plt.legend(loc=2, bbox_to_anchor=(0.1,0.98),fontsize=10,ncol=1, fancybox=False, shadow=False)

        plt.savefig('Daily_V_Run' + str(run_i) + '.' + self.formato)
        
        return
    
    def Extremists_Opinion(self, N_Extremists, run_i, title, label, Dayly_Extremists, fig_size = (8,6)):
        '''
        Function to plot the number of infected nodes at each time step
        Input:
                run_i  - (int) number of simulation run
                title  - (str) title of the graph
            
        Output:
            prints a graph with the x-axis and y-axis representing time and current infections over time respectively
        '''
            
        plt.clf()
        plt.figure(1000+run_i,figsize = fig_size)
        color = iter(cm.rainbow(np.linspace(0, 1, len(N_Extremists))))
        contador_df = 0
        for node_e in N_Extremists:
            Opinion_node_e = []
            contador_df +=1
            c = next(color)
            for G in self.G_path:
                Opinion_node_e.append(G.nodes[node_e]['M_scale'])
            Dayly_Extremists['Nodo'+str(contador_df)] = Opinion_node_e
            plt.plot(Opinion_node_e,'o-', color = c, linewidth=self.linewidth, alpha = self.alpha,label=str(node_e))
        
        plt.title(title)
        
        plt.xlabel("Time",fontsize=self.fontsize)
        plt.ylabel(" number of people with extremists opinions",fontsize=self.fontsize)

        plt.tick_params(axis='y',labelsize=self.labelsize)
        plt.tick_params(axis='x',labelsize=self.labelsize)

        plt.legend(loc=2, bbox_to_anchor=(0.1,0.98),fontsize=10,ncol=1, fancybox=False, shadow=False)

        plt.savefig('Daily_Initial_Ex_Run' + str(run_i) + '.' + self.formato)
        
        return
    
    
    def Extremists(self, N_Extremists,antivaccine_moderate, provaccine_moderate, provaccine_extremist, run_i, title, fig_size = (8,6)):
        '''
        Function to plot the number of infected nodes at each time step
        Input:
                run_i  - (int) number of simulation run
                title  - (str) title of the graph
            
        Output:
            prints a graph with the x-axis and y-axis representing time and current infections over time respectively
        '''
            
        plt.clf()
        plt.figure(1000+run_i,figsize = fig_size)
        plt.plot(N_Extremists,'o-', color = 'red', linewidth=self.linewidth, alpha = self.alpha,label='antivaccine extremists')
        plt.plot(antivaccine_moderate,'o-', color = 'orange', linewidth=self.linewidth, alpha = self.alpha,label='antivaccine moderate')
        plt.plot(provaccine_moderate,'o-', color = 'blue', linewidth=self.linewidth, alpha = self.alpha,label='provaccine moderate')
        plt.plot(provaccine_extremist,'o-', color = 'green', linewidth=self.linewidth, alpha = self.alpha,label='provaccine extremists')
        
        plt.title(title)
        
        plt.xlabel("Time",fontsize=self.fontsize)
        plt.ylabel(" number of people with extremists opinions",fontsize=self.fontsize)

        plt.tick_params(axis='y',labelsize=self.labelsize)
        plt.tick_params(axis='x',labelsize=self.labelsize)

        plt.legend(loc=2, bbox_to_anchor=(0.1,0.98),fontsize=10,ncol=1, fancybox=False, shadow=False)

        plt.savefig('Daily_Ex_Run' + str(run_i) + '.' + self.formato)
        
        return

class Save_Info(object):
    """ docstring for Save_Info"""
    
    def __init__(self):
        return

    def MakeDir(self,path_name):

        if not os.path.exists(path_name):
            os.makedirs(path_name)
        else:
            shutil.rmtree(path_name)
            os.makedirs(path_name)

        return

    def SaveFiles(self,path_name,make_dir=False,file_format='*.hdf5'):
    
        """
          Routine to Copy all the files with a given format into the folder path_name

          path_name: The destination folder
          format:    the files format that will be copy to path_name

        """
    
        if make_dir == True:
            self.MakeDir(path_name)
        else: pass

        #formato1 = '*.hdf5'
        #files1 = glob.iglob(os.path.join(file_format))
        #for ifile in files1: shutil.move(ifile,path_name)

        # Copy  the files in format file_format to path_name
        for ifile in glob.iglob(os.path.join(file_format)):
            # Si el archivo for formato 'file_format' ya existe, lo elimina y guarda el nuevo
            if os.path.isfile(ifile):
                try:
                    shutil.move(ifile,path_name)
                except:
                    os.remove(path_name+'/'+ifile)
                    shutil.move(ifile,path_name)

        return

# class Menu():

#     def __init__(self):  
#         return

#     def main_menu():

#         os.system("clear")
        
#         print('\n \t\t\t\t =================[ MAIN MENU ]=================')
#         #print('\t\t\t\t\t\t'+'-*'*20)

#         print()
#         print('\t\t\t\t\t Step 1:  Initial configuarion: \n')
#         Number_of_nodes    = int(input(""" Enter the number of nodes: """))
        
#         print()
        
#         Number_of_infected = int(input(""" Enter the number of inital infected nodes: """))

#         while Number_of_infected > Number_of_nodes:
#             print(''' The number of initial infected nodes must be less than the total number of nodes ''')
#             Number_of_infected = int(input(""" Enter the number of inital infected nodes: """))

#         print()

#         print('\t\t\t\t\t Step 2:  Select the graph generator: \n')
#         choice_graph =  input("""

#                           1: Regular graph
#                           2: Random graph
#                           3: Small-World graph
#                           4: Scale-free graph

#                           Please enter your choice: """
#                         )

#         if choice_graph == '1' :

#             print('\n \t\t\t\t ===========[ Regular graph menu ]=============')

#             choice_regular_graph = input("""

#                           Step 3:  Select the regular graph generator: 
#                           1: Complete graph
#                           2: Nearest neighbor graph
#                           3: Cycle graph
#                           4: Star graph

#                           Please enter your choice: 
#                           """)

#             if  choice_regular_graph == '1':
#                 Go = nx.complete_graph(Number_of_nodes)
#             elif choice_regular_graph == '2':
#                 k = int( input("Please indicate the number of neighbors: ") )
#                 Go = nx.watts_strogatz_graph(Number_of_nodes,k,0.)
#             elif choice_regular_graph == '3':
#                 Go = nx.cycle_graph(Number_of_nodes)
#             elif choice_regular_graph == '4':
#                 Go = nx.star_graph(Number_of_nodes)

#         elif choice_graph == '2' :

#             print('\n \t\t\t\t ===========[ Random graph menu ]=============')

#             choice_random_graph = input("""

#                           Step 3:  Select the random graph generator: 
#                           1: Binomial graph
#                           2: Erdos-Renyi graph

#                           Please enter your choice: 
#                           """)

#             print()
#             prob = float(input(""" Enter the probability for edge creation: """))

#             while (prob < 0.) or  (prob > 1.):
#                 print('Please choose a value from 0 to 1 \n')
#                 prob = float(input(""" Enter the probability for edge creation: """))

#             # Creación del grafo
#             if choice_random_graph == '1':
#                 Go = nx.binomial_graph(Number_of_nodes,prob) 
#             elif choice_random_graph == '2':
#                 Go = nx.erdos_renyi_graph(Number_of_nodes,prob) 

#         elif choice_graph == '3' :

#             print('\n \t\t\t\t ===========[ Small-Worldgraph menu ]=============')

#             choice_sw_graph = input("""

#                           Step 3:  Select the Small-World graph generator: 
#                           1: Watts-Strogatz
#                           2: Newman-Watts-Strogatz

#                           Please enter your choice: 
#                           """)
            
#             print()
#             k = int( input("Please indicate the number of neighbors: ") )
#             print()
#             prob = float(input(""" Enter the probability for edge creation: """))

#             while (prob < 0.) or  (prob > 1.):
#                 print('Please choose a value from 0 to 1 \n')
#                 prob = float(input(""" Enter the probability for edge creation: """))

#             # Creación del grafo
#             if choice_sw_graph == '1':
#                 Go = nx.watts_strogatz_graph(Number_of_nodes,k,prob)
#             elif choice_sw_graph == '2':
#                 Go = nx.newman_watts_strogatz_graph(Number_of_nodes,k,prob)

#         elif choice_graph == '4' :

#             print('\n \t\t\t\t ===========[ Scale-Free graph menu ]=============')

#             choice_sf_graph = input("""

#                           Step 3:  Select the Scale-Free graph generator: 
#                           1: Barabasi-Albert graph

#                           Please enter your choice: 
#                           """)

#             print()
#             m = int( input("Please indicate the number of edges to attach from a new node to existing nodes: ") )
#             print()
#             if choice_sf_graph == '1':
#                 Go = nx.barabasi_albert_graph(Number_of_nodes,m)
#                 plt.figure(1)
#                 nx.draw_circular(Go)
#                 plt.show()

#         return 

# Obtain the list of suceptible neighbors of each infected node by intersecting the list of susceptible nodes and the list of the neighbors of an infected node
#susceptible_neighbors = list(  set( self.List_of_Nodes_in_State('S',Gc) ) & set( Gc.neighbors(infected_node_i) ) )

# Exposed nodes with t_period = incubation_period, change its state to I
#        for exposed_i in self.Count_Days_in_State('E',self.incubation_period,G_next):
#            G_next.nodes[exposed_i]['state'] = 'I'
#            self.add_or_reset_state_days([exposed_i], G_next, 'I', 'reset')

#        # Infected nodes with t_period = infectious_period, change its state to R
#        for infected_i in self.Count_Days_in_State('I',self.infectious_period,G_next):
#            G_next.nodes[infected_i]['state'] = 'R'
#            self.add_or_reset_state_days([infected_i], G_next, 'R', 'reset')

#        # Recovered nodes, with t_period = inmune_period, change its state to S
#        for recovered_i in self.Count_Days_in_State('R',self.inmmunity_period,G_next):
#            G_next.nodes[recovered_i]['state'] = 'S'
#            self.add_or_reset_state_days([recovered_i], G_next, 'R', 'reset')

#    def WI(self,vj,Gc):
        '''
        Function to assess the sum of weights of node vj in the neighborhood of suceptible nodes.
        Input:
            vj - (int) the label of the node
            GC - (networkx graph) the current graph
        Output:
            (float) - the sum of weights of node vj in the neighborhood of suceptible nodes.
        '''
#        return sum( [ Gc[vj][vk]["weight"][0]  for vk in  Gc.neighbors(vj) ]  ) # self.Neighbors_in_State('S',vj,Gc)

#    def Danger(self,Gc):
        '''
        Function to asses the danger of each node in the graph Gc.

        Input:
            Gc - (networkx graph). The current graph to be updated in the atributed Danger.
        Output:
             G_next - (networkx graph). The updated graph with the current value of Danger.

        '''
#        G_next = Gc.copy() # Create a copy of the current graph

        # For susceptibles nodes in the current network Gc
#        for node_i in self.List_of_Nodes_in_State('S',Gc):
        	#suma = sum([ Gc[node_i][node_j]["weight"][0]/self.WI(node_j,Gc) for node_j in self.Neighbors_in_State('I',node_i,Gc) ])
        	#print('-- Danger node ' + str(node_i) + '= ', suma,'--',self.beta*suma,'===',self.beta)
#        	G_next.nodes[node_i]['Danger'] = self.beta*sum([ Gc[node_i][node_j]["weight"][0]/self.WI(node_j,Gc) for node_j in self.Neighbors_in_State('I',node_i,Gc) ])  #/ Gc.degree(node_i)
        	
        #for node_i in G_next.nodes():

        #	# If the node is surrounded just by infected neighbors, the probability of been ifected (the danger) is the maximum value, i.e. one
        #	#if len( self.Neighbors_in_State('S',node_j,Gc) ) != 0  #len(self.Neighbors_in_State('I',node_i,Gc)) ==  Gc.degree(node_i):
        #	#	G_next.nodes[node_i]['Danger'] = 1.
        #	#else:
        #	#	G_next.nodes[node_i]['Danger'] = self.beta*sum([ Gc[node_i][node_j]["weight"]/self.WI(node_j,Gc) for node_j in self.Neighbors_in_State('I',node_i,Gc) ])

        #    suma = 0
        #    for node_j in self.Neighbors_in_State('I',node_i,Gc):
        #        if len( self.Neighbors_in_State('S',node_j,Gc) ) != 0:
        #            suma += Gc[node_i][node_j]["weight"][0]/self.WI(node_j,Gc)
        #        else:
        #            suma = 1

        #    #G_next.nodes[node_i]['Danger'] = self.beta*suma

        #G_next.nodes[node_i]['Danger'] = self.beta*sum([ Gc[node_j][node_i]["weight"]/self.WI(node_j,Gc) for node_j in self.Neighbors_in_State('I',node_i,Gc) ])
#        return G_next

#    def Risk(self,Gc,J=1,alpha=1):
        '''
        Function to asses the risk of each node in the graph Gc.

        Input:
            Gc - (networkx graph). The current graph to be updated in the atributed Risk.
        Output:
             G_next - (networkx graph). The updated graph with the current value of Risk.

        '''
#        G_next = Gc.copy() # Create a copy of the current graph

#        for node_i in self.List_of_Nodes_in_State('S',Gc): #G_next.nodes():
#            G_next.nodes[node_i]['Risk']     = G_next.nodes[node_i]['Danger']*G_next.nodes[node_i]['Exposed']
#            G_next.nodes[node_i]['Risk_exp'] = np.exp(-J*pow( len( self.Neighbors_in_State('I',node_i,Gc)) /  Gc.degree(node_i),alpha) )

#        return G_next

#def Second_Infected_Neighborhood(self,vi,Gc):
'''
    		Given the label of the node vi, this function returns H_i = sum_{j in N^{S}_{i}} |N^{I}_{j}| / k_{j}   - Danger

    			Inputs:
    				vi    - (int)    the node lable
    				Gc    - (networkx graph) the current graph

    			Output:
    				List of nodes that are neighbors of node i and are in state 'S','E', 'I' or 'R'
'''
#return sum( [ len(self.Neighbors_in_State('I',vj,Gc)) / Gc.degree(vj)  for vj in self.Neighbors_in_State('S',vi,Gc)  ])

'''
Function to asses the danger of each node in the graph Gc.

Input:
            Gc - (networkx graph). The current graph to be updated in the atributed Danger.
        Output:
             G_next - (networkx graph). The updated graph with the current value of Danger.
'''

                                                                                                                                                                                                                                                 