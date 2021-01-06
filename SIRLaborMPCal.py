"""
@author: Constanza Fosco
--------------------------------------------------------------------------------------
SIRLabor_MP- Application to RM Chile

Date (initial date) March 25 2020.
Last version: November 2020.

Labor dynamics: confinement rules (partial or full lockdowns, activity closure, confinement by age)
+teleworkability and essentiality of jobs
Epidemics: SIR metapopulation model with minimal characteristic population size

Project: Fosco & Zurita, 
"Assessing the short-run effects of lockdown policies on economic activity, 
with an application to the Santiago metropolitan area"

All the processes and subprocesses are in this file.

Version for calibration - S0 based, runs a set of configurations {B,pA,Ks}
each one 5 realizations.

Before running, go to the end and check the parameters list.
By default, it is prepared to run as an example the configurations sim_s from 0 to 31
(others configurations setups of the calibration process, in the
file: Calibration_All.xls)

Input files must be located in the same folder.
'Data1_MP.csv', 'Data2_MP.csv','RealDRM.csv','RealDCom.csv'

-------------------------------------------------------------------------------------
Includes all: classes (Agent, SpatialUnit, SystemRM); special functions (readMyfileRM,
pDest, NSE); simulation function (LaborEpiRM).

Last part of this file: setup for simulation, for our implementation.
  

"""
import numpy as np
import random
import csv
import time


"""
General ad-hoc function
--------------------------------------------------------------------

"""
def readMyfileRM(data):
    """
    This function reads files with csv format, delimiter ",", first row
    with variable names. 
    The input is a csv file or txt, with csv format. To simplify, locate this
    file in the same directory.
    Call: x,y=readMyfileRM("myfile.txt") or readMyfileRM("myfile.csv") --depending on
    the case.
    The output yield two lists. The first one (x), a list of names. The second (y), a
    list of lists, each for each observation.
        
    """
    x=[]
    names=[]
    with open(data) as file:
        file1 = csv.reader(file, delimiter=',')
        ncol=len(next(file1))
        file.seek(0)
        count=0
        for line in file1:
            if count == 0:
                for j in range(ncol):
                    names.append(str(line[j]))
                    count+=1
            else:
                y=[]
                for j in range(ncol):
                    if names[j]== 'name' or names[j]== 'CUTh' or names[j]=='CUTw' or names[j]=='idrph':
                        y.append(str(line[j]))
                    else:
                        y.append(float(line[j]))
                    count+=1
                x.append(y)
            
    file.close()
    del file
    
    return names,x


def pDest( x, y ):
    """    
    This function regress y=pDe*x
    x : simulated data, total cummulative infected (array)
    y : real data, detected infected cummulative cases (array)
    and returns
    pDe : slope, i.e. detection fraction
    (used to estimate the average detection)
    
    """
    pDe = 0
    pDe1 = np.sum(x*y)
    pDe2 = np.sum(x**2)
    pDe = pDe1/pDe2
    return pDe


def NSE( simdata, realdata ):
    """
    This function computes the NSE (Nash-Sutcliffe Efficiency Coefficient) 
    = 1- (sum(Oi-Si)^2/sum(Oi-mean(O))^2). 
    simdata: matrix (simulated Si)
    realdata: matrix (observed, real data, Oi)
    Returns the coefficient (-inf,1)
    """
    realdataMean = realdata.mean()
    NSEc = 1 - ((((realdata-simdata)**2).sum())/ ((realdata-realdataMean)**2).sum())
    return NSEc 


"""
Classes: Agent, SpatialUnit, SystemRM
-------------------------------------------------------------------------------
Agent: residents and replicas (of commuters)
SpatialUnit: comuna (municipality)
SystemRM: the system of municipalities (the metropolitean region)

"""       
class Agent:
    def __init__( self, su=None, home=1, status=0, comm=0, activ=0, on=1 ):
        self.su = su            #Order CUT Spatial Unit (municipality) where the agent/replica is located (either on or off)
        self.home = home        #=1 residence, =0 workplace (=0 signals a replica)
        self.status = status    # =0 (S) susceptible, =1 (I) infectious, # =2 (R/D) Removed
        self.comm = comm        #=0 (default) non-commuter, =1 intracity commuter, =2 intercitycommuter, =3 commutes outside RM
        self.activ = activ      #=0 children [0,14], =1 employed , =2 unemployed, =3 inactive
        self.on = on            #=1 (default) staying at current SU, 0 when the agent is not at the corresponding location
        self.conf = 0           #=0 (default) not confined, =1 confined (by territory), =2 (by rama), =21 (by rama, then by comuna), =3 by age
        self.day = 0            #counter for the periods between two status.
        self.isol = 0           # =0 default, =1 isolated when I, =2 not moving when confined (THIS IS NOT PROPER ISOLATION) 
        self.pConf = None       #None (default), probability of moving/working in case that self is confined
        self.pConf7 = None      #same as pConf but applies only for malls closing (we need it when a comuna is not under lockdown but the closure of malls is on)
        self.replica = None     #replica of this Agent when comm in [1,2] (when it is a RM intracity or intercity commuter)
        self.ident = None       #fictional number, one to one relation with idrhp, identifies the group in data2
        self.age = None         #group of age: =1 [0,4];=2 [5,9];=3 [10,14];=4 [15,19];=5 [20,24];=6 [25,29];=7 [30,34];=8 [35,39];
                                #=9 [40,44];=10 [45,49];=11 [50,54];=12 [55,59];=13 [60,64];=14 [65,69];=15 [70,74];=16 [75,79];=17 [80,+]
        self.educ = None        #educ level: =99 unknown; =1 very_low; =2 low; =3 medium; =4 high; =5 very_high
        self.jornada = None     #jornada: =1 full; =2 partial; =99 non-workers
        self.jobcat = None      #jobcat:  from 0 to 7 (categoria_ocupacion), 0 for non-workers
        self.rama = None        #rama: from 1 to 21 (r_p_rev4cl_caenes, rama actividad) =99 for non-workers
        self.telew = None       #=0 (default) , =1 can telework if confined
        self.sector = None      #sector (employer): =1 formal, =2 informal, =3 other household
        self.work = 0           #=0 (default) doesn't work (today or never, depends on activ value), =1 working, but not teleworking, =2 teleworking
        
    
                     
    def creates_Replica ( self, SUwork ):
        #creates a replica in SU=SUwork (id order CUT municipality where the agent works)
        self.replica = Agent( su=SUwork, home=0, status=0, comm=self.comm, activ=self.activ, on=0 )
        self.replica.replica = self
          
    def update_Status ( self, newstatus, qR ):
        #updates the epidemics status. New status may be 1 (I) or 2 (R).
        self.status = newstatus
        if self.replica != None:
            self.replica.status = newstatus
        if newstatus == 1: #I
            self.day = 1
            if self.replica != None:
                self.replica.day = 1
        else: #R
            self.day = 0
            if self.replica != None:
                self.replica.day = 0
            pr_i = random.random()
            if pr_i <= qR:
                self.on = 1
                self.isol = 0
                if self.replica != None:
                    self.replica.on = 0
                    self.replica.isol = 0
            else: #(D)
                self.on = 0
                self.isol = 1
                if self.replica != None:
                    self.replica.on = 0
                    self.replica.isol = 1
        
    def commute_toW ( self ):
        #moves residents to their workplaces (home==1, comm>=1, activ==1)
        #returns Wr=0 if the agent is not working, Wr=1 if the agent is face-to-face working, Wr=2 teleworking
        
        #Special probability of face-to-face working for rama 7 (depending on the type of confinement)
        
        pConf_t = 0
        if self.rama == 7 and self.comm <= 2 and self.sector != 3 and (self.jobcat not in [5,6]) and self.conf == 2 and self.replica.conf == 2:
            pConf_t = self.pConf7
        else:
            pConf_t = self.pConf
        Wr = 0
        if self.comm == 1 or self.comm == 2:
            if self.conf == 0 and self.replica.conf == 0:
                self.on = 0
                self.replica.on = 1
                Wr = 1
            else: #self.conf ==0, but self.replica.conf!=0; or self.conf==1,2,21,3
                if self.telew == 0 and pConf_t == 0:
                    Wr = 0
                elif self.telew == 0 and pConf_t != 0: #not teleworking, but still with some probability of activity
                    pr_ic = random.random()
                    if pr_ic <= pConf_t:
                        self.on = 0
                        self.replica.on = 1
                        Wr = 1
                    else:
                        Wr = 0
                else:
                    Wr = 2
        else:
            if self.conf == 0:
                self.on = 0
                Wr = 1
            else:
                if self.telew == 0 and pConf_t == 0:
                    Wr = 0
                elif self.telew == 0 and pConf_t != 0:
                    pr_ic = random.random()
                    if pr_ic <= pConf_t:
                        self.on = 0
                        Wr = 1
                    else:
                        Wr = 0
                else:
                    Wr = 2
        return Wr
    
    
    def does_Work ( self ):
        #checks if self is working (applied to: commuters who do not move on weekends, and non-commuters)
        #special rules for rama 7
        
        pConf_t = 0
        if self.rama == 7 and self.sector != 3 and (self.jobcat not in [5,6]) and (( self.comm == 0 and self.conf == 2 ) or (self.comm in [1,2] and self.conf == 2 and self.replica.conf == 2)):
            pConf_t = self.pConf7
        else:
            pConf_t = self.pConf
        Wr = 0
        if self.comm == 0 or self.comm == 3:
            if self.conf == 0:
                Wr = 1
            else:
                if self.telew == 0 and pConf_t == 0:
                    Wr = 0
                elif self.telew == 0 and pConf_t != 0:
                    pr_ic = random.random()
                    if pr_ic <= pConf_t:
                        Wr = 1
                    else:
                        Wr = 0
                else:
                    Wr = 2
        else: # self.comm == 1 or self.comm == 2:
            if self.conf == 0 and self.replica.conf == 0:
                Wr = 1
            else: #self.conf ==0, but self.replica.conf!=0; or self.conf==1,2,21,3
                if self.telew == 0 and pConf_t == 0:
                    Wr = 0
                elif self.telew == 0 and pConf_t != 0: #not teleworking, but still with some probability of activity
                    pr_ic = random.random()
                    if pr_ic <= pConf_t:
                        Wr = 1
                    else:
                        Wr = 0
                else:
                    Wr = 2
        return Wr
            
    
    def back_Home ( self ):
        #moves workers to their home - call if you know that the commuter has previously moved
        self.on = 1
        if self.replica != None:
            self.replica.on = 0
               
    def on_Confinement ( self, kindConf ):
        #sets the agent in some sort of confinement (either by rama, by comuna, by age) 
        #(kindConf=1 comuna; =2 rama; =21 rama&comuna; =3 age)
        
        if self.conf == 0:
            self.conf = kindConf
        elif self.conf == 1: #already confined by comuna
            if kindConf == 3:
                self.conf = kindConf
        elif self.conf == 2: #already confined by rama
            if kindConf == 1:
                self.conf = 21
            elif kindConf == 3:
                self.conf = kindConf
            else:
                pass
        else:
            pass
    
    def order_CUTh ( self ):
        #returns the order in ListSU of the comuna (municipality) of residence of any agent/replica
        ind = 0
        if self.home == 1:
            ind = self.su
        else:
            ind = self.replica.su
        return ind
    
    def order_CUTw ( self ):
        #returns the order in ListSU of the comuna (municipality) of workplace of any resident, comm<=2
        ind = 0
        if self.home == 1 and self.comm in [0,1]:
            ind = self.su
        elif self.home == 1 and self.comm == 2:
            ind = self.replica.su
        else:
            print("Error: not a resident or not comm<=2")
        return ind
                         

class SpatialUnit:
    def __init__( self, cut ):
        self.cut = cut          #id (order CUT) of comuna 
        self.agents = [ ]       #list of agents and replicas (all, either on or off) in this SU
        self.confinrules =[]    #ordered list of (21) probabilities of activity level by "rama" (ISIC sector) used when SU=SUwork
                                #these are assigned to workers who work in self (pConf)
        self.confin7 = None     #special rule for rama=7 under agent.conf=2 and replica.conf=2
        
       
        
    def add_Agent ( self, agent ):
        #adds an Agent instance to the list
        self.agents.append( agent )
       
    def start_Confinement( self, partial = 1 ):
        #starts the confinement by comuna (municipality) in the spatial unit (territory dependent, i.e. conf=1).
        #partial is the fraction of population that is locked down (default =1, o.w <1).
        ConfinedAgents = []
        if partial != 1:
            R = [ i for i in self.agents if i.home == 1]
            W = [ i for i in self.agents if i.home == 0]
            RNr = int( partial*len( R ) )
            WNr = int( partial*len( W ) )
            Residents = random.sample( R, RNr ) 
            Workplaces = random.sample( W, WNr ) 
            ConfinedAgents = Residents + Workplaces
        else:
            ConfinedAgents = self.agents
        if len( ConfinedAgents ) > 0:
            for i in ConfinedAgents:
                i.on_Confinement( 1 )
                    
       
    def end_Confinement( self ):
        #ends confinement by comuna (municipality)
        for j in self.agents:
            if j.conf == 1:
                j.conf = 0
            elif j.conf == 21:
                j.conf = 2
            else:
                pass
    
    def get_S_agents( self ):
        #returns the list of all S --susceptible-- agents or replicas
        S = [i for i in self.agents if i.status == 0 ]
        return S
    
       
    def get_N( self ):
        #returns the effective number of agents (agents/replicas who are currently in self, on=1)
        N = 0
        N = len( [i for i in self.agents if i.on == 1 ])
        return N
        
    
    def get_I( self ):
        #returns the number of effective infective agents
        I = 0
        I= len( [i for i in self.agents if i.on == 1 and i.status == 1 ] )
        return I
       

class SystemRM:
    def __init__( self, data1=None, data2=None ):
        self.data1 = data1
        self.data2 = data2
        self.ListSU = []
        #data1: CUTh=w[0], rama3=w[1], rama4=w[2], rama7=w[3], rama8=w[4], rama9=w[5], rama10=w[6], 
        #rama13=w[7], rama14=w[8], rama17=w[9], rama18=w[10], rama19=w[11], rama7Conf2=w[12]
        #data2: nragt=w[0], CUTh=w[1], comm=w[2], activ=w[3], ident=w[4], age=w[5], educ=w[6], jornada=w[7],
        #jobcat=w[8], rama=w[9], telew=w[10], sector=w[11], CUTw=w[12], idrph=w[13], income=w[14],forml=w[15]
    
    
    def InitialSystem( self ):
        # we create the subpopulations, with agents and replicas 
        X,Y=readMyfileRM( self.data1 )
        
        for w in Y:
            z = SpatialUnit( int(w[0]) )
            z.confinrules = [ 1,1,w[1],w[2],1,0,w[3],w[4],w[5],w[6],1,0,w[7],w[8],1,0,w[9],w[10],w[11],0,0 ]
            z.confin7 = w[12]
            self.ListSU.append(z)
    
        #at this point, all the spatial units are created 
        #now we use the second file (see its format)
    
        W,Z=readMyfileRM( self.data2 )
            
        for w in Z: #for each line
            nAgt = int(w[0])  #number of agents to be created with the same set of characteristics
            SUr = self.ListSU[ int(w[1]) ] #the spatial unit of residence
            SUw = None
        
            for i in range( nAgt ): #for each new set of identical agents that must be created (residents)
                newAg = Agent( su=int(w[1]), home=1, status=0, comm=int(w[2]), activ=int(w[3]), on=1 )
                newAg.ident = int(w[4])
                newAg.age = int(w[5])
                newAg.educ = int(w[6])
                
                SUr.add_Agent( newAg )
            
                if int(w[3]) == 1: #activ==1
                    newAg.jornada = int(w[7])
                    newAg.jobcat = int(w[8])
                    newAg.rama = int(w[9])
                    newAg.telew = int(w[10])
                    newAg.sector = int(w[11])
                    
                    if int(w[2]) in [1,2]: #creates the replica (without characteristics)
                        SUw = self.ListSU[ int(w[12])]
                        newAg.creates_Replica( SUwork=int(w[12]) )
                        SUw.add_Agent( newAg.replica )
                    if int(w[9]) == 7 and int(w[2]) <= 2 and int(w[11]) != 3 and (int(w[8]) not in [5,6]):
                        newAg.pConf7 = self.ListSU[int(w[12])].confin7
   
    def get_Agents( self ):
        #retrieves all the agents in all the SU in a single list
        Agents = []
        for x in self.ListSU:
            Agents = Agents + x.agents 
        return Agents
    
    
    def get_S_AgtsSU( self ):
        #returns a list of lists of S agents (each list corresponds to each SU)
        SL = []
        for x in self.ListSU:
            SL.append( x.get_S_agents() )
        return SL

    def reset_Realization( self ):
        #restores the initial conditions
        
        for x in range( len( self.ListSU ) ):
            for i in self.ListSU[x].agents:
                i.status = 0
                i.conf = 0
                i.work = 0
                i.isol = 0
                i.day = 0
                if i.home == 1:
                    i.on = 1 
                else:
                    i.on = 0
  

def LaborEpiRM( sim, a , b, tmax, B, Nm, Ks, Kns, pD, pA, q, tau0, tau1, SystRM, situation ):
    """ Need to create the system and to initialize it as input """   
    #sim: code number of simulation (described in file Codigo)
    #a and b: range for realizations (a<b). For instance: a=0, b=2, will run 2 realizations, starting form rea=0
    #a=2, b=4, two realizations (starting from rea=2)
    #tmax: number of MC steps (days) for each realization
    #B: transmission rate (combination of average contacts and infection rate between a pair S,I)
    #Nm: size of the characteristic population (minimum size) for periphery small comunas (municipalities)
    #Ks: strict confinement factor (0,1) that reduces the interaction time between agents when confined
    #Kns: non-strict confinement factor (0,1) that reduces the interaction time between agents when confined 
    #pD: probability of detection (if I)
    #pA: probability of isolation (if I)
    #q: probability of recovering (1-q = dead), (if I)
    #tau0: fraction of interaction time in the two first rounds of contagion
    #tau1: fraction of interaction time in the last round of contagion (eventually changes with curfew)
    #SystRM: initial system of comunas and agents
    #situation: =0 "real" case; =1 without any confinement; =2 with full confinement; >3 without any measure
    
      
    
    
    DeadRM = np.zeros( ( (b-a),1 )) #vector to keep track of 
    
    AllAgents = SystRM.get_Agents( )
    
    Detected_RM_Calib = np.zeros( ( 22, 51 ) )
    
    for rea in range( a , b ):  #for each realization
                   
        
        tau1 = 6.0/24.0 #initial time fraction of last round of contagion (it will change due to the curfew)
        
        random.shuffle( AllAgents )
        
        #Initial pConf (assignment of initial pConf to each worker, some of them will change in time)
        for i in AllAgents:
            if i.home == 1 and i.activ == 1:
                if i.comm == 0:
                    i.pConf = SystRM.ListSU[ i.su ].confinrules[ int(i.rama)-1 ]
                elif i.comm == 3:
                    if i.sector == 1:
                        i.pConf = 1
                    else:
                        i.pConf = 0
                else:
                    if i.jobcat == 5:
                        i.pConf = 0
                    elif i.jobcat == 6:
                        i.pConf = 1
                    else:
                        if i.sector == 3:
                            i.pConf = 0
                        else:
                            i.pConf = SystRM.ListSU[ i.replica.su ].confinrules[ int(i.rama)-1 ]
        
        
                    
        #Lists for tracking agents for updating
        
        I_rea = []
        R_rea = []
        S_rea = SystRM.get_S_AgtsSU() #list of lists, includes agents and replicas
        
        #Matrix Mobility keeps track of time spent in workplaces whenever these
        #workplaces are inside RM (i.e. excludes commuting outside the RM)
        #This matrix allows us to compare mobility with Google Analytics Reports
        Mobility = np.zeros( ( tmax,1 ) )
        
        #Matrices that keep track of infected in their detection day. Only for showing calibration results        
        Detected_RM = np.zeros( ( tmax, 51 ) )
        Detected_RM_Cum = np.zeros( ( 1, 51 ) ) 
        
        #To keep track the total number of dead at t=151 (i.e. July 30)
        Fall = 0
        
        
        
        #Initial cases (March 1) (week 9 informe epid, pD 0.05)
        #LAST UPDATE: INFORME EPIDEMIOLOGICO NOVEMBER 11, 2020
        #Initial infected at t=0, a list including all the spatial units
        I_Init =  [1,1,1,0,0,1,1,0,0,1,0,1,0,1,0,1,0,0,0,1,0,1,0,1,0,1,0,0,0,1,0,0,3,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1]
        
        for x in range( 51 ):
            if I_Init[x] > 0:
                S_Su = [i for i in S_rea[x] if i.home == 1 ]
                Inx = random.sample( S_Su, I_Init[x] )
                for j in Inx:
                    j.update_Status( 1, q )
                    Detected_RM_Cum[0][j.order_CUTh()] += 15
                    S_rea[x].remove( j )
                    if j.replica != None:
                        S_rea[ j.replica.su ].remove( j.replica )
                    I_rea.append( j ) 
                del( S_Su )
                
        #List of potential commuters. The distribution of commuters who work in either shift (randomly chosen) changes each realization
        
        T0 = [ i for i in AllAgents if i.home == 1 and i.activ == 1 and i.comm == 0 and i.jobcat != 6 ]
        T1 = [ i for i in AllAgents if i.home == 1 and i.activ == 1 and i.jobcat != 6 and ( i.comm == 3 or ( i.comm in [1,2] and i.jornada == 1 ) ) ]
        T2 = [ i for i in AllAgents if i.home == 1 and i.activ == 1 and i.jobcat != 6 and i.jornada == 2 and i.comm in [1,2] ]
        T2_morning = random.sample( T2, int( 0.5*len( T2 ) ) )
        T2_afternoon = [ i for i in T2 if i not in T2_morning ]
        T4 = [ i for i in AllAgents if i.home == 1 and i.jobcat == 6 and i.comm in [1,2] ]
        
        
        
        
        t = 0 #assumed to be March 1, Sunday
        
        d = 7
       
                
        while t < tmax:
            
            
            if situation == 0: #"real case"
            
                
                #Measures and special events (confinement measures, teleworking, etc.)
            
                # 1) Closing of rama = 16 (schools, universities, etc.)- March 15, effective from March 16. i.conf=2
                if t == 15:
                    Agents16 = [ i for i in AllAgents if i.home == 1 and i.rama == 16 ]
                    for i in Agents16:
                        i.on_Confinement( 2 )
                        if i.replica != None:
                            i.replica.on_Confinement( 2 )
                                          
                # 2) Teleworking in Administracion Publica (more or less implemented March 19). i.conf=2.
                # 3)Malls closing - affects rama 7
                if t == 18:
                    Agents15 = [ i for i in AllAgents if i.home == 1 and i.rama == 15 ]
                    for i in Agents15:
                        i.on_Confinement( 2 )
                        if i.replica != None:
                            i.replica.on_Confinement( 2 )
                    Agents7 = [ i for i in AllAgents if i.home == 1 and i.rama == 7 ]
                    for i in Agents7:
                        i.on_Confinement( 2 )
                        if i.replica != None:
                            i.replica.on_Confinement( 2 )
                                 
                # 4) Closing of cinemas, theaters, restaurants, pubs, and sport facilities, March 21.
                # This measure affects part of rama=9 and rama=18. i.conf = 2
                if t == 20:
                    Agents9_18 = [ i for i in AllAgents if i.home == 1 and i.rama in [9,18] ]
                    for i in Agents9_18:
                        i.on_Confinement( 2 )
                        if i.replica != None:
                            i.replica.on_Confinement( 2 )
                            
                # 5) Curfew (toque de queda). March 22. Affects the last round of contagion.
                if t >= 21:
                    tau1 = 5.0/24.0
                
                # 6) Mandatory confinementent of old people (>= 80 years) - March 24 - 
                # Self-confinement: (assumption) People with [50,80) years, educ in 4,5, and not working or teleworking
                
                if t == 23:
                    Agents_age80 = [ i for i in AllAgents if i.home == 1 and i.age == 17 ]
                    for i in Agents_age80:
                        i.on_Confinement( 3 )
                        if i.replica != None:
                            i.replica.on_Confinement( 3 )
                        if i.isol == 0:
                            i.isol = 2
                            if i.replica != None:
                                i.replica.isol = 2
                   
                    Agents_age50_high = [ i for i in AllAgents if i.home == 1 and i.age in [11,12,13,14,15,16] and i.educ in [4,5] and ( i.activ != 1 or ( i.activ == 1 and i.telew == 1 and i.rama != 17 ) ) ]
                    for i in Agents_age50_high:
                        i.on_Confinement( 3 )
                        if i.replica != None:
                            i.replica.on_Confinement( 3 )    
                            
                            
                # 7) Confinement (March 26 night --> March 27 morning): Lo Barnechea (partial 97.3%), Vitacura, Las Condes, Providencia, Santiago, Ñuñoa, Independencia.
                #Confinement rules ramas 9 and 18 according to Instructivo 1.
                            
                if t == 26:
                    SystRM.ListSU[14].start_Confinement( partial = 0.973 )
                    ComConf = [ 0,7,13,19,22,31 ]
                    for x in ComConf:
                        SystRM.ListSU[x].start_Confinement( partial = 1 )
                        
                    Agts_rama918 = [ i for i in AllAgents if i.home == 1 and i.activ == 1 and (i.rama in [9,18]) and i.comm <=2 and (i.jobcat not in [5,6]) and i.sector !=3 ]
                    rules_9 = [0.391976187,0.511875512,0.55,0.633676093,0.629807692,0.558896313,0.388353414,0.81517094,0.433873497,0.473962571,0.427672956,0.400974026,0.457234363,0.484398724,0.452434998,0.555555556,0.486778846,0.521727973,0.566666667,0.473073202,0.439093484,0.462666145,0.536728566,0.694219538,0.88252149,0.520616642,0.449798721,0.803463203,0.403250774,0.531933899,0.480620155,0.416213655,0.541821561,0.397420147,0.343283582,0.326189726,0.674634794,0.841296928,0.477348777,0.465317919,0.488372093,0.473282443,0.544971893,0.38356974,0.239130435,0.459016393,0.507857143,0.435233161,0.264214047,0.428819444,0.458072591]
                    rules_18 = [0.001344011,0.017857143,0,0,0,0.001587302,0,0.000346741,0,0.01715439,0.033333333,0,0.001929012,0.00152391,0.001868207,0,0,0.001583531,0.000597372,0.000125282,0.013513514,0.004065041,0.001447078,0,0,0,0,0,0,0.002257336,0,0.002480022,0.001078749,0,0,0,0,0,0,0.000979432,0.025,0,0,0.019230769,0,0,0,0,0,0,0.003030303]
                    for j in Agts_rama918:
                        if j.rama == 9:
                            j.pConf = rules_9[ j.order_CUTw() ]
                        else:
                            j.pConf = rules_18[ j.order_CUTw() ]
                
                # 8) New confinement rules (Instructivo 2) - April 2
                
                if t == 32:
                    Agts_rama68131418 = [ i for i in AllAgents if i.home == 1 and i.activ == 1 and (i.rama in [6,8,13,14,18]) and i.comm <=2 and (i.jobcat not in [5,6]) and i.sector !=3 ]
                    rules_6 = [0.035758243,0.000249906,0.009450473,0.00287234,0.000401445,0.004246815,0.035407433,0,0.000134898,0.012295299,0,0.000319285,0.000870133,0.053334693,0.027056633,0.000468604,0.000472813,0.000121743,0.006339982,0.006661749,0.036572248,0.000723327,0.004784538,0.004072609,0.002893947,0.010225231,0.001117545,0.053275662,0.091791553,0.099163059,0,0.002688807,0.019400786,0,0,0.01843318,0.002733598,0,0.018754423,0.001283285,0,0.165266106,0.005059631,0,0.006711409,0,0,0.00128041,0,0,0.000274499]
                    rules_8 = [0.784726596,0.954000436,0.56043956,0.916374562,0.888663968,0.490662438,0.632779161,0.738021638,0.758106022,0.800846177,0.910922587,0.903013699,0.701140065,0.724905382,0.848421053,0.942982456,0.776386404,0.950196592,0.874547312,0.802397149,0.762886598,0.874945151,0.81243997,0.86889332,0.974326402,0.519360902,0.868006993,0.887513751,0.968149646,0.760928962,0.9125,0.109670638,0.84676354,0.952714536,0.723076923,0.891231286,0.882196466,0.89456869,0.812801285,0.873200443,0.942105263,0.895746888,0.846153846,0.904051173,0.846153846,0.823529412,0.827642276,0.872222222,0.878331402,0.926169591,0.728547154]
                    rules_13 =[0.02762702,0.067650677,0.18487395,0.12329932,0.076205288,0.067354699,0.003307607,0.023462783,0.034839204,0.057528343,0.086956522,0.032258065,0.046686511,0.005067366,0.044543984,0.245,0.162094763,0.028865164,0.144761397,0.020182374,0.021526419,0.054545455,0.005312832,0.025935532,0.032979639,0.030917553,0.021047479,0.004213327,0.03816047,0.032494197,0.025531915,0.011162066,0.14720986,0.028225806,0.014792899,0.075689784,0.070619587,0.116883117,0.082922014,0.141821112,0.033980583,0.012931034,0.086474501,0.022082019,0.111111111,0.375,0.203669725,0.007692308,0.064516129,0.071428571,0.108956602]
                    rules_14 = [0.409129886,0.193347193,0.440316206,0.564774656,0.336471551,0.420006517,0.359611559,0.670524412,0.210947931,0.434208638,0.054764513,0.293494705,0.27464367,0.129914829,0.14563591,0.177897574,0.322580645,0.445812266,0.225176568,0.215839575,0.426487093,0.725606963,0.19005309,0.242094017,0.161824295,0.121266428,0.113072766,0.103007878,0.358832225,0.473580643,0.411483254,0.134729294,0.380033685,0.323076923,0.010273973,0.610806306,0.081185567,0.448113208,0.349690804,0.550053438,0.32238193,0.365327381,0.218444968,0.338461538,0.034482759,0.06,0.4017991,0.125,0.030744337,0.119897959,0.540950455]
                    rules_18b = [0.125217002,0.017857143,0,0,0,0.007936508,0,0.723300971,0,0.025227043,0.033333333,0,0.001929012,0.004098791,0.001868207,0.022222222,0,0.001583531,0.347072879,0.000375846,0.189189189,0.01300813,0.004754686,0.040816327,0,0,0,0.007246377,0.002293578,0.002257336,0,0.002480022,0.005393743,0,0,0.006635071,0.094488189,0.242424242,0.11751663,0.001958864,0.025,0,0,0.019230769,0,0,0,0,0,0.024390244,0.006060606]
                    for j in Agts_rama68131418:
                        if j.rama == 6:
                            j.pConf = rules_6[ j.order_CUTw() ]
                        elif j.rama == 8:
                            j.pConf = rules_8[ j.order_CUTw() ]
                        elif j.rama == 13:
                            j.pConf = rules_13[ j.order_CUTw() ]
                        elif j.rama == 14:
                            j.pConf = rules_14[ j.order_CUTw() ]
                        else:
                            j.pConf = rules_18b[ j.order_CUTw() ]
                            
                # 9) Ending confinement: Independencia (April 2, night -> April 3 morning)
                
                if t == 33:
                    SystRM.ListSU[7].end_Confinement()
                    
                
                # 10) Partial confinement: Puente Alto (49.9%) (April 9 night --> April 10)
                if t == 40:
                    SystRM.ListSU[32].start_Confinement( partial = 0.499 )
                
                # 11) Ending confinement, total: Lo Barnechea, Vitacura, Providencia. Partial: Ñuñoa, Santiago
                # Partial confinement: Santiago (76.8%), Ñuñoa (76.8%). (April 13)
                if t == 43:
                    ComEnd = [ 0,19,14,22,31 ]
                    for x in ComEnd:
                        SystRM.ListSU[x].end_Confinement( )
                    SystRM.ListSU[0].start_Confinement( partial = 0.768 )
                    SystRM.ListSU[19].start_Confinement( partial = 0.768 )
                            
                # 12) Ending confinement: Las Condes. Confinement: El Bosque. Partial confinement: 
                # San Bernardo (42.1%) (April 17)
                # New confinement rules (Instructivo 3)
                if t == 47:
                    SystRM.ListSU[13].end_Confinement()
                    SystRM.ListSU[4].start_Confinement( partial = 1 )
                    SystRM.ListSU[38].start_Confinement( partial = 0.421 )
                    
                    Agts_rama8 = [ i for i in AllAgents if i.home == 1 and i.activ == 1 and i.rama == 8 and i.comm <=2 and (i.jobcat not in [5,6]) and i.sector !=3 ]
                    rules_8b = [0.839078886,0.962502725,0.611158073,0.932398598,0.904453441,0.521212121,0.649937767,0.773570325,0.782810087,0.828649139,0.928950159,0.910684932,0.708469055,0.798267121,0.878596491,0.96125731,0.828264758,0.957011796,0.906397994,0.850016197,0.769072165,0.883282141,0.877256318,0.955051512,0.98562546,0.528759398,0.921328671,0.892051705,0.984327604,0.862021858,0.920833333,0.129000234,0.858872743,0.954465849,0.738461538,0.916895814,0.886946608,0.897763578,0.95661489,0.896456257,0.943157895,0.897302905,0.895604396,0.946695096,0.846153846,0.904411765,0.828455285,0.877777778,0.908458864,0.975146199,0.83942226]
                    for j in Agts_rama8:
                        j.pConf = rules_8b[ j.order_CUTw() ]
            
                # 13) Confinement: Pedro Aguirre Cerda, Quinta Normal. Partial confinement: Independencia (29.7%)
                # (April 23 night --> April 24 morning)
                if t == 54:
                    SystRM.ListSU[20].start_Confinement( partial = 1 )
                    SystRM.ListSU[25].start_Confinement( partial = 1 )
                    SystRM.ListSU[7].start_Confinement( partial = 0.297 )
                        
                #14) Partial Confinement: La Pintana (46.7%), San Ramon (56.5%). 
                # Confinement: Estacion Central, Independencia (rest). (April 30 night --> May 1 morning).
                
                if t == 61:
                    SystRM.ListSU[5].start_Confinement( partial = 1 )
                    SystRM.ListSU[7].start_Confinement( partial = 1 )
                    SystRM.ListSU[11].start_Confinement( partial = 0.467 )
                    SystRM.ListSU[30].start_Confinement( partial = 0.565 )
                    
                    
                # 15) Confinement: Cerrillos, Recoleta, Santiago (rest).
                # Partial confinement: Quilicura (99.4%). (May 5 night -> May 6 morning) 
                if t == 66:
                    SystRM.ListSU[0].start_Confinement( partial = 1 )
                    SystRM.ListSU[1].start_Confinement( partial = 1 )
                    SystRM.ListSU[24].start_Confinement( partial = 0.994 )
                    SystRM.ListSU[26].start_Confinement( partial = 1 )
            
                # 16) Ending confinement Ñuñoa (May 7 night --> May 8 morning)
                if t == 68:
                    SystRM.ListSU[19].end_Confinement( )                
                
                #17) Confinement: San Ramon (rest), La Pintana (rest), San Miguel, San Joaquín, Renca,
                # Cerro Navia, Conchali, La Cisterna, La Florida, La Granja, Lo Espejo, Lo Prado, Macul, Peñalolen
                # Partial confinement: Puente Alto (increase until 97.7%), San Bernardo (increase until 62.7%)
                # (May 8 night --> May 9 morning) 
                if t == 69:
                    ComConf = [ 3,27,2,16,29,28,17,21,15,8,30,10,9,11 ]
                    for x in ComConf:
                        SystRM.ListSU[x].start_Confinement( partial = 1 )
                    SystRM.ListSU[32].start_Confinement( partial = 0.977 )
                    SystRM.ListSU[38].start_Confinement( partial = 0.627 )
            
                #18) Confinement of people with 75 or more years (i.e until 79, because >=80 are already confined).
                # Confinement: Lampa, Colina, Quilicura (rest), Huechuraba, Pudahuel, Providencia, Vitacura, Lo Barnechea, 
                # Las Condes, Maipú, Padre Hurtado, Ñuñoa, Puente Alto (rest), San Bernardo (rest), Buin, La Reina.
                # (May 15 night --> May 16 morning).
                if t == 76:
                    Agents_75 = [ i for i in AllAgents if i.home == 1 and i.age == 16 ]
                    for i in Agents_75:
                        i.on_Confinement( 3 )
                        if i.replica != None:
                            i.replica.on_Confinement( 3 )
                        if i.isol == 0:
                            i.isol = 2
                            if i.replica != None:
                                i.replica.isol = 2
                        
                                
                    ComConf = [ 36,35,24,6,23,22,31,14,13,18,49,19,32,38,39,12 ]
                    for x in ComConf:
                        SystRM.ListSU[x].start_Confinement( partial = 1 )
                
                
                #19) New confinement rules (May 27). Ramas 10 and 19 (Instructivo 6)
                if t == 87:
                    rules_10 = [0.363608563,0.979029605,0.913978495,0.835185185,0.681818182,0.260135135,0.759002338,0.691943128,0.950504125,0.210626186,0.865853659,0.834782609,0.115076014,0.572575546,0.665594855,0.988764045,0.666666667,0.916751269,0.50309119,0.798590131,0.6,0.612648221,0.520888993,0.65325285,0.590772317,0.63880289,0.893125671,0.687272727,0.749094671,0.583993661,0.744186047,0.743822076,0.574529667,0.784313725,0.878787879,0.445283019,0.531147541,0.204545455,0.474332649,0.696864111,0.966216216,0.770114943,0.426829268,0.394736842,0.92,0.833333333,0.736434109,0.842105263,0.386363636,0.740112994,0.758490566]
                    rules_19 = [0.00183531,0.000656599,0.016722408,0.001011122,0.00267666,0.003370614,0.002629602,0.059984896,0.029272899,0.02179676,0.027104137,0.004864489,0.00147232,0.008024586,0.000930665,0.012669683,0.005263158,0.000740741,0.004814765,0.001929571,0.00097229,0.002251472,0.092587216,0.004060456,0.070404172,0.00497822,0.002403021,0.007168459,0.010692178,0.015388097,0.030534351,0.003627428,0.013106525,0,0,0.00295858,0.005263158,0.008810573,0.018510158,0.006242906,0,0.009950249,0.009615385,0.017326733,0,0,0.006648936,0.012269939,0.001342282,0.00310559,0.02510917]
                    Agts_rama1019 = [ i for i in AllAgents if i.home == 1 and i.activ == 1 and i.rama in [10,19] and i.comm <=2 and (i.jobcat not in [5,6]) and i.sector !=3 ]
                    for j in Agts_rama1019:
                        if j.rama == 10:
                            j.pConf = rules_10[ j.order_CUTw() ]
                        else:
                            j.pConf = rules_19[ j.order_CUTw() ]
                
                #20) Confinement of Melipilla (58.4%), Curacavi (57.4%) , Tiltil (27.7%), San José de Maipo (34.3%) 
                # and Peñaflor. (June 12 night --> June 13 morning). 
                
                if t == 104:
                    SystRM.ListSU[42].start_Confinement( partial = 0.584 )
                    SystRM.ListSU[43].start_Confinement( partial = 0.574 )
                    SystRM.ListSU[37].start_Confinement( partial = 0.277 )
                    SystRM.ListSU[34].start_Confinement( partial = 0.343 )
                    SystRM.ListSU[50].start_Confinement( partial = 1 )
                                 
                #21) Confinement Calera de Tango, El Monte, Talagante, conf=1. ( June 27 )
                if t == 118:
                    ComConf = [ 40, 46, 47 ]
                    for x in ComConf:
                        SystRM.ListSU[x].start_Confinement( partial = 1 )
                
                #22) New confinement rules (July 9) (Instructivo 9)
                if t == 130:
                    rules_6b = [0.570463821,0.17693365,0.879243962,0.226702128,0.732838218,0.645515863,0.410841427,0.415950093,0.646297046,0.573533309,0.781087634,0.727650064,0.296606482,0.127764298,0.166940519,0.939081537,0.845390071,0.43450207,0.720908731,0.739454094,0.673359627,0.671850512,0.244017749,0.63870142,0.145650256,0.532126572,0.433912425,0.507437493,0.779291553,0.706378066,0.85824123,0.153852972,0.815692534,0.628571429,0.542168675,0.47281106,0.476267396,0.83125,0.632814343,0.561116458,0.614876033,0.560690943,0.745934225,0.776274714,0.355704698,0.602678571,0.738831615,0.774647887,0.722998729,0.668604651,0.555037057]
                    rules_13b = [0.02762702,0.067650677,0.193277311,0.12329932,0.076205288,0.067354699,0.003365636,0.023462783,0.034839204,0.057528343,0.086956522,0.189964158,0.046686511,0.006084101,0.044543984,0.245,0.16957606,0.028865164,0.144761397,0.020182374,0.021526419,0.054545455,0.011007632,0.025935532,0.032979639,0.030917553,0.021047479,0.056436412,0.175146771,0.032494197,0.029787234,0.011162066,0.155083875,0.028225806,0.014792899,0.07606264,0.070619587,0.116883117,0.09970385,0.141821112,0.033980583,0.012931034,0.086474501,0.022082019,0.111111111,0.375,0.203669725,0.007692308,0.064516129,0.071428571,0.108956602]
                    Agts_rama613 = [ i for i in AllAgents if i.home == 1 and i.activ == 1 and i.rama in [6,13] and i.comm <=2 and (i.jobcat not in [5,6]) and i.sector !=3 ]
                    for j in Agts_rama613:
                        if j.rama == 6:
                            j.pConf = rules_6b[ j.order_CUTw() ]
                        else:
                            j.pConf = rules_13b[ j.order_CUTw() ]
                
                #23) Confinement ending in: Colina, La Reina, Las Condes, Lo Barnechea, Ñuñoa, Vitacura, Tiltil
                # Confinement: Isla de Maipo (July 28)
                if t == 149:
                    SystRM.ListSU[48].start_Confinement( partial = 1 )
                    ComEnd = [ 35,12,13,14,19,31,37 ]
                    for x in ComEnd:
                        SystRM.ListSU[x].end_Confinement( )
                    
                        
                            
                    
            elif situation == 1: #without any comuna/rama confinement. Just voluntary confinement of people >= 65, not working. From March 24.
                if t == 23:
                    Agents_65 =[ i for i in AllAgents if i.home == 1 and i.activ != 1 and i.age in [14,15,16,17] ]
                    for i in Agents_65:
                        i.on_Confinement( 3 )                       
            
                                   
            elif situation == 2:  #follows sit=0 until t=25. On t=26, instead of locking down a selective number of comunas, they are all locked down
                # 1) Closing of rama = 16 (schools, universities, etc.)- March 15, effective from March 16. i.conf=2
                if t == 15:
                    Agents16 = [ i for i in AllAgents if i.home == 1 and i.rama == 16 ]
                    for i in Agents16:
                        i.on_Confinement( 2 )
                        if i.replica != None:
                            i.replica.on_Confinement( 2 )
                                          
                # 2) Teleworking in Administracion Publica (more or less implemented March 19). i.conf=2.
                # 3)Malls closing - affects rama 7
                if t == 18:
                    Agents15 = [ i for i in AllAgents if i.home == 1 and i.rama == 15 ]
                    for i in Agents15:
                        i.on_Confinement( 2 )
                        if i.replica != None:
                            i.replica.on_Confinement( 2 )
                    Agents7 = [ i for i in AllAgents if i.home == 1 and i.rama == 7 ]
                    for i in Agents7:
                        i.on_Confinement( 2 )
                        if i.replica != None:
                            i.replica.on_Confinement( 2 )
                                 
                # 4) Closing of cinemas, theaters, restaurants, pubs, and sport facilities, March 21.
                # This measure affects part of rama=9 and rama=18. i.conf = 2
                if t == 20:
                    Agents9_18 = [ i for i in AllAgents if i.home == 1 and i.rama in [9,18] ]
                    for i in Agents9_18:
                        i.on_Confinement( 2 )
                        if i.replica != None:
                            i.replica.on_Confinement( 2 )
                            
                # 5) Curfew (toque de queda). March 22. Affects the last round of contagion.
                if t >= 21:
                    tau1 = 5.0/24.0
                
                # 6) Mandatory confinementent of old people (>= 80 years) - March 24 - 
                # Self-confinement: (assumption) People with [50,80) years, educ in 4,5, and not working or teleworking
                
                if t == 23:
                    Agents_age80 = [ i for i in AllAgents if i.home == 1 and i.age == 17 ]
                    for i in Agents_age80:
                        i.on_Confinement( 3 )
                        if i.replica != None:
                            i.replica.on_Confinement( 3 )
                        if i.isol == 0:
                            i.isol = 2
                            if i.replica != None:
                                i.replica.isol = 2
                   
                    Agents_age50_high = [ i for i in AllAgents if i.home == 1 and i.age in [11,12,13,14,15,16] and i.educ in [4,5] and ( i.activ != 1 or ( i.activ == 1 and i.telew == 1 and i.rama != 17 ) ) ]
                    for i in Agents_age50_high:
                        i.on_Confinement( 3 )
                        if i.replica != None:
                            i.replica.on_Confinement( 3 )    
                            
                            
                # 7) Confinement (March 26 night --> March 27 morning): All the comunas
                #Confinement rules ramas 9 and 18 according to Instructivo 1.
                            
                if t == 26:
                    for x in range( 51 ):
                        SystRM.ListSU[x].start_Confinement( partial = 1 )
                        
                    Agts_rama918 = [ i for i in AllAgents if i.home == 1 and i.activ == 1 and (i.rama in [9,18]) and i.comm <=2 and (i.jobcat not in [5,6]) and i.sector !=3 ]
                    rules_9 = [0.391976187,0.511875512,0.55,0.633676093,0.629807692,0.558896313,0.388353414,0.81517094,0.433873497,0.473962571,0.427672956,0.400974026,0.457234363,0.484398724,0.452434998,0.555555556,0.486778846,0.521727973,0.566666667,0.473073202,0.439093484,0.462666145,0.536728566,0.694219538,0.88252149,0.520616642,0.449798721,0.803463203,0.403250774,0.531933899,0.480620155,0.416213655,0.541821561,0.397420147,0.343283582,0.326189726,0.674634794,0.841296928,0.477348777,0.465317919,0.488372093,0.473282443,0.544971893,0.38356974,0.239130435,0.459016393,0.507857143,0.435233161,0.264214047,0.428819444,0.458072591]
                    rules_18 = [0.001344011,0.017857143,0,0,0,0.001587302,0,0.000346741,0,0.01715439,0.033333333,0,0.001929012,0.00152391,0.001868207,0,0,0.001583531,0.000597372,0.000125282,0.013513514,0.004065041,0.001447078,0,0,0,0,0,0,0.002257336,0,0.002480022,0.001078749,0,0,0,0,0,0,0.000979432,0.025,0,0,0.019230769,0,0,0,0,0,0,0.003030303]
                    for j in Agts_rama918:
                        if j.rama == 9:
                            j.pConf = rules_9[ j.order_CUTw() ]
                        else:
                            j.pConf = rules_18[ j.order_CUTw() ]

            else:
                pass
                
            
            #A day begins
            #Detecting, isolating I
            
            if len( I_rea ) > 0: #if there are I (notice that in this list there are only agents, there are not replicas)
                for i in I_rea:
                    if i.day == 6: # Possible detection/isolation
                        i.day = i.day + 1
                        pr_i = random.random()
                        if pr_i <= pA*pD: # I isolated and detected.
                            i.on = 0 
                            i.isol = 1
                            if i.replica != None:
                                i.replica.on = 0
                                i.replica.isol = 1
                            Detected_RM_Cum[0][i.order_CUTh()] += 15
                            
                        elif  (pA*pD) < pr_i <= ((pA*pD)+(pA*(1-pD))): #Isolated, not detected
                            i.on = 0
                            i.isol = 1
                            if i.replica != None:
                                i.replica.on = 0
                                i.replica.isol = 1
                        elif  ((pA*pD)+(pA*(1-pD))) < pr_i <= ((pA*pD)+(pA*(1-pD))+((1-pA)*pD)): #Not isolated, detected
                            Detected_RM_Cum[0][i.order_CUTh()] += 15
                            
                        else: #Not isolated, not detected
                            pass   
                    elif i.day >= 13: #Recovered
                        i.update_Status ( 2, q ) #I--> R/D
                        R_rea.append( i )
                        I_rea.remove( i )
                        if i.isol == 1:
                            Fall += 15
                    else:
                        i.day = i.day + 1
                                    
            #Move people jobcat==6 (servicio doméstico puertas adentro) (we assume that in case of confinement the employee stays with the employer)
            
            for i in T4: 
                if i.isol == 1:
                    i.work = 0
                    i.on = 0
                    i.replica.on = 0
                elif i.isol == 2:
                    i.work = 1
                    i.on = 0
                    i.replica.on = 1
                else:
                    if d >= 1 and d <= 6 and (t not in [40,61,81,120,137]):
                        i.on = 0
                        i.replica.on = 1
                        i.work = 1
                        if i.comm == 1 or i.comm == 2:
                            Mobility[ t ][0] += 1
                    else: # d == 7
                        i.work = 1
                        if i.conf == 0:
                            i.on = 1
                            i.replica.on = 0
                        else:
                            i.on = 0
                            i.replica.on = 1
                            
                    
            #Checking who is working at home (from those workers that always work at home)
            
            for i in T0:
                if i.isol == 1:
                    i.work = 0
                else:
                    w_i = i.does_Work()
                    i.work = int( w_i )
                
            #Moving commuters
            
            CommRMT1=[]  
            CommRMT2=[]
            
            
            if len( T1 ) > 0:
                for i in T1:
                    if i.isol == 0:
                        if d <= 5 and (t not in [40,61,81,120,137]):
                            c_i = i.commute_toW()
                            if i.on == 0:
                                CommRMT1.append( i )
                                if i.comm == 1 or i.comm == 2:
                                    Mobility[ t ][0] += 1
                                i.work = 1
                            else: 
                                i.work = int(c_i)
                        elif d == 6 or (t in [40,120,137]):
                            if i.rama not in [11,15,16,21]:
                                c_i = i.commute_toW()
                                if i.on == 0:
                                    CommRMT1.append( i )
                                    if i.comm == 1 or i.comm == 2:
                                        Mobility[ t ][0] += 1
                                    i.work = 1
                                else: 
                                    i.work = int(c_i)
                            else: #We check if the agent is working 
                                w_i = i.does_Work()
                                i.work = int(w_i)
                        elif d == 7 or (t in [61,81]):
                            if i.rama not in [3,6,10,11,12,13,14,15,16,19,20,21]:
                                c_i = i.commute_toW()
                                if i.on == 0:
                                    CommRMT1.append( i )
                                    if i.comm == 1 or i.comm == 2:
                                        Mobility[ t ][0] += 1
                                    i.work = 1
                                else: 
                                    i.work = int(c_i)
                            else: #We check if the agent is working
                                w_i = i.does_Work()
                                i.work = int(w_i)
                        else:
                            pass
                    elif i.isol == 2:
                        w_i = i.does_Work()
                        i.work = int(w_i)
                    else:    
                        i.work = 0
                                                               
                            
            if len( T2_morning ) > 0:
                for i in T2_morning:
                    if i.isol == 0:
                        if d <= 5 and (t not in [40,61,81,120,137]):
                            c_i = i.commute_toW()
                            if i.on == 0:
                                CommRMT2.append( i )
                                if i.comm == 1 or i.comm == 2:
                                    Mobility[ t ][0] += 0.5
                                i.work = 1
                            else: 
                                i.work = int(c_i)
                        elif d == 6 or (t in [40,120,137]):
                            if i.rama not in [11,15,16,21]:
                                c_i = i.commute_toW()
                                if i.on == 0:
                                    CommRMT2.append( i )
                                    if i.comm == 1 or i.comm == 2:
                                        Mobility[ t ][0] += 0.5
                                    i.work = 1
                                else: 
                                    i.work = int(c_i)
                            else: #We check if the agent is working 
                                w_i = i.does_Work()
                                i.work = int(w_i)
                        elif d == 7 or (t in [61,81]):
                            if i.rama not in [3,6,10,11,12,13,14,15,16,19,20,21]:
                                c_i = i.commute_toW()
                                if i.on == 0:
                                    CommRMT2.append( i )
                                    if i.comm == 1 or i.comm == 2:
                                        Mobility[ t ][0] += 0.5
                                    i.work = 1
                                else: 
                                    i.work = int(c_i)
                            else: #We check if the agent is working
                                w_i = i.does_Work()
                                i.work = int(w_i)
                        else:
                            pass
                    elif i.isol == 2:
                        w_i = i.does_Work()
                        i.work = int(w_i)
                    else:
                        i.work = 0
                        
                                                    
            #First round of contagion
            
            agents_to_update = []
            
            for x in range( 51 ):
                if len( S_rea[x] ) > 0:
                    Nx = 0
                    if x >= 32:
                        Nx = max( int( SystRM.ListSU[x].get_N() ), Nm )
                    else:
                        Nx = int( SystRM.ListSU[x].get_N() )
                    Ix = int( SystRM.ListSU[x].get_I() )
                    if Ix > 0:
                        for i in S_rea[x]:
                            if i.on == 1:
                                pr_i = random.random()
                                pc = 0
                                if i.home == 1 and (i.conf in [1,21,3]):
                                    if i.conf == 3:
                                        pc = 1-( (1-((B/Nx)*tau0*Ks))**Ix )
                                    else:
                                        pc = 1-( (1-((B/Nx)*tau0*Kns))**Ix )
                                else:
                                    pc = 1-( (1-((B/Nx)*tau0))**Ix )
                                if pr_i <= pc:
                                    agents_to_update.append( i )
            
                            
            #T2_morning return
            
            for i in CommRMT2:
                i.back_Home()
                
            
            #T2_afternoon commute
                
            CommRMT2_after = []
                                   
            if len( T2_afternoon ) > 0:
                for i in T2_afternoon:
                    if i.isol == 0:
                        if d <= 5 and (t not in [40,61,81,120,137]):
                            c_i = i.commute_toW()
                            if i.on == 0:
                                CommRMT2_after.append( i )
                                if i.comm == 1 or i.comm == 2:
                                    Mobility[ t ][0] += 0.5
                                i.work = 1
                            else: 
                                i.work = int(c_i)
                        elif d == 6 or (t in [40,120,137]):
                            if i.rama not in [11,15,16,21]:
                                c_i = i.commute_toW()
                                if i.on == 0:
                                    CommRMT2_after.append( i )
                                    if i.comm == 1 or i.comm == 2:
                                        Mobility[ t ][0] += 0.5
                                    i.work = 1
                                else: 
                                    i.work = int(c_i)
                            else: #We check if the agent is working
                                w_i = i.does_Work()
                                i.work = int(w_i)
                        elif d == 7 or (t in [61,81]):
                            if i.rama not in [3,6,10,11,12,13,14,15,16,19,20,21]:
                                c_i = i.commute_toW()
                                if i.on == 0:
                                    CommRMT2_after.append( i )
                                    if i.comm == 1 or i.comm == 2:
                                        Mobility[ t ][0] += 0.5
                                    i.work = 1
                                else: 
                                    i.work = int(c_i)
                            else: #We check if the agent is working
                                w_i = i.does_Work()
                                i.work = int(w_i)
                        else:
                            pass
                    elif i.isol == 2:
                        w_i = i.does_Work()
                        i.work = int(w_i)
                    else:
                        i.work = 0
                            
            
                            
            #Second round of contagion
               
            for x in range( 51 ):
                if len( S_rea[x] ) > 0:
                    Nx = 0
                    if x >= 32:
                        Nx = max( int( SystRM.ListSU[x].get_N() ), Nm )
                    else:
                        Nx = int( SystRM.ListSU[x].get_N())
                    Ix = int( SystRM.ListSU[x].get_I() )
                    if Ix > 0:
                        for i in S_rea[x]:
                            if i.on == 1:
                                pr_i = random.random()
                                pc = 0
                                if i.home == 1 and (i.conf in [1,21,3]):
                                    if i.conf == 3:
                                        pc = 1-( (1-((B/Nx)*tau0*Ks))**Ix )
                                    else:
                                        pc = 1-( (1-((B/Nx)*tau0*Kns))**Ix )
                                else:
                                    pc = 1-( (1-((B/Nx)*tau0))**Ix )
                                if pr_i <= pc:
                                    agents_to_update.append( i )
                                                        
                        
            #All the commuters return
           
            for i in CommRMT1 + CommRMT2_after:
                i.back_Home()
            
            
            #Third round of contagion
                                        
            for x in range( 51 ):
                if len( S_rea[x] ) > 0:
                    Nx = 0
                    if x >= 32:
                        Nx = max( int( SystRM.ListSU[x].get_N() ), Nm )
                    else:
                        Nx = int( SystRM.ListSU[x].get_N())
                    Ix = int( SystRM.ListSU[x].get_I() )
                    if Ix > 0:
                        for i in S_rea[x]:
                            if i.on == 1:
                                pr_i = random.random()
                                pc = 0
                                if i.home == 1 and (i.conf in [1,21,3]):
                                    if i.conf == 3:
                                        pc = 1-( (1-((B/Nx)*tau1*Ks))**Ix )
                                    else:
                                        pc = 1-( (1-((B/Nx)*tau1*Kns))**Ix )
                                else:
                                    pc = 1-( (1-((B/Nx)*tau1))**Ix )
                                if pr_i <= pc:
                                    agents_to_update.append( i )
                                
            
            for i in agents_to_update:
                if i.status == 0:
                    i.update_Status( 1, q )
                    if i.home == 1:
                        I_rea.append( i )
                        
                        S_rea[i.su].remove( i )
                        if i.replica != None:
                            S_rea[ i.replica.su ].remove( i.replica )
                    else:
                        I_rea.append( i.replica )
                        
                        S_rea[i.replica.su].remove( i.replica )
                        S_rea[i.su].remove( i )
           
                        
            #Writing information rea, t, to files (if calib 1=1)
            #Output: for rea, day for each ident type (row) the distribution of 
            #agents among the compartiments {S,I,R}x{work=0,work=1,work=2}
                        
            #DistDay = np.zeros( ( 19584,9 ) )
            #S_extended = []
            #for x in range( 51 ):
            #    S_extended += S_rea[x]
            #UpdateDay = S_extended + I_rea + R_rea
            #for i in UpdateDay:
            #    if i.home == 1:
            #        DistDay[i.ident][i.status + int(3*i.work)] += 1
            #np.savetxt("S"+str(sim)+"_rea_"+str(rea)+"_day_"+str(t)+".csv",DistDay,delimiter=",",fmt="%s")
            #del( DistDay )
            
            
            Detected_RM [t] = Detected_RM_Cum [0]
            
            if t == 6:
                Detected_RM_Calib [0] += Detected_RM_Cum [0]
            elif t == 13:
                Detected_RM_Calib [1] += Detected_RM_Cum [0]
            elif t == 20:
                Detected_RM_Calib [2] += Detected_RM_Cum [0]
            elif t == 27:
                Detected_RM_Calib [3] += Detected_RM_Cum [0]
            elif t == 34:
                Detected_RM_Calib [4] += Detected_RM_Cum [0]
            elif t == 41:
                Detected_RM_Calib [5] += Detected_RM_Cum [0]
            elif t == 48:
                Detected_RM_Calib [6] += Detected_RM_Cum [0]
            elif t == 55:
                Detected_RM_Calib [7] += Detected_RM_Cum [0]                
            elif t == 62:
                Detected_RM_Calib [8] += Detected_RM_Cum [0]
            elif t == 69:
                Detected_RM_Calib [9] += Detected_RM_Cum [0]
            elif t == 76:
                Detected_RM_Calib [10] += Detected_RM_Cum [0]
            elif t == 83:
                Detected_RM_Calib [11] += Detected_RM_Cum [0]
            elif t == 90:
                Detected_RM_Calib [12] += Detected_RM_Cum [0]
            elif t == 97:
                Detected_RM_Calib [13] += Detected_RM_Cum [0]
            elif t == 104:
                Detected_RM_Calib [14] += Detected_RM_Cum [0]
            elif t == 111:
                Detected_RM_Calib [15] += Detected_RM_Cum [0]
            elif t == 118:
                Detected_RM_Calib [16] += Detected_RM_Cum [0]
            elif t == 125:
                Detected_RM_Calib [17] += Detected_RM_Cum [0]
            elif t == 132:
                Detected_RM_Calib [18] += Detected_RM_Cum [0]
            elif t == 139:
                Detected_RM_Calib [19] += Detected_RM_Cum [0]
            elif t == 146:
                Detected_RM_Calib [20] += Detected_RM_Cum [0]
            elif t == 153:
                Detected_RM_Calib [21] += Detected_RM_Cum [0]
            else:
                pass
            
            
            if t == 151:
                DeadRM[rea-a][0] = Fall
                                
            t += 1
            if d == 7:
                d = 1
            else:
                d += 1
                
        
                 
        
        #np.savetxt("S"+str(sim)+"_Mob_Tot_rea_"+str(rea)+".csv",Mobility,delimiter=",",fmt="%s")
        #np.savetxt("S"+str(sim)+"_Detected_rea_"+str(rea)+".csv",Detected_RM,delimiter=",",fmt="%s")
        
        
        del(Mobility)
        
        SystRM.reset_Realization( )
    
        

    np.savetxt("S"+str(sim)+"_Dead_a_"+str(a)+"_b_"+str(b)+".csv",DeadRM,delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_Calib_data_a_"+str(a)+"_b_"+str(b)+".csv",Detected_RM_Calib/(b-a),delimiter=",",fmt="%s")



## Parameters setup for calibration


data1_s='Data1_MP.csv'
data2_s='Data2_MP.csv'

a_s = 0 #first number of realization
b_s = 5 #total number of realizations 
tmax_s = 154
Kns_s = 0.56
pD_s = 1 
tau0_s = 6.0/24.0
tau1_s = 6.0/24.0 #initial value
q_s = 0.996
Nm_s = 10530 #this affects only small villages outside Great Santiago, ad-hoc to this implementation.
             #to remove this assumption, set Nm_s=0.

situation_s = 0 #always =0 for calibrating


               

t1=time.time()

print("Initiating simulation")

#Load and prepare real data for estimating pD and computing NSE

RealRM = np.loadtxt( 'RealDRM.csv', delimiter=",") #real data RM
RealRM_last = RealRM[-1] #last observation

IncDistRM = np.zeros((1,21)) #incremental RM
RMt = RealRM[-1]-RealRM[0]
for j in range(21):
    IncDistRM[0,j]=(RealRM[j+1]-RealRM[j])/RMt

RealCom = np.loadtxt( 'RealDCom.csv', delimiter=",") #real data cumulative by comuna

RealComInc = np.zeros( (21,51) ) #incremental by comuna real data
for x1 in range( 21 ):
    RealComInc[x1] = RealCom[x1+1]-RealCom[x1]


#Create initial system

RM = SystemRM( data1=data1_s, data2=data2_s )
RM.InitialSystem()

#set random seed
rseed_s = 3  #valid for replicating this example

random.seed( rseed_s )

#simulation

sim_s = 0 #first number of this particular set of configurations

for B_s in [0.24,0.25,0.26,0.27]:
    for pA_s in [0.13,0.15,0.17,0.19]:
        for Ks_s in [0.1, 0.15]:
                        
            outres = np.zeros( (1,9) )
            outres[:,0] = sim_s
            outres[:,1] = B_s
            outres[:,2] = pA_s
            outres[:,3] = Ks_s
            outres[:,4] = q_s
            
            LaborEpiRM( sim=sim_s, a=a_s , b=b_s, tmax=tmax_s, B=B_s, Nm=Nm_s, Ks=Ks_s, Kns=Kns_s, pD=pD_s, pA=pA_s, q=q_s, tau0=tau0_s, tau1=tau1_s, SystRM=RM, situation=situation_s ) 
            RM.reset_Realization( )
                                           
            Cal = np.loadtxt( "S"+str(sim_s)+"_Calib_data_a_"+str(a_s)+"_b_"+str(b_s)+".csv", delimiter=",")
            RMCum = np.sum( Cal, axis=1 )
            RMCum_last = RMCum[-1]
            
            Cal1 = np.zeros( ( 22,51 ) )
            
            if RMCum_last > RealRM_last: #if simulated cumulative are at least the observed detected
                pDsim = pDest( RMCum, RealRM ) #detection probability estimation
                outres[:,5] = pDsim
                Cal1 = pDsim*Cal  #detected cumulative data by comuna, simulation
            else: #When detected simulated are less than real cases, pD is loaded as 999 (as signal of not valid), but assumed to be 1 for computing
                outres[:,5] = 999
                Cal1 = 1*Cal  #detected cumulative data by comuna, simulation    
                    
            NSE_AcumCom = NSE( Cal1, RealCom ) #compute NSE coefficient cummulative, by comuna
            outres[:,6] = NSE_AcumCom
                    
            #compute incremental, detected simulated, by comuna
            CalInc = np.zeros( (21,51) ) 
            for y in range( 21 ):
                CalInc[y] = Cal1[y+1]-Cal1[y]
        
            NSE_IncCom = NSE( CalInc, RealComInc ) #compute NSE coefficient, incremental, by comuna
            outres[:,7] = NSE_IncCom
                    
            #compute distribution of increments SMR
            IncDistRMCal = np.zeros((1,21)) 
            CalRMt = RMCum[-1]-RMCum[0]
            for i in range(21):
                IncDistRMCal[0,i] = (RMCum[i+1]-RMCum[i])/CalRMt
                    
            NSE_slopeRM = NSE( IncDistRMCal, IncDistRM) #compute NSE coefficient, incremental, distribution SMR
            outres[ :, 8 ] = NSE_slopeRM
                    
            np.savetxt("S"+str(sim_s)+"_Outcome_a_"+str(a_s)+"_b_"+str(b_s)+".csv",outres,delimiter=",",fmt="%s")    
            
            print ("fin"+str(sim_s))
            
            sim_s += 1


print ( "End", "processing time seconds",time.time()-t1 )



