# SIRLabor_SMR_Chile_Covid19
Contain all the codes, file inputs, and information related to Fosco &amp; Zurita (2021) - Assessing the short-run effects of lockdown policies on economic activity, with an application to the Santiago metropolitan area.

It is an ad-hoc simulation model for the Santiago metropolitan region, Chile. It models the coevolution of covid-19 and the effects on the labor market of the mitigation measures adopted between March 1st and August 1st (scenario S0). It also contains two other scenarios (S1 and S2). S1 simulates the counterfactual without any measure, and S2 follows S0 until March 26. On March 27 (morning) when the first targeted lockdown took place in the region (seven comunas), we assume a full lockdown.

SIRLaborMPSim.py  (Python) : to replicate the three scenarios S0, S1, S2 of our paper.
This program contains all the classes, processes, etc. and at the end of it, the parameters setting (prepared for replicating S0 by default).
It uses two input files, Data1_MP.csv and Data2_MP.csv. Please, locate these files and the program in the same folder.
Warning: this program delivers outcomes files for each day and realization. 
The post-processing of the raw data can be done with OutcomeProcessSIRLabor.py. It uses the the file Data2_MP.csv 

Data1_MP.csv contains the estimated probabilities by municipality (comuna) and economic (sector) of working in an essential activity (the initial values).
Data2_Mp.csv contains 19854 types of agents, the number of each type, and characteristics. The description of each can be found in the paper. Data elaborated from Encuesta Nacional de Empleo dic. 2019 (https://www.ine.cl/docs/default-source/ocupacion-y-desocupacion/bbdd).

SIRLaborMPCal.py (Python): to replicate the simulation of several configurations sequentially (used in the first step of the manual calibration).
It uses four input files: Data1_MP.csv, Data2_MP.csv, RealDCom.csv, RealDRM.csv. The last two files contain real data by comuna (excluding Alhué) 
of cumulative cases, 22 epidemiological weeks, and cumulative cases for SRM (excluding Alhué) (from Retrieved from https://github.com/MinCiencia/Datos-COVID19, Nov 03 2020).

