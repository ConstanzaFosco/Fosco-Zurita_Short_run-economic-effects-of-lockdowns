# SIRLabor_SMR_Chile_Covid19

# Constanza Fosco &amp; Felipe Zurita (2021) - Assessing the short-run effects of lockdown policies on economic activity, with an application to the Santiago Metropolitan Region, Chile.

**Programs, input files, and results.**

This is a simulation model applied to the Santiago metropolitan region, Chile. It models the coevolution of covid-19 and the effects on the labor market of the containment measures adopted between March 1st and August 1st (scenario S0). The model is data-driven, with a metapopulation spatial structure, and agent-based.
Two other scenarios (S1 and S2) can be simulated. S1: a counterfactual without any measure. S2: follows S0 until March 26; on March 27 (morning) when the first targeted lockdown took place in the region (seven comunas), it is assumed a full lockdown instead.
The information about "cuarentenas" (lockdown) was retrieved from https://github.com/MinCiencia/Datos-COVID19, Product 29.

**SIRLaborMP.py** for simulating the three scenarios S0, S1, S2 of the paper.
This program contains all the classes, processes, etc. and at the end of it, the parameters setting.
Requires two input files, Data1_MP.csv and Data2_MP.csv. Please, locate these files and the program in the same folder.
The program delivers one outcome file for each day and realization. Each file "SX_rea_u_day_v.csv" (X=scenario, u=number of realization, v=day simulated)
contains a matrix of dimension (19584,9). Each row represent (in order) a type of agent (see the description of Data2_MP.csv below). Columns are the number
of agents in the compartiments {not working, working on-site, teleworking}x{susceptible, infected, removed}.
The post-processing of the raw data can be done with **OutcomeProcessSIRLabor.py**. It requires the file Data2_MP.csv and links the raw outcome to the full set
of characteristics.

**Data1_MP.csv** contains the estimated probabilities by municipality (comuna) and economic sector of working in an essential activity - own elaboration, based on the official definitions of Chilean authorities (Instructivo Cuarentena) and firms statistics by municipality (https://www.sii.cl/sobre_el_sii/estadisticas_de_empresas.html).

**Data2_MP.csv** includes 19584 types of agents, the number of each type, and characteristics. The description of each can be found in the paper. Data elaborated based on the Encuesta Nacional de Empleo, INE, dic. 2019 (https://www.ine.cl/docs/default-source/ocupacion-y-desocupacion/bbdd), Encuesta Encuesta Suplementaria de Ingresos, INE, 2018 (https://www.ine.cl/estadisticas/sociales/ingresos-y-gastos/encuesta-suplementaria-de-ingresos), Nominal remuneration index (base 2016=100), National according to economic section (CIIU4.CL 2012), monthly, INE (https://stat.ine.cl), Proyecciones de Población, INE (https://www.ine.cl/estadisticas/sociales/demografia-y-vitales/proyecciones-de-poblacion), Census data 2017, INE (https://www.ine.cl/estadisticas/sociales/censos-de-poblacion-y-vivienda/poblacion-y-vivienda).

**VariablesData2.csv**: brief description of variables included in Data2_MP.csv.

**Municipalities.csv**: the list of municipalities (comunas) and their id as spatial units.

**Scenarios S0_S1_S2 Outcomes.zip**: processed outcomes of each scenario. 

**Calibration Outcomes.zip**: explanation and outcomes of the two calibration steps (not automatic).

**Rt_CaseReproductive Outcomes.zip**: outcomes of effective (case) reproductive number, includes a variant of the main code that allows for the counting of secondary cases.

Simulation software: **Python 3.7.6**, within the Anaconda open-source distribution package (https://docs.anaconda.com/). 
