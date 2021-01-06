# -*- coding: utf-8 -*-
"""
Created on Thu May 28 10:32:15 2020

@author: Constanza

This program takes the output files of SIRLabourRM and computes the results
"""
import numpy as np
import time


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

def Mean_Day( sim, realizations, days ):
    #mean of all realizations for each kind of clone in data2
    #sim = number of simulation
    #each file represent a day, within each file: each row represent a row in data 2, columns:
    #S_N, I_N, R_N; S_P, I_P, R_P; S_T, I_T, R_T.
    #N: non working, P: presentially working, T: teleworking
    #Useful for averages of linear functions
    for day in range( days ):
        Matrix_day = np.zeros( ( 19584,9 ) )
        for rea in range( realizations ):
            m_rea = np.loadtxt( "S"+str(sim)+"_rea_"+str(rea)+"_day_"+str(day)+".csv", delimiter=",")
            Matrix_day = Matrix_day + m_rea
        Mean_day = Matrix_day/float( realizations )
        if day <=9:
            np.savetxt("S"+str(sim)+"_P_00"+str(day)+".csv",Mean_day,delimiter=",",fmt="%s")
        elif day>=10 and day<=99:
            np.savetxt("S"+str(sim)+"_P_0"+str(day)+".csv",Mean_day,delimiter=",",fmt="%s")
        else:
            np.savetxt("S"+str(sim)+"_P_"+str(day)+".csv",Mean_day,delimiter=",",fmt="%s")

def SeriesRM( sim, realizations, days, data2 ):
    #returns time series, mean of time series and percentiles of time series
    #of number of agents in each health status
    #files SX_MH_uv for X=sim, u=0,1,2,3 (S,I,R,I+R), v=0 all the realizations, v=1 mean and std
    
    Y = np.loadtxt( data2, delimiter=",", skiprows=1 ) #file with data for each group of clones identified with i.ident
    
    RM_S = np.zeros( ( days, realizations ) )
    RM_I = np.zeros( ( days, realizations ) )
    RM_R = np.zeros( ( days, realizations ) )
    RM_Cum = np.zeros( ( days, realizations ) )
    RM_S_summary = np.zeros( ( days, 3 ) )
    RM_I_summary = np.zeros( ( days, 3 ) )
    RM_R_summary = np.zeros( ( days, 3 ) )
    RM_Cum_summary = np.zeros( ( days, 3 ) )
    
    RM_NR = np.zeros( ( days, realizations ) ) #residentes no trabajan
    RM_PTR = np.zeros( ( days, realizations ) ) #residentes trabajan
    RM_PR = np.zeros( ( days, realizations ) ) #residentes trabajan presencialmente
    RM_TR = np.zeros( ( days, realizations ) ) #residentes teletrabajan
    RM_PMR = np.zeros( ( days, realizations ) ) #residentes trabajan presencialmente y se movilizan
    RM_NRR = np.zeros( ( days, realizations ) ) #residentes no trabajan y están en riesgo de no percibir ingreso
    RM_PTWP = np.zeros( ( days, realizations ) ) #trabajan, con lugar de trabajo en la RM (excluye comm=3)
    RM_NWP = np.zeros( ( days, realizations ) ) #no trabajan, con lugar de trabajo en la RM (excluye comm=3)
    
    WRM_NR = np.zeros( ( days, realizations ) ) #W residentes no trabajan
    WRM_PTR = np.zeros( ( days, realizations ) ) #W residentes trabajan
    WRM_PR = np.zeros( ( days, realizations ) ) #W residentes trabajan presencialmente
    WRM_TR = np.zeros( ( days, realizations ) ) #W residentes teletrabajan
    WRM_NRR = np.zeros( ( days, realizations ) ) #W residentes no trabajan y están en riesgo de no percibir ingreso
    WRM_PTWP = np.zeros( ( days, realizations ) ) #W trabajan, con lugar de trabajo en la RM (excluye comm=3)
    WRM_NWP = np.zeros( ( days, realizations ) ) #W no trabajan, con lugar de trabajo en la RM (excluye comm=3)
    
    Activos = Y[:,19]
    W = Y[:,17]
    Risk = Y[:,16]
    WorkInRM = Y[:,20]
    Commuter = Y[:,21]
    TotalLRes = sum(Y[:,0]*Activos*15)
    TotalWRes = sum(Y[:,0]*Activos*W*15)
    TotalLWP = sum(Y[:,0]*Activos*WorkInRM*15)
    TotalWWP = sum(Y[:,0]*Activos*WorkInRM*W*15)
                         
    for day in range( days ):
        for rea in range( realizations ):
            m_rea = np.loadtxt( "S"+str(sim)+"_rea_"+str(rea)+"_day_"+str(day)+".csv", delimiter=",")
            m_rea_sum = np.sum( m_rea, axis = 0 )
            S_dr = (m_rea_sum[0]+m_rea_sum[3]+m_rea_sum[6])*15
            I_dr = (m_rea_sum[1]+m_rea_sum[4]+m_rea_sum[7])*15
            R_dr = (m_rea_sum[2]+m_rea_sum[5]+m_rea_sum[8])*15
            RM_S[day][rea] += S_dr
            RM_I[day][rea] += I_dr
            RM_R[day][rea] += R_dr
            RM_Cum[day][rea] += I_dr+R_dr
            
            N = np.sum(m_rea[:,0:3],axis=1)*Activos*15
            P = np.sum(m_rea[:,3:6],axis=1)*Activos*15
            T = np.sum(m_rea[:,6:9],axis=1)*Activos*15
            PT = P + T
            PM = P*Commuter
            NR = N*Risk
            PTWP = PT*WorkInRM
            NWP = N*WorkInRM
            
            
            
            N_dr = sum(N)
            P_dr = sum(P)
            T_dr = sum(T)
            PT_dr = sum(PT)
            
            PM_dr = sum(PM)
            NR_dr = sum(NR)
            PTWP_dr = sum(PTWP)
            NWP_dr = sum(NWP)
            
            WN_dr = sum(N*W)
            WP_dr = sum(P*W)
            WT_dr = sum(T*W)
            WPT_dr = WP_dr + WT_dr
            WNR_dr = sum(NR*W)
            WPTWP_dr = sum( PTWP*W )
            WNWP_dr = sum( NWP*W )
            
            RM_NR[day][rea] += N_dr
            RM_PR[day][rea] += P_dr
            RM_TR[day][rea] += T_dr
            RM_PTR[day][rea] += PT_dr
            RM_PMR[day][rea] += PM_dr
            RM_NRR[day][rea] += NR_dr
            RM_PTWP[day][rea] += PTWP_dr
            RM_NWP[day][rea] += NWP_dr
            
            WRM_NR[day][rea] += WN_dr
            WRM_PTR[day][rea] += WPT_dr
            WRM_PR[day][rea] += WP_dr
            WRM_TR[day][rea] += WT_dr
            WRM_NRR[day][rea] += WNR_dr
            WRM_PTWP[day][rea] += WPTWP_dr
            WRM_NWP[day][rea] += WNWP_dr
            
      
    RM_S_summary[:,0] = np.mean( RM_S, axis= 1 )
    RM_S_summary[:,1] = np.percentile( RM_S, 5, axis= 1 )
    RM_S_summary[:,2] = np.percentile( RM_S, 95, axis= 1 )
    RM_I_summary[:,0] = np.mean( RM_I, axis= 1 )
    RM_I_summary[:,1] = np.percentile( RM_I, 5, axis= 1 )
    RM_I_summary[:,2] = np.percentile( RM_I, 95, axis= 1 )
    RM_R_summary[:,0] = np.mean( RM_R, axis= 1 )
    RM_R_summary[:,1] = np.percentile( RM_R, 5, axis= 1 )
    RM_R_summary[:,2] = np.percentile( RM_R, 95, axis= 1 )
    RM_Cum_summary[:,0] = np.mean( RM_Cum, axis= 1 )
    RM_Cum_summary[:,1] = np.percentile( RM_Cum, 5, axis= 1 )
    RM_Cum_summary[:,2] = np.percentile( RM_Cum, 95, axis= 1 )
    
    
    np.savetxt("S"+str(sim)+"_MH_00.csv",RM_S,delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_MH_10.csv",RM_I,delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_MH_20.csv",RM_R,delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_MH_30.csv",RM_Cum,delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_MH_01.csv",RM_S_summary,delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_MH_11.csv",RM_I_summary,delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_MH_21.csv",RM_R_summary,delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_MH_31.csv",RM_Cum_summary,delimiter=",",fmt="%s")
    
    np.savetxt("S"+str(sim)+"_ML_"+str(0)+str(1)+str(0)+".csv", RM_NR, delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_ML_"+str(0)+str(0)+str(0)+".csv", RM_PTR, delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_ML_"+str(0)+str(2)+str(0)+".csv", RM_PR, delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_ML_"+str(0)+str(3)+str(0)+".csv", RM_TR, delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_ML_"+str(0)+str(4)+str(0)+".csv", RM_PMR, delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_ML_"+str(0)+str(5)+str(0)+".csv", RM_NRR, delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_ML_"+str(0)+str(6)+str(0)+".csv", RM_PTWP, delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_ML_"+str(0)+str(7)+str(0)+".csv", RM_NWP, delimiter=",",fmt="%s")
    
    np.savetxt("S"+str(sim)+"_ML_"+str(1)+str(1)+str(0)+".csv", WRM_NR, delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_ML_"+str(1)+str(0)+str(0)+".csv", WRM_PTR, delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_ML_"+str(1)+str(2)+str(0)+".csv", WRM_PR, delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_ML_"+str(1)+str(3)+str(0)+".csv", WRM_TR, delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_ML_"+str(1)+str(5)+str(0)+".csv", WRM_NRR, delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_ML_"+str(1)+str(6)+str(0)+".csv", WRM_PTWP, delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_ML_"+str(1)+str(7)+str(0)+".csv", WRM_NWP, delimiter=",",fmt="%s")
    
    np.savetxt("S"+str(sim)+"_ML_"+str(0)+str(1)+str(1)+".csv", (RM_NR/float(TotalLRes))*100, delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_ML_"+str(0)+str(0)+str(1)+".csv", (RM_PTR/float(TotalLRes))*100, delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_ML_"+str(0)+str(2)+str(1)+".csv", (RM_PR/float(TotalLRes))*100, delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_ML_"+str(0)+str(3)+str(1)+".csv", (RM_TR/float(TotalLRes))*100, delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_ML_"+str(0)+str(4)+str(1)+".csv", (RM_PMR/float(TotalLRes))*100, delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_ML_"+str(0)+str(5)+str(1)+".csv", (RM_NRR/float(TotalLRes))*100, delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_ML_"+str(0)+str(6)+str(1)+".csv", (RM_PTWP/float(TotalLWP))*100, delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_ML_"+str(0)+str(7)+str(1)+".csv", (RM_NWP/float(TotalLWP))*100, delimiter=",",fmt="%s")
    
    np.savetxt("S"+str(sim)+"_ML_"+str(1)+str(1)+str(1)+".csv", (WRM_NR/float(TotalWRes))*100, delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_ML_"+str(1)+str(0)+str(1)+".csv", (WRM_PTR/float(TotalWRes))*100, delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_ML_"+str(1)+str(2)+str(1)+".csv", (WRM_PR/float(TotalWRes))*100, delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_ML_"+str(1)+str(3)+str(1)+".csv", (WRM_TR/float(TotalWRes))*100, delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_ML_"+str(1)+str(5)+str(1)+".csv", (WRM_NRR/float(TotalWRes))*100, delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_ML_"+str(1)+str(6)+str(1)+".csv", (WRM_PTWP/float(TotalWWP))*100, delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_ML_"+str(1)+str(7)+str(1)+".csv", (WRM_NWP/float(TotalWWP))*100, delimiter=",",fmt="%s")

    del(Y,RM_S,RM_I,RM_R,RM_Cum,RM_S_summary,RM_I_summary,RM_R_summary,RM_Cum_summary,
        RM_NR,RM_PTR,RM_PR,RM_TR,RM_PMR,RM_NRR,RM_PTWP,RM_NWP,
        WRM_NR,WRM_PTR,WRM_PR,WRM_TR,WRM_NRR,WRM_PTWP,WRM_NWP) 


def Detected_Series( sim, realizations, days, RMdata ):
    #returns: (i) time series of all realizations of DETECTED cumulative cases of RM
    #(ii) mean, percentiles 5% and 95%
    #sim: simulation code (ALWAYS 001)
    #realizations: number of realizations
    #days: number of days simulated
    #RMdata: observed (real) weekly serie of cumulative detected cases RM
    #pD estimated as the slope of RealDetected=pD*Mean(DetectedSimulated,pD=1)
    #(USING RM WEEKLY CUMULATIVE RMdata)
    ReasDet = np.zeros( ( days,realizations ) )
    DetRMsumm = np.zeros( ( days, 3) )
    DetecComunas = np.zeros( ( days, 51 ))
    for rea in range( realizations ):
        m_rea = np.loadtxt( "S"+str(sim)+"_Detected_rea_"+str(rea)+".csv", delimiter=",")
        DetecComunas = DetecComunas + m_rea
        ReasDet[:,rea] = np.sum(m_rea, axis = 1 )
    RMMeanpD1 = np.mean(ReasDet, axis = 1 )
    RMSimpD1 = np.zeros(22)
    
    w = 6
    z = 0
    while w < days:
        RMSimpD1[z] = RMMeanpD1[w]
        z += 1
        w += 7
    RMReal = np.loadtxt( RMdata, delimiter=",")
    
    pD_estim = pDest( RMSimpD1, RMReal )
    print("pD",pD_estim)
    
    ReasDet_pD = pD_estim * ReasDet
    DetRMsumm[:,0] = np.mean( ReasDet_pD, axis= 1 )
    DetRMsumm[:,1] = np.percentile( ReasDet_pD, 5, axis= 1 )
    DetRMsumm[:,2] = np.percentile( ReasDet_pD, 95, axis= 1 )
    
    DetecComunas1 = pD_estim * (DetecComunas/float(realizations)) #mean detected by comuna
    
    DetecComunas2 = np.zeros( ( 22,51 ) )
    
    w1 = 6
    z1 = 0
    while w1 < days:
        DetecComunas2[z1] = DetecComunas1[w1]
        w1 += 7
        z1 += 1
    
    np.savetxt( "S"+str(sim)+"_M_detected.csv",ReasDet_pD,delimiter=",",fmt="%s" )
    np.savetxt( "S"+str(sim)+"_M_detectedSum.csv",DetRMsumm,delimiter=",",fmt="%s" )
    np.savetxt( "S"+str(sim)+"_CH_detectedFullserie.csv",DetecComunas1,delimiter=",",fmt="%s" )
    np.savetxt( "S"+str(sim)+"_CH_detectedPlanilla.csv",DetecComunas2,delimiter=",",fmt="%s" )
    del(ReasDet, ReasDet_pD, DetRMsumm, DetecComunas,DetecComunas1,DetecComunas2)


def Detected_Series2( sim, realizations, days, RMdata ):
    #returns: (i) time series of all realizations of DETECTED cumulative cases of RM
    #(ii) mean, percentiles 5% and 95%
    #sim: simulation code (ALWAYS 001)
    #realizations: number of realizations
    #days: number of days simulated
    #RMdata: observed (real) weekly serie of cumulative detected cases RM
    #pD estimated as the slope of RealDetected=pD*Mean(DetectedSimulated,pD=1)
    #(USING RM WEEKLY CUMULATIVE RMdata)
    pdEstim=np.zeros(realizations)
    RMReal = np.loadtxt( RMdata, delimiter=",")
    for rea in range( realizations ):
        m_rea = np.loadtxt( "S"+str(sim)+"_Detected_rea_"+str(rea)+".csv", delimiter=",")
        RMDetrea = np.sum(m_rea, axis = 1 )
        RMWeekDetrea=np.zeros(22)
        w=6
        z=0
        while w < days:
            RMWeekDetrea[z]=RMDetrea[w]
            z += 1
            w += 7
        pD_estimrea = pDest(RMWeekDetrea,RMReal)
        pdEstim[rea]=pD_estimrea
    np.savetxt( "S"+str(sim)+"pDestimados.csv",pdEstim,delimiter=",",fmt="%s" )
    

def incrementComuna( sim, realdatcom ):
    Simdata = np.loadtxt("S"+str(sim)+"_CH_detectedPlanilla.csv", delimiter=",")  
    Realdata = np.loadtxt(realdatcom, delimiter=",")  
    IncrSim = np.zeros( (21,51 ) )
    IncrReal = np.zeros( (21,51 ) )
    for x in range( 21 ):
        IncrSim[x,:] = Simdata[x+1,:]-Simdata[x,:]
        IncrReal[x,:] = Realdata[x+1,:]-Realdata[x,:]
    BothIncr = np.zeros((21,153))
    y=1
    z=2
    w=0
    for x in range( 51 ):
        BothIncr[:,y] =IncrSim[:,w]
        BothIncr[:,z] =IncrReal[:,w]
        w += 1
        y += 3
        z += 3
        
    
    np.savetxt( "AABIncrementosCom.csv",BothIncr,delimiter=",",fmt="%s" )
        
        
        
def Comunas_HealthSeries( sim, data2, days ):
    #Time series for each comuna of the "stock" of individuals within each health status at each t, and RM
    #Call after Mean_Day
    Y = np.loadtxt( data2, delimiter=",", skiprows=1 ) #file with data for each group of clones identified with i.ident
    rowsY, colsY = Y.shape
    S_evo = np.zeros( ( days, 51 ) )
    I_evo = np.zeros( ( days, 51 ) )
    R_evo = np.zeros( ( days, 51 ) )
    C_evo = np.zeros( (days, 51 ) )
    MS_evo = np.zeros( ( days, 1 ) )
    MI_evo = np.zeros( ( days, 1 ) )
    MR_evo = np.zeros( ( days, 1 ) )
    MC_evo = np.zeros( ( days, 1 ) )
    
    for day in range( days ):
        P = None
        if day <= 9:
            P = np.loadtxt("S"+str(sim)+"_P_00"+str(day)+".csv",delimiter=",")
        elif day >=10 and day <=99:
            P = np.loadtxt("S"+str(sim)+"_P_0"+str(day)+".csv",delimiter=",") 
        else:
            P = np.loadtxt("S"+str(sim)+"_P_"+str(day)+".csv",delimiter=",") 
        for x in range( rowsY ):
            CutInd = int(Y[x][1])
            
            Sx = (P[x][0]+P[x][3]+P[x][6])*15
            Ix = (P[x][1]+P[x][4]+P[x][7])*15
            Rx = (P[x][2]+P[x][5]+P[x][8])*15
           
            S_evo[day][CutInd] += Sx
            I_evo[day][CutInd] += Ix
            R_evo[day][CutInd] += Rx
            C_evo[day][CutInd] += Ix + Rx 
            
            MS_evo[day][0] += Sx
            MI_evo[day][0] += Ix
            MR_evo[day][0] += Rx
            MC_evo[day][0] += Ix + Rx
        
    np.savetxt("S"+str(sim)+"_CH_"+str(0)+".csv", S_evo, delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_CH_"+str(1)+".csv", I_evo, delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_CH_"+str(2)+".csv", R_evo, delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_CH_"+str(3)+".csv", C_evo, delimiter=",",fmt="%s")
        
    np.savetxt("S"+str(sim)+"_MH_"+str(0)+".csv", MS_evo, delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_MH_"+str(1)+".csv", MI_evo, delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_MH_"+str(2)+".csv", MR_evo, delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_MH_"+str(3)+".csv", MC_evo, delimiter=",",fmt="%s")
            
    del( Y,S_evo,I_evo,R_evo,MS_evo,MI_evo,MR_evo,C_evo, MC_evo)



def Percent_Calculus( matr, vect):
    #matr: matrix of data, vect: array with correspondent column total (it is not a sum, but the total of that case)
    rowsM, colsM = matr.shape
    New_matr = np.zeros( (rowsM, colsM) )
    for x in range( colsM ):
        for y in range( rowsM ):
            New_matr[y][x] = ( matr[y][x]/float(vect[x]) )*100
    return New_matr
        
    

def LabourSeriesComuna( sim, data2, days ):
    #Time series by Comuna of all the relevant labour variables, with the exception of Atkinson Index
    #Call after Mean_Day
    Y = np.loadtxt( data2, delimiter=",", skiprows=1 ) #file with data for each group of clones identified with i.ident
    rowsY, colsY = Y.shape
    
    #The following correspond to comuna of residence_then should be used for welfare considerations
    #People Comunas Residence
    N_ComR = np.zeros(( days,  51 ) )         #Employed, not working
    PT_ComR = np.zeros( ( days,  51 ) )       #Employed, working
    P_ComR = np.zeros( ( days,  51 ) )        #Employed, presentially working
    T_ComR = np.zeros( ( days,  51 ) )        #Employed, teleworking
    PM_ComR = np.zeros( ( days,  51 ) )       #Employed, presentially working and moving
    NR_ComR = np.zeros( ( days,  51 ) )       #Employed, not working, and in "risk": jobcat=1,2,7; or 3,4,5,6 informal sector/home sector
    #Wages Comunas Residence
    WN_ComR = np.zeros(( days,  51 ) )         #Employed, not working
    WPT_ComR = np.zeros( ( days,  51 ) )       #Employed, working
    WP_ComR = np.zeros( ( days,  51 ) )        #Employed, presentially working
    WT_ComR = np.zeros( ( days,  51 ) )        #Employed, teleworking
    WNR_ComR = np.zeros( ( days,  51 ) )       #Employed, not working, and in "risk": jobcat=1,2,7; or 3,4,5,6 informal sector/home sector
    
    #The following correspond to comuna of WORKPLACE_then should be used AS A PROXY OF PRODUCTION
    #EXCLUDING PEOPLE WORKING OUTSIDE THE RM (COMM==3)
    #People Comunas WORKPLACE
    N_ComWP = np.zeros(( days,  51 ) )         #Employed, not working
    PT_ComWP = np.zeros( ( days,  51 ) )       #Employed, working
    
    #Wages Comunas WORKPLACE
    WN_ComWP = np.zeros(( days,  51 ) )         #Employed, not working
    WPT_ComWP = np.zeros( ( days,  51 ) )       #Employed, working
    
    for day in range( days ):
        P = None
        if day <= 9:
            P = np.loadtxt("S"+str(sim)+"_P_00"+str(day)+".csv",delimiter=",")
        elif day >=10 and day <=99:
            P = np.loadtxt("S"+str(sim)+"_P_0"+str(day)+".csv",delimiter=",") 
        else:
            P = np.loadtxt("S"+str(sim)+"_P_"+str(day)+".csv",delimiter=",") 
        for x in range( rowsY ):
            if int(Y[x][3]) == 1: #activ
                CutHome = int(Y[x][1])
                CutWork = int(Y[x][12])
                Wx = Y[x][17]
                Nx = ( np.sum( P[x,[0,1,2]]))*15
                PTx = ( np.sum( P[x,[3,4,5,6,7,8]]))*15
                Px = ( np.sum( P[x,[3,4,5]]))*15
                Tx = ( np.sum( P[x,[6,7,8]]))*15
                #People Comunas Residence
                N_ComR[day][CutHome] += Nx
                PT_ComR[day][CutHome] += PTx
                P_ComR[day][CutHome] += Px
                T_ComR[day][CutHome] += Tx
                                
                #Wages Comunas Residence
                WN_ComR[day][CutHome] += Nx * Wx
                WPT_ComR[day][CutHome] += PTx * Wx
                WP_ComR[day][CutHome] += Px * Wx
                WT_ComR[day][CutHome] += Tx * Wx
                
                #People and income, non working in risk (comunas residence)
                if (int(Y[x][16])) == 1: #worker in risk in the case that she didn't work
                    NR_ComR[day][CutHome] += Nx
                    WNR_ComR[day][CutHome] += Nx * Wx
                   
                #Mobility comunas
                if int(Y[x][2]) in [1,2,3] and int(Y[x][8]) != 6: #excludes jobcat6
                    PM_ComR[day][CutHome] += Px
                
                #As a proxy of production
                if int(Y[x][2]) in [0,1,2]: #works in the RM
                    #People
                    N_ComWP[day][CutWork] += Nx
                    PT_ComWP[day][CutWork] += PTx
                    #Wages
                    WN_ComWP[day][CutWork] += Nx * Wx
                    WPT_ComWP[day][CutWork] += PTx * Wx
                
    LComunaR = N_ComR[0]+PT_ComR[0] #Total workers by comuna of residence
    WComunaR = WN_ComR[0]+WPT_ComR[0] #Total daily wages by comuna of residence  
    LComunaWP = N_ComWP[0]+PT_ComWP[0] #Total workers by workplaces' comuna (excluding comm=3)
    WComunaWP = WN_ComWP[0]+WPT_ComWP[0] #Total daily wages by workplaces' comuna (excluding comm=3)
    
    #Comunas (people, wages, percent people, percent wages)        
    np.savetxt("S"+str(sim)+"_CL_"+str(0)+str(1)+str(0)+".csv", N_ComR, delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_CL_"+str(0)+str(0)+str(0)+".csv", PT_ComR, delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_CL_"+str(0)+str(2)+str(0)+".csv", P_ComR, delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_CL_"+str(0)+str(3)+str(0)+".csv", T_ComR, delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_CL_"+str(0)+str(4)+str(0)+".csv", PM_ComR, delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_CL_"+str(0)+str(5)+str(0)+".csv", NR_ComR, delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_CL_"+str(0)+str(6)+str(0)+".csv", PT_ComWP, delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_CL_"+str(0)+str(7)+str(0)+".csv", N_ComWP, delimiter=",",fmt="%s")
    
    np.savetxt("S"+str(sim)+"_CL_"+str(1)+str(1)+str(0)+".csv", WN_ComR, delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_CL_"+str(1)+str(0)+str(0)+".csv", WPT_ComR, delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_CL_"+str(1)+str(2)+str(0)+".csv", WP_ComR, delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_CL_"+str(1)+str(3)+str(0)+".csv", WT_ComR, delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_CL_"+str(1)+str(5)+str(0)+".csv", WNR_ComR, delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_CL_"+str(1)+str(6)+str(0)+".csv", WPT_ComWP, delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_CL_"+str(1)+str(7)+str(0)+".csv", WN_ComWP, delimiter=",",fmt="%s")
    
    np.savetxt("S"+str(sim)+"_CL_"+str(0)+str(1)+str(1)+".csv", Percent_Calculus( N_ComR,LComunaR ), delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_CL_"+str(0)+str(0)+str(1)+".csv", Percent_Calculus( PT_ComR,LComunaR ), delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_CL_"+str(0)+str(2)+str(1)+".csv", Percent_Calculus( P_ComR,LComunaR ), delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_CL_"+str(0)+str(3)+str(1)+".csv", Percent_Calculus( T_ComR,LComunaR ), delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_CL_"+str(0)+str(4)+str(1)+".csv", Percent_Calculus( PM_ComR,LComunaR ), delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_CL_"+str(0)+str(5)+str(1)+".csv", Percent_Calculus( NR_ComR,LComunaR ), delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_CL_"+str(0)+str(6)+str(1)+".csv", Percent_Calculus( PT_ComWP,LComunaWP ), delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_CL_"+str(0)+str(7)+str(1)+".csv", Percent_Calculus( N_ComWP,LComunaWP ), delimiter=",",fmt="%s")
    
    np.savetxt("S"+str(sim)+"_CL_"+str(1)+str(1)+str(1)+".csv", Percent_Calculus( WN_ComR,WComunaR ), delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_CL_"+str(1)+str(0)+str(1)+".csv", Percent_Calculus( WPT_ComR,WComunaR ), delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_CL_"+str(1)+str(2)+str(1)+".csv", Percent_Calculus( WP_ComR,WComunaR ), delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_CL_"+str(1)+str(3)+str(1)+".csv", Percent_Calculus( WT_ComR,WComunaR ), delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_CL_"+str(1)+str(5)+str(1)+".csv", Percent_Calculus( WNR_ComR,WComunaR ), delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_CL_"+str(1)+str(6)+str(1)+".csv", Percent_Calculus( WPT_ComWP,WComunaWP ), delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_CL_"+str(1)+str(7)+str(1)+".csv", Percent_Calculus( WN_ComWP,WComunaWP ), delimiter=",",fmt="%s")
      
    del( Y, N_ComR,PT_ComR,P_ComR,T_ComR, PM_ComR, NR_ComR, WN_ComR,WPT_ComR,WP_ComR,WT_ComR,WNR_ComR,
        N_ComWP,WN_ComWP,PT_ComWP,WPT_ComWP)

        
def OD_RM_Day( sim, data2, days ):
    
    #Call after Mean_Day
    Y = np.loadtxt( data2, delimiter=",", skiprows=1 ) #file with data for each group of clones identified with i.ident
    rowsY, colsY = Y.shape
    
    for day in range( days ):
        OD_Day = np.zeros( ( 51,51 ) )  #OD matrix per day from RM to RM
        P = None
        if day <= 9:
            P = np.loadtxt("S"+str(sim)+"_P_00"+str(day)+".csv",delimiter=",")
        elif day >=10 and day <=99:
            P = np.loadtxt("S"+str(sim)+"_P_0"+str(day)+".csv",delimiter=",") 
        else:
            P = np.loadtxt("S"+str(sim)+"_P_"+str(day)+".csv",delimiter=",") 
        
        for x in range( rowsY ):
            if int(Y[x][3]) == 1 and int(Y[x][2]) in [ 2 ] : #activ and commuter inside RM
                CUThInd = int(Y[x][1])
                CUTwInd = int(Y[x][12])
                Px = ( np.sum( P[x,[3,4,5]]))*15
                OD_Day[CUThInd][CUTwInd] += Px
        if day <= 9:
            np.savetxt("S"+str(sim)+"_OD_00"+str(day)+".csv", OD_Day, delimiter=",",fmt="%s")
        elif day >=10 and day <=99:
            np.savetxt("S"+str(sim)+"_OD_0"+str(day)+".csv", OD_Day, delimiter=",",fmt="%s")
        else:
            np.savetxt("S"+str(sim)+"_OD_"+str(day)+".csv", OD_Day, delimiter=",",fmt="%s")
    del(OD_Day)

def Strenght_comunas(sim, realizations, days):
    InStrength = np.zeros((days,51))
    OutStrength = np.zeros((days,51))
    Od = None
    for day in range( days ):
        if day <= 9:
            Od=np.loadtxt("S"+str(sim)+"_OD_00"+str(day)+".csv", delimiter=",")
        elif day >=10 and day <=99:
            Od=np.loadtxt("S"+str(sim)+"_OD_0"+str(day)+".csv", delimiter=",")
        else:
            Od=np.loadtxt("S"+str(sim)+"_OD_"+str(day)+".csv", delimiter=",")
        InStrength[day,:] = np.sum(Od,axis=0)
        OutStrength[day,:] = np.sum(Od,axis=1)
    np.savetxt("S"+str(sim)+"_InStrengthComunas.csv", InStrength, delimiter=",",fmt="%s")
    np.savetxt("S"+str(sim)+"_OutStrengthComunas.csv", OutStrength, delimiter=",",fmt="%s")
    
    
def Mobility_Google(sim, realizations, days):
    #time series of all realizations of "mobility" following Google Analytics logic: i.e. time spend in workplace per day
    #Only RM
    ReasMob = np.zeros( ( days,realizations ) )
    RMbasic = np.loadtxt( "RM_basicMob.csv", delimiter=",")
    
    for rea in range( realizations ):
        mob_rea = np.loadtxt( "S"+str(sim)+"_Mob_Tot_rea_"+str(rea)+".csv", delimiter=",")
        ReasMob[:,rea] = ((mob_rea - RMbasic)/RMbasic)*100
    
    np.savetxt( "S"+str(sim)+"_MobilityRM.csv",ReasMob,delimiter=",",fmt="%s" )
    del(ReasMob,RMbasic)
    
def Atkinson( sim, realizations, days, data2 ):
    #Returns XXX files. Three for RM (Atkinson per day, rea epsilon=0.25,0.5,0.75)
    #and three for the daily average of Atkinson epsilon 0.25,0.5,0.75 by comuna
    
    RM_Atk025 = np.zeros( ( days,realizations ) )
    RM_Atk050 = np.zeros( ( days,realizations ) )
    RM_Atk075 = np.zeros( ( days,realizations ) )
    RM_WEDE025 = np.zeros( ( days,realizations ) )
    RM_WEDE050 = np.zeros( ( days,realizations ) )
    RM_WEDE075 = np.zeros( ( days,realizations ) )
    RM_WEDE025Perc = np.zeros( ( days,realizations ) ) #% Wede_t/Wede_0
    RM_WEDE050Perc = np.zeros( ( days,realizations ) )
    RM_WEDE075Perc = np.zeros( ( days,realizations ) )
    RM_Ut025 = np.zeros( ( days,realizations ) ) #% Ut_t/Ut_0
    RM_Ut050 = np.zeros( ( days,realizations ) )
    RM_Ut075 = np.zeros( ( days,realizations ) )
        
    Com_Atk025 = np.zeros( ( days,51 ) )
    Com_Atk050 = np.zeros( ( days,51 ) )
    Com_Atk075 = np.zeros( ( days,51 ) )
    Com_WEDE025 = np.zeros( ( days,51 ) )
    Com_WEDE050 = np.zeros( ( days,51 ) )
    Com_WEDE075 = np.zeros( ( days,51 ) )
    Com_WEDE025Perc = np.zeros( ( days,51 ) )
    Com_WEDE050Perc = np.zeros( ( days,51 ) )
    Com_WEDE075Perc = np.zeros( ( days,51 ) )
    Com_Ut025 = np.zeros( ( days,51 ) ) #% Ut_t/Ut_0
    Com_Ut050 = np.zeros( ( days,51 ) )
    Com_Ut075 = np.zeros( ( days,51 ) )
    
    
    Y = np.loadtxt( data2, delimiter=",", skiprows=1 ) #file with data for each group of clones identified with i.ident

    for rea in range( realizations ):
        for day in range( days ):
            m_dr = np.zeros( ( 52, 8) )
            out_dr = np.loadtxt( "S"+str(sim)+"_rea_"+str(rea)+"_day_"+str(day)+".csv", delimiter=",")
            for x in range( 19584 ):
                if int(Y[x][3]) == 1: #activ
                    Cuthx = int(Y[x][1])
                    Wx = Y[x][17]
                    WAt025x = Wx**(1-0.25)
                    WAt050x = Wx**(1-0.5)
                    WAt075x = Wx**(1-0.75)
                    NbrTL = (int(Y[x][0])*15)
                    NotWx = (out_dr[x][0]+out_dr[x][1]+out_dr[x][2])*15*int(Y[x][16]) #Number of clones who probably do not perceive income
                    NbrWx = NbrTL-NotWx #number of clones who perceive the daily income
                    m_dr[Cuthx][0] += NbrTL
                    m_dr[51][0] += NbrTL
                    m_dr[Cuthx][1] += NbrWx*Wx
                    m_dr[51][1] += NbrWx*Wx
                    m_dr[Cuthx][2] += NbrWx*WAt025x
                    m_dr[51][2] += NbrWx*WAt025x
                    m_dr[Cuthx][3] += NbrWx*WAt050x
                    m_dr[51][3] += NbrWx*WAt050x
                    m_dr[Cuthx][4] += NbrWx*WAt075x
                    m_dr[51][4] += NbrWx*WAt075x
                    m_dr[Cuthx][5] += NbrTL*WAt025x
                    m_dr[51][5] += NbrTL*WAt025x
                    m_dr[Cuthx][6] += NbrTL*WAt050x
                    m_dr[51][6] += NbrTL*WAt050x
                    m_dr[Cuthx][7] += NbrTL*WAt075x
                    m_dr[51][7] += NbrTL*WAt075x
            
            for y in range( 51 ): #compute comunas WEDE & atkinson index
                WEDE025y = ((m_dr[y][2]/m_dr[y][0])**(1.0/(1-0.25)))
                WEDE050y = ((m_dr[y][3]/m_dr[y][0])**(1.0/(1-0.50)))
                WEDE075y = ((m_dr[y][4]/m_dr[y][0])**(1.0/(1-0.75)))
                WEDETot025y = ((m_dr[y][5]/m_dr[y][0])**(1.0/(1-0.25)))
                WEDETot050y = ((m_dr[y][6]/m_dr[y][0])**(1.0/(1-0.50)))
                WEDETot075y = ((m_dr[y][7]/m_dr[y][0])**(1.0/(1-0.75)))
                WEDE025Percy = (WEDE025y/WEDETot025y)*100
                WEDE050Percy = (WEDE050y/WEDETot050y)*100
                WEDE075Percy = (WEDE075y/WEDETot075y)*100
                At025y = 1-(WEDE025y/(m_dr[y][1]/m_dr[y][0]))
                At050y = 1-(WEDE050y/(m_dr[y][1]/m_dr[y][0]))
                At075y = 1-(WEDE075y/(m_dr[y][1]/m_dr[y][0]))
                Ut025y = (m_dr[y][2]/m_dr[y][5])*100
                Ut050y = (m_dr[y][3]/m_dr[y][6])*100
                Ut075y = (m_dr[y][4]/m_dr[y][7])*100
                Com_Atk025[day][y] +=  At025y
                Com_Atk050[day][y] +=  At050y
                Com_Atk075[day][y] +=  At075y
                Com_WEDE025[day][y] += WEDE025y
                Com_WEDE050[day][y] += WEDE050y
                Com_WEDE075[day][y] += WEDE075y
                Com_WEDE025Perc[day][y] += WEDE025Percy
                Com_WEDE050Perc[day][y] += WEDE050Percy
                Com_WEDE075Perc[day][y] += WEDE075Percy    
                Com_Ut025[day][y] += Ut025y
                Com_Ut050[day][y] += Ut050y
                Com_Ut075[day][y] += Ut075y
            
            RM_WEDE025dr = ((m_dr[51][2]/m_dr[51][0])**(1.0/(1-0.25)))
            RM_WEDE050dr = ((m_dr[51][3]/m_dr[51][0])**(1.0/(1-0.50)))
            RM_WEDE075dr = ((m_dr[51][4]/m_dr[51][0])**(1.0/(1-0.75)))
            RM_WEDE025T = ((m_dr[51][5]/m_dr[51][0])**(1.0/(1-0.25)))
            RM_WEDE050T = ((m_dr[51][6]/m_dr[51][0])**(1.0/(1-0.50)))
            RM_WEDE075T = ((m_dr[51][7]/m_dr[51][0])**(1.0/(1-0.75)))
            RM_Atk025[day][rea] = 1-( RM_WEDE025dr/(m_dr[51][1]/m_dr[51][0]))        
            RM_Atk050[day][rea] = 1-( RM_WEDE050dr/(m_dr[51][1]/m_dr[51][0]))
            RM_Atk075[day][rea] = 1-( RM_WEDE075dr/(m_dr[51][1]/m_dr[51][0]))
            RM_Ut025dr = (m_dr[51][2]/m_dr[51][5])*100
            RM_Ut050dr = (m_dr[51][3]/m_dr[51][6])*100
            RM_Ut075dr = (m_dr[51][4]/m_dr[51][7])*100
            RM_WEDE025[day][rea] = RM_WEDE025dr
            RM_WEDE050[day][rea] = RM_WEDE050dr
            RM_WEDE075[day][rea] = RM_WEDE075dr
            RM_WEDE025Perc[day][rea] = ( RM_WEDE025dr/RM_WEDE025T )*100
            RM_WEDE050Perc[day][rea] = ( RM_WEDE050dr/RM_WEDE050T )*100
            RM_WEDE075Perc[day][rea] = ( RM_WEDE075dr/RM_WEDE075T )*100
            RM_Ut025[day][rea] = RM_Ut025dr 
            RM_Ut050[day][rea] = RM_Ut050dr
            RM_Ut075[day][rea] = RM_Ut075dr
            
    Com_Atk025 = Com_Atk025/float(realizations)
    Com_Atk050 = Com_Atk050/float(realizations)
    Com_Atk075 = Com_Atk075/float(realizations)
    Com_WEDE025 = Com_WEDE025/float(realizations)
    Com_WEDE050 = Com_WEDE050/float(realizations)
    Com_WEDE075 = Com_WEDE075/float(realizations)
    Com_WEDE025Perc = Com_WEDE025Perc/float(realizations)
    Com_WEDE050Perc = Com_WEDE050Perc/float(realizations)
    Com_WEDE075Perc = Com_WEDE075Perc/float(realizations)
    Com_Ut025 = Com_Ut025/float(realizations)
    Com_Ut050 = Com_Ut050/float(realizations)
    Com_Ut075 = Com_Ut075/float(realizations) 
    
    np.savetxt( "S"+str(sim)+"_CAt_025.csv",Com_Atk025,delimiter=",",fmt="%s" )
    np.savetxt( "S"+str(sim)+"_CAt_050.csv",Com_Atk050,delimiter=",",fmt="%s" )
    np.savetxt( "S"+str(sim)+"_CAt_075.csv",Com_Atk075,delimiter=",",fmt="%s" )
    np.savetxt( "S"+str(sim)+"_CWede_025.csv",Com_WEDE025,delimiter=",",fmt="%s" )
    np.savetxt( "S"+str(sim)+"_CWede_050.csv",Com_WEDE050,delimiter=",",fmt="%s" )
    np.savetxt( "S"+str(sim)+"_CWede_075.csv",Com_WEDE075,delimiter=",",fmt="%s" )
    np.savetxt( "S"+str(sim)+"_CWedeP_025.csv",Com_WEDE025Perc,delimiter=",",fmt="%s" )
    np.savetxt( "S"+str(sim)+"_CWedeP_050.csv",Com_WEDE050Perc,delimiter=",",fmt="%s" )
    np.savetxt( "S"+str(sim)+"_CWedeP_075.csv",Com_WEDE075Perc,delimiter=",",fmt="%s" )
    np.savetxt( "S"+str(sim)+"_CUt_025.csv",Com_Ut025,delimiter=",",fmt="%s" )
    np.savetxt( "S"+str(sim)+"_CUt_050.csv",Com_Ut050,delimiter=",",fmt="%s" )
    np.savetxt( "S"+str(sim)+"_CUt_075.csv",Com_Ut075,delimiter=",",fmt="%s" )
    np.savetxt( "S"+str(sim)+"_MAt_025.csv",RM_Atk025,delimiter=",",fmt="%s" )
    np.savetxt( "S"+str(sim)+"_MAt_050.csv",RM_Atk050,delimiter=",",fmt="%s" )
    np.savetxt( "S"+str(sim)+"_MAt_075.csv",RM_Atk075,delimiter=",",fmt="%s" )
    np.savetxt( "S"+str(sim)+"_MWede_025.csv",RM_WEDE025,delimiter=",",fmt="%s" )
    np.savetxt( "S"+str(sim)+"_MWede_050.csv",RM_WEDE050,delimiter=",",fmt="%s" )
    np.savetxt( "S"+str(sim)+"_MWede_075.csv",RM_WEDE075,delimiter=",",fmt="%s" )
    np.savetxt( "S"+str(sim)+"_MWedeP_025.csv",RM_WEDE025Perc,delimiter=",",fmt="%s" )
    np.savetxt( "S"+str(sim)+"_MWedeP_050.csv",RM_WEDE050Perc,delimiter=",",fmt="%s" )
    np.savetxt( "S"+str(sim)+"_MWedeP_075.csv",RM_WEDE075Perc,delimiter=",",fmt="%s" )
    np.savetxt( "S"+str(sim)+"_MUt_025.csv",RM_Ut025,delimiter=",",fmt="%s" )
    np.savetxt( "S"+str(sim)+"_MUt_050.csv",RM_Ut050,delimiter=",",fmt="%s" )
    np.savetxt( "S"+str(sim)+"_MUt_075.csv",RM_Ut075,delimiter=",",fmt="%s" )
    
    del(Y,RM_Atk025,RM_Atk050,RM_Atk075,RM_WEDE025,RM_WEDE050,RM_WEDE075,RM_WEDE025Perc,
        RM_WEDE050Perc,RM_WEDE075Perc,Com_Atk025,Com_Atk050,Com_Atk075,Com_WEDE025,
        Com_WEDE050,Com_WEDE075,Com_WEDE025Perc,Com_WEDE050Perc,Com_WEDE075Perc)

def Production( sim, realizations, days, data2 ):
    #Workforce, wages as a proxy of production. Includes only workers working in RM (excludes comm=3)
    
    Y = np.loadtxt( data2, delimiter=",", skiprows=1 ) #file with data for each group of clones identified with i.ident
    
    #People Comunas WORKPLACE
    N_ComWP = np.zeros(( days,  51 ) )         #Employed, not working
    PT_ComWP = np.zeros( ( days,  51 ) )       #Employed, working
    
    #Wages Comunas WORKPLACE
    WN_ComWP = np.zeros(( days,  51 ) )         #Employed, not working
    WPT_ComWP = np.zeros( ( days,  51 ) )       #Employed, working
    
    #People RM WORKPLACE
    RM_NWP = np.zeros( ( days, realizations ) ) #Employed, not working (excludes comm=3)
    RM_PTWP = np.zeros( ( days, realizations ) ) #trabajan, con lugar de trabajo en la RM (excluye comm=3)
    
    #Wages RM WORKPLACE
    WRM_NWP = np.zeros( ( days, realizations ) ) #Employed, not working (excludes comm=3)
    WRM_PTWP = np.zeros( ( days, realizations ) ) #trabajan, con lugar de trabajo en la RM (excluye comm=3)
    
    dayweek = 7
    
    for day in range( days ):
        for rea in range( realizations ):
            m_rea = np.loadtxt( "S"+str(sim)+"_rea_"+str(rea)+"_day_"+str(day)+".csv", delimiter=",")
            for w in range( 19584 ):
                activox = int(Y[w][19])
                workRMx = int(Y[w][20])
                Wagex = Y[w][17]
                jobcatx = int(Y[w][8])
                ramax = int(Y[w][9])
                Cutwx = int(Y[w][12])
                if activox == 1 and workRMx == 1: #if agent is a worker laboring in RM
                    if jobcatx != 6: #if agent is not an "indoor service worker"
                        if dayweek <= 5 and day not in [40,61,81,120,137]:
                            Nx = (m_rea[w][0]+ m_rea[w][1]+m_rea[w][2])*15
                            Tx = (m_rea[w][3]+ m_rea[w][4]+m_rea[w][5]+m_rea[w][6]+ m_rea[w][7]+m_rea[w][8])*15
                            WNx = Wagex * Nx
                            WTx = Wagex * Tx
                            RM_NWP[day][rea] += Nx
                            RM_PTWP[day][rea] += Tx
                            WRM_NWP[day][rea] += WNx
                            WRM_PTWP[day][rea] += WTx
                            N_ComWP[day][Cutwx] += Nx
                            PT_ComWP[day][Cutwx] += Tx
                            WN_ComWP[day][Cutwx] += WNx
                            WPT_ComWP[day][Cutwx] += WTx
                        elif dayweek == 6 or ( day in [40,120,137]):
                            if ramax in [11,15,16,21]: #do not work
                                Nx = (m_rea[w][0]+ m_rea[w][1]+m_rea[w][2]+m_rea[w][3]+ m_rea[w][4]+m_rea[w][5]+m_rea[w][6]+ m_rea[w][7]+m_rea[w][8])*15
                                WNx = Wagex * Nx
                                RM_NWP[day][rea] += Nx
                                WRM_NWP[day][rea] += WNx
                                N_ComWP[day][Cutwx] += Nx
                                WN_ComWP[day][Cutwx] += WNx
                            else:
                                Nx = (m_rea[w][0]+ m_rea[w][1]+m_rea[w][2])*15
                                Tx = (m_rea[w][3]+ m_rea[w][4]+m_rea[w][5]+m_rea[w][6]+ m_rea[w][7]+m_rea[w][8])*15
                                WNx = Wagex * Nx
                                WTx = Wagex * Tx
                                RM_NWP[day][rea] += Nx
                                RM_PTWP[day][rea] += Tx
                                WRM_NWP[day][rea] += WNx
                                WRM_PTWP[day][rea] += WTx
                                N_ComWP[day][Cutwx] += Nx
                                PT_ComWP[day][Cutwx] += Tx
                                WN_ComWP[day][Cutwx] += WNx
                                WPT_ComWP[day][Cutwx] += WTx
                        elif dayweek == 7 or ( day in [61,81] ):
                            if ramax in [3,6,10,11,12,13,14,15,16,19,20,21]: #do not work
                                Nx = (m_rea[w][0]+ m_rea[w][1]+m_rea[w][2]+m_rea[w][3]+ m_rea[w][4]+m_rea[w][5]+m_rea[w][6]+ m_rea[w][7]+m_rea[w][8])*15
                                WNx = Wagex * Nx
                                RM_NWP[day][rea] += Nx
                                WRM_NWP[day][rea] += WNx
                                N_ComWP[day][Cutwx] += Nx
                                WN_ComWP[day][Cutwx] += WNx
                            else:
                                Nx = (m_rea[w][0]+ m_rea[w][1]+m_rea[w][2])*15
                                Tx = (m_rea[w][3]+ m_rea[w][4]+m_rea[w][5]+m_rea[w][6]+ m_rea[w][7]+m_rea[w][8])*15
                                WNx = Wagex * Nx
                                WTx = Wagex * Tx
                                RM_NWP[day][rea] += Nx
                                RM_PTWP[day][rea] += Tx
                                WRM_NWP[day][rea] += WNx
                                WRM_PTWP[day][rea] += WTx
                                N_ComWP[day][Cutwx] += Nx
                                PT_ComWP[day][Cutwx] += Tx
                                WN_ComWP[day][Cutwx] += WNx
                                WPT_ComWP[day][Cutwx] += WTx
                        else:
                            pass
                    else: #jobcat=6
                        if dayweek == 7: #do not work
                            Nx = (m_rea[w][0]+ m_rea[w][1]+m_rea[w][2]+m_rea[w][3]+ m_rea[w][4]+m_rea[w][5]+m_rea[w][6]+ m_rea[w][7]+m_rea[w][8])*15
                            WNx = Wagex * Nx
                            RM_NWP[day][rea] += Nx
                            WRM_NWP[day][rea] += WNx
                            N_ComWP[day][Cutwx] += Nx
                            WN_ComWP[day][Cutwx] += WNx
                        else:
                            Nx = (m_rea[w][0]+ m_rea[w][1]+m_rea[w][2])*15
                            Tx = (m_rea[w][3]+ m_rea[w][4]+m_rea[w][5]+m_rea[w][6]+ m_rea[w][7]+m_rea[w][8])*15
                            WNx = Wagex * Nx
                            WTx = Wagex * Tx
                            RM_NWP[day][rea] += Nx
                            RM_PTWP[day][rea] += Tx
                            WRM_NWP[day][rea] += WNx
                            WRM_PTWP[day][rea] += WTx
                            N_ComWP[day][Cutwx] += Nx
                            PT_ComWP[day][Cutwx] += Tx
                            WN_ComWP[day][Cutwx] += WNx
                            WPT_ComWP[day][Cutwx] += WTx
                                
        if dayweek == 7:
            dayweek = 1
        else:
            dayweek += 1
    #Day averages (over 100 realizations) for comunas
    N_ComWP = N_ComWP/float(realizations)
    PT_ComWP = PT_ComWP/float(realizations)
    WN_ComWP = WN_ComWP/float(realizations)
    WPT_ComWP = WPT_ComWP/float(realizations)
    #Totals by comuna and RM
    TotLCom =  N_ComWP + PT_ComWP
    TotWCom = WN_ComWP + WPT_ComWP
    TotLRM = RM_NWP + RM_PTWP
    TotWRM = WRM_NWP + WRM_PTWP 
    
    np.savetxt( "S"+str(sim)+"_CL_NoProdNumberPeople.csv", N_ComWP,delimiter=",",fmt="%s" )
    np.savetxt( "S"+str(sim)+"_CL_ProdNumberPeople.csv", PT_ComWP,delimiter=",",fmt="%s" )
    np.savetxt( "S"+str(sim)+"_ML_NoProdNumberPeople.csv", RM_NWP,delimiter=",",fmt="%s" )
    np.savetxt( "S"+str(sim)+"_ML_ProdNumberPeople.csv", RM_PTWP,delimiter=",",fmt="%s" )
    
    np.savetxt( "S"+str(sim)+"_CL_NoProdPesosWages.csv", WN_ComWP,delimiter=",",fmt="%s" )
    np.savetxt( "S"+str(sim)+"_CL_ProdPesosWages.csv", WPT_ComWP,delimiter=",",fmt="%s" )
    np.savetxt( "S"+str(sim)+"_ML_NoProdPesosWages.csv", WRM_NWP,delimiter=",",fmt="%s" )
    np.savetxt( "S"+str(sim)+"_ML_ProdPesosWages.csv", WRM_PTWP,delimiter=",",fmt="%s" )
    
    
    np.savetxt( "S"+str(sim)+"_CL_NoProdFracPeople.csv", (N_ComWP/TotLCom)*100,delimiter=",",fmt="%s" )
    np.savetxt( "S"+str(sim)+"_CL_ProdFracPeople.csv", (PT_ComWP/TotLCom)*100,delimiter=",",fmt="%s" )
    np.savetxt( "S"+str(sim)+"_ML_NoProdFracPeople.csv", (RM_NWP/TotLRM)*100,delimiter=",",fmt="%s" )
    np.savetxt( "S"+str(sim)+"_ML_ProdFracPeople.csv", (RM_PTWP/TotLRM)*100,delimiter=",",fmt="%s" )
    
    np.savetxt( "S"+str(sim)+"_CL_NoProdFracWages.csv", (WN_ComWP/TotWCom)*100,delimiter=",",fmt="%s" )
    np.savetxt( "S"+str(sim)+"_CL_ProdFracWages.csv", (WPT_ComWP/TotWCom)*100,delimiter=",",fmt="%s" )
    np.savetxt( "S"+str(sim)+"_ML_NoProdFracWages.csv", (WRM_NWP/TotWRM)*100,delimiter=",",fmt="%s" )
    np.savetxt( "S"+str(sim)+"_ML_ProdFracWages.csv", (WRM_PTWP/TotWRM)*100,delimiter=",",fmt="%s" )
    
    np.savetxt( "S"+str(sim)+"_CL_CheckComunaPeople.csv", TotLCom,delimiter=",",fmt="%s" )
    np.savetxt( "S"+str(sim)+"_CL_CheckComunaWages.csv", TotWCom,delimiter=",",fmt="%s" )
    np.savetxt( "S"+str(sim)+"_ML_CheckRMPeople.csv", TotLRM,delimiter=",",fmt="%s" )
    np.savetxt( "S"+str(sim)+"_ML_CheckRMWages.csv", TotWRM,delimiter=",",fmt="%s" )
    
    del(Y, N_ComWP, PT_ComWP, WN_ComWP, WPT_ComWP, RM_NWP, RM_PTWP, WRM_NWP, WRM_PTWP)
    
#print ("inicio", time.ctime())
#t1=time.time()
#Mean_Day( sim="001", realizations=100, days=154 )
#print ("finMean", time.time()-t1)

#t11=time.time()
#Detected_Series( sim="001", realizations=100, days=154, RMdata="RealDRM.csv" )
#print ("finDetectedSeries", time.time()-t11)

#t2=time.time()
#LabourSeriesComuna( sim="001", data2="Data2_MP.csv", days=154 )
#print ("finLabourSeries", time.time()-t2)

#t3=time.time()
#Comunas_HealthSeries( sim="001", data2="Data2_MP.csv", days=154 )
#print ("finHealthSeries", time.time()-t3)

#t4=time.time()
#OD_RM_Day( sim="001", data2="Data2_MP.csv", days=154 )
#print ("finOD", time.time()-t4)

#t5=time.time()
#SeriesRM ( sim="001", realizations=100, days=154, data2="Data2_MP.csv" )
#print ("finSeriesRM", time.time()-t5)

#t6=time.time()
#Atkinson( sim="001", realizations=100, days=154, data2="Data2_MP.csv" )
#print("finAtkinson",time.time()-t6)

#incrementComuna( sim="001", realdatcom="RealDCom.csv" )

#Mobility_Google(sim="001", realizations=100, days=154)

#Strenght_comunas(sim="001", realizations=100, days=154)

#Production( sim="001", realizations=100, days=154, data2="Data2_MP.csv" )

#Detected_Series2( sim="001", realizations=100, days=154, RMdata="RealDRM.csv" )

#print("FIN PROCESO",time.time()-t1)
     
