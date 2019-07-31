# -*- coding: utf-8 -*-
"""
Created on Tue Jul 17 11:06:44 2018

@author: gjoterorodrigu
"""

## START: Initializing Packages 
import sys                         # for sys.exit(-1)
import time                        # for timing the code
import matplotlib.pyplot as plt    # for plotting routines
import os                          # Check Existance of File
import pdb                         # Debugging Module
import numpy as np
import math as calc
from src.ISEN_RELATIONS import *
from src.SHOCKWAVES import *
from src.FANNO_RAYLEIGH import *
#from src.EoS_RG import calcRG_gamma
from copy import deepcopy
import CoolProp
import CoolProp.CoolProp as CP
fluid = 'Toluene' #'R245FA'
## END: Initializing Packages 

      
#Oblique shock example
if 0>1:
    theta1=5*calc.pi/180
    theta2=7*calc.pi/180
    theta3=theta1+theta2
    
    gam =1.4
    M1=8.5
    
    p1= 12000
    T1= 220
    R=287
    rho1=p1/R/T1
    p01 = p1*(1+(gam-1)*0.5*M1*M1)**(gam/(gam-1))
    
    #solving shock 1
    beta1, M2, pR, TR, p0R = ObliqueShock(gam,M1,theta1)
    p2=p1*pR
    T2=T1*TR
    p02=p01*p0R
    Mn1= M1*np.sin(beta1)
    Mn2= M2*np.sin(beta1-theta1)
    rho2=p2/R/T2
    
    #solving shock 2
    beta2, M3, pR, TR, p0R  = ObliqueShock(gam,M2,theta2)
    p3=p2*pR
    T3=T2*TR
    p03=p02*p0R
    Mn2= M2*np.sin(beta2)
    Mn3= M3*np.sin(beta2-theta2)
    rho2=p2/R/T2
    
    #solving shock 3
    M3= 6
    p3=23
    T3=600
    beta3, M4, pR, TR, p0R   = ObliqueShock(gam,M3,theta3)
    p4=p3*pR
    T4=T3*TR
    p04=p03*p0R
    Mn3= M3*np.sin(beta3)
    Mn4= M3*np.sin(beta3-theta3)
    rho3=p3/R/T3




if 0>1:
    
    #Fanno flow example
    f=0.002
    D=0.4
    L=0.5
    gam=1.4
    Min=3.5
    
    fanno_in= FannoFunction(gam,Min)
    fanno_out= fanno_in-4*f*L/D
       
    
    
    Minn = MachFromFannoFunction(gam,fanno_in)
    Mout = MachFromFannoFunction(gam,fanno_out)
    
    Lstar=ChockedLengthFanno(gam,Min,f,D)
    #fanno_star= fanno_in-4*f*Lstar/D   
    #Mstar = MachFromFannoFunction(gam,fanno_star)
    
    P0P0star_in  = CriticalTotalPressureRatioFanno(gam,Min)   
    P0P0star_out = CriticalTotalPressureRatioFanno(gam,Mout)    
    
    print("Error fanno out and Min", (FannoFunction(gam,Mout)-fanno_out)*100/fanno_out, (Min-Minn)*100/Min)
    print("Results: M_out, Missing legnth for Chocked and p0out/p0in", Mout, Lstar-L, P0P0star_out/P0P0star_in)
    
    #Fanno flow excercise
    D=0.4
    L=0.5
    gam=1.4
    Min=2.5
    
    fanno_in   = FannoFunction(gam,Min)
    fanno_star = FannoFunction(gam,1)
    f = D/L/4*(fanno_in-fanno_star)
    
    print("f ", f)
    
    f_new = f/2
    
    fanno_out= fanno_in-4*f_new*L/D
    Mout= MachFromFannoFunction(gam,fanno_out)
    print("Mout ", Mout)
    
    P0P0star_in  = CriticalTotalPressureRatioFanno(gam,Min)   
    P0P0star_out = CriticalTotalPressureRatioFanno(gam,Mout)    
    print("P0_exit/P0_in ", P0P0star_out/P0P0star_in)
    

if 0>1:
    
    #Rayleigh flow example
    Min=3.5
    Tin=600
    dh0=500*(10**3)   #[J/kg]
    gam=1.4
    cp=1.005*(10**3)
    
    # Calculating the stagnation temperature
    T0in=Tin*TotalTemperatureToStaticTemperatureRatio(gam,Min,gam)
    print("T0in ", T0in)
    #energy balance
    T0out=dh0/cp + T0in
    print("T0out ", T0out)
    #calculating the Rayleigh function value at the inlet and outlet
    ray_in  = RayleighFunction(gam,Min)
    ray_out = T0out/T0in*ray_in
    
    
        
    print("T0in/T0star ", CriticalTotalTemperatureRatioRayleigh(gam,Min)) 
    print("Tin/Tstar ", CriticalTemperatureRatioRayleigh(gam,Min)) 
    print("P0in/P0star ", CriticalTotalPressureRatioRayleigh(gam,Min)) 
    
    Mout= MachFromRayleighFunction(gam,ray_out)
    print("Result 1, Mout ", Mout)
    
    Tout=T0out/TotalTemperatureToStaticTemperatureRatio(gam,Mout,gam)
    print("Result 2, Tout ", Tout)
    
    p02_p01=CriticalTotalPressureRatioRayleigh(gam,Mout)/CriticalTotalPressureRatioRayleigh(gam,Min)
    print("Result 3, p02 to p01 ", p02_p01)
    
    #Rayleigh flow excercise
    D=0.4         #diameter of the pipe
    Min=3.2
    T0in=700
    dh0=0.2*(10**6)    #[J/kg]
    gam=1.4
    cp=1.005*(10**3)   #[J/kgK]
    
    #energy balance
    T0out=dh0/cp + T0in
    print("T0out ", T0out)
    
    ray_in  = RayleighFunction(gam,Min)
    ray_out = T0out/T0in*ray_in
    
    Mout= MachFromRayleighFunction(gam,ray_out)
    print("Result excercise 1, Mout ", Mout)
    p02_p01=CriticalTotalPressureRatioRayleigh(gam,Mout)/CriticalTotalPressureRatioRayleigh(gam,Min)
    print("Result excercise 2, p02 to p01 ", p02_p01)
    
    ray_outChocked  = RayleighFunction(gam,1.0000000)
    T0outChocked =ray_outChocked/ray_in*T0in
    dh0Chocked=cp*(T0outChocked-T0in)
    print("Result excercise 3, chocked q in [MJ/kg] ", dh0Chocked/(10**6))
    
    p02_p01Chocked=CriticalTotalPressureRatioRayleigh(gam,1.0)/CriticalTotalPressureRatioRayleigh(gam,Min)
    print("Result excercise 4, T0out-T0in/T0in Chocked ", (T0outChocked-T0in)/T0in )    
    print("Result excercise 5, p02 to p01 Chocked ", p02_p01Chocked)
    
if 0>1:
    
    #Combustion problem
    #C3H8
    x=3
    y=8
    LHV=46.35*(10**6) 
    MWair=28.85
    MWfuel=44
    eta=0.9
    phi=0.8
    
    mair=1
  
    
    a=x+y/4
    print("a  ", a)
    n_st=a*4.76
    print("stoichiometric molar air to fuel ratio ", n_st)
    m_st=n_st*MWair/MWfuel
    print("stoichiometric mass air to fuel ratio ", m_st)
    m_real=m_st/phi
    print("actual mass air to fuel ratio ", m_real)
    mfuel=mair/m_real
    Qcomb=eta*mfuel*LHV
    print("Amount of heat released in [MW]", Qcomb/(10**6))
    print("Amount of heat released in [MJ/kg]", Qcomb/(10**6)/(mair+mfuel))
    
    #Hydrogen now
    LHVh2=119.96*(10**6) 
    mfuelh2= Qcomb/eta/LHVh2
    print("Amount of h2 ", mfuelh2)
    x=0
    y=2
    a=x+y/4
    MWh2=2
    n_st=a*4.76
    m_st=n_st*MWair/MWh2
    print("stoichiometric mass air to fuel ratio for hydrogen ", m_st)
    m_real=mair/mfuelh2
    phih2=m_st/m_real
    print("equivalence ratio for hydrogen-air ", phih2)
    
if 0>1:
    #Assignment
    gam=1.4 
    M1=4.5
    T1=1200
    D=0.4
    f=2*(10**(-3))
    Lmin=0.5
    Lmax=1.0
    
    
    #PART A: Fanno part
    b=0.75 #0.25
    L=Lmin+b*(Lmax-Lmin)
    fanno_1= FannoFunction(gam,M1)
    fanno_2= fanno_1-4*f*L/D
    M1n = MachFromFannoFunction(gam,fanno_1)
    M2  = MachFromFannoFunction(gam,fanno_2)
    print("PART A: Value of Minterim: ",M2,b)
    P0P0star_1  = CriticalTotalPressureRatioFanno(gam,M1)   
    P0P0star_2 = CriticalTotalPressureRatioFanno(gam,M2)    
    print("PART A: P0_interim/P0_inlet ", P0P0star_2/P0P0star_1,b)
    
    #PART B: Combustion part
    LHV=119.96*(10**6) #hydrogen
    phi=0.7
    MWair=28.85
    MWh2=2
    x=0
    y=2
    a=x+y/4
    
    n_st=a*4.76
    m_st=n_st*MWair/MWh2
    n_real=n_st/phi
    m_real=m_st/phi
    print("PART B: actual mass air to fuel ratio ", m_real)
    print("PART B: actual molar air to fuel ratio ", n_real)
    
    x=0*0.5+1*0.5
    y=2*0.5+4*0.5
    a=x+y/4
    print("PART B: number of moles of oxygen ", a)
    
    #PART C: Rayleigh flow
    eta_min=0.73
    eta_max=0.98
    
    eta_comb=eta_min+(L-Lmin)*(eta_max-eta_min)/(Lmax-Lmin)
    cp=1.005*(10**3)
    
    print("PART C: combustion efficiency ", eta_comb, b)
    dh0= eta_comb*LHV/(m_real+1)
    T01=T1*TotalTemperatureToStaticTemperatureRatio(gam,M1,gam)
    print("PART C: heat added per unit mass in MJ/kg", dh0/(10**6), b)
    T03=dh0/cp+T01
    print("PART C: temperature ratio (T03-T01)/T01", (T03-T01)/T01, b)
    ray_2  = RayleighFunction(gam,M2)
    ray_3 = T03/T01*ray_2  #T02=T01
    
    M3= MachFromRayleighFunction(gam,ray_3)
    ray_3_  = RayleighFunction(gam,M3)
    
    P020Pstar=CriticalTotalPressureRatioRayleigh(gam,M2)
    P030Pstar=CriticalTotalPressureRatioRayleigh(gam,M3)

    # They key of this solution is that p* from fanno flow is not the same as
    # p* from Rayleigh
    print("PART C: P0_exit/P0_inlet ", P0P0star_2/P0P0star_1*P030Pstar/P020Pstar  ,b)
    
    
    
if 0>1:
    gam=1.4
    R=287
    
    U=5000
    theta=8*calc.pi/180
    Tinf= 250
    sos=(gam*R*Tinf)**0.5    
    
    Up=5000*np.tan(theta)
    
    Us=HypersonicEquivalencePrinciple(gam,Up,Tinf,R)
    sys.exit()
    a_in=330
    T_in=a_in**2/R/gam
    up=250
    Us=HypersonicEquivalencePrinciple(gam,Up,Tinf,R)
    print("Shockwave velocity",Us)
    
    #given by the question
    u_piston=250    
    u_sw=500    
    
    theta=6*calc.pi/180
    Mach=7.2
    Velocity=Mach*a_in
    Upiston=Velocity*np.tan(theta) 
    print("Upiston velocity",Upiston, u_piston)
    
    beta=calc.atan(u_sw/Velocity)
    print("Beta angle: ",beta*180/calc.pi)
    
    
    u_piston=100    
    u_sw=392   
    beta=9.5*calc.pi/180
    u_shock=Velocity*np.tan(beta)
    print("Shcokwave velocity",u_shock, u_sw)  
    theta=calc.atan(u_piston/Velocity)
    print("Theta angle",theta*180/calc.pi)
    
        #given by the question
    u_piston=250    
    u_sw=500    
    
    theta=2*calc.pi/180
    Mach=21.7
    Velocity=Mach*a_in
    Upiston=Velocity*np.tan(theta) 
    print("Upiston velocity",Upiston, u_piston)
    
    beta=calc.atan(u_sw/Velocity)
    print("Beta angle: ",beta*180/calc.pi)    
    
    u_piston=100    
    u_sw=392   
    beta=3.2*calc.pi/180
    u_shock=Velocity*np.tan(beta)
    print("Shcokwave velocity",u_shock, u_sw)  
    theta=calc.atan(u_piston/Velocity)
    print("Theta angle",theta*180/calc.pi)
    
    
    
    mass=5500
    cd=1.3
    ventry=11.2*1000
    area=12
    theta=7*calc.pi/180
    e=calc.exp(1)
    y0=6620
    rho0=1.22
    
    gmax=ventry*ventry*np.sin(theta)/2/e/y0
    y= y0*calc.log(cd*rho0*area*y0/mass/np.sin(theta))
    V= calc.exp(-0.5)*ventry
    
    print("max deacceleration, y and vmin ", gmax, y ,V)

##############################################################
##############################################################
### ### ### ###  ###  ###   Project   ### ### ### ### ### ### 
##############################################################
##############################################################
if 2>1:
    #Location 2: after first shock
    # Need to calculate the Drag: Use the F2=p2*Asurface2 (hypotenusas of the triangle)
    #                         D2=F2_x=F2_x*sin or cos of theta1
        #second compression shock
    #Location 3: after second shock
        #third compression shock
        # Need to calculate the Drag: Use the F3=p3*Asurface3 (hypotenusas of the triangle)
        #                         D3=F2_x=F3_x*sin or cos of theta2
        # Inlet Drag = D2+D3
    #Location 4: after third shock/inlet of the combustor
        #fuel is added, therefore heat is introduced
    #Location 5: end of the combustor/inlet of the isentropic nozzle
        #expansion in the isentropic nozzle
    #Location 6: end of the nozzle/exit of the engine
        #Need to Calculate the Nozzle Thrust: Control volume over the nozzle
        # All needs to be balance by the nozzle thrust: Thrust_nozzle=p5*A5*(1+gam*M5)-p6*A6*(1+gam*M6)
    
    
    # net Thrust= Thrust_nozzle - Inlet Drag
    # Specific impulse Isp= net thrust/(mfuel*g)
    
    #MY CODE!!   30113647100
    gam=1.4
    R = 287    #[J/kg/K]
    cp = 1003  #[J/kg/K]    
    w = 1
    hc = 1
    A1 = w*hc
    

    M1= 11.36
    p1= 1.1969*1000
    T1= 226.5
    
    ht1= 0.3279
    theta1 = 4 *calc.pi/180
        
    ht2= 0.6006
    theta2 = 7 *calc.pi/180
    
    theta3 = theta1+theta2
    ht3 = 0.0715    
    
    ht4=ht3
    ht6=1.5
    phi = 1.0
    MWair=28.85
    MWh2=2
    LHVh2= 120*(10**6)
    etaComb=100
    g=9.81
    
    
    
    #Location 1: Inlet of the engine
        #first compression shock
    rho1=p1/R/T1
    p01 = p1*TotalPressureToStaticPressureRatio(gam,M1)
    T01 = T1*TotalTemperatureToStaticTemperatureRatio(gam,M1,gam)
    v1 = M1*((gam*R*T1)**0.5)
    print("state 1: v, rho, p0, T0" , v1, rho1, p01/1000, T01)

    beta1, M2, pR, TR, p0R = ObliqueShock(gam,M1,theta1)
    
    p2=p1*pR
    T2=T1*TR
    rho2=p2/R/T2
    p02=p01*p0R
    T02 = T2*TotalTemperatureToStaticTemperatureRatio(gam,M2,gam)
    v2 = M2*((gam*R*T2)**0.5)
    print("state 2: beta, p, T, M" , beta1*180/calc.pi, p2/1000, T2, M2)
    
    beta2, M3, pR, TR, p0R  = ObliqueShock(gam,M2,theta2)
    p3=p2*pR
    T3=T2*TR
    T03 = T3*TotalTemperatureToStaticTemperatureRatio(gam,M3,gam)
    p03=p02*p0R
    Mn2= M2*np.sin(beta2)
    Mn3= M3*np.sin(beta2-theta2)
    
    print("state 3: beta, p, T, M" , beta2*180/calc.pi, p3/1000, T3, M3)
    
    beta3, M4, pR, TR, p0R   = ObliqueShock(gam,M3,theta3)
    p4=p3*pR
    T4=T3*TR
    T04 = T4*TotalTemperatureToStaticTemperatureRatio(gam,M4,gam)
    p04=p03*p0R
    Mn3= M3*np.sin(beta3)
    Mn4= M3*np.sin(beta3-theta3)
    rho3=p3/R/T3
    
    print("state 4: beta, p, T, M" , beta3*180/calc.pi, p4/1000, T4, M4)
    
    #calculating drag
    D1= p2*w*ht1
    D2= p3*w*ht2
    D =  D1+D2
    print("Drag force in the inlet strucutre" , D)
    
    rho4=p4/R/T4
    v4 = M4*((gam*R*T4)**0.5)
    
    x=0
    y=2
    a=x+y/4
    n_st=a*4.76
    m_st=n_st*MWair/MWh2
    
    mair=  rho4*v4*w*ht3
    m_real=m_st/phi
    mfuel=mair/m_real
    print("Mass flows of air and fuel" , mair,mfuel)
    
    dh0_combustor= mfuel*LHVh2/mair
    T05=dh0_combustor/cp + T04
    print("Heat released per unit mas of air" , dh0_combustor)    
    
    ray_in  = RayleighFunction(gam,M4)
    ray_out = T05/T04*ray_in   
    M5= MachFromRayleighFunction(gam,ray_out)
    T5 = T05/TotalTemperatureToStaticTemperatureRatio(gam,M5,gam)
    p05=p04*CriticalTotalPressureRatioRayleigh(gam,M5)/CriticalTotalPressureRatioRayleigh(gam,M4)    
    p5 = p05/TotalPressureToStaticPressureRatio(gam,M5)
    v5 = M5*((gam*R*T5)**0.5)
    print("State 5, T0, M, T, p" , T05, M5, T5, p5/1000) 
    
    Astar = ht3*w/CriticalAreaRatio(gam,M5)
    A6byAstar= ht6*w/Astar
    
    #Using a newton raphson
    tol=1e-06
    error=100
    n=0
    dM=0.001
    M0=M5
    relaxfact=0.85  
    while error>tol and n<1000:
        derivative = ((A6byAstar-CriticalAreaRatio(gam,M0+dM))-(A6byAstar-CriticalAreaRatio(gam,M0)))/dM           
        M6 = M0 - relaxfact*(A6byAstar-CriticalAreaRatio(gam,M0))/derivative
        error=abs(M6-M0)
        M0=M6              
        n=n+1    
        #print("errror ", error, CriticalAreaRatio(gam,M0), A6byAstar, M0)
    
    T06=T05
    T6 = T06/TotalTemperatureToStaticTemperatureRatio(gam,M6,gam)
    pstar = p05*CriticalPressureRatio(gam,M5)
    Tstar = T05*CriticalTemperatureRatio(gam,M5,gam)    
    p6 = pstar/(CriticalPressureRatio(gam,M5)*TotalPressureToStaticPressureRatio(gam,M6))
    T6s = Tstar/(CriticalTemperatureRatio(gam,M5,gam)*TotalTemperatureToStaticTemperatureRatio(gam,M6,gam))
    v6 = M6*((gam*R*T6)**0.5)
    print("State 6, T0, M, T, p" , T06, M6, T6, p6/1000, T6s) 
    
    Thrust_nozzle=mair*v6+p6*ht6*w-mair*v5-p5*ht3*w
    print("Thrust produce by the nozzle" , Thrust_nozzle) 
    
    Net_Thrust=Thrust_nozzle-D
    Isp = Net_Thrust/(mfuel*g)
    print("Net Thrust of the engine" , Net_Thrust)    
    print("Specific impulse " ,Isp)   
    print("flight Mach" , M1)   
sys.exit()
