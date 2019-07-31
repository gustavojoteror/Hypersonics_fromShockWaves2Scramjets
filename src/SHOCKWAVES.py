# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 15:16:55 2018

@author: gjoterorodrigu
"""
import sys
import os        
import math as calc
import numpy as np
###########################################
##For normal shockwaves ideal gases
def MachDownstreamShock(gam,Min):
    return ((1+(gam-1)*0.5*Min*Min)/(gam*Min*Min-(gam-1)*0.5))**(0.5)
    
def PressureRatioShock(gam,Min):  #p2/p1
    return (2*gam*Min*Min-gam+1)/(gam+1)
    
def TotalPressureRatioShock(gam,Min):  #p02/p01
    return ((gam+1)*Min*Min/((gam-1)*Min*Min+2))**(gam/(gam-1)) * ((gam+1)/(Min*Min*gam*2-gam+1))**(1/(gam-1))
    
def TemperatureRatioShock(gam,Min):  #T2/T1 
    return ((2*gam*Min*Min-gam+1)*((gam-1)*Min*Min+2))/((gam+1)**2*Min*Min)
    
def TotalTemperatureRatioShock(gam,Min):  #T02/T01
    return 1

def DensityRatioShock(gam, Min): #rho2/rho1
    return (gam+1)*(Min)**2/(2+(gam-1)*(Min)**2)
    
def DensityRatioObliqueShock(gam,Min,beta,theta,switch):
    if switch==1:
        return np.tan(beta*(calc.pi/180))/np.tan((beta-theta)*(calc.pi/180)) 
    else:
        return DensityRatioShock(gam, Min*np.sin(beta*(calc.pi/180))) #(gam+1)*(mach*np.sin(beta*(calc.pi/180)))**2/(2+(gam-1)*(mach*np.sin(beta*(calc.pi/180)))**2)
     
def ObliqueShockFunction(gam,Min,theta,beta):
    return np.tan(beta)/np.tan(beta-theta)-(gam+1)*(Min*np.sin(beta))**2/(2+(gam-1)*(Min*np.sin(beta))**2)

def WaveAngleObliqueShock(gam,Min,theta,waveType):
    tol=1e-06
    error=100
    n=0
    dbeta = 0.01*calc.pi/180
    if waveType==1:
        beta0=90*calc.pi/180 # strong shock wave
    else:
        beta0=(theta+1*calc.pi/180) #(theta+1*calc.pi/180)  # weak shock wave

    #Using a newton raphson
    while error>tol and n<1000:
        derivative = (ObliqueShockFunction(gam,Min,theta,beta0+dbeta)-ObliqueShockFunction(gam,Min,theta,beta0))/dbeta
        beta = beta0 - ObliqueShockFunction(gam,Min,theta,beta0)/derivative
        error=abs(beta-beta0)
        beta0=beta
        n=n+1

    return beta
    
def TurningAngleObliqueShock(gam,Min,beta):
    tol=1e-06
    error=100
    n=0
    dtheta = 0.01*calc.pi/180
    theta0 =beta    

    #Using a newton raphson
    while error>tol and n<1000:
        derivative = (ObliqueShockFunction(gam,Min,theta0+dtheta,beta)-ObliqueShockFunction(gam,Min,theta0,beta))/dtheta
        theta = theta0 - ObliqueShockFunction(gam,Min,theta0,beta)/derivative
        error=abs(theta-theta0)
        theta0=theta
        n=n+1

    return beta
 
def ObliqueShock(gam,M1,theta): 
    beta1 = WaveAngleObliqueShock(gam,M1,theta,0)
    Mn1= M1*np.sin(beta1)
    TR=TemperatureRatioShock(gam,Mn1)   
    pR=PressureRatioShock(gam,Mn1) 
    p0R=TotalPressureRatioShock(gam,Mn1)
    Mn2=MachDownstreamShock(gam,Mn1) 
    M2 = Mn2/np.sin(beta1-theta)
    
    return beta1, M2, pR, TR, p0R
