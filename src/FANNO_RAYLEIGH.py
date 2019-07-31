# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 15:18:18 2018

@author: gjoterorodrigu
"""
import sys
import os        
import math as calc
import numpy as np
##############################################################
## Fanno flow (friction flow in a pipe with compressible flow)

def FannoFunction(gam,M):
    return (1/(gam*M*M))+(gam+1)*0.5/gam*np.log(M*M/(1+0.5*M*M*(gam-1)))
    
def ChockedLengthFanno(gam,M,f,D):
    return D/f/4*(FannoFunction(gam,M)-FannoFunction(gam,1))

def MachFromFannoFunction(gam,Fanno,SupersonicFlag=1):
    tol=1e-06
    error=100
    n=0
    dM= 0.005         #for the derivative
    if SupersonicFlag>0:  
        M0 =1.01           #first guess   
    else:
        M0 =0.99           #first guess     
    relaxfact = 0.15  #relaxation factor (derivative are to big)
    
    #Using a newton raphson
    while error>tol and n<1000:
        derivative = ((Fanno-FannoFunction(gam,M0+dM))-(Fanno-FannoFunction(gam,M0)))/dM
        #print("calculating Mach in rayleigh: ", Fanno,FannoFunction(gam,M0),M0)
        M = M0 - relaxfact*(Fanno-FannoFunction(gam,M0))/derivative
        error=abs(M-M0)
        M0=M
        n=n+1
        
    return abs(M)
    
def CriticalPressureRatioFanno(gam,M):  #p/p*
## function to calculate the critical pressure ratio 
#
#  Inputs:
#  gam: the insentropic exponent in ideal gas
#  M: Mach number = vel/a=vel/sqrt(gamPv*P*v)
#  (v: specific volume, p:pressure, vel:velocity magnitude, T:temperature, a:sound speed)
#  Outputs:
#  p/p*, where p* (critical pressure) is the pressure when M=1
    return 1/M*calc.sqrt((gam+1)/(2+(gam-1)*M*M))   

def CriticalTotalPressureRatioFanno(gam,M):  #p0/p0*
## function to calculate the critical pressure ratio 
#
#  Inputs:
#  gam: the insentropic exponent in ideal gas
#  M: Mach number = vel/a=vel/sqrt(gamPv*P*v)
#  (v: specific volume, p:pressure, vel:velocity magnitude, T:temperature, a:sound speed)
#  Outputs:
#  p0/p0*, where p0* (critical total pressure) is the pressure when M=1
    return ((1+(gam-1)*0.5*M*M)/(1+(gam-1)*0.5))**(gam/(gam-1))*CriticalPressureRatioFanno(gam,M)  
    
def RayleighFunction(gam,M):
    return (1+(gam-1)*0.5*M*M)*M*M/((1+gam*M*M)**2)

def MachFromRayleighFunction(gam,Rayleigh,SupersonicFlag=1):
    tol=1e-06
    error=100
    n=0
    dM= 0.005         #for the derivative
    if SupersonicFlag>0:  
        M0 =1.01           #first guess   
    else:
        M0 =0.99           #first guess     
    
    relaxfact = 0.15  #relaxation factor (derivative are to big)
    
    #Using a newton raphson
    while error>tol and n<1000:
        derivative = ((Rayleigh-RayleighFunction(gam,M0+dM))-(Rayleigh-RayleighFunction(gam,M0)))/dM
        #print("calculating Mach in rayleigh: ", Rayleigh,RayleighFunction(gam,M0),M0)
        M = M0 - relaxfact*(Rayleigh-RayleighFunction(gam,M0))/derivative
        error=abs(M-M0)
        M0=M
        n=n+1
        
    return abs(M)
   
def CriticalTotalTemperatureRatioRayleigh(gam,M):  #T0/T0*
## function to calculate the critical pressure ratio 
#
#  Inputs:
#  gam: the insentropic exponent in ideal gas
#  M: Mach number = vel/a=vel/sqrt(gamPv*P*v)
#  (v: specific volume, p:pressure, vel:velocity magnitude, T:temperature, a:sound speed)
#  Outputs:
#  T0/T0*, where p0* (critical total pressure) is the pressure when M=1
    return RayleighFunction(gam,M)/RayleighFunction(gam,1) 
    
def CriticalTemperatureRatioRayleigh(gam,M):  #T/T*
## function to calculate the critical pressure ratio 
#
#  Inputs:
#  gam: the insentropic exponent in ideal gas
#  M: Mach number = vel/a=vel/sqrt(gamPv*P*v)
#  (v: specific volume, p:pressure, vel:velocity magnitude, T:temperature, a:sound speed)
#  Outputs:
#  T/T*, where p0* (critical total pressure) is the pressure when M=1
    return (1/(1+(gam-1)*0.5*M*M))/(1/(1+(gam-1)*0.5))*CriticalTotalTemperatureRatioRayleigh(gam,M) 
    
def CriticalPressureRatioRayleigh(gam,M):  #p/p*
## function to calculate the critical pressure ratio 
#
#  Inputs:
#  gam: the insentropic exponent in ideal gas
#  M: Mach number = vel/a=vel/sqrt(gamPv*P*v)
#  (v: specific volume, p:pressure, vel:velocity magnitude, T:temperature, a:sound speed)
#  Outputs:
#  p/p*, where p0* (critical total pressure) is the pressure when M=1
    return (gam+1)/(1+gam*M*M) 
    
def CriticalTotalPressureRatioRayleigh(gam,M):  #p0/p0*
## function to calculate the critical pressure ratio 
#
#  Inputs:
#  gam: the insentropic exponent in ideal gas
#  M: Mach number = vel/a=vel/sqrt(gamPv*P*v)
#  (v: specific volume, p:pressure, vel:velocity magnitude, T:temperature, a:sound speed)
#  Outputs:
#  p0/p0*, where p0* (critical total pressure) is the pressure when M=1
    return (((2/(gam+1))*(1+(gam-1)*0.5*M*M))**(gam/(gam-1)))*CriticalPressureRatioRayleigh(gam,M)
    
def CriticalDensityRatioRayleigh(gam,M):  #rho/rho*
## function to calculate the critical pressure ratio 
#
#  Inputs:
#  gam: the insentropic exponent in ideal gas
#  M: Mach number = vel/a=vel/sqrt(gamPv*P*v)
#  (v: specific volume, p:pressure, vel:velocity magnitude, T:temperature, a:sound speed)
#  Outputs:
#  rho/rho*, where rho* (critical total pressure) is the pressure when M=1
    return CriticalPressureRatioRayleigh(gam,M)/CriticalTemperatureRatioRayleigh(gam,M)