# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 15:13:35 2018

@author: gjoterorodrigu
"""
import sys
import os        
import math as calc
import numpy as np

###########################################
## For real gases  
def TotalPressureToStaticPressureRatio(gamPv,M):  #p0/p
## function to calculate the total to static pressure ratio 
#
#  Inputs:
#  gamPv: -cp/cv v/P dp/dv|T (the insentropic exponent in ideal gas)
#  M: Mach number = vel/a=vel/sqrt(gamPv*P*v)
#  (v: specific volume, p:pressure, vel:velocity magnitude, T:temperature, a:sound speed)
#  Outputs:
#  p0/p
    return (1.0+(gamPv-1)*0.5*M*M)**(gamPv/(gamPv-1))

def TotalTemperatureToStaticTemperatureRatio(gamPv,M,gamTv):  #T0/T
## function to calculate the total to static temperature ratio 
#
#  Inputs:
#  gamPv: -cp/cv v/P dp/dv|T (the insentropic exponent in ideal gas)
#  M: Mach number = vel/a=vel/sqrt(gamPv*P*v)
#  gamTv: 1+v/cv dp/dT|v (the insentropic exponent in ideal gas)
#  (v: specific volume, p:pressure, vel:velocity magnitude, T:temperature, a:sound speed)
#  Outputs:
#  T0/T
    return (1.0+(gamPv-1)*0.5*M*M)**((gamTv-1)/(gamPv-1))
    
def TotalDensityToStaticDensityRatio(gamPv,M):  #rho0/rho
## function to calculate the total to static density ratio 
#
#  Inputs:
#  gamPv: -cp/cv v/P dp/dv|T (the insentropic exponent in ideal gas)
#  M: Mach number = vel/a=vel/sqrt(gamPv*P*v)
#  (v: specific volume, p:pressure, vel:velocity magnitude, T:temperature, a:sound speed)
#  Outputs:
#  rho0/rho
    return (1.0+(gamPv-1)*0.5*M*M)**(1/(gamPv-1))

def TotalCompressFactToStaticCompressFactRatio(gamPv,M,gamTv):  #z0/z
## function to calculate the total to static compressibility factor ratio 
#
#  Inputs:
#  gamPv: -cp/cv v/P dp/dv|T (the insentropic exponent in ideal gas)
#  M: Mach number = vel/a=vel/sqrt(gamPv*P*v)
#  (v: specific volume, p:pressure, vel:velocity magnitude, T:temperature, a:sound speed)
#  
#  Outputs:
#  z0/z  (z: compressibility factor as in: Pv=zRT)
    return (1.0+(gamPv-1)*0.5*M*M)**((gamPv-gamTv)/(gamPv-1))
    
def TotalSpeedSoundToStaticSpeedSoundRatio(gamPv,M):  #a0/a
## function to calculate the total to static sound speed ratio 
#
#  Inputs:
#  gamPv: -cp/cv v/P dp/dv|T (the insentropic exponent in ideal gas)
#  M: Mach number = vel/a=vel/sqrt(gamPv*P*v)
#  (v: specific volume, p:pressure, vel:velocity magnitude, T:temperature, a:sound speed)
#  Outputs:
#  a0/a
    return (1.0+(gamPv-1)*0.5*M*M)
    
def CriticalMach(gamPv,M):  #Mach*
## function to calculate the critical Mach number 
#
#  Inputs:
#  gamPv: -cp/cv v/P dp/dv|T (the insentropic exponent in ideal gas)
#  M: Mach number = vel/a=vel/sqrt(gamPv*P*v)
#  (v: specific volume, p:pressure, vel:velocity magnitude, T:temperature, a:sound speed)
#  Outputs:
#  M*=vel/a*, where a* is the critical sound speed when M=1
    return (((gamPv+1)*M*M)/(2+(gamPv-1)*M*M))**(0.5)

def CriticalAreaRatio(gamPv,M):  #A/A*
## function to calculate the critical area ratio 
#
#  Inputs:
#  gamPv: -cp/cv v/P dp/dv|T (the insentropic exponent in ideal gas)
#  M: Mach number = vel/a=vel/sqrt(gamPv*P*v)
#  (v: specific volume, p:pressure, vel:velocity magnitude, T:temperature, a:sound speed)
#  Outputs:
#  A/A*, where A* (critical area) stands for the throat area where M=1
    return (1/M)*((2+(gamPv-1)*M*M)/(gamPv+1))**((gamPv+1)/(2*(gamPv-1)))
    
def CriticalPressureRatio(gamPv,M):  #p*/p0
## function to calculate the critical pressure ratio 
#
#  Inputs:
#  gamPv: -cp/cv v/P dp/dv|T (the insentropic exponent in ideal gas)
#  M: Mach number = vel/a=vel/sqrt(gamPv*P*v)
#  (v: specific volume, p:pressure, vel:velocity magnitude, T:temperature, a:sound speed)
#  Outputs:
#  p/p*, where p* (critical pressure) is the pressure when M=1
    return (2/(gamPv+1))**(gamPv/(gamPv-1))

def CriticalDensityRatio(gamPv,M):  #rho*/rho0
## function to calculate the critical density ratio 
#
#  Inputs:
#  gamPv: -cp/cv v/P dp/dv|T (the insentropic exponent in ideal gas)
#  M: Mach number = vel/a=vel/sqrt(gamPv*P*v)
#  (v: specific volume, p:pressure, vel:velocity magnitude, T:temperature, a:sound speed)
#  Outputs:
#  rho/rho*, where rho* (critical density) is the density when M=1
    return (2/(gamPv+1))**(1/(gamPv-1))
    
def CriticalTemperatureRatio(gamPv,M,gamTv):  #T*/T0
## function to calculate the critical temperature ratio 
#
#  Inputs:
#  gamPv: -cp/cv v/P dp/dv|T (the insentropic exponent in ideal gas)
#  M: Mach number = vel/a=vel/sqrt(gamPv*P*v)
#  gamTv: 1+v/cv dp/dT|v (the insentropic exponent in ideal gas)
#  (v: specific volume, p:pressure, vel:velocity magnitude, T:temperature, a:sound speed)
#  Outputs:
#  T/T*, where T* (critical temperature) is the temperature when M=1
    return (2/(gamPv+1))**((gamTv-1)/(gamPv-1))
    
def CriticalCompressFactRatio(gamPv,M,gamTv):  #z/z*
## function to calculate the critical temperature ratio 
#
#  Inputs:
#  gamPv: -cp/cv v/P dp/dv|T (the insentropic exponent in ideal gas)
#  M: Mach number = vel/a=vel/sqrt(gamPv*P*v)
#  gamTv: 1+v/cv dp/dT|v (the insentropic exponent in ideal gas)
#  (v: specific volume, p:pressure, vel:velocity magnitude, T:temperature, a:sound speed)
#  Outputs:
#  z/z*, where z* (critical compressibility factor) is the compressibility when M=1
#  where pv=ZRT
    return (2/(gamPv+1))**((gamPv-gamTv)/(gamPv-1))

def CriticalMassFlow(Astar,rho0,p0,gamPv,M):  #m*
## function to calculate the critical mass flow 
#
#  Inputs:
#  Astar: area at the throat (minimal area)
#  rho0: stagnation density
#  p0: stagnation pressure
#  gamPv: -cp/cv v/P dp/dv|T (the insentropic exponent in ideal gas)
#  M: Mach number = vel/a=vel/sqrt(gamPv*P*v)
#  (v: specific volume, p:pressure, vel:velocity magnitude, T:temperature, a:sound speed)
#  Outputs:
#  m*, where m* (critical mass flow) is the chocked mass flow for the given pressure, density and area
    return (Astar*(gamPv*rho0*p0)**0.5)*(2/(gamPv+1))**((gamPv+1)/(2*(gamPv-1))) 