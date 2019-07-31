# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 15:26:30 2018

@author: gjoterorodrigu
"""
import sys
import os        
import math as calc
import numpy as np
def HypersonicEquivalencePrinciple(gam,Up,T,R):  
    
    # Procedure to calculate the "Piston" Velocity
    def function(M,gam,T1,R):
        sos1=(gam*R*T1)**0.5
        Us=sos1*M
        T2=T1*TemperatureRatioShock(gam,M)
        sos2=(gam*R*T2)**0.5
        V2=MachDownstreamShock(gam,M)*sos2
        return Us-V2
                
    
    #Using a newton raphson
    tol=1e-06
    error=100
    n=0
    dM=0.001
    M0=1.1
    relaxfact=0.85  
    while error>tol and n<1000:
        derivative = ((Up-function(M0+dM,gam,T,R))-(Up-function(M0,gam,T,R)))/dM           
        M = M0 - relaxfact*(Up-function(M0,gam,T,R))/derivative
        error=abs(M-M0)
        M0=M              
        n=n+1    
        #print("errror ", error, function(M0,gam,T,R), Up, M0*(gam*R*T)**0.5)
        
    sos1=(gam*R*T)**0.5
    Ushock= sos1*M
    return Ushock