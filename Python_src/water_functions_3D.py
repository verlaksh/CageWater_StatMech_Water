import os ; import sys ; import matplotlib.pyplot as plt ; from math import *; import numpy as np
from scipy.special import erf ; from scipy.special import expi ; #from scipy.integrate import quad 
from default_parameters import * ; import mpmath as mp 


#Note: This code contains the auxilary functions called in watermain_2D.py
#      para.py imported here contains the contstant parameter necessary

# r_HB and eps_HB real units
real_r	   = r_HB_real	# r_HB in nm 
#real_e     = para.eps_HB_real	# HB energy in kJ/mol
#T_real     = para.T_real        # T in K
#P_real     = para.P_real        # P in bar

# Constants
#NA         = 6.022*10**23	# Avogadro number in per mole
#kB         = 1.381*10**-23      # Boltzmann constant in J/K

# Ideal gases and nonpolar molecules
def ss(gas):                                              #####              
   if   gas == 'He' :  sig_s = 0.727999999                #####   HELIUM     
   elif gas == 'Ne' :  sig_s = 0.804363636                #####   NEON       
   elif gas == 'Ar' :  sig_s = 0.987636364                #####   ARGON      
   elif gas == 'Kr' :  sig_s = 1.053818182                #####   KRYPTON    
   elif gas == 'Xe' :  sig_s = 1.160727273                #####   XENON      
   elif gas == 'Rn' :  sig_s = 1.221818181                #####   REDON      
   elif gas == 'Me' :  sig_s = 1.048727273                #####   METHANE    
   elif gas == 'bn' :  sig_s = 1.710545454                #####   BENZENE    
   elif gas == 'np' :  sig_s = 2.271563636                #####   NAPHTHALENE
   elif gas == 'C60':  sig_s = 2.548                      #####   C60        
   elif gas == 'Li' :  sig_s = 0.24                       #####   Lithium
   elif gas == 'Na' :  sig_s = 0.37                       #####   Sodium
   elif gas == 'K'  :  sig_s = 0.52                       #####   Postassium
   elif gas == 'Rb' :  sig_s = 0.58                       #####   Rubidium
   elif gas == 'Cs' :  sig_s = 0.66                       #####   Cesium
   elif gas == 'F'  :  sig_s = 0.53                       #####   Florine
   elif gas == 'Cl' :  sig_s = 0.71                       #####   Chlorine
   elif gas == 'Br' :  sig_s = 0.77                       #####   Bromine
   elif gas == 'I'  :  sig_s = 0.85                       #####   Iodine
   else:               #sig_s = float(sys.argv()[1])       #####                                 #####
       exit("gas name should be one of 'He','Ne','Ar','Kr','Xe','Rn','Me','bn','np','C60' ")
   return sig_s 

#a,eps_LJ,eps_c,KT_HB,KT_LJ,Kx,Kth,r_d,r_mS,r_mHB,r_mLJ,xHB,xLJ = opt(a,eps_LJ,eps_c,KT_HB,KT_LJ,Kx,Kth,real_e,a,eps_LJ,eps_c,KT_HB,KT_LJ,Kx,Kth,r_d,r_mS,r_mHB,r_mLJ,xHB,xLJ)

## volume per water molecule in bulk in different state 
def bulk_wat(p_0,p,T,real_e,a,eps_LJ,eps_c,KT_HB,KT_LJ,Kx,Kth,r_d,r_mS,r_mHB,r_mLJ,xHB,xLJ):  
   # volume factors 
   A_S        = (8.0*(sqrt(3.0))/9.0)
   A_HB       = A_S*xHB
   A_LJ       = A_S*xLJ
   A_O        = A_S*1700
   
   v_s        = A_S*sigma_LJ**3.0                # solid cage
   v_HB       = A_HB*sigma_LJ**3.0
   v_LJ       = A_LJ*(sigma_LJ*2**(1/6.0))**3.0
   v_d        = (sigma_LJ*(1-r_d))**3.0

   v_s,v_HB,v_LJ,v_O = volume(T,p,A_HB,A_S,A_LJ,A_O,r_d,r_mS,r_mHB,r_mLJ)

   # High Density case

   delta_HB =del_HB(p,T,A_HB,sigma_LJ*(1-r_d),sigma_LJ*(1+r_mHB),KT_HB,Kx,Kth,eps_HB)
   delta_LJ =del_LJ(p,T,A_LJ,sigma_LJ*(1-r_d),sigma_LJ*(1+r_mLJ),KT_LJ,Kx,eps_LJ)
   delta_O  =del_O(p,T,A_O) #,sigma_LJ*(1-r_d))
   delta_s  =del_S(p,T,A_S,sigma_LJ*(1-r_d),sigma_LJ*(1+r_mS),KT_HB,Kx,Kth,eps_HB)

   # Total partition function of single 12 membered cage
   Q1 = ( ( delta_HB + delta_LJ + delta_O )**12.0 )  -  ( ( delta_HB )**12.0 )  +  exp(-12*eps_c/T) * delta_s**12.0  #####                                      

   # fractional populations 
   fHB = ( delta_HB *(( delta_HB + delta_LJ + delta_O )**11.0) - (delta_HB**12.0) ) / Q1  # fractional populations of HB state                      
   fLJ = ( delta_LJ *(( delta_HB + delta_LJ + delta_O )**11.0) )/ Q1                     # fractional populations of LJ state                      
   fO  = ( delta_O  *(( delta_HB + delta_LJ + delta_O )**11.0) ) / Q1                    # fractional populations of OPEN(noninetracting) state    
   fS  = (exp(-12*eps_c/T) * delta_s**12.0 )/ Q1                                          # fractional populations of SOLID state                   

   # Molar volume
   v =  (fS * v_s ) + (fHB * v_HB) + (fLJ * v_LJ) + (fO * v_O)

   # Average energy in bulk and partition function
   av_uHB,av_uS,av_uLJ,av_uO = avgE(T,r_d,r_mS,r_mHB,r_mLJ,KT_HB,KT_LJ,Kx,Kth)
   avg_ene_bulk_water,q_b = Bulk_EQ(p,T,av_uHB,av_uS,av_uLJ,av_uO,fHB,fLJ,fO,fS,v) 

   return v,q_b,Q1,fHB,fS,fLJ,fO,v_HB,v_s,v_LJ,v_O,av_uHB,av_uS,av_uLJ,av_uO,avg_ene_bulk_water

   ################


# Partition function contribution from different states in the bulk

# HB state 
def del_HB(p,T,A,rl,rh,KT_HB,Kx,Kth,eps_HB):
    tran_term  = 4.0*pi*T/(3.0*A*p) * (exp(-A*p*(rl**3.0)/T)-exp(-A*p*(rh**3.0)/T))    
    rota_term  = 4.0*pi*pi*sqrt(pi*T/(8.0*Kth))*erf((2.0-sqrt(3.0))*sqrt((Kth)/(2.0*T) ) )                   
    const_term = exp(2.0*eps_HB/T-KT_HB)
    rl = rl -sigma_LJ; rh = rh- sigma_LJ
    Flex_term  = T*(exp(-Kx*rl**2/T) - exp(-Kx*rh**2/T))**2/(pi*(erf(sqrt(Kx/T)*rh) - erf(sqrt(Kx/T)*rl))**2)
    Flex_term  = exp(-Flex_term/T)
    delta_HB = const_term*tran_term*rota_term*Flex_term                                   
    return delta_HB

# LJ state 
def del_LJ(p,T,A,rl,rh,KT_LJ,Kx,eps_LJ):
    tran_term = 4.0*pi*T/(3.0*A*p) * (exp(-A*p*(rl**3.0)/T)-exp(-A*p*(rh**3.0)/T)) 
    rota_term = 8.0*pi*pi
    const_term = exp(2.0*eps_LJ/T-KT_LJ)
    rl = rl -sigma_LJ; rh = rh- sigma_LJ
    Flex_term  = T*(exp(-Kx*rl**2/T) - exp(-Kx*rh**2/T))**2/(pi*(erf(sqrt(Kx/T)*rh) - erf(sqrt(Kx/T)*rl))**2)
    Flex_term  = exp(-Flex_term/T)
    delta_LJ = const_term*tran_term*rota_term*Flex_term        
    return delta_LJ                          #####                 

# Open state
def del_O(p,T,A):
    tran_term = 4.0*pi*T/(3.0*A*p)#*exp(-A*p*(rl**3.0)/T) 
    rota_term = 8.0*pi*pi
    const_term = 1.0              
    delta_O = const_term*tran_term*rota_term
    return delta_O

# Solid state
def del_S(p,T,A,rl,rh,KT_HB,Kx,Kth,eps_HB):
    tran_term = 4.0*pi*T/(3.0*A*p) * (exp(-A*p*(rl**3.0)/T)-exp(-A*p*(rh**3.0)/T))
    rota_term = 4.0*pi*pi*sqrt(pi*T/(8.0*Kth))*erf((2.0-sqrt(3.0))*sqrt((Kth)/(2.0*T) ) ) 
    const_term = exp(2.0*eps_HB/T-KT_HB)
    rl = rl -sigma_LJ; rh = rh- sigma_LJ
    Flex_term  = T*(exp(-Kx*rl**2/T) - exp(-Kx*rh**2/T))**2/(pi*(erf(sqrt(Kx/T)*rh) - erf(sqrt(Kx/T)*rl))**2)
    Flex_term  = exp(-Flex_term/T)
    delta_S = const_term*tran_term*rota_term*Flex_term
    return delta_S    

# Partition function and population
def population(delta_HB,delta_LJ,delta_O,delta_s,T):
    #    # Total partition function of single hexagon
    Q1 = ( ( delta_HB + delta_LJ + delta_O )**12.0 )  -  ( ( delta_HB )**12.0 )  +  exp(-12*eps_c/T) * delta_s**12.0  #####                                      
 
    # fractional populations 
    fHB = ( delta_HB *(( delta_HB + delta_LJ + delta_O )**11.0) - (delta_HB**12.0) ) / Q1  # fractional populations of HB state                      
    fLJ = ( delta_LJ *(( delta_HB + delta_LJ + delta_O )**11.0) )/ Q1                     # fractional populations of LJ state                      
    fO  = ( delta_O  *(( delta_HB + delta_LJ + delta_O )**11.0) ) / Q1                    # fractional populations of OPEN(noninetracting) state    
    fS  = (exp(-12*eps_c/T) * delta_s**12.0 )/ Q1                                          # fractional populations of SOLID state                   
    return Q1,fHB,fLJ,fO,fS


def Bulk_EQ(p,T,av_uHB,av_uS,av_uLJ,av_uO,fHB,fLJ,fO,fS,v):
    # Average energy in bulk state
    avg_ene_bulk = ( av_uHB*fHB +  av_uLJ*fLJ + av_uS*fS + av_uO*fO  )
    #partition function of bulk water
    q_b = v*( 8.0*pi*pi) * mp.exp( -(avg_ene_bulk + p*v)/T)                              #####
    return avg_ene_bulk,q_b



# Molar  Volume  per  Water  Molecule  in  the  First  Solvation  Shell  of  the  Solute
def v_shell(v,sigma_LJ,sigma_s):
    t_a =((3.0/4.0)*(v/pi))**(1.0/3.0)                       # radius of bulk H2O
    t_b = (2.0**(1.0/6.0))*sigma_s/2.0                       # radius of solute
    t_d = (2.0**(1.0/6.0))*(sigma_LJ + sigma_s)/2.0          # distance between H2O and sollute
    if t_a<=t_b+t_d:
       dv = pi*((t_a+t_b-t_d)**2)*(-3*t_a*t_a+2*t_a*t_d+2*t_b*t_d-3*t_b*t_b+6*t_a*t_b+t_d*t_d)/(12*t_d)
    else:
        dv = dv = (4./3.)*pi*(t_b**3.)
    v_h = v - dv                                 #molar volume of H2O in solvation shell
    return v_h

# average zeta - factor accounting for H-bond occulation due to solute size
def avg_zeta(sigma_s):                                                                    
    alpha = (pi-acos(1./3.))                                                ### TETRAHEDRAL ANGLE 
    phi_c = critical_angle(sigma_s )                                                         
    ct0 = cos(phi_c) ; st0 = sin(phi_c)                                                           
    if (phi_c <= alpha/2. ):                                                                      
       V3 = 0. ; V2 = 0. ; V1 = 4.*v1(ct0,st0) ; V0 = (1.-V1-V2-V3)                               
    elif alpha/2. < phi_c and phi_c <= pi-alpha :                                                 
       V3 = 0. ; V2 = 6.*v2(ct0,st0) ; V1 = 4.*v1(ct0,st0) - 12.*v2(ct0,st0) ; V0 = (1.-V1-V2-V3) 
    else :                                                                                        
       V3 = 4.*v3(ct0,st0) ; V2 = 6.*v2(ct0,st0) - 12.*v3(ct0,st0)                                
       V1 = 4.*v1(ct0,st0) - 12.*v2(ct0,st0) + 12*v3(ct0,st0) ; V0 = (1.-V1-V2-V3)                
    avg_n_HB = 4.*V0 + 3.*V1 + 2.*V2 + 1.*V3                                  ### AVG nHB per H2O 
    return avg_n_HB                                                                               

def v1(ct0,st0): return (1.-ct0)/2.                                                                 # 1 ar_m occultation ##
def v2(ct0,st0):                                                                                    #                   ##
    return (.50/pi)*acos((1./3+ct0**2)/st0**2) - (1.0/pi)*ct0*(acos(sqrt(2.)*ct0/st0))              # 2 ar_m occultation ##
def v3(ct0,st0):                                                                                                     
    return (.75/pi)*acos((1./3+ct0**2)/st0**2) - (1.5/pi)*ct0*(acos(sqrt(2.)*ct0/st0)-pi/3.) - 0.25 # 3 ar_m occultation ##



# AVG. SOLV ENG IN 1ST HYDRATION SHELL
def avg_ene_solvation_wat(zeta,avg_en_HB,fHB_h,avg_en_S,fS_h,avg_en_LJ,fLJ_h,eps_SW,avg_elec):
    eps_h = (1.0/2.0)*(zeta*(avg_en_HB*fHB_h + avg_en_S*fS_h) + 4.*avg_en_LJ*fLJ_h - eps_SW+avg_elec)                                                     
    return eps_h                     

def q_shell(p_h,T,v_h,avg_en_shell):
    q_h  = 8.*pi*pi*exp((-(avg_en_shell + p_h*v_h))/T)
    return q_h


# Meanfield approximation

#low precision
def meanfield_pressure(p_0,p,T,v_0,real_e,a,eps_LJ,eps_c,KT_HB,KT_LJ,Kx,Kth,r_d,r_mS,r_mHB,r_mLJ,xHB,xLJ):
    p1=p_0
    p2 =p_0+(a/(v_0*v_0))
    i = 0
    while abs((p2-p1)/(p1+p2))>10**-5 and  i<=10**6:
        if i>20000 and abs((p2-p1)/(p1+p2))>10**-6 and p_0<0:
            #p1 = 0.1*(p1 +p2)/2. 
            p1 = p1 +0.5*(p2-p1)
        elif i>100000 and abs((p2-p1)/(p1+p2))>10**-7 and p_0>0:
            p1 = (p1 +p2)/2.
        else:
            p1 = p2 #+ 0.5*abs(p2-p1)
        v1 = ((bulk_wat(p_0,p1,T,real_e,a,eps_LJ,eps_c,KT_HB,KT_LJ,Kx,Kth,r_d,r_mS,r_mHB,r_mLJ,xHB,xLJ))[0])
        p2 = p_0+(a/v1**2)                  # checking condition p_0=P+a/v_0^2
        i = i+1
        if i==10**5:
              p2 = 0.5*(p1+p2)
              break;
    print("i=",i)
    return p2

def kap(p_0,p,T,real_e,a,eps_LJ,eps_c,KT_HB,KT_LJ,Kx,Kth,r_d,r_mS,r_mHB,r_mLJ,xHB,xLJ):
    dp = 0.0001
    pp = p*(1+dp)
    pm = p*(1-dp)
    v,q_b,Q1,fHBpm,fSpm,fLJpm,fOpm,v_HBpm,v_spm,v_LJpm,v_Opm,av_uHB,av_uS,av_uLJ,av_uO,avg_ene_bulk_water = bulk_wat(p_0,pm,T,real_e,a,eps_LJ,eps_c,KT_HB,KT_LJ,Kx,Kth,r_d,r_mS,r_mHB,r_mLJ,xHB,xLJ)
    v,q_b,Q1,fHBpp,fSpp,fLJpp,fOpp,v_HBpp,v_spp,v_LJpp,v_Opp,av_uHB,av_uS,av_uLJ,av_uO,avg_ene_bulk_water = bulk_wat(p_0,pp,T,real_e,a,eps_LJ,eps_c,KT_HB,KT_LJ,Kx,Kth,r_d,r_mS,r_mHB,r_mLJ,xHB,xLJ)
    vmol,q_b,Q1,fHB,fS,fLJ,fO,v_HB,v_s,v_LJ,v_O,av_uHB,av_uS,av_uLJ,av_uO,avg_ene_bulk_water = bulk_wat(p_0,p,T,real_e,a,eps_LJ,eps_c,KT_HB,KT_LJ,Kx,Kth,r_d,r_mS,r_mHB,r_mLJ,xHB,xLJ)

    del_p = 2*p*dp
    kaHB =-(1.0/vmol)*(v_HB*(fHBpp-fHBpm) +fHB*(v_HBpp-v_HBpm)) /del_p
    kaS  =-(1.0/vmol)*(v_s*(fSpp-fSpm)    +fS*(v_spp-v_spm))    /del_p
    kaLJ =-(1.0/vmol)*(v_LJ*(fLJpp-fLJpm) +fLJ*(v_LJpp-v_LJpm)) /del_p
    kaO  =-(1.0/vmol)*(v_O*(fOpp-fOpm)    +fO*(v_Opp-v_Opm))    /del_p

    ka = kaHB + kaS + kaLJ +kaO
    return ka,kaHB,kaS,kaLJ,kaO

def alp(p_0,p,T,real_e,a,eps_LJ,eps_c,KT_HB,KT_LJ,Kx,Kth,r_d,r_mS,r_mHB,r_mLJ,xHB,xLJ):
    dt = 0.001
    tp = T*(1+dt)
    tm = T*(1-dt)
    vtm,q_b,Q1,fHBtm,fStm,fLJtm,fOtm,v_HBtm,v_stm,v_LJtm,v_Otm,av_uHB,av_uS,av_uLJ,av_uO,avg_ene_bulk_water = bulk_wat(p_0,p,tm,real_e,a,eps_LJ,eps_c,KT_HB,KT_LJ,Kx,Kth,r_d,r_mS,r_mHB,r_mLJ,xHB,xLJ)
    vtp,q_b,Q1,fHBtp,fStp,fLJtp,fOtp,v_HBtp,v_stp,v_LJtp,v_Otp,av_uHB,av_uS,av_uLJ,av_uO,avg_ene_bulk_water = bulk_wat(p_0,p,tp,real_e,a,eps_LJ,eps_c,KT_HB,KT_LJ,Kx,Kth,r_d,r_mS,r_mHB,r_mLJ,xHB,xLJ)
    vmol,q_b,Q1,fHB,fS,fLJ,fO,v_HB,v_s,v_LJ,v_O,av_uHB,av_uS,av_uLJ,av_uO,avg_ene_bulk_water = bulk_wat(p_0,p,T,real_e,a,eps_LJ,eps_c,KT_HB,KT_LJ,Kx,Kth,r_d,r_mS,r_mHB,r_mLJ,xHB,xLJ)

    del_T = 2.0*dt*T
    alHB =(1.0/vmol)*(v_HB*(fHBtp-fHBtm) +fHB*(v_HBtp-v_HBtm)) /del_T
    alS  =(1.0/vmol)*(v_s*(fStp-fStm)    +fS*(v_stp-v_stm))    /del_T
    alLJ =(1.0/vmol)*(v_LJ*(fLJtp-fLJtm) +fLJ*(v_LJtp-v_LJtm)) /del_T
    alO  =(1.0/vmol)*(v_O*(fOtp-fOtm)    +fO*(v_Otp-v_Otm))    /del_T

    al = alHB + alS + alLJ +alO
    return al,alHB,alS,alLJ,alO

def avgE(T,r_d,r_mS,r_mHB,r_mLJ,KT_HB,KT_LJ,Kx,Kth):
    av_uHB1 = T/2.0 - 2.0*eps_HB + T*KT_HB
    av_uHB2 = (2.0-sqrt(3))*sqrt(T*Kth)*exp((-7+4.0*sqrt(3.0))*Kth/(2.0*T))
    av_uHB2 = av_uHB2/(sqrt(2.0*pi)*erf((-2+sqrt(3))*sqrt(Kth/(2.0*T))))
    av_uHB = av_uHB1 + av_uHB2 
    av_uS  = av_uHB + 1*eps_c + flex_bond(T,r_d,r_mS,Kx)[0] #+0.5*T
    av_uHB = av_uHB + flex_bond(T,r_d,r_mHB,Kx)[0]
    av_uLJ = -2.0*(eps_LJ)+ flex_bond(T,r_d,r_mLJ,Kx)[0] + T*KT_LJ
    av_uO  = 0
    return av_uHB,av_uS,av_uLJ,av_uO
    
def dudt(T,r_d,r_mS,r_mHB,r_mLJ,KT_HB,KT_LJ,Kx,Kth):
    dudt_HB1 = (-7.0+4.0*sqrt(3.0))*Kth*exp((-7+4.0*sqrt(3.0))*Kth/(2.0*T))
    dudt_HB1 = dudt_HB1/(2.0*pi*T*(erf((-2+sqrt(3))*sqrt(Kth/(2.0*T))))**2)
    dudt_HB2 = sqrt(Kth)*((26-15*sqrt(3))*Kth + (2-sqrt(3))*T)*exp((-7+4.0*sqrt(3.0))*Kth/(2.0*T))
    dudt_HB2 = dudt_HB2/(sqrt(8.0*pi*T)*T*erf((-2+sqrt(3))*sqrt(Kth/(2.0*T))))
    dudt_HB  = 1/2.0 + dudt_HB1 + dudt_HB2 + KT_HB 
    dudt_S   = dudt_HB + flex_bond(T,r_d,r_mS,Kx)[1] #+0.5
    dudt_LJ  = flex_bond(T,r_d,r_mLJ,Kx)[1] + KT_LJ 
    dudt_HB  = dudt_HB  + flex_bond(T,r_d,r_mHB,Kx)[1]
    dudt_O   = 0
    return dudt_HB,dudt_S,dudt_LJ,dudt_O

def cp(p_0,p,T,real_e,a,eps_LJ,eps_c,KT_HB,KT_LJ,Kx,Kth,r_d,r_mS,r_mHB,r_mLJ,xHB,xLJ,al0):
    dt = 0.001
    tp = T*(1+dt)
    tm = T*(1-dt)
    
    v,q_b,Q1,fHB,fS,fLJ,fO,v_HB,v_s,v_LJ,v_O,av_uHB,av_uS,av_uLJ,av_uO,avg_ene_bulk_water = bulk_wat(p_0,p,T,real_e,a,eps_LJ,eps_c,KT_HB,KT_LJ,Kx,Kth,r_d,r_mS,r_mHB,r_mLJ,xHB,xLJ)
    vtm,q_b,Q1,fHBtm,fStm,fLJtm,fOtm,v_HBtm,v_stm,v_LJtm,v_Otm,av_uHB,av_uS,av_uLJ,av_uO,avg_ene_bulk_water = bulk_wat(p_0,p,tm,real_e,a,eps_LJ,eps_c,KT_HB,KT_LJ,Kx,Kth,r_d,r_mS,r_mHB,r_mLJ,xHB,xLJ)
    vtp,q_b,Q1,fHBtp,fStp,fLJtp,fOtp,v_HBtp,v_stp,v_LJtp,v_Otp,av_uHB,av_uS,av_uLJ,av_uO,avg_ene_bulk_water = bulk_wat(p_0,p,tp,real_e,a,eps_LJ,eps_c,KT_HB,KT_LJ,Kx,Kth,r_d,r_mS,r_mHB,r_mLJ,xHB,xLJ)

    # sum(dfihi)
    av_uHB,av_uS,av_uLJ,av_uO = avgE(T,r_d,r_mS,r_mHB,r_mLJ,KT_HB,KT_LJ,Kx,Kth)

    uHB = av_uHB *fHB
    uS  = av_uS  *fS
    uLJ = av_uLJ *fLJ
    uO  = av_uO  *fO
    u   = uHB+uS+uLJ+uO

    hHB = (av_uHB + p*v_HB)*fHB
    hS  = (av_uS  + p*v_s )*fS
    hLJ = (av_uLJ + p*v_LJ)*fLJ
    hO  = (av_uO  + p*v_O )*fO
    h   = hHB+hS+hLJ+hO

    dudt_HB,dudt_S,dudt_LJ,dudt_O = dudt(T,r_d,r_mS,r_mHB,r_mLJ,KT_HB,KT_LJ,Kx,Kth)
    del_T = 2.0*dt*T

    dufHB = av_uHB *(fHBtp-fHBtm)/del_T + dudt_HB*fHB
    dufS  = av_uS  *(fStp-fStm)/del_T   + dudt_S*fS
    dufLJ = av_uLJ *(fLJtp-fLJtm)/del_T + dudt_LJ*fLJ
    dufO  = av_uO  *(fOtp-fOtm)/del_T   + dudt_O*fO

    dvfHB  = (v_HB*(fHBtp-fHBtm) +fHB*(v_HBtp-v_HBtm)) /del_T
    dvfS   = (v_s *(fStp-fStm)   +fS *(v_stp-v_stm))   /del_T
    dvfLJ  = (v_LJ*(fLJtp-fLJtm) +fLJ*(v_LJtp-v_LJtm)) /del_T
    dvfO   = (v_O *(fOtp-fOtm)   +fO *(v_Otp-v_Otm))   /del_T

    dfdt  = ((fHBtp-fHBtm) + (fStp-fStm) + (fLJtp-fLJtm) + (fOtp-fOtm))/del_T

    cpHB = dufHB + p*dvfHB
    cpS  = dufS  + p*dvfS
    cpLJ = dufLJ + p*dvfLJ
    cpO  = dufO  + p*dvfO
    cpm  = 2*a*(al0-1*dfdt)/v
    Cp   = cpHB + cpS + cpLJ + cpO

    return h,hHB,hO,hLJ,hO,Cp,cpHB,cpS,cpLJ,cpO,cpm,dufHB,dufS,dufLJ,dufO,dvfHB,dvfS,dvfLJ,dvfO

def cv(p_0,p,T,real_e,a,eps_LJ,eps_c,KT_HB,KT_LJ,Kx,Kth,r_d,r_mS,r_mHB,r_mLJ,xHB,xLJ,al0):
    dt = 0.001
    tp = T*(1+dt)
    tm = T*(1-dt)

    v,q_b,Q1,fHB,fS,fLJ,fO,v_HB,v_s,v_LJ,v_O,av_uHB,av_uS,av_uLJ,av_uO,avg_ene_bulk_water = bulk_wat(p_0,p,T,real_e,a,eps_LJ,eps_c,KT_HB,KT_LJ,Kx,Kth,r_d,r_mS,r_mHB,r_mLJ,xHB,xLJ)
    vtm,q_b,Q1,fHBtm,fStm,fLJtm,fOtm,v_HBtm,v_stm,v_LJtm,v_Otm,av_uHB,av_uS,av_uLJ,av_uO,avg_ene_bulk_water = bulk_wat(p_0,p,tm,real_e,a,eps_LJ,eps_c,KT_HB,KT_LJ,Kx,Kth,r_d,r_mS,r_mHB,r_mLJ,xHB,xLJ)
    vtp,q_b,Q1,fHBtp,fStp,fLJtp,fOtp,v_HBtp,v_stp,v_LJtp,v_Otp,av_uHB,av_uS,av_uLJ,av_uO,avg_ene_bulk_water = bulk_wat(p_0,p,tp,real_e,a,eps_LJ,eps_c,KT_HB,KT_LJ,Kx,Kth,r_d,r_mS,r_mHB,r_mLJ,xHB,xLJ)

    # sum(dfihi)
    av_uHB,av_uS,av_uLJ,av_uO = avgE(T,r_d,r_mS,r_mHB,r_mLJ,KT_HB,KT_LJ,Kx,Kth)

    uHB = av_uHB *fHB
    uS  = av_uS  *fS
    uLJ = av_uLJ *fLJ
    uO  = av_uO  *fO
    u   = uHB+uS+uLJ+uO

    dudt_HB,dudt_S,dudt_LJ,dudt_O = dudt(T,r_d,r_mS,r_mHB,r_mLJ,KT_HB,KT_LJ,Kx,Kth)
    del_T = 2.0*dt*T

    dufHB = av_uHB *(fHBtp-fHBtm)/del_T + dudt_HB*fHB
    duf_S = av_uS  *(fStp-fStm)/del_T   + dudt_S*fS
    dufLJ = av_uLJ *(fLJtp-fLJtm)/del_T + dudt_LJ*fLJ
    dufO  = av_uO  *(fOtp-fOtm)/del_T   + dudt_O*fO

    dfdt  = ((fHBtp-fHBtm) + (fStp-fStm) + (fLJtp-fLJtm) + (fOtp-fOtm))/del_T

    cvHB = dufHB
    cvS  = duf_S
    cvLJ = dufLJ
    cvO  = dufO
    cvm  = a*(al0-dfdt)/v
    Cv   = cvHB + cvS + cvLJ + cvO
    
    return u,uHB,uO,uLJ,uO,Cv,cvHB,cvS,cvLJ,cvO,cvm

def flex_bond(T,rl,rh,Kx):
    rl = -rl*sigma_LJ; rh = rh*sigma_LJ
    avgU = T*(exp(-Kx*rl**2/T) - exp(-Kx*rh**2/T))**2/(pi*(erf(sqrt(Kx/T)*rh) - erf(sqrt(Kx/T)*rl))**2)
    dudt1 = (-2*sqrt(T*Kx)*(exp(Kx*rh**2/T) - exp(Kx*rl**2/T))*(exp(Kx*rh**2/T)*rl - exp(-Kx*rl**2/T)*rh))
    dudt2 = sqrt(pi)*(exp(Kx*(rh**2+2*rl**2)/T)*(2*Kx*rh**2 + T)-exp(Kx*(2*rh**2+rl**2)/T)*(2*Kx*rl**2 + T))
    dudt2 = dudt2*(-erf(sqrt(Kx/T)*rh) + erf(sqrt(Kx/T)* rl))
    dudt = exp(-2*Kx*(rh**2-rl**2)/T)*(dudt1+dudt2)
    dudt = dudt/(2*T*sqrt(Kx*T)*pi*(erf(sqrt(Kx/T)*rh) - erf(sqrt(Kx/T)*rl))**2)
    return avgU,dudt

def volume(T,p,AHB,AS,ALJ,AO,r_d,r_mS,r_mHB,r_mLJ):
    r_d   = sigma_LJ*(1-r_d)
    r_mS  = sigma_LJ*(1+r_mS)
    r_mHB = sigma_LJ*(1+r_mHB)
    r_mLJ = sigma_LJ*(1+r_mLJ)

    v_s  =  (-exp(-AS*p*(r_mS**3)/T)*(T + AS*p*(r_mS**3)) + exp(-AS*p*(r_d**3)/T)*(T + AS*p*(r_d**3)))
    v_s  =  v_s /(p*(-exp(-AS*p*(r_mS**3)/T) + exp(-AS*p*(r_d**3)/T)))

    v_HB =  (-exp(-AHB*p*(r_mHB**3)/T)*(T + AHB*p*(r_mHB**3)) + exp(-AHB*p*(r_d**3)/T)*(T + AHB*p*(r_d**3)))
    v_HB =  v_HB /(p*(-exp(-AHB*p*(r_mHB**3)/T) + exp(-AHB*p*(r_d**3)/T)))

    v_LJ =  (-exp(-ALJ*p*(r_mLJ**3)/T)*(T + ALJ*p*(r_mLJ**3)) + exp(-ALJ*p*(r_d**3)/T)*(T + ALJ*p*(r_d**3)))
    v_LJ =  v_LJ /(p*(-exp(-ALJ*p*(r_mLJ**3)/T) + exp(-ALJ*p*(r_d**3)/T)))

    v_O  = v_LJ + T/p
    return v_s,v_HB,v_LJ,v_O



def vol(T,p,rl,rh,Aj):
    vj =  T/p + Aj*(rh^3.0-rl^3.0)
    return vj

