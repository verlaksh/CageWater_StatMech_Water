#Copyright (c) 2026 Lakshmanji Verma
# This is run script for calculating bulk properties of water as described in 
# Lakshmanji Verma and Ken A. Dill, Statistical mechanical theory of liquid water, JCTC, 2025
# It requires "water_functions_3D.py" and "default_parameters.py"
# To run : bulk_water_script.py 

import os ; import sys ; from math import *; import numpy as np
from scipy.special import erf ; from scipy.special import expi  
from water_functions_3D import *  
from default_parameters import *
import mpmath as mp


# INPUTS PARAMETERS
real_r     = r_HB_real 
real_e     = eps_HB_real
sigma_s    = 0.0   #dummy sol

### Parameters in reduced units to be read from command line if needed
# look at the default_parameters.py for parameter definitions
#to run: python bulk_water_script.py 14 parameter values

#real_e     = float(sys.argv[1])
#a          = float(sys.argv[2])  
#eps_LJ     = float(sys.argv[3])  
#eps_c      = float(sys.argv[4])  
#KT_HB      = float(sys.argv[5])
#KT_LJ      = float(sys.argv[6])
#Kx         = float(sys.argv[7])
#Kth        = float(sys.argv[8])  
#r_d        = float(sys.argv[9])  
#r_mS       = float(sys.argv[10]) 
#r_mHB      = float(sys.argv[11])
#r_mLJ      = float(sys.argv[12])
#xHB        = float(sys.argv[13])
#xLJ        = float(sys.argv[14])



# Thermodynamic parameters
T_factor   = (kB*NA)/(real_e*1000.0)
T_0        = T_real*T_factor                       # initial reduced temperature                                                                    
P_factor   = NA*(real_r*nm_to_m)**3/real_e/1000.0
p_0        = P_real*P_factor                       # initial reduced pressure
d_factor   = 18*0.001/(NA*(real_r*nm_to_m)**3) # Density factor

# Log file
flog= open("./logfile",'w')
flog.write ('Default inputs in reduced units\n')
flog.write ('Pure water and nonpolar parameters:\n')
flog.write ('T\t\t= %s\nP\t\t= %s\nNA\t\t= %s\nkB\t\t= %s\nnmtom\t\t= %s\n'%(T_0,p_0,NA,kB,nm_to_m))
flog.write ('vdW_a\t\t= %s\neps_HB\t\t= %s\neps_LJ\t\t= %s\n'%(a,eps_HB,eps_LJ))
flog.write ('eps_c\t\t= %s\nk_HB\t\t= %s\n'%(eps_c,Kth)) 
flog.write ('r_HB\t\t= %s\nsigma_LJ\t\t= %s\nr_d\t\t= %s\n'%(r_HB,sigma_LJ,r_d))
flog.write ('\nElectrostatics Parameters:\n')


flog.write ('\n\nRun parameters:\n')
P = [1]
#P = [4000, 3800, 3600, 3400, 3200, 3000, 2900, 2800, 2700, 2600, 2500, 2400, 2300, 2200, 2100, 2000, 1900, 1800, 1700, 1600, 1500, 1400, 1200, 1000, 800, 600, 500, 400, 220, 200, 150, 100, 50, 20, 10, 1,  -200, -400, -600, -800, -1000]


for p in P:
   fvarT= open("../Press/real_units_%sbar"%(p), 'w')
   fpopT= open("../Press/population_real_units_%sbar"%(p), 'w')
   fvarT.write('%10s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n'%("#sigma_s","Temp(C)","Press(Mpa)","p_meanfield","density (kg/m^3)","kappa","alpha","Cp(J/molK)","Cv(J/molK)","mu(J/mol)"))
   fpopT.write('%10s %10s %10s %10s %11s %10s %10s %10s %10s \n' %("#solute size.","Temp(C)","Press(Mpa)","fHB" ,"fS","fLJ","fO","xHB","nHB"))
   #for p in range(1,601):
   if p ==0:
      continue
   elif p == 1:
      p_0 = (p+0.01325)*(10**5)*P_factor
      # Creating files to save data
      f12    = open("../Atm_P/qh_qb",'w')
      f13    = open("../Atm_P/f_HB_S_LJ_O",'w')
      f14    = open("../Atm_P/delta_HB_S_LJ_O",'w')
      f1     = open("../Atm_P/molv_kappa", 'w')
      f2     = open("../Atm_P/G_H_S_Cp", 'w')
      f3     = open("../Atm_P/avg_nHB", 'w')
      fene   = open("../Atm_P/ene_contri", 'w')
      ftrans = open("../Atm_P/trans_contri", 'w')
      fpv    = open("../Atm_P/ph_vh", 'w')
      fcp    = open("../Atm_P/cp", 'w')
      fdudv  = open("../Atm_P/dudv", 'w')
      fcv    = open("../Atm_P/cv", 'w')
      fka    = open("../Atm_P/kappa", 'w')
      fal    = open("../Atm_P/alpha", 'w')
      fvar   = open("../Atm_P/real_units", 'w')
      fpop   = open("../Atm_P/population_real_units", 'w')
      fv     = open("../Atm_P/vol", 'w')
      favgv  = open("../Atm_P/avgv", 'w')

      f1.write( '%10s %10s %15s %15s \n'  %("#solute size.","Temp","Molar Volume","kappa"))
      f2.write( '%10s %10s %10s %10s %10s %10s %10s\n' %("#solute size.","Temp(C)","P(MPa)", "G", "H", "TdS ", "Cp"))
      f3.write( '%10s %10s %15s \n'  %("#solute size.","Temp", "avg_nHB"))
      f12.write('%10s %10s %10s %10s %10s %10s %10s\n'  %("#solute size.","Temp","q_h","q_b","avg_wat_shell","ph","v_h" ))
      f13.write('%10s %10s %6s  %11s %10s %10s \n' %("#solute size.","Temp","fHB" ,"fS","fLJ","fO"  ))
      f14.write('%10s %10s %10s %10s %10s %10s \n' %("#solute size.","Temp"," Delta_HB"," Delta_S"," Delta_LJ"," Delta_O"  ))
      fene.write('%10s %10s %10s  %10s %10s %10s %10s\n' %("#solute size.","Temp(C)","P(MPa)","av_uHB","av_uS","av_uLJ","avg_ene_bulk"))
      ftrans.write('%10s %10s %10s  %10s %10s %10s %10s\n' %("#solute size.","Temp(C)","P(MPa)","HB","S","LJ","avg_trans"))
      fpv.write('%10s %10s %10s %10s %10s %10s %10s %10s %10s %10s \n' %("#solute size","Temp","p_0","p_meanfield","v_h","kappa","alpha","g", "h", "Cp"))
      fcp.write('%10s %10s %10s %10s %10s %10s %10s %10s %10s %10s \n' %("#solute size","Temp","p_0","p_meanfield","Cp","meanC","cpHB","cpS", "cpLJ", "cpO"))
      fdudv.write('%10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s \n' %("#solute size","Temp","p_0","p_meanfield","duHB","duS", "duLJ", "duO", "dvHB","dvS", "dvLJ", "dvO"))
      fcv.write('%10s %10s %10s %10s %10s %10s %10s %10s %10s %10s \n' %("#solute size","Temp","p_0","p_meanfield","Cv","meanC","cvHB","cvS", "cvLJ", "cvO"))
      fka.write('%10s %10s %10s %10s %10s %10s %10s %10s %10s \n' %("#solute size","Temp","p_0","p_meanfield","ka","kaHB","kaS", "kaLJ", "kaO"))
      fal.write('%10s %10s %10s %10s %10s %10s %10s %10s %10s \n' %("#solute size","Temp","p_0","p_meanfield","al","alHB","alS", "alLJ", "alO"))
      fvar.write('%10s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n'%("#sigma_s","Temp(C)","Press(Mpa)","p_meanfield","density (kg/m^3)","kappa","alpha","Cp(J/molK)","Cv(J/molK)","mu(J/mol)"))
      fpop.write('%10s %10s %10s %10s  %11s %10s %10s %10s %10s \n' %("#solute size.","Temp(C)","Press(Mpa)","fHB" ,"fS","fLJ","fO","xHB","nHB"))
      fv.write('%10s %10s %10s %10s %10s %10s %10s %10s %10s %10s \n' %("#solute size","Temp","p_0","p_meanfield","v","vT(l/kg)","vHB","vS", "vLJ", "vO"))
    
   else:
      p_0 = p*(10**5)*P_factor
   
   ## Default loop for quick calculaitons 
   for t in range(100,400):
      print( "Temp=%s Press=%sMPa"%(t,0.1*p))
      T_0=(t+0.15)*R/real_e/1000.0

   ## for high definition data to find LLCP one need finely placed data
   #for t in range(1000,4001):
   #   print( "Temp=%s Press=%sMPa"%(0.1*t,0.1*p))
   #   T_0=(t*0.1)*R/real_e/1000.0

      flog.write('T_real\t\t= %s\tin K\nT*\t\t= %s\n'%(t,T_0))

      # Initial Calculation at P_0
      if p_0>0:
         p_i = p_0    #initial guess
      else :
         p_i = -0.1
      
      v_0,q_b,Q1,fHB,fS,fLJ,fO,v_HB,v_s,v_LJ,v_O,av_uHB,av_uS,av_uLJ,av_uO,avg_ene_bulk_water = bulk_wat(p_0,p_i,T_0,real_e,a,eps_LJ,eps_c,KT_HB,KT_LJ,Kx,Kth,r_d,r_mS,r_mHB,r_mLJ,xHB,xLJ)
      # Meanfield pressure estimation 
      p_meanfield=meanfield_pressure(p_0,p_i,T_0,v_0,real_e,a,eps_LJ,eps_c,KT_HB,KT_LJ,Kx,Kth,r_d,r_mS,r_mHB,r_mLJ,xHB,xLJ)

      # Calculations at meanfield p
      v,q_b,Q1,fHB,fS,fLJ,fO,v_HB,v_s,v_LJ,v_O,av_uHB,av_uS,av_uLJ,av_uO,avg_ene_bulk_water = bulk_wat(p_0,p_meanfield,T_0,real_e,a,eps_LJ,eps_c,KT_HB,KT_LJ,Kx,Kth,r_d,r_mS,r_mHB,r_mLJ,xHB,xLJ)

      # Compressibility, Exnpansivity, and heat capacity of bulk water
      ka0,kaHB,kaS,kaLJ,kaO = kap(p_0,p_meanfield,T_0,real_e,a,eps_LJ,eps_c,KT_HB,KT_LJ,Kx,Kth,r_d,r_mS,r_mHB,r_mLJ,xHB,xLJ)
      al0,alHB,alS,alLJ,alO = alp(p_0,p_meanfield,T_0,real_e,a,eps_LJ,eps_c,KT_HB,KT_LJ,Kx,Kth,r_d,r_mS,r_mHB,r_mLJ,xHB,xLJ)
      h,hHB,hO,hLJ,hO,Cp0,cpHB,cpS,cpLJ,cpO,cpm ,dufHB,dufS,dufLJ,dufO,dvfHB,dvfS,dvfLJ,dvfO = cp(p_0,p_meanfield,T_0,real_e,a,eps_LJ,eps_c,KT_HB,KT_LJ,Kx,Kth,r_d,r_mS,r_mHB,r_mLJ,xHB,xLJ,al0)
      u,uHB,uO,uLJ,uO,Cv,cvHB,cvS,cvLJ,cvO,cvm = cv(p_0,p_meanfield,T_0,real_e,a,eps_LJ,eps_c,KT_HB,KT_LJ,Kx,Kth,r_d,r_mS,r_mHB,r_mLJ,xHB,xLJ,al0)

      # Meanfield Correction and uniconversions
      # Compressibility
      ka   = ka0 /(1.0-2.0*a*ka0/v**2)*P_factor*10**5
      kaHB = kaHB/(1.0-2.0*a*ka0/v**2)*P_factor*10**5
      kaS  = kaS /(1.0-2.0*a*ka0/v**2)*P_factor*10**5
      kaLJ = kaLJ/(1.0-2.0*a*ka0/v**2)*P_factor*10**5
      kaO  = kaO /(1.0-2.0*a*ka0/v**2)*P_factor*10**5

      # Exnpansivity
      al    = al0 /(1.0-2.0*a*ka0/v**2)*T_factor
      alHB  = alHB/(1.0-2.0*a*ka0/v**2)*T_factor
      alS   = alS /(1.0-2.0*a*ka0/v**2)*T_factor
      alLJ  = alLJ/(1.0-2.0*a*ka0/v**2)*T_factor
      alO   = alO /(1.0-2.0*a*ka0/v**2)*T_factor

      # heat capacity
      Cp = Cp0 + 2*a*al0*ka0*(p_0-1*a/v**2)/v/(1-2.0*a*ka0/v**2)
      cpm = 2*a*al0*ka0*(p_0-1*a/v**2)/v/(1-2.0*a*ka0/v**2)
      meanC = cpm*real_e*T_factor*1000/(1-2.0*a*ka0/v**2)
      cpO = cpO*real_e*T_factor*1000/(1-2.0*a*ka0/v**2)
      cpLJ = cpLJ*real_e*T_factor*1000/(1-2.0*a*ka0/v**2)
      cpS = cpS*real_e*T_factor*1000/(1-2.0*a*ka0/v**2)
      cpHB = cpHB*real_e*T_factor*1000/(1-2.0*a*ka0/v**2)
      Cp = Cp*real_e*T_factor*1000
      
      # const v
      Cv = Cv #+ cvm)/(1-2.0*a*ka0/v**2)
      meanC = cvm*real_e*T_factor*1000/(1-2.0*a*ka0/v**2) 
      cvO  = cpO*real_e*T_factor*1000/(1-2.0*a*ka0/v**2)
      cvLJ = cpLJ*real_e*T_factor*1000/(1-2.0*a*ka0/v**2)
      cvS  = cpS*real_e*T_factor*1000/(1-2.0*a*ka0/v**2)
      cvHB = cpHB*real_e*T_factor*1000/(1-2.0*a*ka0/v**2)
      Cv   = Cv*real_e*T_factor*1000

      # Energetics of bulk water
      u = (u-a/v) *real_e
      h = (h-2*a/v) *real_e             
      mu= (-1*T_0*mp.log(q_b)) - 2*a/v#/T_factor
      g = mu *real_e

      av_uHB             = (av_uHB -a/v) *real_e
      av_uS              = (av_uS  -a/v) *real_e 
      av_uLJ             = (av_uLJ -a/v) *real_e
      avg_ene_bulk_water = (avg_ene_bulk_water -a/v) *real_e

      #volume in litre/kg
      vS  = 1000*(v_s*fS)/d_factor
      vHB = 1000*(v_HB*fHB)/d_factor
      vLJ = 1000*(v_LJ*fLJ)/d_factor
      vO  = 1000*(v_O*fO)/d_factor
      vT  = vS+vHB+vLJ+vO
      #v   = 1000*NA*((real_r*nm_to_m)**3)*v

      # Energy contributions

      trans_E_cage =  flex_bond(T_0,r_d,r_mS,Kx)[0]  *real_e 
      trans_E_HB   =  flex_bond(T_0,r_d,r_mHB,Kx)[0] *real_e 
      trans_E_LJ   =  flex_bond(T_0,r_d,r_mLJ,Kx)[0] *real_e 
      avg_trans = trans_E_cage*fS + trans_E_HB*fHB+trans_E_LJ*fLJ

      fvarT.write('%15.10f %10.3f %15f %10f %15.4f %15.10f %15.10f %15.10f %15.10f %15.10f \n'%(sigma_s/2.0,T_0/T_factor-273.15,p_0/P_factor/1000000,p_meanfield,d_factor/v,ka,al,Cp,Cv,mu))
      fpopT.write('%15.10f %10.3f %10.3f %10f  %10f %10f %10f %10f %10f \n' %(sigma_s/2.,T_0/T_factor-273.15,p_0/P_factor/1000000,fHB,fS,fLJ,fO,fHB+fS,4*(fHB+fS) )) 

      if p == 1:

         fka.write('%15.10f %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f \n' %(sigma_s,T_0/T_factor-273.15,p_0,p_meanfield,ka,kaHB,kaS,kaLJ,kaO))       
         fal.write('%15.10f %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f \n' %(sigma_s,T_0/T_factor-273.15,p_0,p_meanfield,al,alHB,alS,alLJ,alO))
         favgv.write('%10f %10f %10f %10f %10f %10f %10f %10f %10f %10f \n' %(sigma_s,T_0/T_factor-273.15,p_0,p_meanfield,v,v_HB+v_s+v_LJ+v_O,v_HB,v_s,v_LJ,v_O))
         fv.write('%10f %10f %10f %10f %10f %10f %10f %10f %10f %10f \n' %(sigma_s,T_0/T_factor-273.15,p_0,p_meanfield,v,vT,vHB,vS,vLJ,vO))
         fcp.write('%15.10f %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f \n' %(sigma_s,T_0/T_factor-273.15,p_0,p_meanfield,Cp,meanC,cpHB,cpS,cpLJ,cpO))
         fdudv.write('%15.10f %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f \n' %(sigma_s,T_0/T_factor-273.15,p_0,p_meanfield,dufHB,dufS,dufLJ,dufO,dvfHB,dvfS,dvfLJ,dvfO))
         fcv.write('%15.10f %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f \n' %(sigma_s,T_0/T_factor-273.15,p_0,p_meanfield,Cv,meanC,cvHB,cvS,cvLJ,cvO))
         fvar.write('%15.10f %10.3f %10.3f %10f %15.4f %15.10f %15.10f %15.10f %15.10f  %15.10f \n'%(sigma_s,T_0/T_factor-273.15,p_0/P_factor/1000000,p_meanfield,d_factor/v,ka,al,Cp,Cv,mu))
         fene.write('%15.10f %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f\n' %(sigma_s/2.,T_0/T_factor-273.15,p_0/P_factor/1000000,av_uHB,av_uS,av_uLJ,avg_ene_bulk_water))
         ftrans.write('%15.10f %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f\n' %(sigma_s/2.,T_0/T_factor-273.15,p_0/P_factor/1000000,trans_E_HB,trans_E_cage,trans_E_LJ,avg_trans))
         fpv.write('%15.10f %15.10f %s10 %10s %10s %10s %10s %10s %10s %10s\n' %(sigma_s/2.,T_0/T_factor-273.15,p_0,p_meanfield,v,ka,al,g,h,Cp))


         f14.write('%15.10f %15.10f %10f  %10f %10f %10f  \n' %(sigma_s/2.,T_0,fHB,fS,fLJ,fO ))
         f13.write('%15.10f %15.10f %10f  %10f %10f %10f \n' %(sigma_s/2.,T_0,fHB,fS,fLJ,fO ))
         fpop.write('%15.10f %10.3f %10.3f %10f  %10f %10f %10f %10f %10f \n' %(sigma_s/2.,T_0/T_factor-273.15,p_0/P_factor,fHB,fS,fLJ,fO,fHB+fS,4*(fHB+fS)))
         
         f1.write('%15.10f %15.10f %15f %15f\n' %(sigma_s/2.,T_0/T_factor-273.15,v,ka))
         f2.write('%10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n' %(sigma_s/2.,T_0/T_factor-273.15,p_0/P_factor/1000000, g, h, -1*(g-h), Cp))
   if p == 1:
      f14.close()
      f13.close()
      f12.close()
      f3.close()
      f2.close()
      f1.close()
      fpv.close()
      fka.close()
      fal.close()
      fcp.close()
      fene.close()
      fvar.close()
      fv.close() 
   fvarT.close()
   fpopT.close()

flog.close()

print("Done")
