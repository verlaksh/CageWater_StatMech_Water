import scipy
import numpy as np
from pylab import *
from scipy.optimize import curve_fit
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)



SMALL_SIZE = 10
MEDIUM_SIZE = 16
BIGGER_SIZE = 22
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize

#extracting exmodel from colormap
cmap = cm.get_cmap('plasma', 25)    # PiYG
plasma=["" for x in range(25)]

for i in range(cmap.N):
        rgb = cmap(i)[:3] # will return rgba, we take only first 3 so we get rgb
        plasma[i]=matplotlib.colors.rgb2hex(rgb)

cmap = cm.get_cmap('Reds', 25)    # PiYG
Reds=["" for x in range(25)]

for i in range(cmap.N):
        rgb = cmap(i)[:3] # will return rgba, we take only first 3 so we get rgb
        Reds[i]=matplotlib.colors.rgb2hex(rgb)

cmap = cm.get_cmap('Blues', 25)    # PiYG
Blues=["" for x in range(25)]

for i in range(cmap.N):
        rgb = cmap(i)[:3] # will return rgba, we take only first 3 so we get rgb
        Blues[i]=matplotlib.colors.rgb2hex(rgb)


f = open("./PT_data",'w')
#fslope = open("./LLPT",'w') 
f.write("# max min Temp(K) of thermophysical properties vs pressure\n")
for p in range(-1400, -1000, 10):
    #if p == -1460:
    #    continue
    #else:
    #p = p*100
        print(p,"bar")
        model = np.loadtxt('./Press/real_units_%sbar'%p, comments='#')
        max_rho = model[np.argmax(model[:2100,4]),1]+273.15
        dtmin = 0
        min_rho = model[np.argmin(model[dtmin:1900,4])+dtmin,1]+273.15
        print(np.argmin(model[930:2330,4]))
        print(min_rho)
        max_kap = model[np.argmax(model[1570:2050,5])+1570,1]+273.15
        min_kap = model[np.argmin(model[930:2300,5])+930,1]+273.15
        max_alp = model[np.argmax(model[930:2330,6])+930,1]+273.15
        min_alp = model[np.argmin(model[1000:2330,6])+1000,1]+273.15
        max_Cp  = model[np.argmax(model[:2284,7]),1]+273.15
        min_Cp  = model[np.argmin(model[930:2200,7])+930,1]+273.15
        
        model = np.loadtxt('./Press/population_real_units_%sbar'%p, comments='#')
        fHB = np.round(model[:2230,3]/0.5,decimals=2)
        fHB = abs(model[:2230,3]-0.5)
        eq_K    = model[np.argmin(fHB),1]+273.15
        f.write("%7.2f  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f %7.2f\n"%(p*0.1,max_rho, min_rho, max_kap, min_kap, max_Cp, min_Cp, max_alp, min_alp, eq_K))

for p in range(-1000, 2600, 10):
    if p == 0:
        p=p+1
    #else:
    #p = p*100
    print(p,"bar")
    model = np.loadtxt('./Press/real_units_%sbar'%p, comments='#')
    max_rho = model[np.argmax(model[:2330,4]),1]+273.15
    dtmin = 0
    min_rho = model[np.argmin(model[dtmin:1900,4])+dtmin,1]+273.15
    print(np.argmin(model[930:2330,4]))
    print(min_rho)
    max_kap = model[np.argmax(model[:2050,5]),1]+273.15
    min_kap = model[np.argmin(model[930:2300,5])+930,1]+273.15
    max_alp = model[np.argmax(model[930:2330,6])+930,1]+273.15
    min_alp = model[np.argmin(model[:2330,6]),1]+273.15
    max_Cp  = model[np.argmax(model[:2284,7]),1]+273.15
    min_Cp  = model[np.argmin(model[930:2200,7])+930,1]+273.15

    model = np.loadtxt('./Press/population_real_units_%sbar'%p, comments='#')
    fHB = np.round(model[:2230,3]/0.5,decimals=2)
    fHB = abs(model[:2230,3]-0.5)
    eq_K    = model[np.argmin(fHB),1]+273.15
    f.write("%7.2f  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f %7.2f\n"%(p*0.1,max_rho, min_rho, max_kap, min_kap, max_Cp, min_Cp, max_alp, min_alp, eq_K))


fslope = open("./LLPT",'w')
fslope.write("# Slope of density and LLPT\n")
dt1=0
dt2=2
for p in range(1500, 2641, 10):
    print(p,"bar")
    model = np.loadtxt('./Press/real_units_%sbar'%p, comments='#')
    drdT    = (model[dt1+dt2:2100+dt2,4]-model[dt1:2100,4])
    #drdT    = (model[dt1+dt2:918+dt2,4]-model[dt1:918,4])
    drdT_ind= np.argmax(drdT)
    min_rho = model[drdT_ind+dt1,4]; print(min_rho)
    max_rho = model[drdT_ind+dt1+dt2,4]
    Temp_max= model[drdT_ind+dt1+dt2,1]+273.15
    Temp_min= model[drdT_ind+dt1,1]+273.15
    fslope.write("%7.2f  %7.2f  %7.2f  %7.2f %7.2f  %7.2f\n"%(p*0.1,Temp_min, min_rho, Temp_max, max_rho, (min_rho+max_rho)/2.))

fslope.close()
print("Done")
