import scipy
import numpy as np
import matplotlib.pyplot as plt; #plt.switch_backend('agg')
from pylab import *
from scipy.optimize import curve_fit
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)



SMALL_SIZE = 14
MEDIUM_SIZE = 16
BIGGER_SIZE = 22
#plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
#plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

#extracting exmodel from colormap
cmap = cm.get_cmap('plasma',11)    # PiYG
#cmap = plt.colormaps['plasma']

print (cmap)
hexmodel=["" for x in range(11)]
for i in range(cmap.N):
        rgb = cmap(i)[:3] # will return rgba, we take only first 3 so we get rgb
        hexmodel[i]=matplotlib.colors.rgb2hex(rgb)

cmap = cm.get_cmap('Reds', 11)    # PiYG
Reds=["" for x in range(11)]
for i in range(cmap.N):
        rgb = cmap(i)[:3] # will return rgba, we take only first 3 so we get rgb
        Reds[i]=matplotlib.colors.rgb2hex(rgb)

cmap = cm.get_cmap('Blues', 11)    # PiYG
Blues=["" for x in range(11)]
for i in range(cmap.N):
        rgb = cmap(i)[:3] # will return rgba, we take only first 3 so we get rgb
        Blues[i]=matplotlib.colors.rgb2hex(rgb)


alldata = np.loadtxt('./rho_all', comments='#')

sort = alldata[np.lexsort((alldata[:,0],alldata[:,1]))]

print(sort)

f = open("sorted_Mishima", "w")
f.write("#T(K) P(MPa) rho(g/cc) uncertainity\n")
for i in range(len(sort[:,0])):
    f.write("%10f %10f %10f %10f\n"%(sort[i,0],sort[i,1],sort[i,2],sort[i,3]))


