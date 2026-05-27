import scipy
import numpy as np
from pylab import *
from scipy.optimize import curve_fit
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)



SMALL_SIZE = 11
MEDIUM_SIZE = 14
BIGGER_SIZE = 22
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=12)    # legend fontsize

#extracting exmodel from colormap
cmap = cm.get_cmap('plasma', 11)    # PiYG
Plasma=["" for x in range(11)]

for i in range(cmap.N):
        rgb = cmap(i)[:3] # will return rgba, we take only first 3 so we get rgb
        Plasma[i]=matplotlib.colors.rgb2hex(rgb)

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


factor = 17180.0/8.314



#Separate plots 
fig = plt.figure(figsize=(4,4))

ax =fig.add_axes([0.1, 0.14, 0.8, 0.78], frameon=False, xticklabels=[],yticklabels=[])
ax.tick_params(axis='both',direction = 'in',bottom=False, top=False, left=False, right=False)
ax1 = fig.add_axes([0.508, 0.69, 0.485, 0.305], xticklabels=[])#, ylim=(0, 0.6))
ax2 = fig.add_axes([0.18, 0.12, 0.813, 0.505])#, xticklabels=[])#, ylim=(0, 0.6))

model = np.loadtxt('./Atm_P/real_units', comments='#')
#Volume
num = range(0,121)
num = np.array(num)*0.01
xdat = np.zeros(len(num))

exp = np.loadtxt('density', comments='#')
tip3p=np.loadtxt('./MD_Model_data/tip3p_rho', comments='#')
ax1.plot(exp[::2,0],1/exp[::2,2],color='darkorange',linestyle='', marker='^', label=r'$\mathrm{Exp.}$',fillstyle='none',mew = 1.5)
ax1.plot(model[:,1],1000/model[:,4] ,color='k',linestyle='-',label=r'$\mathrm{Theory}$' , linewidth="2")
ax1.plot(xdat,num,color='k',linestyle=':',label=r'$\mathrm{0^oC}$' )
ax1.set_ylabel(r'$\mathrm{\mathbf{v}\/[cm^3/g]}$',{'fontsize':MEDIUM_SIZE})
ax1.yaxis.set_label_coords(-0.9, -0.2)
ax1.set_yticks(np.arange(1.0, 1.047, 0.02))
ax1.set_ylim([0.99,1.047])
ax1.set_xlim([-50, 100])
ax1.set_xticks(np.arange(-50, 100, 50))
ax1.tick_params(axis='both',direction = 'in',bottom=1, top=1, left=1, right=1)
ax1.legend( frameon=False, fontsize = "9", ncol =1, loc = "center", bbox_to_anchor=(-0.55, 0.7))

#micro
model = np.loadtxt('./Atm_P/vol', comments='#')
ax2.plot(model[:,1],model[:,5],color='k',linestyle='-',label=r'$\mathrm{Total}$' , linewidth="2")
ax2.plot(model[:,1],model[:,6], linewidth="2",color=Reds[10],linestyle='-',label=r'$\mathrm{pHB}$' )
ax2.plot(model[:,1],model[:,7], linewidth="2",color=Reds[8],linestyle='-',label=r'$\mathrm{Cage}$' )
ax2.plot(model[:,1],model[:,8], linewidth="2",color=Reds[5],linestyle='-',label=r'$\mathrm{vdW}$' )
ax2.plot(model[:,1],model[:,9], linewidth="2",color=Reds[3],linestyle='-',label=r'$\mathrm{O}$' )
ax2.plot(xdat,num,color='k',linestyle=':')#,label=r'$\mathrm{0^oC}$' )
ax2.yaxis.set_label_coords(-0.18, -0.1)
ax2.set_xlabel(r'$\mathrm{\mathbf{T}\ [{^\circ} C]}$',{'fontsize':MEDIUM_SIZE})
ax2.set_ylim([0,1.2])
ax2.set_yticks(np.arange(0,1.2,0.5))
ax2.set_xlim([-150, 100])
ax2.set_xticks(np.arange(-50, 100, 50))
ax2.tick_params(axis='both',direction = 'in', which= 'both',bottom=1, top=1, left=1,right=1)
ax2.legend( frameon=False, fontsize = "9", ncol=1, loc = "center", bbox_to_anchor=(0.15, 0.5))

plt.savefig("Microscopic_volume.png", dpi =600)
plt.show()
plt.close()

fig = plt.figure(figsize=(3,4))

ax =fig.add_axes([0.1, 0.14, 0.8, 0.78], frameon=False, xticklabels=[],yticklabels=[])
ax.tick_params(axis='both',direction = 'in',bottom=False, top=False, left=False, right=False)


ax1 = fig.add_axes([0.19, 0.55, 0.803, 0.445], xticklabels=[])#, ylim=(0, 0.6))
ax2 = fig.add_axes([0.19, 0.11, 0.803, 0.435])#, xticklabels=[])#, ylim=(0, 0.6))

model = np.loadtxt('./Atm_P/real_units', comments='#')
#Volume
num = range(0,121)
num = np.array(num)*0.01
xdat = np.zeros(len(num))+4

exp = np.loadtxt('density', comments='#')
tip4p=np.loadtxt('./MD_Model_data/tip4p2005_rho', comments='#')
mbpol=np.loadtxt('./Expt_Data/Palos_2024_JCTC_MBpol', comments='#')
mbdnn=np.loadtxt('./Expt_Data/Palos_2024_JCTC_DNN', comments='#')

print(tip4p)
ax1.plot(exp[::2,0],1/exp[::2,1],color='darkorange',linestyle='', marker='^', label=r'$\mathrm{Exp.}$',fillstyle='none',mew = 1.5)
ax1.plot(tip4p[:,0]-273.15,tip4p[:,1],color='green',linestyle='', marker='o', label=r'$\mathrm{TIP4P/2005}$',fillstyle='none',mew = 1.5)
ax1.plot(mbpol[:19,0]-273.15,mbpol[:19,1],color='b',linestyle='', marker='o', label=r'$\mathrm{MB-pol}$',fillstyle='none',mew = 1.5)
ax1.plot(model[:,1],model[:,4]*0.001,color='k',linestyle='-',label=r'$\mathrm{Theory}$' , linewidth="2")
ax1.plot(xdat,num,color='k',linestyle=':')#,label=r'$\mathrm{4^oC}$' )
ax1.set_ylabel(r'$\mathrm{\mathbf{\rho}\/[g cm^3]}$',{'fontsize':MEDIUM_SIZE})
ax1.yaxis.set_label_coords(-0.13, 0.5)
ax1.set_yticks(np.arange(0.9, 1.1, 0.1))
ax1.set_ylim([0.9,1.047])
ax1.set_xlim([-70, 100])
ax1.set_xticks(np.arange(-50, 100, 50))
ax1.tick_params(axis='both',direction = 'in',bottom=1, top=1, left=1, right=1)
ax1.legend( frameon=False, fontsize = "10", ncol =1, labelspacing=0.05, loc = "center", bbox_to_anchor=(0.55, 0.3))

model = np.loadtxt('./Atm_P/population_real_units', comments='#')

ax2.plot(model[:,1],model[:,4],color=Reds[10],linestyle='-',linewidth="2",label=r'$\mathrm{Cage}$' )
ax2.plot(model[:,1],model[:,3],color=Reds[8],linestyle='-',linewidth="2",label=r'$\mathrm{pHB}$' )
ax2.plot(model[:,1],model[:,5],color=Reds[5],linestyle='-',linewidth="2",label=r'$\mathrm{vdW}$' )
ax2.plot(xdat,num,color='k',linestyle=':')
ax2.set_ylabel(r'$\mathrm{\mathbf{f_{j}}}$',{'fontsize':MEDIUM_SIZE})
ax2.yaxis.set_label_coords(-0.13, 0.5)
ax2.set_xlabel(r'$\mathrm{\mathbf{T}\ [{^\circ} C]}$',{'fontsize':MEDIUM_SIZE})
ax2.set_yticks(np.arange(0, 1, 0.5))
ax2.set_ylim([0, 1])
ax2.set_xlim([-70, 100])
ax2.set_xticks(np.arange(-50, 100, 50))
ax2.tick_params(axis='both',direction = 'in', which= 'both',bottom=1, top=1, left=1,right=1)
ax2.legend( frameon=False,fontsize = "10", labelspacing=0.05)

plt.savefig("Microscopic_Density.png", dpi =600)
plt.show()
plt.close()


#Alpha

fig = plt.figure(figsize=(3,4))
ax =fig.add_axes([0.1, 0.14, 0.8, 0.78], frameon=False, xticklabels=[],yticklabels=[])
ax.tick_params(axis='both',direction = 'in',bottom=False, top=False, left=False, right=False)


ax1 = fig.add_axes([0.19, 0.55, 0.803, 0.445], xticklabels=[])#, ylim=(0, 0.6))
ax2 = fig.add_axes([0.19, 0.11, 0.803, 0.435])#, xticklabels=[])#, ylim=(0, 0.6))
ax22 = fig.add_axes([0.69, 0.16, 0.30, 0.185])

num = range(-200,200)
num = np.array(num)
xdat = np.zeros(len(num))+4

model = np.loadtxt('./Atm_P/real_units', comments='#')
tip4p=np.loadtxt('./MD_Model_data/tip4p2005_alp', comments='#')
reddy = np.loadtxt('./Expt_Data/Reddy_JCP_2016', comments='#')
ax1.plot(exp[::2,0],exp[::2,3]*0.001,color='darkorange',linestyle='', marker='^', label=r'$\mathrm{Exp.}$',fillstyle='none',mew = 1.5)
ax1.plot(tip4p[:,0]-273.15,tip4p[:,1]*0.01,color='green',linestyle='', marker='o', label=r'$\mathrm{TIP4P/2005}$',fillstyle='none',mew = 1.5)
ax1.plot(reddy[:,0]-273.15,reddy[:,1]*0.1,color='b',linestyle='', marker='o', label=r'$\mathrm{MB-pol}$',fillstyle='none',mew = 1.5)
ax1.plot(model[:,1],model[:,6]*10**3,color='k',linestyle='-',label=r'$\mathrm{Theory}$' , linewidth="2")
ax1.plot(model[:,1],model[:,1]*0,color='indigo',linestyle='--')#,label=r'$\mathrm{0^oC}$' )
ax1.plot(xdat,num,color='k',linestyle=':')#,label=r'$\mathrm{0^oC}$' )
ax1.set_ylabel(r'$\mathrm{\mathbf{\alpha_P} \times 10^3\/[K^{-1}]}$',{'fontsize':MEDIUM_SIZE})
ax1.yaxis.set_label_coords(-0.13, -0.05)
ax1.set_yticks(np.arange(-2, 2, 1))
ax1.set_ylim([-3,1.2])
ax1.set_xlim([-70, 100])
ax1.set_xticks(np.arange(-50, 100, 50))
ax1.tick_params(axis='both',direction = 'in',bottom=1, top=1, left=1, right=1)
ax1.legend( frameon=False, fontsize = "10", ncol =1, labelspacing=0.05)#, loc = "center", bbox_to_anchor=(-0.55, 0.7))

#micro
model = np.loadtxt('./Atm_P/alpha', comments='#')

ax2.plot(model[:,1],model[:,6]*10**3, linewidth="2",color=Reds[10],linestyle='-',label=r'$\mathrm{{Cage}}$' )
ax2.plot(model[:,1],model[:,5]*10**3, linewidth="2",color=Reds[8],linestyle='-',label=r'$\mathrm{{pHB}}$' )
ax2.plot(model[:,1],model[:,7]*10**3, linewidth="2",color=Reds[5],linestyle='-',label=r'$\mathrm{{vdW}}$' )
ax2.plot(xdat,num,color='k',linestyle=':')
ax2.set_xlabel(r'$\mathrm{\mathbf{T}\ [{^\circ} C]}$',{'fontsize':MEDIUM_SIZE})
ax2.set_ylim([-30,20])
ax2.set_yticks(np.arange(-20, 15, 10))
ax2.set_xlim([-70, 100])
ax2.set_xticks(np.arange(-50, 100, 50))
ax2.tick_params(axis='both',direction = 'in', which= 'both',bottom=1, top=1, left=1,right=1)
ax2.legend( frameon=False, fontsize = "10", ncol=1, loc = "upper right", labelspacing=0.05)#, bbox_to_anchor=(0.23, 0.4))

ax22.plot(model[:,1],model[:,6]*10**3, linewidth="1.5",color=Reds[10],linestyle='-',label=r'$\mathrm{{Cage}}$' )
ax22.plot(model[:,1],model[:,5]*10**3, linewidth="1.5",color=Reds[8],linestyle='-',label=r'$\mathrm{{pHB}}$' )
ax22.plot(model[:,1],model[:,7]*10**3, linewidth="1.5",color=Reds[5],linestyle='-',label=r'$\mathrm{{vdW}}$' )
ax22.yaxis.set_label_coords(-0.18, -0.1)
ax22.set_xlim([20, 75])
ax22.set_xticks(np.arange(25, 70, 20))
ax22.set_ylim([-2,2])
ax22.set_yticks(np.arange(-1, 2, 2))
ax22.tick_params(axis='both',direction = 'in', which= 'both',bottom=1, top=1, left=1,right=1)

plt.savefig("Microscopic_alpha.png", dpi =600)
plt.show()
plt.close()


#kappa
fig = plt.figure(figsize=(3,4))
ax =fig.add_axes([0.1, 0.14, 0.8, 0.78], frameon=False, xticklabels=[],yticklabels=[])
ax.tick_params(axis='both',direction = 'in',bottom=False, top=False, left=False, right=False)

ax1 = fig.add_axes([0.19, 0.55, 0.803, 0.445], xticklabels=[])#, ylim=(0, 0.6))
ax2 = fig.add_axes([0.19, 0.11, 0.803, 0.435])#, xticklabels=[])#, ylim=(0, 0.6))
ax22 = fig.add_axes([0.69, 0.37, 0.30, 0.175])

num = range(-60,150)
num = 0.1* np.array(num)
xdat = np.zeros(len(num))

model = np.loadtxt('./Atm_P/real_units', comments='#')
tip4p=np.loadtxt('./MD_Model_data/tip4p2005_kap', comments='#')
mbpol=np.loadtxt('./Expt_Data/Palos_2024_JCTC_MBpol', comments='#')
mbdnn=np.loadtxt('./Expt_Data/Palos_2024_JCTC_DNN', comments='#')
kim = np.loadtxt('./Expt_Data/Kim_Science_2017', comments='#')

ax1.plot(exp[::2,0],exp[::2,4]*0.1,color='darkorange',linestyle='', marker='^', label=r'$\mathrm{Exp.}$',fillstyle='none',mew = 1.5)
ax1.plot(kim[:,0]-273.15,kim[:,1]*0.1,color='darkorange',linestyle='', marker='^',fillstyle='none',mew = 1.5)
ax1.plot(tip4p[:,0]-273.15,tip4p[:,1]*10**5,color='green',linestyle='', marker='o', label=r'$\mathrm{TIP4P/2005}$',fillstyle='none',mew = 1.5)
ax1.plot(mbpol[20:,0]-273.15,mbpol[20:,1]*0.1,color='b',linestyle='', marker='o', label=r'$\mathrm{MB-pol}$',fillstyle='none',mew = 1.5)
ax1.plot(model[:,1],model[:,5]*10**5,color='k',linestyle='-',label=r'$\mathrm{Theory}$' , linewidth="2")
ax1.plot(xdat,num,color='k',linestyle=':')#,label=r'$\mathrm{0^oC}$' )
ax1.set_ylabel(r'$\mathrm{\mathbf{\kappa_T} \times 10^5\/[bar^{-1}]}$',{'fontsize':MEDIUM_SIZE})
ax1.yaxis.set_label_coords(-0.13, -0.05)
ax1.set_yticks(np.arange(4, 12.5, 4))
ax1.set_ylim([3,12.5])
ax1.set_xlim([-70, 100])
ax1.set_xticks(np.arange(-50, 100, 50))
ax1.tick_params(axis='both',direction = 'in',bottom=1, top=1, left=1, right=1)
ax1.legend( frameon=False, fontsize = "10", ncol =1, labelspacing=0.05)#, loc = "center", bbox_to_anchor=(0.83, 0.23))

#micro
num = range(-60,150)
num = np.array(num)
xdat = np.zeros(len(num))

model = np.loadtxt('./Atm_P/kappa', comments='#')
ax2.plot(model[:,1],model[:,6]*10**5, linewidth="2",color=Reds[10],linestyle='-',label=r'$\mathrm{{Cage}}$' )
ax2.plot(model[:,1],model[:,5]*10**5, linewidth="2",color=Reds[8],linestyle='-',label=r'$\mathrm{{pHB}}$' )
ax2.plot(model[:,1],model[:,7]*10**5, linewidth="2",color=Reds[5],linestyle='-',label=r'$\mathrm{{vdW}}$' )
ax2.plot(xdat,num,color='k',linestyle=':')
ax2.set_xlabel(r'$\mathrm{\mathbf{T}\ [{^\circ} C]}$',{'fontsize':MEDIUM_SIZE})
ax2.set_xlim([-70, 100])
ax2.set_xticks(np.arange(-50, 100, 50))
ax2.set_ylim([-50,60])
ax2.set_yticks(np.arange(-40, 60, 40))
ax2.tick_params(axis='both',direction = 'in', which= 'both',bottom=1, top=1, left=1,right=1)
ax2.legend( frameon=False, fontsize = "10", ncol=1, labelspacing=0.05, loc = "lower right")#, bbox_to_anchor=(0.83, 0.23))

ax22.plot(model[:,1],model[:,6]*10**5, linewidth="1.5",color=Reds[10],linestyle='-',label=r'$\mathrm{{Cage}}$' )
ax22.plot(model[:,1],model[:,5]*10**5, linewidth="1.5",color=Reds[8],linestyle='-',label=r'$\mathrm{{pHB}}$' )
ax22.plot(model[:,1],model[:,7]*10**5, linewidth="1.5",color=Reds[5],linestyle='-',label=r'$\mathrm{{vdW}}$' )
ax22.yaxis.set_label_coords(-0.18, -0.1)
ax22.set_xlim([30, 60])
ax22.set_xticks(np.arange(35, 60, 10))
ax22.set_ylim([-3,6])
ax22.set_yticks(np.arange(0, 10, 5))
ax22.tick_params(axis='both',direction = 'in', which= 'both',bottom=1, top=1, left=1,right=1)

plt.savefig("Microscopic_kappa.png", dpi =600)
plt.show()
plt.close()


# Heat Capacity Cp

fig = plt.figure(figsize=(3,4))
ax =fig.add_axes([0.1, 0.14, 0.8, 0.78], frameon=False, xticklabels=[],yticklabels=[])
ax.tick_params(axis='both',direction = 'in',bottom=False, top=False, left=False, right=False)

ax1 = fig.add_axes([0.19, 0.55, 0.803, 0.445], xticklabels=[])#, ylim=(0, 0.6))
ax2 = fig.add_axes([0.19, 0.11, 0.803, 0.435])#, xticklabels=[])#, ylim=(0, 0.6))
ax22 = fig.add_axes([0.7, 0.16, 0.29, 0.165])#, xticklabels=[])#, ylim=(0, 0.6))

num = range(-500,500)
num = np.array(num)
xdat = np.zeros(len(num))+4
mbpol = np.loadtxt('./Expt_Data/Gartner_JPCL_2022', comments='#')
model = np.loadtxt('./Atm_P/real_units', comments='#')
tip4p=np.loadtxt('./MD_Model_data/tip4p2005_cp', comments='#')
exp = np.loadtxt('exp_cp', comments='#')
pathak= np.loadtxt('./Expt_Data/Pathak_PNAS_2021', comments='#')
ax1.plot(exp[:,0],exp[:,1]/18,color='darkorange',linestyle='', marker='^', label=r'$\mathrm{Exp.}$',fillstyle='none',mew = 1.5)
ax1.plot(pathak[:,0]-273.15,pathak[:,1]/18,color='darkorange',linestyle='', marker='^',fillstyle='none',mew = 1.5)
ax1.plot(tip4p[:,0]-273.15,tip4p[:,1]*4.18/18,color='green',linestyle='', marker='o', label=r'$\mathrm{TIP4P/2005}$',fillstyle='none',mew = 1.5)
ax1.plot(mbpol[:,0]-273.15,mbpol[:,1]/18,color='b',linestyle='', marker='o', label=r'$\mathrm{MB-pol}$',fillstyle='none',mew = 1.5)
ax1.plot(model[:,1],model[:,7]/18,color='k',linestyle='-',label=r'$\mathrm{Theory}$' , linewidth="2")
ax1.plot(xdat,num,color='k',linestyle=':')#,label=r'$\mathrm{0^oC}$' )
ax1.set_ylabel(r'$\mathrm{\mathbf{C_{p}} \/[J K^{-1} g^{-1}]}$',{'fontsize':MEDIUM_SIZE})
ax1.yaxis.set_label_coords(-0.12, -0.05)
ax1.set_ylim([3,13])
ax1.set_yticks(np.arange(4,15,4))
ax1.set_xlim([-70, 100])
ax1.set_xticks(np.arange(-50, 100, 50))
ax1.tick_params(axis='both',direction = 'in',bottom=1, top=1, left=1, right=1)
ax1.legend( frameon=False, fontsize = "10", ncol =1, labelspacing=0.05)#, loc = "upper right")#, bbox_to_anchor=(0.6, 0.7))

#micro
model = np.loadtxt('./Atm_P/cp', comments='#')
ax2.plot(model[:,1],model[:,7]/18, linewidth="2", color=Reds[10],linestyle='-',label=r'$\mathrm{Cage}$' )
ax2.plot(model[:,1],model[:,6]/18, linewidth="2", color=Reds[8],linestyle='-',label=r'$\mathrm{pHB}$' )
ax2.plot(model[:,1],model[:,8]/18, linewidth="2", color=Reds[5],linestyle='-',label=r'$\mathrm{vdW}$' )
ax2.plot(xdat,num,color='k',linestyle=':')
ax2.set_xlabel(r'$\mathrm{\mathbf{T}\ [{^\circ} C]}$',{'fontsize':MEDIUM_SIZE})
ax2.set_ylim([-30,30])
ax2.set_yticks(np.arange(-25,30,25))
ax2.set_xlim([-70, 100])
ax2.set_xticks(np.arange(-50, 100, 50))
ax2.tick_params(axis='both',direction = 'in', which= 'both',bottom=1, top=1, left=1,right=1)
ax2.legend( frameon=False, fontsize = "10", ncol=1, loc = "upper right",labelspacing=0.05)#, bbox_to_anchor=(0.16, 0.66))

ax22.plot(model[:,1],model[:,7]*0.1, linewidth="1.5",color=Reds[10],linestyle='-',label=r'$\mathrm{{Cage}}$' )
ax22.plot(model[:,1],model[:,6]*0.1, linewidth="1.5",color=Reds[8],linestyle='-',label=r'$\mathrm{{pHB}}$' )
ax22.plot(model[:,1],model[:,8]*0.1, linewidth="1.5",color=Reds[5],linestyle='-',label=r'$\mathrm{{vdW}}$' )
ax22.set_xlim([10, 30])
ax22.set_xticks(np.arange(15, 35, 10))
ax22.set_ylim([0,6])
ax22.set_yticks(np.arange(0, 6, 5))
ax22.tick_params(axis='both',direction = 'in', which= 'both',bottom=1, top=1, left=1,right=1)

plt.savefig("Microscopic_Cp.png", dpi =600)
plt.show()
plt.close()


MEDIUM_SIZE = 12
num = range(-500,500)
num = np.array(num)
xdat = np.zeros(len(num))

fig = plt.figure(figsize=(6,4))
ax1 = fig.add_axes([0.10, 0.704, 0.395, 0.285], xticklabels=[])#, ylim=(0, 0.6))
ax2 = fig.add_axes([0.10, 0.412, 0.395, 0.285], xticklabels=[])#, ylim=(0, 0.6))
ax3 = fig.add_axes([0.10, 0.120, 0.395, 0.285])#, xticklabels=[])#, ylim=(0, 0.6))
ax4 = fig.add_axes([0.60, 0.560, 0.395, 0.43], xticklabels=[])#, ylim=(0, 0.6))
ax5 = fig.add_axes([0.60, 0.120, 0.395, 0.43])#, xticklabels=[])#, ylim=(0, 0.6))


model = np.loadtxt('./Atm_P/real_units', comments='#')
exp = np.loadtxt('exp_cp', comments='#')
ax1.plot(exp[:,0],exp[:,1]/18,color='darkorange',linestyle='', marker='^', label=r'$\mathrm{Bulk}$',fillstyle='none',mew = 1.5)

cnf = np.loadtxt('Expt_Data/cp_data_1', comments='#')
cnf = np.loadtxt('Expt_Data/cp_data_2', comments='#')

pathak= np.loadtxt('./Expt_Data/Pathak_PNAS_2021', comments='#')
ax1.plot(pathak[:,0]-273.15,pathak[:,1]/18,color='darkorange',linestyle='', marker='^',fillstyle='none',mew = 1.5)

ax1.plot(model[:,1],model[:,7]/18,color='k',linestyle='-',label=r'$\mathrm{Theory}$' , linewidth="2.5")
ax1.plot(xdat-44,num/18,color='k',linestyle=':')
ax1.plot(xdat-38,num/18,color='k',linestyle=':')
ax1.set_ylabel(r'$\mathrm{\mathbf{C_{p}}\/[Jg^{-1} K^{-1}]}$',{'fontsize':MEDIUM_SIZE})
ax1.yaxis.set_label_coords(-0.14, 0.5)
ax1.set_ylim([3,6])
ax1.set_yticks(np.arange(4,6,1))
ax1.set_xlim([-70, 50])
ax1.set_xticks(np.arange(-50, 50, 25))
ax1.tick_params(axis='both',direction = 'in',bottom=1, top=1, left=1, right=1)

exp = np.loadtxt('density', comments='#')
kim = np.loadtxt('./Expt_Data/Kim_Science_2017', comments='#')

ax2.plot(kim[:,0]-273.15,kim[:,1]*0.1,color='darkorange',linestyle='', marker='^',fillstyle='none',mew = 1.5)
ax2.plot(exp[::2,0],exp[::2,4]*0.1,color='darkorange',linestyle='', marker='^', label=r'$\mathrm{Exp.}$',fillstyle='none',mew = 1.5)
ax2.plot(model[:,1],model[:,5]*10**5,color='k',linestyle='-',label=r'$\mathrm{Theory}$' , linewidth="2.5")
ax2.plot(xdat-44,num,color='k',linestyle=':')
ax2.plot(xdat-38,num,color='k',linestyle=':')
ax2.set_ylabel(r'$\mathrm{\mathbf{\kappa_T} \ 10^5\/[bar^{-1}]}$',{'fontsize':MEDIUM_SIZE})
ax2.yaxis.set_label_coords(-0.14, 0.5)
ax2.set_yticks(np.arange(4, 12, 4))
ax2.set_ylim([3,11])
ax2.set_xlim([-70, 50])
ax2.set_xticks(np.arange(-50, 50, 25))
ax2.tick_params(axis='both',direction = 'in',bottom=1, top=1, left=1, right=1)

exp = np.loadtxt('density', comments='#')
cnf = np.loadtxt('Expt_Data/confined_alpha', comments='#')
model = np.loadtxt('./Atm_P/real_units', comments='#')
ax3.plot(exp[::2,0],exp[::2,3]*0.001,color='darkorange',linestyle='', marker='^', label=r'$\mathrm{Bulk\ Exp.}$',fillstyle='none',mew = 1.5)
ax3.plot(cnf[:16,0]-273.15,cnf[:16,1],color='darkturquoise',linestyle='', marker='o', label=r'$\mathrm{Conf.\ Exp.}$',fillstyle='none',mew = 1.5)
ax3.plot(model[:,1],model[:,6]*10**3,color='k',linestyle='-',label=r'$\mathrm{Theory}$' , linewidth="2.5")
ax3.plot(xdat-44,num,color='k',linestyle=':')
ax3.plot(xdat-38,num,color='k',linestyle=':')
ax3.set_ylabel(r'$\mathrm{\mathbf{\alpha_p} \ 10^4\/[K^{-1}]}$',{'fontsize':MEDIUM_SIZE})
ax3.set_xlabel(r'$\mathrm{\mathbf{T}\ [{^\circ} C]}$',{'fontsize':MEDIUM_SIZE})
ax3.yaxis.set_label_coords(-0.14, 0.4)
ax3.set_yticks(np.arange(-2, 2, 1))
ax3.set_ylim([-2.5,1.5])
ax3.set_xlim([-70, 50])
ax3.set_xticks(np.arange(-50, 50, 25))
ax3.tick_params(axis='both',direction = 'in',bottom=1, top=1, left=1, right=1)
ax3.legend( frameon=False, fontsize = "9", ncol =1, labelspacing=0.01, loc = "center", bbox_to_anchor=(0.75, 0.25))

model = np.loadtxt('./Atm_P/population_real_units', comments='#')
exp = np.loadtxt('./Expt_Data_plot/Xu_NatPhys_2009', comments='#')
ldl = np.loadtxt('./Expt_Data/Tyburski_MolPhys_2024', comments='#')

ax4.plot(model[:,1],model[:,4]+model[:,5],color=Reds[7],linestyle='-',linewidth="2",label=r'$\mathrm{LDL_{Theory}}$')
ax4.plot(model[:,1],model[:,3],color=Reds[10],linestyle='-',linewidth="2")
ax4.plot(exp[:9,0]-273.15,exp[:9,1],color=Reds[10],linestyle=' ', marker='o',fillstyle='none',mew = 1.5)
ax4.plot(ldl[:,0]-273.15,1-ldl[:,1],color=Reds[10],linestyle=' ', marker='^',fillstyle='none',mew = 1.5)
ax4.plot(ldl[:,0]-273.15,ldl[:,1],color=Reds[7],linestyle=' ', marker='v',fillstyle='none',mew = 1.5, label=r'$\mathrm{LDL_{Exp.}}$')
ax4.plot(exp[:9,2]-273.15,exp[:9,3],color=Reds[7],linestyle=' ', marker='o',fillstyle='none',mew = 1.5,label=r'$\mathrm{LDL_{Conf.\ Exp.}}$')
ax4.plot(xdat-44,num,color='k',linestyle=':')
ax4.plot(xdat-38,num,color='k',linestyle=':')
ax4.set_ylabel(r'$\mathrm{\mathbf{f_{j}}}$',{'fontsize':MEDIUM_SIZE})
ax4.yaxis.set_label_coords(-0.13, 0.5)
ax4.set_yticks(np.arange(0, 1.2, 0.5))
ax4.set_ylim([-0.02, 1.02])
ax4.set_xlim([-70, 50])
ax4.set_xticks(np.arange(-50, 50, 25))
ax4.tick_params(axis='both',direction = 'in', which= 'both',bottom=1, top=1, left=1,right=1)
ax4.legend( frameon=False, fontsize = "9", ncol =1, loc = "center", bbox_to_anchor=(0.77, 0.5),labelspacing=0.05,handlelength=1)


model = np.loadtxt('./Atm_P/population_real_units', comments='#')
R = 1.380649*6.022
HDA = log(model[:,3])
LDA = log(model[:,4]+model[:,5])

ax5.plot(model[:,1],LDA-HDA,color="k",linestyle='-',linewidth="2")#,label=r'$\mathrm{}$')
ax5.plot(xdat-44,num,color='k',linestyle=':')
ax5.plot(xdat-38,num,color='k',linestyle=':')
ax5.plot(model[:,1],model[:,1]*0,color='indigo',linestyle='--',linewidth="2",label=r'$\mathrm{0^o C}$')
ax5.set_ylabel(r'$\mathrm{\mathbf{ln K}}$',{'fontsize':MEDIUM_SIZE})
ax5.yaxis.set_label_coords(-0.13, 0.5)
ax5.set_xlabel(r'$\mathrm{\mathbf{T}\ [{^\circ} C]}$',{'fontsize':MEDIUM_SIZE})
ax5.set_yticks(np.arange(-2, 3, 2))
ax5.set_ylim([-3, 3])
ax5.set_xlim([-70, 50])
ax5.set_xticks(np.arange(-50, 50, 25))
ax5.tick_params(axis='both',direction = 'in', which= 'both',bottom=1, top=1, left=1,right=1)

plt.savefig("LLPT_lnK.png", dpi =600)
plt.show()
plt.close()


fig = plt.figure(figsize=(4,3))
ax =fig.add_axes([0.16, 0.14, 0.835, 0.845]) #, frameon=False, xticklabels=[],yticklabels=[])

model = np.loadtxt('./Atm_P/population_real_units', comments='#')
exp = np.loadtxt('./Expt_Data_plot/Xu_NatPhys_2009', comments='#')

ax.plot(model[:,1],model[:,4]+model[:,5],color='b',linestyle='-',linewidth="2")#,label=r'$\mathrm{Theory\ LDL}$' )
ax.plot(model[:,1],model[:,3],color='g',linestyle='-',linewidth="2")#,label=r'$\mathrm{Theory\ HDL}$' )
ax.plot(exp[:,0]-273.15,exp[:,1],color='g',linestyle=' ', marker='^', label=r'$\mathrm{IR \ HDL}$',fillstyle='none',mew = 1.5)
ax.plot(exp[:,2]-273.15,exp[:,3],color='b',linestyle=' ', marker='v', label=r'$\mathrm{IR \ LDL}$',fillstyle='none',mew = 1.5)
ax.plot(xdat-48,num,color='k',linestyle=':')
ax.plot(xdat-41,num,color='k',linestyle=':')
ax.set_ylabel(r'$\mathrm{\mathbf{x_{populations}}}$',{'fontsize':MEDIUM_SIZE})
ax.set_xlabel(r'$\mathrm{\mathbf{T}\ [{^\circ} C]}$',{'fontsize':MEDIUM_SIZE})
ax.yaxis.set_label_coords(-0.12, 0.5)
ax.set_yticks(np.arange(0, 1.2, 0.5))
ax.set_ylim([-0.02, 1.02])
ax.set_xlim([-100, 80])
ax.set_xticks(np.arange(-50, 100, 50))
ax.tick_params(axis='both',direction = 'in', which= 'both',bottom=1, top=1, left=1,right=1)
ax.legend( frameon=False, fontsize = "12", ncol =1, loc = "center", bbox_to_anchor=(0.82, 0.5),labelspacing=0.1,handlelength=1)

plt.savefig("LDL_HDL_IR.png", dpi =600)
plt.show()
plt.close()



fig = plt.figure(figsize=(4,3))
ax =fig.add_axes([0.16, 0.15, 0.835, 0.845]) #, frameon=False, xticklabels=[],yticklabels=[])

model = np.loadtxt('./Atm_P/population_real_units', comments='#')
R = 1.380649*6.022

HDA = log(model[:,3])
LDA = log(model[:,4]+model[:,5])

ax.plot(model[:,1],LDA-HDA,color="k",linestyle='-',linewidth="2")#,label=r'$\mathrm{}$')
ax.plot(xdat-48,num,color='k',linestyle=':')
ax.plot(xdat-41,num,color='k',linestyle=':')
ax.plot(model[:,1],model[:,1]*0,color='indigo',linestyle='--',linewidth="2",label=r'$\mathrm{0^o C}$')
ax.set_ylabel(r'$\mathrm{\mathbf{ln K}}$',{'fontsize':MEDIUM_SIZE})
ax.yaxis.set_label_coords(-0.12, 0.5)
ax.set_xlabel(r'$\mathrm{\mathbf{T}\ [{^\circ} C]}$',{'fontsize':MEDIUM_SIZE})
ax.set_yticks(np.arange(-4, 5, 2))
ax.set_ylim([-5, 5])
ax.set_xlim([-100, 80])
ax.set_xticks(np.arange(-50, 100, 50))
ax.tick_params(axis='both',direction = 'in', which= 'both',bottom=1, top=1, left=1,right=1)

plt.savefig("lnK.png", dpi =600)
plt.show()
plt.close()


#extracting exmodel from colormap
cmap = cm.get_cmap('plasma', 25)    # PiYG
plasma=["" for x in range(25)]

for i in range(cmap.N):
        rgb = cmap(i)[:3] # will return rgba, we take only first 3 so we get rgb
        plasma[i]=matplotlib.colors.rgb2hex(rgb)

cmap = cm.get_cmap('Reds_r', 20)    # PiYG
Reds=["" for x in range(20)]

for i in range(cmap.N):
        rgb = cmap(i)[:3] # will return rgba, we take only first 3 so we get rgb
        Reds[i]=matplotlib.colors.rgb2hex(rgb)

cmap = cm.get_cmap('Blues', 10)    # PiYG
Blues=["" for x in range(10)]

for i in range(cmap.N):
        rgb = cmap(i)[:3] # will return rgba, we take only first 3 so we get rgb
        Blues[i]=matplotlib.colors.rgb2hex(rgb)

#LLCP

fig = plt.figure(figsize=(4,4))
ax2 =fig.add_axes([0.16, 0.16, 0.835, 0.835])
P = [3000, 2800, 2600, 2400, 2200, 2000, 1800, 1600, 1400]#,  1400, 1300, 1200, 1100, 1000]
for i in range(len(P)):
    model = np.loadtxt('./Press/real_units_%sbar'%(P[i]), comments='#')
    ax2.plot(model[::,1]+273.15,model[::,4]*0.001,linestyle='-', color=Reds[i*2] )

model = np.loadtxt('./Press/real_units_3000bar', comments='#')
ax2.plot(model[:,1]+273.15,model[:,4]*0.001,linestyle='-',label=r'$\mathrm{300MPa}$', color='k',linewidth="2")

model = np.loadtxt('./Press/real_units_1400bar', comments='#')
ax2.plot(model[:,1]+273.15,model[:,4]*0.001,linestyle='--',label=r'$\mathrm{140MPa}$', color='k',linewidth="2")

model = np.loadtxt('LLPT', comments='#')
ax2.plot(model[::2,1],model[::2,2]*0.001,linestyle='',marker='.',mew=0.5,linewidth="2", color='g')
ax2.plot(model[::2,3],model[::2,4]*0.001,linestyle='',marker='.',mew=0.5,label=r'$\mathrm{binodal}$', linewidth="2", color='g' )

ax2.plot(159,1.055,color='b',linestyle='', marker='o',  mew = 3, label=r'$\mathrm{LLCP}$',fillstyle='none')
ax2.set_ylabel(r'$\mathrm{\mathbf{\rho}\ [g/cc]}$',{'fontsize':MEDIUM_SIZE})
ax2.set_xlabel(r'$\mathrm{\mathbf{T}\ [{^\circ} C]}$',{'fontsize':MEDIUM_SIZE})
ax2.set_yticks(np.arange(0.9, 1.3, 0.1))
ax2.set_ylim([0.9,1.3])
ax2.set_xlim([100, 320])
ax2.set_xticks(np.arange(100, 320, 50))
ax2.tick_params(axis='both',direction = 'in', which= 'both',bottom=1, top=1, left=1,right=1)
ax2.legend(frameon=False, fontsize = "9")#, ncol=1,loc = "center", bbox_to_anchor=(1.18,
plt.savefig("multi_density_LLCP.png", dpi =600)
plt.show()
plt.close()


# Vertical T 
fig = plt.figure(figsize=(4,4))
ax2 =fig.add_axes([0.16, 0.16, 0.835, 0.835])
P = [3400, 3200, 3000, 2800, 2600, 2400, 2200, 2000, 1800, 1600, 1400,  1200, 1000]#, 800, 600,]
for i in range(len(P)):
    model = np.loadtxt('./Press/real_units_%sbar'%(P[i]), comments='#')
    ax2.plot(model[::,4]*0.001,273.15+model[::,1],linestyle='-', color=Reds[i] )

model = np.loadtxt('./Press/real_units_3400bar', comments='#')
ax2.plot(model[:,4]*0.001,273.15+model[:,1],linestyle='--',label=r'$\mathrm{340MPa}$', color='k',linewidth="2")

model = np.loadtxt('./Press/real_units_1000bar', comments='#')
ax2.plot(model[:,4]*0.001,273.15+model[:,1],linestyle='--',label=r'$\mathrm{100MPa}$', color='k',linewidth="2")

model = np.loadtxt('LLPT', comments='#')
ax2.plot(model[::2,2]*0.001,model[::2,1],linestyle='',marker='.',mew=0.5,linewidth="2", color='g')
ax2.plot(model[::2,4]*0.001,model[::2,3],linestyle='',marker='.',mew=0.5,label=r'$\mathrm{binodal}$', linewidth="2", color='g' )

ax2.plot(1.055,159,color='b',linestyle='', marker='o',  mew = 3, label=r'$\mathrm{LLCP}$',fillstyle='none')
ax2.set_xlabel(r'$\mathrm{\mathbf{\rho}\ [g/cc]}$',{'fontsize':MEDIUM_SIZE})
ax2.set_ylabel(r'$\mathrm{\mathbf{T}\ [{^\circ} C]}$',{'fontsize':MEDIUM_SIZE})
ax2.set_xticks(np.arange(0.9, 1.3, 0.1))
ax2.set_xlim([0.9,1.3])
ax2.set_ylim([100, 320])
ax2.set_yticks(np.arange(100, 320, 50))
ax2.tick_params(axis='both',direction = 'in', which= 'both',bottom=1, top=1, left=1,right=1)
ax2.legend(frameon=False, fontsize = "9")#, ncol=1,loc = "center", bbox_to_anchor=(1.18,
plt.savefig("multi_density_LLCP_T(rho).png", dpi =600)
plt.show()
plt.close()


fig = plt.figure(figsize=(4,4))
ax =fig.add_axes([0.16, 0.16, 0.835, 0.835]) #, frameon=False, xticklabels=[],yticklabels=[])

# Temp of max density
#model = np.loadtxt('PT_data', comments='#')
tip4p = np.loadtxt('tip4p05_data/max_rho', comments='#')
Palla = np.loadtxt('Expt_Data/Pallares_PCCP_2016', comments='#')
Caldw = np.loadtxt('Expt_Data/Caldwell_DSR_1978', comments='#')
Mishi = np.loadtxt('Expt_Data_plot/Mishima_max_rho', comments='#')
Sotan = np.loadtxt('Expt_Data_plot/Sotani_max_rho', comments='#')


model = np.loadtxt('./PT_data', comments='#')


ax.plot(model[20:,5],model[20:,0],linestyle='-',linewidth="2",color="teal",label=r'$\mathrm{C_{p}^{max}}$' )
ax.plot(model[44:,3],model[44:,0],linestyle='-',linewidth="2",color="crimson",label=r'$\mathrm{\kappa_{T}^{max}}$' )
ax.plot(model[20:,8],model[20:,0],linestyle='-',linewidth="2",color="blueviolet",label=r'$\mathrm{\alpha_{p}^{min}}$' )
ax.plot(model[:,9],model[:,0],color="k",linestyle='--',linewidth="2",label=r'$\mathrm{WL}$')

ax.plot(159,176,color='r',linestyle='', marker='o',  mew = 2.5, label=r'$\mathrm{LLCP}$',fillstyle='none')
ax.set_ylabel(r'$\mathrm{\mathbf{p}\ [MPa]}$',{'fontsize':MEDIUM_SIZE})
ax.yaxis.set_label_coords(-0.12, 0.5)
ax.set_xlabel(r'$\mathrm{\mathbf{T}\ [K]}$',{'fontsize':MEDIUM_SIZE})
ax.set_ylim([-300, 300])
ax.set_yticks(np.arange(-200, 300, 100))
ax.set_xlim([80, 320])
ax.set_xticks(np.arange(100, 320, 50))
ax.tick_params(axis='both',direction = 'in', which= 'both',bottom=1, top=1, left=1,right=1)
ax.legend( frameon=False, fontsize = "9", ncol =3,labelspacing=0.1,handlelength=1, loc = "center", bbox_to_anchor=(0.4, 0.1))

plt.savefig("PT_all.png", dpi =600)
plt.show()
plt.close()

fig = plt.figure(figsize=(4,4))
ax =fig.add_axes([0.16, 0.16, 0.835, 0.835]) #, frameon=False, xticklabels=[],yticklabels=[])

#model = np.loadtxt('PT_data', comments='#')
mbpol = np.loadtxt('Expt_Data/Sciortino_NatPhys_2025', comments='#')
tip4p = np.loadtxt('tip4p05_data/max_rho', comments='#')
Palla = np.loadtxt('Expt_Data/Pallares_PCCP_2016', comments='#')
Caldw = np.loadtxt('Expt_Data/Caldwell_DSR_1978', comments='#')
Mishi = np.loadtxt('Expt_Data_plot/Mishima_max_rho', comments='#')
Sotan = np.loadtxt('Expt_Data_plot/Sotani_max_rho', comments='#')

ax.plot(Palla[93::2,0],Palla[93::2,1],color='darkorange',linestyle='', marker='^', mew = 2,label=r'$\mathrm{\rho^{max}_{Expt.}}$',fillstyle='none')
ax.plot(Caldw[::2,0],Caldw[::2,1],color='darkorange',linestyle='', marker='^', mew = 2,fillstyle='none')
ax.plot(Mishi[:,0],Mishi[:,1],color='darkorange',linestyle='', marker='^', mew = 2,fillstyle='none')
ax.plot(Sotan[:,0],Sotan[:,1],color='darkorange',linestyle='', marker='^', mew = 2,fillstyle='none')
ax.plot(tip4p[:,0],tip4p[:,1],color='g',linestyle='', marker='o',  mew = 2, label=r'$\mathrm{\rho^{max}_{TIP4P/2005}}$',fillstyle='none')
ax.plot(mbpol[:,0],mbpol[:,1]*0.1,color='b',linestyle='', marker='o',  mew = 2, label=r'$\mathrm{\rho^{max}_{DNN@MBpol}}$',fillstyle='none')

model = np.loadtxt('./PT_data', comments='#')


#


ax.plot(model[0:,1],model[0:,0],color="b",linestyle='-',linewidth="2",marker="",fillstyle='none',label=r'$\mathrm{\rho^{max}}$' )

ax.plot(model[:156*2,5],model[:156*2,0],linestyle='-',linewidth="2",color="teal",label=r'$\mathrm{C_{p}^{max}}$' )
ax.plot(model[20:156*2,3],model[20:156*2,0],linestyle='-',linewidth="2",color="crimson",label=r'$\mathrm{\kappa_{T}^{max}}$' )
ax.plot(model[:156*2,8],model[:156*2,0],linestyle='-',linewidth="2",color="blueviolet",label=r'$\mathrm{\alpha_{p}^{min}}$' )
ax.plot(model[:156*2,9],model[:156*2,0],color="k",linestyle='--',linewidth="2",label=r'$\mathrm{WL}$')
ax.plot(model[156*2:,9],model[156*2:,0],color="k",linestyle='-',linewidth="4",label=r'$\mathrm{PT line}$')


ax.plot(159,176,color='r',linestyle='', marker='o',  mew = 2.5, label=r'$\mathrm{LLCP}$',fillstyle='none')
ax.set_ylabel(r'$\mathrm{\mathbf{p}\ [MPa]}$',{'fontsize':MEDIUM_SIZE})
ax.yaxis.set_label_coords(-0.12, 0.5)
ax.set_xlabel(r'$\mathrm{\mathbf{T}\ [K]}$',{'fontsize':MEDIUM_SIZE})
ax.set_ylim([-300, 300])
ax.set_yticks(np.arange(-200, 300, 100))
ax.set_xlim([100, 320])
ax.set_xticks(np.arange(100, 320, 50))
ax.tick_params(axis='both',direction = 'in', which= 'both',bottom=1, top=1, left=1,right=1)
ax.legend( frameon=False, fontsize = "9", ncol =1,labelspacing=0.05,handlelength=1.0)

plt.savefig("PT_all2.png", dpi =600)
plt.show()
plt.close()


fig = plt.figure(figsize=(4,4))
ax =fig.add_axes([0.16, 0.16, 0.835, 0.835]) 

# Temp of max density
#model = np.loadtxt('PT_data', comments='#')
mbpol = np.loadtxt('Expt_Data/Sciortino_NatPhys_2025', comments='#')
tip4p = np.loadtxt('tip4p05_data/max_rho', comments='#')
Palla = np.loadtxt('Expt_Data/Pallares_PCCP_2016', comments='#')
Caldw = np.loadtxt('Expt_Data/Caldwell_DSR_1978', comments='#')
Mishi = np.loadtxt('Expt_Data_plot/Mishima_max_rho', comments='#')
Sotan = np.loadtxt('Expt_Data_plot/Sotani_max_rho', comments='#')

ax.plot(Palla[93::2,0],Palla[93::2,1],color='darkorange',linestyle='', marker='^', mew = 2,label=r'$\mathrm{\rho^{max}_{Expt.}}$',fillstyle='none')
ax.plot(Caldw[::2,0],Caldw[::2,1],color='darkorange',linestyle='', marker='^', mew = 2,fillstyle='none')
ax.plot(Mishi[:,0],Mishi[:,1],color='darkorange',linestyle='', marker='^', mew = 2,fillstyle='none')
ax.plot(Sotan[:,0],Sotan[:,1],color='darkorange',linestyle='', marker='^', mew = 2,fillstyle='none')
ax.plot(tip4p[:,0],tip4p[:,1],color='g',linestyle='', marker='o',  mew = 2, label=r'$\mathrm{\rho^{max}_{TIP4P/2005}}$',fillstyle='none')
ax.plot(mbpol[:,0],mbpol[:,1]*0.1,color='b',linestyle='', marker='o',  mew = 2, label=r'$\mathrm{\rho^{max}_{MB-pol}}$',fillstyle='none')

model = np.loadtxt('./PT_data', comments='#')
ax.plot(model[28:,1],model[28:,0],color="k",linestyle='-',linewidth="2",marker="",fillstyle='none',label=r'$\mathrm{\rho^{max}_{Theory}}$' )
ax.set_ylabel(r'$\mathrm{\mathbf{p}\ [MPa]}$',{'fontsize':MEDIUM_SIZE})
ax.yaxis.set_label_coords(-0.12, 0.5)
ax.set_xlabel(r'$\mathrm{\mathbf{T}\ [K]}$',{'fontsize':MEDIUM_SIZE})
ax.set_ylim([-300, 300])
ax.set_yticks(np.arange(-200, 300, 100))
ax.set_xlim([190, 300])
ax.set_xticks(np.arange(200, 300, 20))
ax.tick_params(axis='both',direction = 'in', which= 'both',bottom=1, top=1, left=1,right=1)
ax.legend( frameon=False, fontsize = "10", ncol =1,labelspacing=0.1,handlelength=1, loc = "lower left")#, bbox_to_anchor=(0.4, 0.15))

plt.savefig("PT_dens_max.png", dpi =600)
plt.show()
plt.close()


fig = plt.figure(figsize=(4,6))
ax1 =fig.add_axes([0.15, 0.544, 0.845, 0.455],xticklabels=[])
ax2 =fig.add_axes([0.15, 0.080, 0.845, 0.455])

#240.00   132.95  1085.41   133.15 1098.85  1092.13

P = [3000, 2900, 2800, 2700, 2600, 2500, 2400, 2300, 2200, 2100, 2000, 1900, 1800, 1700, 1600, 1500]#,  1400, 1300, 1200, 1100, 1000]
P = [3000, 2800, 2600, 2400, 2200, 2000, 1800, 1600, 1400]
for i in range(len(P)):
    model = np.loadtxt('./Press/real_units_%sbar'%(P[i]), comments='#')
    ax2.plot(model[::,1]+273.15,model[::,4]*0.001,linestyle='-', color=Reds[i*2] )

model = np.loadtxt('./Press/real_units_3000bar', comments='#')
ax2.plot(model[:,1]+273.15,model[:,4]*0.001,linestyle='-',label=r'$\mathrm{300MPa}$', color='k',linewidth="2")

model = np.loadtxt('./Press/real_units_1400bar', comments='#')
ax2.plot(model[:,1]+273.15,model[:,4]*0.001,linestyle='--',label=r'$\mathrm{140MPa}$', color='k',linewidth="2")

model = np.loadtxt('LLPT', comments='#')
ax2.plot(model[::2,1],model[::2,2]*0.001,linestyle='',marker='.',mew=0.5,linewidth="2", color='g')
ax2.plot(model[::2,3],model[::2,4]*0.001,linestyle='',marker='.',mew=0.5,label=r'$\mathrm{binodal}$', linewidth="2", color='g' )

ax2.plot(159,1.055,color='b',linestyle='', marker='o',  mew = 3, label=r'$\mathrm{LLCP}$',fillstyle='none')
ax2.set_ylabel(r'$\mathrm{\mathbf{\rho}\ [g/cc]}$',{'fontsize':MEDIUM_SIZE})
ax2.set_xlabel(r'$\mathrm{\mathbf{T}\ [K]}$',{'fontsize':MEDIUM_SIZE})
ax2.set_yticks(np.arange(0.9, 1.3, 0.1))
ax2.set_ylim([0.9,1.3])
ax2.set_xlim([100, 320])
ax2.set_xticks(np.arange(100, 320, 50))
ax2.tick_params(axis='both',direction = 'in', which= 'both',bottom=1, top=1, left=1,right=1)
ax2.legend(frameon=False, fontsize = "9")#, ncol=1,loc = "center", bbox_to_anchor=(1.18, 0.5))


# Temp of max density
#model = np.loadtxt('PT_data', comments='#')
tip4p = np.loadtxt('tip4p05_data/max_rho', comments='#')
mbpol = np.loadtxt('Expt_Data/Sciortino_NatPhys_2025', comments='#')
Palla = np.loadtxt('Expt_Data/Pallares_PCCP_2016', comments='#')
Caldw = np.loadtxt('Expt_Data/Caldwell_DSR_1978', comments='#')
Mishi = np.loadtxt('Expt_Data_plot/Mishima_max_rho', comments='#')
Sotan = np.loadtxt('Expt_Data_plot/Sotani_max_rho', comments='#')

ax1.plot(Palla[93::2,0],Palla[93::2,1],color='darkorange',linestyle='', marker='^', mew = 2,label=r'$\mathrm{\rho^{max}_{Expt.}}$',fillstyle='none')
ax1.plot(Caldw[::2,0],Caldw[::2,1],color='darkorange',linestyle='', marker='^', mew = 2,fillstyle='none')
ax1.plot(Mishi[:,0],Mishi[:,1],color='darkorange',linestyle='', marker='^', mew = 2,fillstyle='none')
ax1.plot(Sotan[:,0],Sotan[:,1],color='darkorange',linestyle='', marker='^', mew = 2,fillstyle='none')
model = np.loadtxt('./PT_data', comments='#')

ax1.plot(model[1:,1],model[1:,0],color="b",linestyle='-',linewidth="2",marker="",fillstyle='none',label=r'$\mathrm{\rho^{max}}$' )


ax1.plot(model[:156*2,5],model[:156*2,0],linestyle='-',linewidth="2",color="teal",label=r'$\mathrm{C_{p}^{max}}$' )
ax1.plot(model[20:156*2,3],model[20:156*2,0],linestyle='-',linewidth="2",color="crimson",label=r'$\mathrm{\kappa_{T}^{max}}$' )
ax1.plot(model[:156*2,8],model[:156*2,0],linestyle='-',linewidth="2",color="blueviolet",label=r'$\mathrm{\alpha_{p}^{min}}$' )
ax1.plot(model[:156*2,9],model[:156*2,0],color="k",linestyle='--',linewidth="2",label=r'$\mathrm{WL}$')
ax1.plot(model[156*2:,9],model[156*2:,0],color="k",linestyle='-',linewidth="4",label=r'$\mathrm{PT line}$')


ax1.plot(159,176,color='r',linestyle='', marker='o',  mew = 2.5, label=r'$\mathrm{LLCP}$',fillstyle='none')
ax1.set_ylabel(r'$\mathrm{\mathbf{p}\ [MPa]}$',{'fontsize':MEDIUM_SIZE})
ax1.yaxis.set_label_coords(-0.12, 0.5)
ax1.set_xlabel(r'$\mathrm{\mathbf{T}\ [K]}$',{'fontsize':MEDIUM_SIZE})
ax1.set_ylim([-300, 300])
ax1.set_yticks(np.arange(-200, 300, 100))
ax1.set_xlim([100, 320])
ax1.set_xticks(np.arange(100, 320, 50))
ax1.tick_params(axis='both',direction = 'in', which= 'both',bottom=1, top=1, left=1,right=1)
ax1.legend( frameon=False, fontsize = "9", ncol =1,labelspacing=0.05,handlelength=1.0)
plt.savefig("PT_phase.png", dpi =600)
plt.show()
plt.close()
