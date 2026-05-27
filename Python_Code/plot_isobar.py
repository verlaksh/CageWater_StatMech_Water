import scipy
import numpy as np
from pylab import *
from scipy.optimize import curve_fit
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)



SMALL_SIZE = 11
MEDIUM_SIZE = 14
BIGGER_SIZE = 21
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=12)    # legend fontsize

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

cmap = cm.get_cmap('Blues', 12)    # PiYG
Blues=["" for x in range(12)]

for i in range(cmap.N):
        rgb = cmap(i)[:3] # will return rgba, we take only first 3 so we get rgb
        Blues[i]=matplotlib.colors.rgb2hex(rgb)


cmap = cm.get_cmap('Reds_r', 20)    # PiYG
colors=["" for x in range(20)]

for i in range(cmap.N):
        rgb = cmap(i)[:3] # will return rgba, we take only first 3 so we get rgb
        colors[i]=matplotlib.colors.rgb2hex(rgb)


fig = plt.figure(figsize=(6.2,6))

ax =fig.add_axes([0.1, 0.14, 0.8, 0.78], frameon=False, xticklabels=[],yticklabels=[])
ax.tick_params(axis='both',direction = 'in',bottom=False, top=False, left=False, right=False)
ax1 = fig.add_axes([0.09, 0.545, 0.403, 0.45], xticklabels=[])#, ylim=(0, 0.6))
ax2 = fig.add_axes([0.09, 0.090, 0.403, 0.45])#, xticklabels=[])#, ylim=(0, 0.6))
ax3 = fig.add_axes([0.59, 0.545, 0.403, 0.45], xticklabels=[])#, ylim=(0, 0.6))
ax4 = fig.add_axes([0.59, 0.090, 0.403, 0.45])#, xticklabels=[])#, ylim=(0, 0.6))

P = [3000, 2400, 2000, 1600, 1200, 800, 400, 1]
for i in range(len(P)):
    print(i*2+4)
    model = np.loadtxt('./Press/real_units_%sbar'%(P[i]), comments='#')
    ax1.plot(model[:,1],model[:,4]*0.001,linestyle='-',label=r'$\mathrm{%sMPa}$'%(P[i]*0.1), color=Reds[i+2])#, linewidth =2 )

P = [200, 400, 600, 800, 1000]#, -1200]#, -1400, -1600, -2000]
for i in range(len(P)):
    print(i*2+3)
    model = np.loadtxt('./Press/real_units_-%sbar'%(P[i]), comments='#')
    ax1.plot(model[:,1],model[:,4]*0.001,linestyle='-',label=r'$\mathrm{%sMPa}$'%(P[i]*0.1), color=Blues[i+3])#, linewidth = 2 )

ax1.set_ylabel(r'$\mathrm{\mathbf{\rho}[g/cc]}$',{'fontsize':MEDIUM_SIZE})
ax1.yaxis.set_label_coords(-0.13, 0.5)
ax1.set_yticks(np.arange(0.9, 1.3, 0.1))
ax1.set_ylim([0.85,1.25])
ax1.set_xlim([-130, 100])
ax1.set_xticks(np.arange(-130, 100, 50))
ax1.tick_params(axis='both',direction = 'in', which= 'both',bottom=1, top=1, left=1,right=1)

#Heat Capacity

P = [3000, 2400, 2000, 1600, 1200, 800, 400, 1]
for i in range(len(P)):
    model = np.loadtxt('./Press/real_units_%sbar'%(P[i]), comments='#')
    ax2.plot(model[:,1],model[:,7]/18,linestyle='-',label=r'$\mathrm{%sMPa}$'%(P[i]*0.1), color=Reds[i+2])#, linewidth =1.5 )

P = [-200, -400, -600, -800, -1000]#, -1200]#, -1400, -1600, -2000]
for i in range(len(P)):
    model = np.loadtxt('./Press/real_units_%sbar'%(P[i]), comments='#')
    ax2.plot(model[:,1],model[:,7]/18,linestyle='-',label=r'$\mathrm{%sMPa}$'%(P[i]*0.1), color=Blues[i+3] )

ax2.set_ylabel(r'$\mathrm{\mathbf{C_{p}}\ [J K^{-1} g^{-1}]}$',{'fontsize':MEDIUM_SIZE})
ax2.yaxis.set_label_coords(-0.11, 0.5)
ax2.set_xlabel(r'$\mathrm{\mathbf{T}[{^\circ} C]}$',{'fontsize':MEDIUM_SIZE})
ax2.set_yticks(np.arange(2, 7, 2))
ax2.set_ylim([2,7])
ax2.set_xlim([-130, 100])
ax2.set_xticks(np.arange(-130, 100, 50))
ax2.tick_params(axis='both',direction = 'in', which= 'both',bottom=1, top=1, left=1,right=1)


#Compressibility

P = [3000, 2400, 2000, 1600, 1200, 800, 400, 1]
for i in range(len(P)):
    model = np.loadtxt('./Press/real_units_%sbar'%(P[i]), comments='#')
    ax3.plot(model[:,1],model[:,5]*10**4,linestyle='-',label=r'$\mathrm{%sMPa}$'%(P[i]*0.1), color=colors[i+2])#, linewidth =2 )


P = [-200, -400, -600, -800, -1000]
for i in range(len(P)):
    model = np.loadtxt('./Press/real_units_%sbar'%(P[i]), comments='#')
    ax3.plot(model[:,1],model[:,5]*10**4,linestyle='-',color=Blues[i+3] )

ax3.set_ylabel(r'$\mathrm{\mathbf{\kappa_T} \times 10^4\/[bar^{-1}]}$',{'fontsize':MEDIUM_SIZE})
ax3.yaxis.set_label_coords(-0.11, 0.5)
ax3.set_yticks(np.arange(1, 6, 2))
ax3.set_ylim([0,6])
ax3.set_xlim([-130, 100])
ax3.set_xticks(np.arange(-130, 100, 50))
ax3.tick_params(axis='both',direction = 'in', which= 'both',bottom=1, top=1, left=1,right=1)
ax3.legend(frameon=False, ncol=1,loc = "best", fontsize = "10",handlelength=1)


#Alpha

P = [3000, 2400, 2000, 1600, 1200, 800, 400, 1]
for i in range(len(P)):
    model = np.loadtxt('./Press/real_units_%sbar'%(int(P[i])), comments='#')
    ax4.plot(model[:,1],model[:,6]*10**3,linestyle='-', color=Reds[i+2] )

P = [-200, -400, -600, -800, -1000]
for i in range(len(P)):
    model = np.loadtxt('./Press/real_units_%sbar'%(P[i]), comments='#')
    ax4.plot(model[:,1],model[:,6]*10**3,linestyle='-',label=r'$\mathrm{%sMPa}$'%(P[i]*0.1), color=Blues[i+3] )
ax4.set_ylabel(r'$\mathrm{\mathbf{\alpha_P} \times 10^3\/[K^{-1}]}$',{'fontsize':MEDIUM_SIZE})
ax4.yaxis.set_label_coords(-0.11, 0.5)
ax4.set_xlabel(r'$\mathrm{\mathbf{T}[{^\circ} C]}$',{'fontsize':MEDIUM_SIZE})
ax4.set_yticks(np.arange(-9, 3, 3))
ax4.set_ylim([-10,3])
ax4.set_xlim([-130, 100])
ax4.set_xticks(np.arange(-130, 100, 50))
ax4.tick_params(axis='both',direction = 'in', which= 'both',bottom=1, top=1, left=1,right=1)
ax4.legend(frameon=False, ncol=1,loc = "best", fontsize = "10",handlelength=1)

plt.savefig("isobars2.png", dpi =600)
plt.show()
plt.close()


cmap = cm.get_cmap('Reds_r', 20)    # PiYG
Reds=["" for x in range(20)]

for i in range(cmap.N):
        rgb = cmap(i)[:3] # will return rgba, we take only first 3 so we get rgb
        Reds[i]=matplotlib.colors.rgb2hex(rgb)

cmap = cm.get_cmap('Blues', 12)    # PiYG
Blues=["" for x in range(12)]

for i in range(cmap.N):
        rgb = cmap(i)[:3] # will return rgba, we take only first 3 so we get rgb
        Blues[i]=matplotlib.colors.rgb2hex(rgb)



fig = plt.figure(figsize=(5,4))
ax =fig.add_axes([0.13, 0.14, 0.65, 0.85])

exp = np.loadtxt('./Expt_Data_plot/JCP_Mishima_2010', comments='#')
ax.plot(exp[:,0]-273.15,exp[:,2]*0.001,color='k',linestyle=' ', marker='.', fillstyle='none', label=r'$\mathrm{Exp.}$', mew = 2)#,fillstyle='none')

exp = np.loadtxt('./Expt_Data_plot/pccp_pallares_2016', comments='#')
ax.plot(exp[:,0]-273.15,exp[:,2]*0.001,color='k',linestyle=' ', marker='.', fillstyle='none', mew = 2)#, label=r'$\mathrm{Exp.}$')


exp = np.loadtxt('density', comments='#')
ax.plot(exp[::3,0],exp[::3,2],color='k',linestyle=' ', marker='.', fillstyle='none', mew = 2)#, label=r'$\mathrm{Exp.}$',fillstyle='none')


P = [4000, 3000, 2400, 2000, 1600, 1200, 1000, 800, 600, 400, 1]

for i in range(len(P)):
    model = np.loadtxt('./Press/real_units_%sbar'%(P[i]), comments='#')
    ax.plot(model[:,1],model[:,4]*0.001,linestyle='-',label=r'$\mathrm{%s}$'%(P[i]*0.1), color=Reds[i], linewidth = 2)
P = [-200, -400, -600, -800, -1000]#, -1200, -1400, -1500, -2000]
for i in range(len(P)):
    model = np.loadtxt('./Press/real_units_%sbar'%(P[i]), comments='#')
    ax.plot(model[:,1],model[:,4]*0.001,linestyle='-',label=r'$\mathrm{%s}$'%(P[i]*0.1), color=Blues[i+3], linewidth = 2 )


ax.set_ylabel(r'$\mathrm{\mathbf{\rho}[g/cc]}$',{'fontsize':MEDIUM_SIZE})
ax.set_xlabel(r'$\mathrm{\mathbf{T}[{^\circ} C]}$',{'fontsize':MEDIUM_SIZE})
ax.set_yticks(np.arange(0.9, 1.2, 0.1))
ax.set_ylim([0.88,1.2])
ax.set_xlim([-130, 100])
ax.set_xticks(np.arange(-50, 100, 50))
ax.tick_params(axis='both',direction = 'in', which= 'both',bottom=1, top=1, left=1,right=1)
plt.legend(frameon=False, ncol=1,loc = "center", bbox_to_anchor=(1.18, 0.45), fontsize = "11",labelspacing=0.2)
plt.savefig("multi_density_nPi_small.png", dpi =600)
plt.show()
plt.close()

fig = plt.figure(figsize=(3,4))
ax =fig.add_axes([0.18, 0.11, 0.815, 0.885])
P = [4000, 3000, 2400, 2000, 1600, 1200, 800, 400, 1]
for i in range(len(P)):
    print(i*2+4)
    model = np.loadtxt('./Press/real_units_%sbar'%(P[i]), comments='#')
    ax.plot(model[:,1],model[:,4]*0.001,linestyle='-',label=r'$\mathrm{%sMPa}$'%(P[i]*0.1), color=Reds[i+2], linewidth =2 )
    exp = np.loadtxt('./Expt_Data_plot/Dens/%sbar'%(P[i]), comments='#')
    ax.plot(exp[:,0]-273.15,exp[:,2]*0.001,color=Reds[i+2],linestyle=' ', marker='.', fillstyle='none', label=r'$\mathrm{Exp.}$', mew =1.5)
    if P[i] == 1:
        print(P[i])
        exp = np.loadtxt('./Expt_Data_plot/Dens/%sbar'%(P[i]), comments='#')
        ax1.plot(exp[::5,0],exp[::5,2],color=Reds[i+2],linestyle=' ', marker='.', fillstyle='none', mew =1.5)

P = [200, 400, 600, 800, 1000]#, -1200]#, -1400, -1600, -2000]
for i in range(len(P)):
    print(i*2+3)
    model = np.loadtxt('./Press/real_units_-%sbar'%(P[i]), comments='#')
    ax.plot(model[:,1],model[:,4]*0.001,linestyle='-',label=r'$\mathrm{%sMPa}$'%(P[i]*0.1), color=Blues[i+5], linewidth = 2 )
    exp = np.loadtxt('./Expt_Data_plot/Dens/n%sbar'%(P[i]), comments='#')
    ax.plot(exp[:,0]-273.15,exp[:,2]*0.001,color=Blues[i+5],linestyle=' ', marker='.', fillstyle='none', mew =1.5)#, label=r'$\mathrm{Exp.}$')

ax.set_ylabel(r'$\mathrm{\mathbf{\rho}[g cm^3]}$',{'fontsize':MEDIUM_SIZE})
ax.yaxis.set_label_coords(-0.115, 0.5)
ax.set_xlabel(r'$\mathrm{\mathbf{T}[{^\circ} C]}$',{'fontsize':MEDIUM_SIZE})
ax.set_yticks(np.arange(0.9, 1.2, 0.1))
ax.set_ylim([0.88,1.2])
ax.set_xlim([-70, 100])
ax.set_xticks(np.arange(-50, 100, 50))
ax.tick_params(axis='both',direction = 'in', which= 'both',bottom=1, top=1, left=1,right=1)
plt.show()
plt.savefig("Rho_isobars.png", dpi =600)
plt.close()

fig = plt.figure(figsize=(4,4))
ax =fig.add_axes([0.13, 0.12, 0.61, 0.875])
P = [4000, 3000, 2400, 2000, 1600, 1200, 800, 400, 1]
for i in range(len(P)):
    if P[i]>1:
        model = np.loadtxt('./Press/real_units_%sbar'%(P[i]), comments='#')
        ax.plot(model[:,1],model[:,4]*0.001,linestyle='-',label=r'$\mathrm{%sMPa}$'%int(P[i]*0.1), color=Reds[i+2], linewidth =2 )
        exp = np.loadtxt('./Expt_Data_plot/Dens/%sbar'%(P[i]), comments='#')
        ax.plot(exp[:,0]-273.15,exp[:,2]*0.001,color=Reds[i+2],linestyle=' ', marker='.', fillstyle='none', mew =1.5)
    if P[i] == 1:
        exp = np.loadtxt('./Expt_Data_plot/Dens/%sbar'%(P[i]), comments='#')
        ax.plot(exp[::3,0],exp[::3,2],color=Reds[i+2],linestyle=' ', marker='.', fillstyle='none', mew =1.5)
        model = np.loadtxt('./Press/real_units_%sbar'%(P[i]), comments='#')
        ax.plot(model[:,1],model[:,4]*0.001,linestyle='-',label=r'$\mathrm{%sMPa}$'%(P[i]*0.1), color=Reds[i+2], linewidth =2 )

P = [200, 400, 600, 800, 1000]#, -1200]#, -1400, -1600, -2000]
for i in range(len(P)):
    print(i*2+3)
    model = np.loadtxt('./Press/real_units_-%sbar'%(P[i]), comments='#')
    ax.plot(model[:,1],model[:,4]*0.001,linestyle='-',label=r'$\mathrm{-%sMPa}$'%int(P[i]*0.1), color=Blues[i+5], linewidth = 2 )
    exp = np.loadtxt('./Expt_Data_plot/Dens/n%sbar'%(P[i]), comments='#')
    ax.plot(exp[:,0]-273.15,exp[:,2]*0.001,color=Blues[i+5],linestyle=' ', marker='.', fillstyle='none', mew =1.5)#, label=r'$\mathrm{Exp.}$')

ax.set_ylabel(r'$\mathrm{\mathbf{\rho}[g cm^3]}$',{'fontsize':MEDIUM_SIZE})
ax.yaxis.set_label_coords(-0.12, 0.52)
ax.set_xlabel(r'$\mathrm{\mathbf{T}[{^\circ} C]}$',{'fontsize':MEDIUM_SIZE})
ax.set_yticks(np.arange(0.9, 1.2, 0.1))
ax.set_ylim([0.88,1.2])
ax.set_xlim([-70, 100])
ax.set_xticks(np.arange(-50, 100, 50))
ax.tick_params(axis='both',direction = 'in', which= 'both',bottom=1, top=1, left=1,right=1)
plt.legend(frameon=False, ncol=1,loc = "center", bbox_to_anchor=(1.23, 0.5), fontsize = "10",handlelength=1)
plt.savefig("Rho_isobars2.png", dpi =600)
plt.show()
plt.close()

fig = plt.figure(figsize=(3.5,4))
ax =fig.add_axes([0.18, 0.11, 0.815, 0.885])
# Temp of max density
#model = np.loadtxt('PT_data', comments='#')
mbpol = np.loadtxt('Expt_Data/Sciortino_NatPhys_2025', comments='#')
tip4p = np.loadtxt('tip4p05_data/max_rho', comments='#')
Palla = np.loadtxt('Expt_Data/Pallares_PCCP_2016', comments='#')
Caldw = np.loadtxt('Expt_Data/Caldwell_DSR_1978', comments='#')
Mishi = np.loadtxt('Expt_Data_plot/Mishima_max_rho', comments='#')
Sotan = np.loadtxt('Expt_Data_plot/Sotani_max_rho', comments='#')

ax.plot(Palla[93::2,0]-273.15,Palla[93::2,1],color='darkorange',linestyle='', marker='^', mew = 2,label=r'$\mathrm{\rho^{max}_{Expt.}}$',fillstyle='none')
ax.plot(Caldw[::2,0]-273.15,Caldw[::2,1],color='darkorange',linestyle='', marker='^', mew = 2,fillstyle='none')
ax.plot(Mishi[:,0]-273.15,Mishi[:,1],color='darkorange',linestyle='', marker='^', mew = 2,fillstyle='none')
ax.plot(Sotan[:,0]-273.15,Sotan[:,1],color='darkorange',linestyle='', marker='^', mew = 2,fillstyle='none')
ax.plot(tip4p[:,0]-273.15,tip4p[:,1],color='g',linestyle='', marker='o',  mew = 2, label=r'$\mathrm{\rho^{max}_{TIP4P/2005}}$',fillstyle='none')
ax.plot(mbpol[:,0]-273.15,mbpol[:,1]*0.1,color='b',linestyle='', marker='o',  mew = 2, label=r'$\mathrm{\rho^{max}_{DNN@MB-pol}}$',fillstyle='none')

model = np.loadtxt('./PT_data', comments='#')
ax.plot(model[20:,1]-273.15,model[20:,0],color="k",linestyle='-',linewidth="2",marker="",fillstyle='none',label=r'$\mathrm{\rho^{max}_{Theory}}$' )
ax.set_ylabel(r'$\mathrm{\mathbf{p}\ [MPa]}$',{'fontsize':MEDIUM_SIZE})
ax.yaxis.set_label_coords(-0.12, 0.5)
ax.set_xlabel(r'$\mathrm{\mathbf{T}\ [^\circ C]}$',{'fontsize':MEDIUM_SIZE})
ax.set_ylim([-200, 210])
ax.set_yticks(np.arange(-150, 200, 150))
ax.set_xlim([-70, 30])
ax.set_xticks(np.arange(-60, 30, 20))
ax.tick_params(axis='both',direction = 'in', which= 'both',bottom=1, top=1, left=1,right=1)
ax.legend( frameon=False, fontsize = "10", ncol =1,labelspacing=0.1,handlelength=1, loc = "lower left")#, bbox_to_anchor=(0.4, 0.15))

plt.savefig("PT_rho_max.png", dpi =600)
plt.show()
plt.close()


fig = plt.figure(figsize=(3,4))
ax =fig.add_axes([0.18, 0.11, 0.815, 0.885])
P = [2000, 1500, 1000, 500, 1]
for i in range(len(P)):
    model = np.loadtxt('./Press/real_units_%sbar'%(P[i]), comments='#')
    ax.plot(model[:,1],model[:,7]/18,linestyle='-',label=r'$\mathrm{%sMPa}$'%(P[i]*0.1), color=Reds[i+2], linewidth =1.5 )
    Lin = np.loadtxt('Expt_Data_plot/CP/%sbar'%(P[i]), comments='#')
    ax.plot(Lin[:,0]-273.15,Lin[:,2],color=Reds[i+2],linestyle='', marker='.',fillstyle='none',mew =1.5)

ax.set_ylabel(r'$\mathrm{\mathbf{C_{p}}\ [J K^{-1} g{-1}]}$',{'fontsize':MEDIUM_SIZE})
ax.yaxis.set_label_coords(-0.10, 0.5)
ax.set_xlabel(r'$\mathrm{\mathbf{T}[{^\circ} C]}$',{'fontsize':MEDIUM_SIZE})
ax.set_yticks(np.arange(2, 7, 2))
ax.set_ylim([2,7])
ax.set_xlim([-70, 100])
ax.set_xticks(np.arange(-50, 100, 50))
ax.tick_params(axis='both',direction = 'in', which= 'both',bottom=1, top=1, left=1,right=1)

plt.savefig("Cp_isobars.png", dpi =600)
plt.show()
plt.close()

fig = plt.figure(figsize=(4,4))
ax =fig.add_axes([0.13, 0.11, 0.61, 0.885])
P = [2000, 1500, 1000, 500, 1]
for i in range(len(P)):
    if P[i]>1:
       model = np.loadtxt('./Press/real_units_%sbar'%(P[i]), comments='#')
       ax.plot(model[:,1],model[:,7]/18,linestyle='-',label=r'$\mathrm{%sMPa}$'%int(P[i]*0.1), color=Reds[i*3], linewidth =1.5 )
       Lin = np.loadtxt('Expt_Data_plot/CP/%sbar'%(P[i]), comments='#')
       ax.plot(Lin[:,0]-273.15,Lin[:,2],color=Reds[i*3],linestyle='', marker='.',fillstyle='none',mew =1.5)
    else:
       model = np.loadtxt('./Press/real_units_%sbar'%(P[i]), comments='#')
       ax.plot(model[:,1],model[:,7]/18,linestyle='-',label=r'$\mathrm{%sMPa}$'%(P[i]*0.1), color=Reds[i*3], linewidth =1.5 )
       Lin = np.loadtxt('Expt_Data_plot/CP/%sbar'%(P[i]), comments='#')
       ax.plot(Lin[:,0]-273.15,Lin[:,2],color=Reds[i*3],linestyle='', marker='.',fillstyle='none',mew =1.5)

ax.set_ylabel(r'$\mathrm{\mathbf{C_{p}}\ [J K^{-1} g{-1}]}$',{'fontsize':MEDIUM_SIZE})
ax.yaxis.set_label_coords(-0.10, 0.5)
ax.set_xlabel(r'$\mathrm{\mathbf{T}[{^\circ} C]}$',{'fontsize':MEDIUM_SIZE})
ax.set_yticks(np.arange(2, 7, 2))
ax.set_ylim([2,7])
ax.set_xlim([-70, 100])
ax.set_xticks(np.arange(-50, 100, 50))
ax.tick_params(axis='both',direction = 'in', which= 'both',bottom=1, top=1, left=1,right=1)
plt.legend(frameon=False, ncol=1,loc = "center", bbox_to_anchor=(1.23, 0.5), fontsize = "10",handlelength=1)
plt.savefig("Cp_isobars2.png", dpi =600)
plt.show()
plt.close()


fig = plt.figure(figsize=(3,4))
ax =fig.add_axes([0.18, 0.11, 0.815, 0.885])
P = [4000, 3000, 2400, 2000, 1600, 1200, 800, 400, 1]
for i in range(len(P)):
    model = np.loadtxt('./Press/real_units_%sbar'%(P[i]), comments='#')
    ax.plot(model[:,1],model[:,5]*10**5,linestyle='-',label=r'$\mathrm{%sMPa}$'%(P[i]*0.1), color=colors[i+2], linewidth =2 )
    Mishima = np.loadtxt('./Expt_Data_plot/Kappa/%sbar'%(P[i]), comments='#')
    ax.plot(Mishima[:,1]-273.15,Mishima[:,2]*0.1,color=colors[i+2],linestyle='', marker='.', mew =1.5, label=r'$\mathrm{Exp.}$',fillstyle='none')


P = [-200, -400, -600, -800, -1000]
for i in range(len(P)):
    model = np.loadtxt('./Press/real_units_%sbar'%(P[i]), comments='#')
    ax.plot(model[:,1],model[:,5]*10**5,linestyle='-',label=r'$\mathrm{%sMPa}$'%(P[i]*0.1), color=Blues[i+5] )

ax.set_ylabel(r'$\mathrm{\mathbf{\kappa_T} \times 10^5\/[bar^{-1}]}$',{'fontsize':MEDIUM_SIZE})
ax.yaxis.set_label_coords(-0.08, 0.5)
ax.set_xlabel(r'$\mathrm{\mathbf{T}[{^\circ} C]}$',{'fontsize':MEDIUM_SIZE})
ax.set_yticks(np.arange(2, 10, 2))
ax.set_ylim([1.5,10])
ax.set_xlim([-70, 100])
ax.set_xticks(np.arange(-50, 100, 50))
ax.tick_params(axis='both',direction = 'in', which= 'both',bottom=1, top=1, left=1,right=1)
plt.show()
plt.savefig("Kappa_isobars.png", dpi =600)
plt.close()

fig = plt.figure(figsize=(4,4))
ax =fig.add_axes([0.13, 0.11, 0.61, 0.885])
P = [4000, 3000, 2400, 2000, 1600, 1200, 800, 400, 1]
P = [4000, 3000, 2000, 800, 1]
for i in range(len(P)):
    if P[i]>1:
       model = np.loadtxt('./Press/real_units_%sbar'%(P[i]), comments='#')
       ax.plot(model[:,1],model[:,5]*10**5,linestyle='-',label=r'$\mathrm{%sMPa}$'%int(P[i]*0.1), color=colors[i*3], linewidth =2 )
       Mishima = np.loadtxt('./Expt_Data_plot/Kappa/%sbar'%(P[i]), comments='#')
       ax.plot(Mishima[:,1]-273.15,Mishima[:,2]*0.1,color=colors[i*3],linestyle='', marker='.', mew =1.5, fillstyle='none')
    else:
       model = np.loadtxt('./Press/real_units_%sbar'%(P[i]), comments='#')
       ax.plot(model[:,1],model[:,5]*10**5,linestyle='-',label=r'$\mathrm{%sMPa}$'%(P[i]*0.1), color=colors[i*3], linewidth =2 )
       Mishima = np.loadtxt('./Expt_Data_plot/Kappa/%sbar'%(P[i]), comments='#')
       ax.plot(Mishima[:,1]-273.15,Mishima[:,2]*0.1,color=colors[i*3],linestyle='', marker='.', mew =1.5, fillstyle='none')



ax.set_ylabel(r'$\mathrm{\mathbf{\kappa_T} \times 10^5\/[bar^{-1}]}$',{'fontsize':MEDIUM_SIZE})
ax.yaxis.set_label_coords(-0.11, 0.5)
ax.set_xlabel(r'$\mathrm{\mathbf{T}[{^\circ} C]}$',{'fontsize':MEDIUM_SIZE})
ax.set_yticks(np.arange(2, 13, 3))
ax.set_ylim([1.5,13])
ax.set_xlim([-70, 100])
ax.set_xticks(np.arange(-50, 100, 50))
ax.tick_params(axis='both',direction = 'in', which= 'both',bottom=1, top=1, left=1,right=1)
plt.legend(frameon=False, ncol=1,loc = "center", bbox_to_anchor=(1.23, 0.5), fontsize = "10",handlelength=1)
plt.savefig("Kappa_isobars2.png", dpi =600)
plt.show()
plt.close()



#Alpha
fig = plt.figure(figsize=(3,4))
ax =fig.add_axes([0.18, 0.11, 0.815, 0.885])
P = [400, 300, 200, 160, 140, 100, 80, 40, 0.1]
for i in range(len(P)):
    exp = np.loadtxt('./Expt_Data_plot/Alpha/alpha_%sMPa'%(P[i]), comments='#')
    ax.plot(exp[:,1]-273.15,exp[:,2]*10,linestyle='',marker='.',label=r'$\mathrm{%sMPa}$'%(P[i]), color=Reds[i+2],fillstyle='none')

    model = np.loadtxt('./Press/real_units_%sbar'%(int(P[i]*10)), comments='#')
    ax.plot(model[:,1],model[:,6]*10**4,linestyle='-',label=r'$\mathrm{%sMPa}$'%(P[i]*0.1), color=Reds[i+2] )
    if P[i] == 0.1:
        exp = np.loadtxt('./Expt_Data_plot/Alpha/alpha_%sMPa'%(P[i]), comments='#')
ax.set_ylabel(r'$\mathrm{\mathbf{\alpha_P} \times 10^4\/[K^{-1}]}$',{'fontsize':MEDIUM_SIZE})
ax.yaxis.set_label_coords(-0.08, 0.5)
ax.set_xlabel(r'$\mathrm{\mathbf{T}[{^\circ} C]}$',{'fontsize':MEDIUM_SIZE})
ax.set_yticks(np.arange(-4, 8, 4))
ax.set_ylim([-4,8])
ax.set_xlim([-70, 100])
ax.set_xticks(np.arange(-50, 100, 50))
ax.tick_params(axis='both',direction = 'in', which= 'both',bottom=1, top=1, left=1,right=1)

plt.savefig("Alpha_isobars.png", dpi =600)
plt.show()
plt.close()

fig = plt.figure(figsize=(4,4))
ax =fig.add_axes([0.13, 0.11, 0.61, 0.885])
P = [160, 100, 80, 60, 40, 0.1]
for i in range(len(P)):
    if P[i] == 0.1:
        exp = np.loadtxt('./Expt_Data_plot/Alpha/alpha_%sMPa'%(P[i]), comments='#')
    elif P[i]<=100:
        exp = np.loadtxt('./Expt_Data_plot/Alpha/alpha_%sMPa'%(P[i]), comments='#')
        ax.plot(exp[:,0],exp[:,2]*0.01,linestyle='',marker='.', color=Reds[i*3],fillstyle='none')
    else: 
        exp = np.loadtxt('./Expt_Data_plot/Alpha/alpha_%sMPa'%(P[i]), comments='#')
        ax.plot(exp[:,1]-273.15,exp[:,2]*10,linestyle='',marker='.', color=Reds[i*3],fillstyle='none')

    model = np.loadtxt('./Press/real_units_%sbar'%(int(P[i]*10)), comments='#')
    ax.plot(model[:,1],model[:,6]*10**4,linestyle='-',label=r'$\mathrm{%sMPa}$'%(P[i]), color=Reds[i*3] )
    #if P[i] == 0.1:
    #    exp = np.loadtxt('./Expt_Data_plot/Alpha/alpha_%sMPa'%(P[i]), comments='#')
ax.set_ylabel(r'$\mathrm{\mathbf{\alpha_P} \times 10^4\/[K^{-1}]}$',{'fontsize':MEDIUM_SIZE})
ax.yaxis.set_label_coords(-0.11, 0.5)
ax.set_xlabel(r'$\mathrm{\mathbf{T}[{^\circ} C]}$',{'fontsize':MEDIUM_SIZE})
ax.set_yticks(np.arange(-4, 8, 4))
ax.set_ylim([-4,8])
ax.set_xlim([-70, 100])
ax.set_xticks(np.arange(-50, 100, 50))
ax.tick_params(axis='both',direction = 'in', which= 'both',bottom=1, top=1, left=1,right=1)
plt.legend(frameon=False, ncol=1,loc = "center", bbox_to_anchor=(1.23, 0.5), fontsize = "10",handlelength=1)
plt.savefig("Alpha_isobars2.png", dpi =600)
plt.show()
plt.close()

