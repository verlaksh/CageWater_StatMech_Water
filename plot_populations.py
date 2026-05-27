import scipy
import numpy as np
from pylab import *
from scipy.optimize import curve_fit
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)



SMALL_SIZE = 12
MEDIUM_SIZE = 14
BIGGER_SIZE = 22
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
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


#with tip4p

model = np.loadtxt('./Atm_P/real_units', comments='#')

fig = plt.figure(figsize=(3.5,8))


ax1 = fig.add_axes([0.22, 0.820, 0.775, 0.177], xticklabels=[])#, ylim=(0, 0.6))
ax2 = fig.add_axes([0.22, 0.640, 0.775, 0.177], xticklabels=[])#, ylim=(0, 0.6))
ax3 = fig.add_axes([0.22, 0.460, 0.775, 0.177], xticklabels=[])#, ylim=(0, 0.6))
ax4 = fig.add_axes([0.22, 0.280, 0.775, 0.177], xticklabels=[])#, ylim=(0, 0.6))
ax5 = fig.add_axes([0.22, 0.070, 0.775, 0.207])#, ylim=(0, 0.6))

exp = np.loadtxt('density', comments='#')
#dens=697.15223/model[:,4]
tip4p=np.loadtxt('./MD_Model_data/tip4p2005_rho', comments='#')
tip3p=np.loadtxt('./MD_Model_data/tip3p_rho', comments='#')

ax1.plot(exp[::5,0],exp[::5,2],color='darkorange',linestyle='', marker='^', label=r'$\mathrm{Exp.}$',fillstyle='none')
ax1.plot(tip4p[:,0]-273.15,tip4p[:,1],color='green',linestyle='', marker='o', label=r'$\mathrm{TIP4P/2005}$',fillstyle='none')
ax1.plot(tip3p[:,0],tip3p[:,1],color='blue',linestyle='', marker='o', label=r'$\mathrm{Tip3p}$',fillstyle='none')
ax1.plot(model[:,1],model[:,4]*0.001 ,color='k',linestyle='-',label=r'$\mathrm{Model}$' )
ax1.set_ylabel(r'$\mathrm{\mathbf{\rho}\/[kg/m^3]}$',{'fontsize':MEDIUM_SIZE})
ax1.yaxis.set_label_coords(-0.15, 0.5)
ax1.set_yticks(np.arange(0.9, 1.05, 0.1))
ax1.set_ylim([0.9,1.1])
ax1.set_xlim([-60, 110])
ax1.set_xticks(np.arange(-50, 110, 50))
ax1.tick_params(axis='both',direction = 'in',bottom=1, top=1, left=1, right=1)

tip4p=np.loadtxt('./MD_Model_data/tip4p2005_alp', comments='#')
tip3p=np.loadtxt('./MD_Model_data/tip3p_alp', comments='#')
ax2.plot(exp[::5,0],exp[::5,3]*0.001,color='darkorange',linestyle='', marker='^', label=r'$\mathrm{Exp.}$',fillstyle='none')
ax2.plot(tip4p[:,0]-273.15,tip4p[:,1]*0.01,color='green',linestyle='', marker='o', label=r'$\mathrm{TIP4P/2005}$',fillstyle='none')
ax2.plot(tip3p[:,0],tip3p[:,1]*0.01,color='blue',linestyle='', marker='o', label=r'$\mathrm{Tip3p}$',fillstyle='none')
ax2.plot(model[:,1],model[:,6]*10**3,color='k',linestyle='-',label=r'$\mathrm{Theory}$' )
ax2.set_ylabel(r'$\mathrm{\mathbf{\alpha_P}10^3\/[K^{-1}]}$',{'fontsize':MEDIUM_SIZE})
ax2.yaxis.set_label_coords(-0.15, 0.5)
ax2.set_yticks(np.arange(-1, 2, 1))
ax2.set_ylim([-2,2])
ax2.set_xlim([-60, 110])
ax2.set_xticks(np.arange(-50, 110, 50))
ax2.tick_params(axis='both',direction = 'in',bottom=1, top=1, left=1, right=1)

tip4p=np.loadtxt('./MD_Model_data/tip4p2005_kap', comments='#')
tip3p=np.loadtxt('./MD_Model_data/tip3p_kap', comments='#')
kim = np.loadtxt('./Expt_Data/Kim_Science_2017', comments='#')

ax3.plot(kim[:,0]-273.15,kim[:,1]*0.1,color='darkorange',linestyle='', marker='^',fillstyle='none',mew = 1.5)
ax3.plot(exp[::5,0],exp[::5,4]*0.1,color='darkorange',linestyle='', marker='^', label=r'$\mathrm{Exp.}$',fillstyle='none')
ax3.plot(tip4p[:,0]-273.15,tip4p[:,1]*10**5,color='green',linestyle='', marker='o', label=r'$\mathrm{TIP4P/2005}$',fillstyle='none')
ax3.plot(tip3p[:,0],tip3p[:,1]*0.1,color='blue',linestyle='', marker='o', label=r'$\mathrm{Tip3p}$',fillstyle='none')
ax3.plot(model[:,1],model[:,5]*10**5,color='k',linestyle='-',label=r'$\mathrm{Theory}$' )
ax3.set_ylabel(r'$\mathrm{\mathbf{\kappa_T}10^5\/[bar^{-1}]}$',{'fontsize':MEDIUM_SIZE})
ax3.yaxis.set_label_coords(-0.15, 0.5)
ax3.set_yticks(np.arange(1, 10, 3))
ax3.set_ylim([0,11])
ax3.set_xlim([-60, 110])
ax3.set_xticks(np.arange(-50, 110, 50))
ax3.tick_params(axis='both',direction = 'in',bottom=1, top=1, left=1, right=1)


tip4p=np.loadtxt('./MD_Model_data/tip4p2005_cp', comments='#')
tip3p=np.loadtxt('./MD_Model_data/tip3p_cp', comments='#')
exp = np.loadtxt('exp_cp', comments='#')
pathak= np.loadtxt('./Expt_Data/Pathak_PNAS_2021', comments='#')

ax4.plot(pathak[:,0]-273.15,pathak[:,1]/18,color='darkorange',linestyle='', marker='^',fillstyle='none',mew = 1.5)
ax4.plot(exp[::,0],exp[::,1]/18,color='darkorange',linestyle='', marker='^', label=r'$\mathrm{Exp.}$',fillstyle='none')
ax4.plot(tip4p[:,0]-273.15,tip4p[:,1]*4.18/18,color='green',linestyle='', marker='o', label=r'$\mathrm{TIP4P/2005}$',fillstyle='none')
ax4.plot(tip3p[:,0],tip3p[:,1]*4.18/18,color='blue',linestyle='', marker='o', label=r'$\mathrm{TIP3P}$',fillstyle='none')
ax4.plot(model[:,1],model[:,7]/18,color='k',linestyle='-',label=r'$\mathrm{Theory}$' )
ax4.set_ylabel(r'$\mathrm{\mathbf{C_{p}}\/[J/ K g]}$',{'fontsize':MEDIUM_SIZE})
ax4.yaxis.set_label_coords(-0.15, 0.5)
ax4.set_xlabel(r'$\mathrm{\mathbf{T}\/[{^\circ} C]}$',{'fontsize':MEDIUM_SIZE})
ax4.set_yticks(np.arange(4,7,2))
ax4.set_ylim([3,7])
ax4.set_xlim([-60, 110])
ax4.set_xticks(np.arange(-50, 110, 50))
ax4.tick_params(axis='both',direction = 'in',bottom=1, top=1, left=1, right=1)
ax4.legend(frameon=False, ncol=2, loc="center", bbox_to_anchor=(0.55, 0.8), fontsize = "9", handlelength=1.5)

model = np.loadtxt('./Atm_P/population_real_units', comments='#')
ax5.plot(model[:,1],model[:,4],color=Plasma[0],linestyle='-',label=r'$\mathrm{Cage}$' )
ax5.plot(model[:,1],model[:,3],color=Plasma[3],linestyle='-',label=r'$\mathrm{pHB}$' )
ax5.plot(model[:,1],model[:,5],color=Plasma[5],linestyle='-',label=r'$\mathrm{vdW}$' )
ax5.plot(model[:,1],model[:,6],color=Plasma[7],linestyle='-',label=r'$\mathrm{NI}$' )
ax5.set_ylabel(r'$\mathrm{\mathbf{f_i}}$',{'fontsize':MEDIUM_SIZE})
ax5.yaxis.set_label_coords(-0.15, 0.5)
ax5.set_xlabel(r'$\mathrm{\mathbf{T}\/[{^\circ} C]}$',{'fontsize':MEDIUM_SIZE})
ax5.set_ylim([0, 1])
ax5.set_yticks(np.arange(0, 1, 0.5))
ax5.set_xlim([-60, 110])
ax5.set_xticks(np.arange(-50, 110, 50))
ax5.tick_params(axis='both',direction = 'in',bottom=1, top=1, left=1, right=1)
ax5.legend(frameon=False, loc="center", bbox_to_anchor=(0.75, 0.5), fontsize = "9", handlelength=1.5)

plt.savefig("resp_funct_temp_pop_tips.png", dpi =600)
plt.show()
plt.close()

fig = plt.figure(figsize=(4,6))

ax =fig.add_axes([0.1, 0.14, 0.8, 0.78],frameon=False,xticklabels=[],yticklabels=[])
ax.tick_params(axis='both',direction = 'in',bottom=False, top=False, left=False, right=False)

ax1 = fig.add_axes([0.15, 0.770, 0.845, 0.225], xticklabels=[])#, ylim=(0, 0.6))
ax2 = fig.add_axes([0.15, 0.540, 0.845, 0.225], xticklabels=[])#, ylim=(0, 0.6))
ax3 = fig.add_axes([0.15, 0.310, 0.845, 0.225], xticklabels=[])#, ylim=(0, 0.6))
ax4 = fig.add_axes([0.15, 0.080, 0.845, 0.225])#, ylim=(0, 0.6))

cmap = cm.get_cmap('Reds_r', 15)    # PiYG
Reds=["" for x in range(15)]

for i in range(cmap.N):
        rgb = cmap(i)[:3] # will return rgba, we take only first 3 so we get rgb
        Reds[i]=matplotlib.colors.rgb2hex(rgb)

cmap = cm.get_cmap('Blues', 12)    # PiYG
Blues=["" for x in range(12)]

for i in range(cmap.N):
        rgb = cmap(i)[:3] # will return rgba, we take only first 3 so we get rgb
        Blues[i]=matplotlib.colors.rgb2hex(rgb)


P = [300, 240, 200, 160, 120, 80, 40, 0.1]
for i in range(len(P)):
    model = np.loadtxt('./Press/population_real_units_%sbar'%(int(P[i]*10)), comments='#')
    ax1.plot(model[:,1],model[:,4],linestyle='-',label=r'$\mathrm{%sMPa}$'%P[i], color = Reds[i+2] )
    ax2.plot(model[:,1],model[:,3],linestyle='-',label=r'$\mathrm{%sMPa}$'%P[i], color = Reds[i+2] )
    ax3.plot(model[:,1],model[:,5],linestyle='-',label=r'$\mathrm{%sMPa}$'%P[i], color = Reds[i+2] )
    ax4.plot(model[:,1],model[:,8],linestyle='-',label=r'$\mathrm{%sMPa}$'%P[i], color = Reds[i+2] )

P = [20, 40, 60, 80, 100]
for i in range(len(P)):
    model = np.loadtxt('./Press/population_real_units_-%sbar'%(P[i]*10), comments='#')
    ax1.plot(model[:,1],model[:,4],linestyle='-',label=r'$\mathrm{-%sMPa}$'%P[i], color = Blues[i+5] )
    ax2.plot(model[:,1],model[:,3],linestyle='-',label=r'$\mathrm{-%sMPa}$'%P[i], color = Blues[i+5] )
    ax3.plot(model[:,1],model[:,5],linestyle='-',label=r'$\mathrm{-%sMPa}$'%P[i], color = Blues[i+5] )
    ax4.plot(model[:,1],model[:,8],linestyle='-',label=r'$\mathrm{-%sMPa}$'%P[i], color = Blues[i+5] )


ax1.set_ylabel(r'$\mathrm{\mathbf{f_{cage}}}$',{'fontsize':MEDIUM_SIZE})
ax1.set_yticks(np.arange(0, 1, 0.5))
ax1.set_ylim([0, 1])
ax1.set_xlim([-200, 100])
ax1.set_xticks(np.arange(-200, 100, 100))
ax1.tick_params(axis='both',direction = 'in',bottom=1, top=1, left=1, right=1)

ax2.set_ylabel(r'$\mathrm{\mathbf{f_{pHB}}}$',{'fontsize':MEDIUM_SIZE})
ax2.set_yticks(np.arange(0, 1, 0.5))
ax2.set_ylim([0, 1])
ax2.set_xlim([-200, 100])
ax2.set_xticks(np.arange(-200, 100, 100))
ax2.tick_params(axis='both',direction = 'in',bottom=1, top=1, left=1, right=1)

ax3.set_ylabel(r'$\mathrm{\mathbf{f_{vdW}}}$',{'fontsize':MEDIUM_SIZE})
ax3.set_yticks(np.arange(0, 0.2, 0.1))
ax3.set_ylim([0, 0.2])
ax3.set_xlim([-200, 100])
ax3.set_xticks(np.arange(-200, 100, 100))
ax3.tick_params(axis='both',direction = 'in',bottom=1, top=1, left=1, right=1)

#tip4p=np.loadtxt('./MD_Model_data/nHB_tip4p2005_100MPa', comments='#')
ax4.set_ylabel(r'$\mathrm{\mathbf{n_{HB}}}$',{'fontsize':MEDIUM_SIZE})
ax4.set_xlabel(r'$\mathrm{\mathbf{T}\/[{^\circ} C]}$',{'fontsize':MEDIUM_SIZE})
ax4.set_ylim([3, 4])
ax4.set_yticks(np.arange(3, 4, 0.5))
ax4.set_xlim([-200, 100])
ax4.set_xticks(np.arange(-200, 100, 100))
ax4.tick_params(axis='both',direction = 'in',bottom=1, top=1, left=1, right=1)
ax4.legend(frameon=False, ncol=3, fontsize ="9", labelspacing=0.1, handlelength=0.5,loc = "center", bbox_to_anchor=(0.41, 0.3))

plt.savefig("multi_population_and_HBf.png", dpi =600)
plt.show()
plt.close()


cmap = cm.get_cmap('Reds_r', 15)    # PiYG
Reds=["" for x in range(15)]

for i in range(cmap.N):
        rgb = cmap(i)[:3] # will return rgba, we take only first 3 so we get rgb
        Reds[i]=matplotlib.colors.rgb2hex(rgb)

cmap = cm.get_cmap('Blues', 15)    # PiYG
Blues=["" for x in range(15)]

for i in range(cmap.N):
        rgb = cmap(i)[:3] # will return rgba, we take only first 3 so we get rgb
        Blues[i]=matplotlib.colors.rgb2hex(rgb)


fig = plt.figure(figsize=(5,4))
ax =fig.add_axes([0.08, 0.14, 0.8, 0.78],frameon=False,xticklabels=[],yticklabels=[])
ax.tick_params(axis='both',direction = 'in',bottom=False, top=False, left=False, right=False)
ax.set_ylabel(r'$\mathrm{\mathbf{f_{j}}}$',{'fontsize':MEDIUM_SIZE})

ax1 = fig.add_axes([0.12, 0.56, 0.875, 0.435], xticklabels=[])#, ylim=(0, 0.6))
ax2 = fig.add_axes([0.12, 0.12, 0.875, 0.435])#, ylim=(0, 0.6))


P = [200, 160, 120, 80, 40, 0.1]
for i in range(len(P)):
    model = np.loadtxt('./Press/population_real_units_%sbar'%(int(P[i]*10)), comments='#')
    ax1.plot(model[:,1],model[:,4],linestyle='-',label=r'$\mathrm{%sMPa}$'%P[i], color = Reds[2*i] )
    ax1.plot(model[:,1],model[:,3],linestyle='--', color = Reds[2*i] )
    ax2.plot(model[:,1],model[:,4],linestyle='-',  color = Reds[2*i] )
    ax2.plot(model[:,1],model[:,5],linestyle=':',  color = Reds[2*i] )
j=i+1
P = [20, 40, 60, 80, 100]
for i in range(len(P)):
    model = np.loadtxt('./Press/population_real_units_-%sbar'%(P[i]*10), comments='#')
    ax1.plot(model[:,1],model[:,4],linestyle='-',label=r'$\mathrm{-%sMPa}$'%P[i], color = Blues[2*(i+2)] )
    ax1.plot(model[:,1],model[:,3],linestyle='--', color = Blues[2*(i+2)] )
    ax2.plot(model[:,1],model[:,4],linestyle='-',  color = Blues[2*(i+2)] )
    ax2.plot(model[:,1],model[:,5],linestyle=':',  color = Blues[2*(i+2)] )

ax1.set_yticks(np.arange(0.4, 0.7, 0.2))
ax1.set_ylim([0.35, 0.65])
ax1.set_xlim([-100, 50])
ax1.set_xticks(np.arange(-100, 50, 50))
ax1.tick_params(axis='both',direction = 'in',bottom=1, top=1, left=1, right=1)
ax1.legend(frameon=False, ncol=1, fontsize ="9", labelspacing=0.1,handlelength=1)#, loc="upper right")

ax2.set_xlabel(r'$\mathrm{\mathbf{T}\/[{^\circ} C]}$',{'fontsize':MEDIUM_SIZE})
ax2.set_yticks(np.arange(0.07, 0.1, 0.02))
ax2.set_ylim([0.07, 0.1])
ax2.set_xlim([-100, 50])
ax2.set_xticks(np.arange(-100, 50, 50))
ax2.tick_params(axis='both',direction = 'in',bottom=1, top=1, left=1, right=1)

plt.savefig("crossovers.png", dpi =600)
plt.show()
plt.close()

