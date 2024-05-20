#%%
"""
VSM 
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import os

def lineal(x,m,n):
    return m*x+n
#%% LAURICO S/ NPM

C = np.loadtxt('Laurico1.txt',skiprows=12)
mL=0.2874 #g
f = C[:,0]
g = C[:,1]/mL

f1=f[np.nonzero(f>=2500)]
g1=g[np.nonzero(f>=2500)]
f2=f[np.nonzero(f<=-2500)]
g2=g[np.nonzero(f<=-2500)]

(m1,n1), pcov= curve_fit(lineal,f1,g1)
(m2,n2), pcov= curve_fit(lineal,f2,g2)


chi_mass_laurico = np.mean([m1,m2])
g_ajustado = lineal(f,chi_mass_laurico,0)

plt.plot(f,g,'.-')
plt.plot(f1,g1,'.-')
plt.plot(f2,g2,'.-')
plt.plot(f,g_ajustado,'-',c='tab:red',label=f'$\chi$ = {chi_mass_laurico:.2e} emu/gG')
plt.legend()
plt.grid()
plt.xlabel('H (G)')
plt.ylabel('m (emu/g)')
plt.title('Capsula ac. laurico')

#%% Muestra RANDOM
files=os.listdir(os.path.join(os.getcwd(),'Random'))
files.sort()
masa_muestraRd=0.20945 #g
masa_npmRd=0.0168 #g
masa_lauricoRd=masa_muestraRd-masa_npmRd 

campos=[]
m_correg_norm=[]
fname=[]
Ms_Rd=[]

H_aux = np.linspace(-50,50,1000)
campo_lineal=[]
m_interp=[]

pendientes_Rd=[]
err_pend_Rd=[]
ordenadas_Rd=[]
Mr=[]
Hc=[]
for f in files:
    B = np.loadtxt(os.path.join(os.getcwd(),'Random',f),skiprows=12)
    H = B[:,0]
    m = B[:,1]
    #calculo la contribucion diamagnetica del ac laurico y se la descuento
    contrib_diamagRd = chi_mass_laurico*masa_lauricoRd*H
    m_correg = (m-contrib_diamagRd)
    m_correg_norm_masa=m_correg/masa_npmRd
    
    interpolador=interp1d(H,m_correg_norm_masa,fill_value='extrapolate')
    m_interp.append(interpolador(H_aux))
    
    m_norm=m_correg_norm_masa/max(m_correg_norm_masa) #Normalizo por valor maximo
    
    m_correg_norm.append(m_norm)
    campos.append(np.array(H))
    fname.append(f.split('_')[-1].split('.')[0])
    Ms_Rd.append(max(m_correg_norm_masa))
    
    (chi,n),pcov=curve_fit(lineal, H[np.nonzero(H<10)],m_norm[np.nonzero(H<10)]) #Ajuste lineal 
    err_m=pcov[0][0]
    pendientes_Rd.append(chi)
    err_pend_Rd.append(err_m)
    ordenadas_Rd.append(n)
    
    H_aux = np.linspace(-10,10,5000)
    H_recortado = H[np.nonzero(abs(H)<10)]
    m_recortado = m_correg_norm_masa[np.nonzero(abs(H)<10)]
    
    interpolador2=interp1d(H_recortado,m_recortado,fill_value='extrapolate')
    #m_new = interpolador2(H_aux)
    m_new = lineal(H_aux,chi,n)
    indx_H = np.nonzero(H_aux>=0)[0][0]
    indx_M = np.nonzero(m_new>0)[0][0]
    Hc.append(-H_aux[indx_M])
    Mr.append(m_new[indx_H])
    print('*'*50)
    print(f.split('_')[0], f.split('_')[-1])
    print('Susceptibilidad =',f'{chi:.3e}','+/-',f'{err_m:.3e}')
    print(f'Mag Remanente = {m_new[indx_H]:.3e}')
    print(f'Campo Coercitivo = {H_aux[indx_M]:.3e} A/m')
    fig,ax=plt.subplots(constrained_layout=True)
    ax.plot(H,m/masa_npmRd,'.-',label='V')
    ax.plot(H,m_correg_norm_masa,'.-',label='V (s/ laurico)')
    ax.legend(ncol=2,loc='lower right')
    ax.grid()
    ax.set_xlabel('H (G)')
    ax.set_ylabel('m (emu/g)')
    ax.set_title(f)
    
    axins = ax.inset_axes([0.3, 0.2, 0.69, 0.55])
    axins.axvline(0,0,1,c='k',lw=0.8)
    axins.axhline(0,0,1,c='k',lw=0.8)
    axins.plot(H,m_norm,'.-',label='m norm')
    axins.plot(H_aux,lineal(H_aux,chi,n),'r-',label=f'$\chi$ = {chi:.3e}\nHc = {H_aux[indx_M]:.3e} A/m')
    axins.set_xlabel('H (G)')
    axins.grid()
    axins.legend(loc='lower right')
    # axins.legend()
    axins.set_xlim(-2,11)
    axins.set_ylim(-0.01,.06)
    ax.indicate_inset_zoom(axins, edgecolor="black")
    plt.savefig('VSM_'+f[:-4]+'.png',dpi=300)
    plt.show()

#%% Vertical Rotado

files=os.listdir(os.path.join(os.getcwd(),'V_rotado'))
files.sort()
masa_muestraV=0.31170 #g
masa_npmV=0.020 #g
masa_lauricoV=masa_muestraV-masa_npmV 

campos=[]
m_correg_norm=[]
fname=[]
Ms_V=[]

H_aux = np.linspace(-50,50,1000)
campo_lineal=[]
m_interp=[]

pendientes_V=[]
err_pend_V=[]
ordenadas_V=[]
Mr=[]
Hc=[]
for f in files:
    B = np.loadtxt(os.path.join(os.getcwd(),'V_rotado',f),skiprows=12)
    H = B[:,0]
    m = B[:,1]
    #calculo la contribucion diamagnetica del ac laurico y se la descuento
    contrib_diamagV = chi_mass_laurico*masa_lauricoV*H
    m_correg = (m-contrib_diamagV)
    m_correg_norm_masa=m_correg/masa_npmV
    
    interpolador=interp1d(H,m_correg_norm_masa,fill_value='extrapolate')
    m_interp.append(interpolador(H_aux))
    
    m_norm=m_correg_norm_masa/max(m_correg_norm_masa) #Normalizo por valor maximo
    
    m_correg_norm.append(m_norm)
    campos.append(np.array(H))
    fname.append(f.split('_')[-1].split('.')[0])
    Ms_V.append(max(m_correg_norm_masa))
    
    (chi,n),pcov=curve_fit(lineal, H[np.nonzero(H<10)],m_norm[np.nonzero(H<10)]) #Ajuste lineal 
    err_m=pcov[0][0]
    pendientes_V.append(chi)
    err_pend_V.append(err_m)
    ordenadas_V.append(n)
    
    H_aux = np.linspace(-10,10,5000)
    H_recortado = H[np.nonzero(abs(H)<10)]
    m_recortado = m_correg_norm_masa[np.nonzero(abs(H)<10)]
    
    interpolador2=interp1d(H_recortado,m_recortado,fill_value='extrapolate')
    #m_new = interpolador2(H_aux)
    m_new = lineal(H_aux,chi,n)
    indx_H = np.nonzero(H_aux>=0)[0][0]
    indx_M = np.nonzero(m_new>0)[0][0]
    Hc.append(-H_aux[indx_M])
    Mr.append(m_new[indx_H])
    print('*'*50)
    print(f.split('_')[0], f.split('_')[-1])
    print('Susceptibilidad =',f'{chi:.3e}','+/-',f'{err_m:.3e}')
    print(f'Mag Remanente = {m_new[indx_H]:.3e}')
    print(f'Campo Coercitivo = {H_aux[indx_M]:.3e} A/m')

    fig,ax=plt.subplots(constrained_layout=True)
    ax.plot(H,m/masa_npmV,'.-',label='V')
    ax.plot(H,m_correg_norm_masa,'.-',label='V (s/ laurico)')
    ax.legend(ncol=2,loc='lower right')
    ax.grid()
    ax.set_xlabel('H (G)')
    ax.set_ylabel('m (emu/g)')
    ax.set_title(f)
    
    axins = ax.inset_axes([0.3, 0.2, 0.69, 0.55])
    axins.axvline(0,0,1,c='k',lw=0.8)
    axins.axhline(0,0,1,c='k',lw=0.8)
    axins.plot(H,m_norm,'.-',label='m norm')
    axins.plot(H_aux,lineal(H_aux,chi,n),'r-',label=f'$\chi$ = {chi:.3e}\nHc = {H_aux[indx_M]:.3e} A/m')
    axins.set_xlabel('H (G)')
    axins.grid()
    axins.legend(loc='lower right')
    # axins.legend()
    axins.set_xlim(-2,11)
    axins.set_ylim(-0.01,.06)
    ax.indicate_inset_zoom(axins, edgecolor="black")
    plt.savefig('VSM_'+f[:-4]+'.png',dpi=300)
    plt.show()

#%% Vertical Rotado 2

files=os.listdir(os.path.join(os.getcwd(),'Vertical_2'))
files.sort()
masa_muestraV2=0.501 #g
masa_npmV2=0.032 #g
masa_lauricoV2=masa_muestraV2-masa_npmV2 

campos=[]
m_correg_norm=[]
fname=[]
Ms_V2=[]

H_aux = np.linspace(-50,50,1000)
campo_lineal=[]
m_interp=[]

pendientes_V2=[]
err_pend_V2=[]
ordenadas_V2=[]
Mr=[]
Hc=[]
for f in files:
    B = np.loadtxt(os.path.join(os.getcwd(),'Vertical_2',f),skiprows=12)
    H = B[:,0]
    m = B[:,1]
    #calculo la contribucion diamagnetica del ac laurico y se la descuento
    contrib_diamagV2 = chi_mass_laurico*masa_lauricoV2*H
    m_correg = (m-contrib_diamagV2)
    m_correg_norm_masa=m_correg/masa_npmV2
    
    interpolador=interp1d(H,m_correg_norm_masa,fill_value='extrapolate')
    m_interp.append(interpolador(H_aux))
    
    m_norm=m_correg_norm_masa/max(m_correg_norm_masa) #Normalizo por valor maximo
    
    m_correg_norm.append(m_norm)
    campos.append(np.array(H))
    fname.append(f.split('_')[-1].split('.')[0])
    Ms_V2.append(max(m_correg_norm_masa))
    
    (chi,n),pcov=curve_fit(lineal, H[np.nonzero(H<10)],m_norm[np.nonzero(H<10)]) #Ajuste lineal 
    err_m=pcov[0][0]
    pendientes_V2.append(chi)
    err_pend_V2.append(err_m)
    ordenadas_V2.append(n)
    
    H_aux = np.linspace(-10,10,5000)
    H_recortado = H[np.nonzero(abs(H)<10)]
    m_recortado = m_correg_norm_masa[np.nonzero(abs(H)<10)]
    
    interpolador2=interp1d(H_recortado,m_recortado,fill_value='extrapolate')
    #m_new = interpolador2(H_aux)
    m_new = lineal(H_aux,chi,n)
    indx_H = np.nonzero(H_aux>=0)[0][0]
    indx_M = np.nonzero(m_new>0)[0][0]
    Hc.append(-H_aux[indx_M])
    Mr.append(m_new[indx_H])
    print('*'*50)
    print(f.split('_')[0], f.split('_')[-1])
    print('Susceptibilidad =',f'{chi:.3e}','+/-',f'{err_m:.3e}')
    print(f'Mag Remanente = {m_new[indx_H]:.3e}')
    print(f'Campo Coercitivo = {H_aux[indx_M]:.3e} A/m')

    fig,ax=plt.subplots(constrained_layout=True)
    ax.plot(H,m/masa_npmV2,'.-',label='V')
    ax.plot(H,m_correg_norm_masa,'.-',label='V (s/ laurico)')
    ax.legend(ncol=2,loc='lower right')
    ax.grid()
    ax.set_xlabel('H (G)')
    ax.set_ylabel('m (emu/g)')
    ax.set_title(f)
    
    axins = ax.inset_axes([0.3, 0.2, 0.69, 0.55])
    axins.axvline(0,0,1,c='k',lw=0.8)
    axins.axhline(0,0,1,c='k',lw=0.8)
    axins.plot(H,m_norm,'.-',label='m norm')
    axins.plot(H_aux,lineal(H_aux,chi,n),'r-',label=f'$\chi$ = {chi:.3e}\nHc = {H_aux[indx_M]:.3e} A/m')
    axins.set_xlabel('H (G)')
    axins.grid()
    axins.legend(loc='lower right')
    # axins.legend()
    axins.set_xlim(-2,11)
    axins.set_ylim(-0.01,.06)
    ax.indicate_inset_zoom(axins, edgecolor="black")
    plt.savefig('VSM_'+f[:-4]+'.png',dpi=300)
    plt.show()

#%%SUCEPTIBILIDAD vs ANGULO VERTICAL 
#plt.close('all')
rot_all_V = np.array([0,15,30,-15,-30,-45,-60,-75,0,0,-75,-90,-105,-120,-135,-150,-165,-180,-195,-210,-225,-240,-255])
rot_all_V2 =  np.arange(0,375,15)
rot_all_V=rot_all_V[0]-rot_all_V + 90
fig,ax=plt.subplots(nrows=1,figsize=(10,4),sharex=True,constrained_layout=True)
ax.errorbar(x=rot_all_V,y=pendientes_V,yerr=err_pend_V, fmt='.', capsize=5,label='$\chi_V$$_1$')
ax.errorbar(x=rot_all_V2,y=pendientes_V2,yerr=err_pend_V2, fmt='.', capsize=5,label='$\chi_V$$_2$')
ax.grid()
ax.legend()
ax.set_ylabel('$\chi$')
ax.set_title('Susceptibilidad vs angulo')

xticks_values = np.arange(0,375,15)
xticks_labels = [str(i) for i in xticks_values]
plt.xticks(xticks_values, xticks_labels,rotation=45)

# ax1.plot(rot_all_V,Ms,'o-',label='M$_s$')
# ax1.set_ylabel('M$_s$')
# ax1.grid(True)
# ax1.legend()
plt.xlabel('Ángulo (º)')
delta_V=(max(pendientes_V)- min(pendientes_V))/max(pendientes_V)
print(delta_V)
delta_V2=(max(pendientes_V2)- min(pendientes_V2))/max(pendientes_V2)
print(delta_V2)
#plt.savefig('Susceptibilidad_vs_angulo.png',dpi=300,facecolor='w')

#%% HORIZONTAL ROTADO 
# %matplotlib 
# %matplotlib inline
files=os.listdir(os.path.join(os.getcwd(),'H_rotado'))
files.sort()
masa_muestraH=0.32994 #g
masa_npmH=0.020 #g
masa_lauricoH=masa_muestraH-masa_npmH 

campos=[]
m_correg_norm=[]
fname=[]
Ms_H=[]

H_aux = np.linspace(-50,50,1000)
campo_lineal=[]
m_interp=[]

pendientes_H=[]
err_pend_H=[]
ordenadas_H=[]
Mr=[]
Hc=[]
#%
for f in files:
    B = np.loadtxt(os.path.join(os.getcwd(),'H_rotado',f),skiprows=12)
    H = B[:,0]
    m = B[:,1]
    #calculo la contribucion diamagnetica del ac laurico y se la descuento
    contrib_diamagV = chi_mass_laurico*masa_lauricoV*H
    m_correg = (m-contrib_diamagV)
    m_correg_norm_masa=m_correg/masa_npmV
    
    interpolador=interp1d(H,m_correg_norm_masa,fill_value='extrapolate')
    m_interp.append(interpolador(H_aux))
    
    m_norm=m_correg_norm_masa/max(m_correg_norm_masa) #Normalizo por valor maximo
    
    m_correg_norm.append(m_norm)
    campos.append(np.array(H))
    fname.append(f.split('_')[-1].split('.')[0])
    Ms_H.append(max(m_correg_norm_masa))
    
    (chi,n),pcov=curve_fit(lineal, H[np.nonzero(H<10)],m_norm[np.nonzero(H<10)]) #Ajuste lineal 
    err_m=pcov[0][0]
    pendientes_H.append(chi)
    err_pend_H.append(err_m)
    ordenadas_H.append(n)
    
    H_aux = np.linspace(-10,10,5000)
    H_recortado = H[np.nonzero(abs(H)<10)]
    m_recortado = m_correg_norm_masa[np.nonzero(abs(H)<10)]
    
    interpolador2=interp1d(H_recortado,m_recortado,fill_value='extrapolate')
    #m_new = interpolador2(H_aux)
    m_new = lineal(H_aux,chi,n)
    indx_H = np.nonzero(H_aux>=0)[0][0]
    indx_M = np.nonzero(m_new>0)[0][0]
    Hc.append(-H_aux[indx_M])
    Mr.append(m_new[indx_H])
    print('*'*50)
    print(f.split('_')[0], f.split('_')[-1])
    print('Susceptibilidad =',f'{chi:.3e}','+/-',f'{err_m:.3e}')
    print(f'Mag Remanente = {m_new[indx_H]:.3e}')
    print(f'Campo Coercitivo = {H_aux[indx_M]:.3e} A/m')

    fig,ax=plt.subplots(constrained_layout=True)
    ax.plot(H,m/masa_npmV,'.-',label='V')
    ax.plot(H,m_correg_norm_masa,'.-',label='V (s/ laurico)')
    ax.legend(ncol=2,loc='lower right')
    ax.grid()
    ax.set_xlabel('H (G)')
    ax.set_ylabel('m (emu/g)')
    ax.set_title(f)
    
    axins = ax.inset_axes([0.3, 0.2, 0.69, 0.55])
    axins.axvline(0,0,1,c='k',lw=0.8)
    axins.axhline(0,0,1,c='k',lw=0.8)
    axins.plot(H,m_norm,'.-',label='m norm')
    axins.plot(H_aux,lineal(H_aux,chi,n),'r-',label=f'$\chi$ = {chi:.3e}\nHc = {H_aux[indx_M]:.3e} A/m')
    axins.set_xlabel('H (G)')
    axins.grid()
    axins.legend(loc='lower right')
    # axins.legend()
    axins.set_xlim(-2,11)
    axins.set_ylim(-0.01,.06)
    ax.indicate_inset_zoom(axins, edgecolor="black")
    #plt.savefig('VSM_'+f[:-4]+'.png',dpi=300)
    plt.show()
#%% SUCEPTIBILIDAD vs ANGULO HORIZONTAL 
#   HORIZONTAL

rot_all_H= np.array([00,15,30,45,60,75,90,105,120,135,150,165,180,195,210,225,240,255,270,285,300,315,330,345,360])

fig,ax=plt.subplots(nrows=1,figsize=(9,4),sharex=True,constrained_layout=True)
ax.errorbar(x=rot_all_H,y=pendientes_H,yerr=err_pend_H, fmt='.', capsize=5,label='$\chi$')
ax.grid()
ax.legend()
ax.set_ylabel('$\chi$')
ax.set_title('Susceptibilidad vs angulo')

xticks_values = [00,15,30,45,60,75,90,105,120,135,150,165,180,195,210,225,240,255,270,285,300,315,330,345,360]
xticks_labels = ['00','15','30','45','60','75','90','105','120','135','150','165','180','195','210','225','240','255','270','285','300','315','330','345','360']
plt.xticks(xticks_values, xticks_labels)

# ax1.plot(rot_all_H,Ms,'o-',label='M$_s$')
# ax1.set_ylabel('M$_s$')
# ax1.grid(True)
# ax1.legend()
plt.xlabel('Ángulo (º)')
delta_H=(max(pendientes_H)- min(pendientes_H))/max(pendientes_H)
print(delta_H)
plt.savefig('Susceptibilidad_vs_angulo.png',dpi=300,facecolor='w')

# %%COMPARATIVA ALL
rot_all_V=rot_all_V-rot_all_V[0]
fig,ax=plt.subplots(nrows=2,figsize=(12,6),sharex=True,constrained_layout=True)
ax[0].errorbar(x=rot_all_H,y=pendientes_H,yerr=err_pend_H, fmt='.-', capsize=5,label=f'$\chi_H$ - $\Delta$ = {delta_H:.2f}')
ax[0].errorbar(x=rot_all_V,y=pendientes_V,yerr=err_pend_V, fmt='.-', capsize=5,label=f'$\chi_V$ - $\Delta$ = {delta_V:.2f}')
ax[0].errorbar(x=rot_all_V2,y=pendientes_V2,yerr=err_pend_V2, fmt='.-', capsize=5,label=f'$\chi_V$$_2$ - $\Delta$ = {delta_V2:.2f}')
ax[0].axhline(pendientes_Rd,0,1,c='k',label='$\chi_R$$_d$')

ax[0].grid()
ax[0].legend()
ax[0].set_ylabel('$\chi$')
ax[0].set_title('Susceptibilidad vs angulo')

ax[1].plot(rot_all_H,Ms_H,'o-',label='M$_s$ - H')
ax[1].plot(rot_all_V,Ms_V,'o-',label='M$_s$ - V')
ax[1].plot(rot_all_V2,Ms_V2,'o-',label='M$_s$ - V2')
# ax[1].axhline(Ms_Rd,0,1,label='M$_s$ - Rd')

ax[1].set_ylabel('M$_s$')
ax[1].grid(True)
ax[1].legend()
ax[1].set_title('Momento de saturación vs angulo')
ax[1].set_ylabel('m (emu/g)')
xticks_values =  np.arange(-30,375,15)
xticks_labels = [str(f) for f in xticks_values]
plt.xticks(xticks_values, xticks_labels,rotation=45)
plt.savefig('chi_Ms_vs_angulo.png',dpi=300)
plt.show()
# %%
