import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import seaborn as sns
from scipy.linalg import norm
import funciones as fc
import scipy as sc
import pandas as pd
from random import shuffle
from random import seed
#==============================================================================
# Cargo los datos
#==============================================================================
path = 'C:\\Users\\Usuario\\Documents\\Python Scripts\\scripts Ari\\Datos.pkl'
tabla = pd.read_pickle(path)
tabla = tabla.reset_index()
tabla.set_index(['Experiment','Well', 'Field', 'Cell'], inplace=True)
for i in range(len(tabla.index)):
    if tabla.index[i][1]=='W3':
        tabla.loc[tabla.index[i],'Dif']=True
    else:
        tabla.loc[tabla.index[i],'Dif']=False
    tabla.loc[tabla.index[i],'Area_media']=np.mean(tabla.Area[i])
    tabla.loc[tabla.index[i],'Area_inicial']=tabla.Area[i][0]
#==============================================================================
# Defino CC y G1 corregido
#==============================================================================            
CC_dif=tabla.Cycle_Length[(tabla.Is_Complete) & (tabla.Validated) & (tabla.Monster==False) & (tabla.Apoptosis==False) & (tabla.Dif)]
CC_UNdif=tabla.Cycle_Length[(tabla.Is_Complete) & (tabla.Validated) & (tabla.Monster==False) & (tabla.Apoptosis==False) & (tabla.Dif==False)]
Gen_dif=tabla.Generation[(tabla.Is_Complete) & (tabla.Validated) & (tabla.Monster==False) & (tabla.Apoptosis==False) & (tabla.Dif)]
Gen_UNdif=tabla.Generation[(tabla.Is_Complete) & (tabla.Validated) & (tabla.Monster==False) & (tabla.Apoptosis==False) & (tabla.Dif==False)]
G1_dif=tabla.G1_Length[(np.isnan(tabla.G1_Length)==False) & (tabla.Validated) & (tabla.Monster==False) & (tabla.Apoptosis==False) & (tabla.Dif)]
G1_UNdif=tabla.G1_Length[(np.isnan(tabla.G1_Length)==False) & (tabla.Validated) & (tabla.Monster==False) & (tabla.Apoptosis==False) & (tabla.Dif==False)]
Area_dif_media=tabla.Area_media[(tabla.Is_Complete) & (tabla.Validated) & (tabla.Monster==False) & (tabla.Apoptosis==False) & (tabla.Dif)]
Area_UNdif_media=tabla.Area_media[(tabla.Is_Complete) & (tabla.Validated) & (tabla.Monster==False) & (tabla.Apoptosis==False) & (tabla.Dif==False)]
Area_dif_inicial=tabla.Area_inicial[(tabla.Generation !=1) & (tabla.Validated) & (tabla.Monster==False) & (tabla.Apoptosis==False) & (tabla.Dif)]
Area_UNdif_inicial=tabla.Area_inicial[(tabla.Generation !=1) & (tabla.Validated) & (tabla.Monster==False) & (tabla.Apoptosis==False) & (tabla.Dif==False)]
## ESTOS SON PARA CC
#Val_dif=[(tabla.Is_Complete) & (tabla.Validated) & (tabla.Monster==False) & (tabla.Apoptosis==False) & (tabla.Dif)]
#Val_UNdif=[(tabla.Is_Complete) & (tabla.Validated) & (tabla.Monster==False) & (tabla.Apoptosis==False) & (tabla.Dif==False)]
### ESTOS SON PARA G1
#Val_dif=[(np.isnan(tabla.G1_Length)==False) & (tabla.Validated) & (tabla.Monster==False) & (tabla.Apoptosis==False) & (tabla.Dif)]
#Val_UNdif=[(np.isnan(tabla.G1_Length)==False) & (tabla.Validated) & (tabla.Monster==False) & (tabla.Apoptosis==False) & (tabla.Dif==False)]
# ESTOS SON PARA G2
G1_dif=tabla.G1_Length[(((tabla.Generation)==1) & (tabla.Lineage_in_G1==True) & ((np.isnan(tabla.G1_Length)==False) & (tabla.Validated) & (tabla.Monster==False) & (tabla.Apoptosis==False) & (tabla.Dif))) | ((np.isnan(tabla.G1_Length)==False) & ((tabla.Generation)!=1) & (tabla.Validated) & (tabla.Monster==False) & (tabla.Apoptosis==False) & (tabla.Dif))]
G1_UNdif=tabla.G1_Length[(((tabla.Generation)==1) & (tabla.Lineage_in_G1==True) & ((np.isnan(tabla.G1_Length)==False) & (tabla.Validated) & (tabla.Monster==False) & (tabla.Apoptosis==False) & (tabla.Dif==False))) | ((np.isnan(tabla.G1_Length)==False) & ((tabla.Generation)!=1) & (tabla.Validated) & (tabla.Monster==False) & (tabla.Apoptosis==False) & (tabla.Dif==False))]
##ESTOS SON PARA Área inicial
#Val_dif=[(tabla.Generation !=1) & (tabla.Validated) & (tabla.Monster==False) & (tabla.Apoptosis==False) & (tabla.Dif)]
#Val_UNdif=[(tabla.Generation !=1) & (tabla.Validated) & (tabla.Monster==False) & (tabla.Apoptosis==False) & (tabla.Dif==False)]
for i in G1_dif.index:
    if tabla.G1_Length_Manual[i] is not None:
        G1_dif[i]=tabla.G1_Length_Manual[i]
for i in G1_UNdif.index:
    if tabla.G1_Length_Manual[i] is not None:
        G1_UNdif[i]=tabla.G1_Length_Manual[i]
Val_dif=[(((tabla.Generation)==1) & (tabla.Lineage_in_G1==True) & ((np.isnan(tabla.G1_Length)==False) & (tabla.Validated) & (tabla.Monster==False) & (tabla.Apoptosis==False) & (tabla.Dif))) | ((np.isnan(tabla.G1_Length)==False) & ((tabla.Generation)!=1) & (tabla.Validated) & (tabla.Monster==False) & (tabla.Apoptosis==False) & (tabla.Dif))]
Val_UNdif=[(((tabla.Generation)==1) & (tabla.Lineage_in_G1==True) & ((np.isnan(tabla.G1_Length)==False) & (tabla.Validated) & (tabla.Monster==False) & (tabla.Apoptosis==False) & (tabla.Dif==False))) | ((np.isnan(tabla.G1_Length)==False) & ((tabla.Generation)!=1) & (tabla.Validated) & (tabla.Monster==False) & (tabla.Apoptosis==False) & (tabla.Dif==False))]
G2_dif=tabla.Cycle_Length[(((tabla.Generation)==1) & (tabla.Lineage_in_G1==True) & ((np.isnan(tabla.G1_Length)==False) & (tabla.Validated) & (tabla.Monster==False) & (tabla.Apoptosis==False) & (tabla.Dif))) | ((np.isnan(tabla.G1_Length)==False) & ((tabla.Generation)!=1) & (tabla.Validated) & (tabla.Monster==False) & (tabla.Apoptosis==False) & (tabla.Dif))]-G1_dif
G2_UNdif=tabla.Cycle_Length[(((tabla.Generation)==1) & (tabla.Lineage_in_G1==True) & ((np.isnan(tabla.G1_Length)==False) & (tabla.Validated) & (tabla.Monster==False) & (tabla.Apoptosis==False) & (tabla.Dif==False))) | ((np.isnan(tabla.G1_Length)==False) & ((tabla.Generation)!=1) & (tabla.Validated) & (tabla.Monster==False) & (tabla.Apoptosis==False) & (tabla.Dif==False))]-G1_UNdif
#==============================================================================
## Armo los vectores T con la duración de los CC
## V-A-M-H es T; A-M-H es T2; M-H es T1;
#==============================================================================
from Grass_proc_func_var import Grass_proc_variable
NNN=40;CCC=[0.75,1.05,1.55];D=[2,3,4]
#==============================================================================
## Lista de los radios máximos (Undif - Dif):
# CC-L _ [0.55,0.55,0.525] _ [1.2765,1.5,1.2765]
# G1-L _  [1.3,1.3,2.25] _ [2.25,2.25,2.25]
# S/G2/M-L _ [0.87,1.2,1.3] _ [0.75,1.05,1.55]
# Área (control) _ [1.5,2,1.6] _ [1.7,2.25,2.1]
#==============================================================================
Pendiente_control=np.zeros((NNN,3));Pendiente=np.zeros(3)
for dd in range(len(D)):
    CC=G2_dif*1
    Val=Val_dif*1
    Pendiente[dd]=Grass_proc_variable(CC,Val,tabla.Generation,tabla.Daughters,CCC[dd],D[dd])
    for jjj in range(NNN):
        seed(jjj)
        CC=G2_dif*1
        shuffle(CC)
        seed(jjj)
        Val=Val_dif*1
        shuffle(Val)
        Pendiente_control[jjj,dd]=Grass_proc_variable(CC,Val,tabla.Generation,tabla.Daughters,CCC[dd],D[dd])
        print([D[dd],jjj])
Pendiente_control_mean=np.zeros(3)
Pendiente_control_sd=np.zeros(3)
Pendiente_control_median=np.zeros(3)
Pendiente_control_97=np.zeros(3)
Pendiente_control_2=np.zeros(3)
Pendiente_control_75=np.zeros(3)
Pendiente_control_25=np.zeros(3)
for dd in range(len(D)):
    Pendiente_control_mean[dd]=np.mean(Pendiente_control[:,dd])
    Pendiente_control_sd[dd]=np.std(Pendiente_control[:,dd])
    Pendiente_control_median[dd]=np.percentile(Pendiente_control[:,dd],50)
    Pendiente_control_97[dd]=np.percentile(Pendiente_control[:,dd],97.5)
    Pendiente_control_2[dd]=np.percentile(Pendiente_control[:,dd],2.5)    
    Pendiente_control_75[dd]=np.percentile(Pendiente_control[:,dd],75)
    Pendiente_control_25[dd]=np.percentile(Pendiente_control[:,dd],25)
#==============================================================================
# Ploteo d vs De
#==============================================================================
# De=[1,2,3,4]
#plt.figure();plt.hold(True)        
#fig5=plt.fill_between([2,3,4], Pendiente_control_2, Pendiente_control_97,color=(0,1,0,1),label='Shuffle 95%')
#fig4=plt.fill_between([2,3,4], Pendiente_control_25, Pendiente_control_75,color=(0,0.5,0,1),label='Shuffle 50%')
#fig3, =plt.plot([2,3,4], Pendiente_control_median,color=(0,1,0,1),label='Shuffle Median') 
#fig2, =plt.plot([2,3,4],[2,3,4],lw=3,color='r',label='Random')
#fig1, =plt.plot([2,3,4],Pendiente,lw=3,color='b',label='Datos')
#plt.xlabel('De');plt.ylabel('d');plt.title('CC-L - En Diferenciación');plt.legend(handles=[fig1,fig2,fig3,fig4,fig5],loc='upper left');
          
def gpplot(fplip=1):
     plt.figure();ax=plt.subplot();plt.hold(True)         
     fig4=plt.fill_between([2,3,4], Pendiente_control_2, Pendiente_control_97,color=(0,1,0,0.4),label='Shuffle 95%')
     fig3=plt.fill_between([2,3,4], Pendiente_control_25, Pendiente_control_75,color=(0,0.5,0,0.4),label='Shuffle 50%')
     fig2, =plt.plot([2,3,4],[2,3,4],color=(0.8,0.1,0.1,1),label='Random')
     fig1, =plt.plot([2,3,4],Pendiente,color=(0.1,0.1,0.8,1),label='Data')
     plt.xlabel(r'D$_E$', fontsize=40);plt.ylabel('d', fontsize=40);plt.title('$S/G2/M-L$ - Dif', fontweight='bold', fontsize=40);plt.legend(handles=[fig1,fig2,fig3,fig4],loc='upper left');
     ax.set_xticks([2,3,4]);ax.set_yticks([1,2,3,4,5,6])
sns.set_context("poster", font_scale=1.5, rc={"lines.linewidth": 10})
gpplot()    
