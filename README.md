# Fede
Probando Grassberguer-Proccatia

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from scipy.linalg import norm
import funciones as fc
import scipy as sc
#==============================================================================
## Armo los vectores T con la duración de los CC
## V-A-M-H es T; A-M-H es T2; M-H es T1;
#==============================================================================
T,T1,G,G1,T2,G2=[],[],[],[],[],[]
for i in range(len(Daughters.index)):
    if Is_Complete[(Daughters.index[i])]==True and (Daughters.index[i][1])=='W3': #Filtro las madres completas Dif
        if Is_Complete[(Daughters[i][0])]==True: #Filtro hija 0 completa
            G1.append(Generation[(Daughters.index[i])])
            T1.append((Cycle_Length[(Daughters.index[i])],Cycle_Length[(Daughters[i][0])]))
            if Is_Complete[(Daughters[(Daughters[i][0])][0])]==True: #Filtro nieta 0 0 completa
                G2.append(Generation[(Daughters.index[i])])
                T2.append((Cycle_Length[(Daughters.index[i])],Cycle_Length[(Daughters[i][0])],Cycle_Length[(Daughters[(Daughters[i][0])][0])]))
                if Is_Complete[(Daughters[(Daughters[(Daughters[i][0])][0])][0])]==True:
                    T.append((Cycle_Length[(Daughters.index[i])],Cycle_Length[(Daughters[i][0])],Cycle_Length[(Daughters[(Daughters[i][0])][0])],Cycle_Length[(Daughters[(Daughters[(Daughters[i][0])][0])][0])]))
                    G.append(Generation[(Daughters.index[i])])
                if Is_Complete[(Daughters[(Daughters[(Daughters[i][0])][0])][1])]==True:
                    T.append((Cycle_Length[(Daughters.index[i])],Cycle_Length[(Daughters[i][0])],Cycle_Length[(Daughters[(Daughters[i][0])][0])],Cycle_Length[(Daughters[(Daughters[(Daughters[i][0])][0])][1])]))
                    G.append(Generation[(Daughters.index[i])])
            if Is_Complete[(Daughters[(Daughters[i][0])][1])]==True: #Filtro nieta 0 1 completa
                G2.append(Generation[(Daughters.index[i])])
                T2.append((Cycle_Length[(Daughters.index[i])],Cycle_Length[(Daughters[i][0])],Cycle_Length[(Daughters[(Daughters[i][0])][1])]))
                if Is_Complete[(Daughters[(Daughters[(Daughters[i][0])][1])][0])]==True:
                    T.append((Cycle_Length[(Daughters.index[i])],Cycle_Length[(Daughters[i][0])],Cycle_Length[(Daughters[(Daughters[i][0])][1])],Cycle_Length[(Daughters[(Daughters[(Daughters[i][0])][1])][0])]))
                    G.append(Generation[(Daughters.index[i])])
                if Is_Complete[(Daughters[(Daughters[(Daughters[i][0])][1])][1])]==True:
                    T.append((Cycle_Length[(Daughters.index[i])],Cycle_Length[(Daughters[i][0])],Cycle_Length[(Daughters[(Daughters[i][0])][1])],Cycle_Length[(Daughters[(Daughters[(Daughters[i][0])][1])][1])]))
                    G.append(Generation[(Daughters.index[i])])
        if Is_Complete[(Daughters[i][1])]==True: #Filtro hija 1 completa
            G1.append(Generation[(Daughters.index[i])])
            T1.append((Cycle_Length[(Daughters.index[i])],Cycle_Length[(Daughters[i][1])]))
            if Is_Complete[(Daughters[(Daughters[i][1])][0])]==True: #Filtro nieta 1 0 completa
                G2.append(Generation[(Daughters.index[i])])
                T2.append((Cycle_Length[(Daughters.index[i])],Cycle_Length[(Daughters[i][1])],Cycle_Length[(Daughters[(Daughters[i][1])][0])]))
                if Is_Complete[(Daughters[(Daughters[(Daughters[i][1])][0])][0])]==True:
                    T.append((Cycle_Length[(Daughters.index[i])],Cycle_Length[(Daughters[i][1])],Cycle_Length[(Daughters[(Daughters[i][1])][0])],Cycle_Length[(Daughters[(Daughters[(Daughters[i][1])][0])][0])]))
                    G.append(Generation[(Daughters.index[i])])
                if Is_Complete[(Daughters[(Daughters[(Daughters[i][1])][0])][1])]==True:
                    T.append((Cycle_Length[(Daughters.index[i])],Cycle_Length[(Daughters[i][1])],Cycle_Length[(Daughters[(Daughters[i][1])][0])],Cycle_Length[(Daughters[(Daughters[(Daughters[i][1])][0])][1])]))
                    G.append(Generation[(Daughters.index[i])])
            if Is_Complete[(Daughters[(Daughters[i][1])][1])]==True: #Filtro nieta 1 1 completa
                G2.append(Generation[(Daughters.index[i])])
                T2.append((Cycle_Length[(Daughters.index[i])],Cycle_Length[(Daughters[i][1])],Cycle_Length[(Daughters[(Daughters[i][1])][1])]))
                if Is_Complete[(Daughters[(Daughters[(Daughters[i][1])][1])][0])]==True:
                    T.append((Cycle_Length[(Daughters.index[i])],Cycle_Length[(Daughters[i][1])],Cycle_Length[(Daughters[(Daughters[i][1])][1])],Cycle_Length[(Daughters[(Daughters[(Daughters[i][1])][1])][0])]))
                    G.append(Generation[(Daughters.index[i])])
                if Is_Complete[(Daughters[(Daughters[(Daughters[i][1])][1])][1])]==True:
                    T.append((Cycle_Length[(Daughters.index[i])],Cycle_Length[(Daughters[i][1])],Cycle_Length[(Daughters[(Daughters[i][1])][1])],Cycle_Length[(Daughters[(Daughters[(Daughters[i][1])][1])][1])]))
                    G.append(Generation[(Daughters.index[i])])
                    
#==============================================================================
## Redefino T1 y T2 para tomar la misma cantidad de puntos que T
#==============================================================================
T11, G11=T1, G1;T1, G1=[],[];T1.append(T11[1]);G1.append(G11[1])
for k in range(1,len(T11)):
    if T11[k]!=T11[k-1]:
            T1.append(T11[k]);G1.append(G11[k])
T22, G22=T2, G2;T2,G2=[],[];T2.append(T22[1]);G2.append(G22[1])
for k in range(1,len(T22)):
    if T22[k]!=T22[k-1]:
            T2.append(T22[k]);G2.append(G22[k])
#==============================================================================
## Calculo el punto medio para centrar la esfera en el punto más centrado
#==============================================================================
## T1
T1m=(0,0);R=[0,0];T1mean1,T1mean1=[],[];dT1=np.zeros(len(T))
for i in range(len(T)):
    T1mean1=T1[i][0];T1mean2=T1[i][1];T1m=(T1m[0]+T1mean1,T1m[1]+T1mean2)
T1m=(T1m[0]/len(T),T1m[1]/len(T))
## T2
T2m=(0,0,0);R=[0,0,0];T2mean1,T2mean2,T2mean3=[],[],[];dT2=np.zeros(len(T))
for i in range(len(T)):
    T2mean1=T2[i][0]
    T2mean2=T2[i][1]
    T2mean3=T2[i][2]
    T2m=(T2m[0]+T2mean1,T2m[1]+T2mean2,T2m[2]+T2mean3)
T2m=(T2m[0]/len(T),T2m[1]/len(T),T2m[2]/len(T))
## T
Tm=(0,0,0,0);R=[0,0,0,0];Tmean1,Tmean2,Tmean3,Tmean4=[],[],[],[];dT=np.zeros(len(T))
for i in range(len(T)):
    Tmean1=T[i][0]
    Tmean2=T[i][1]
    Tmean3=T[i][2]
    Tmean4=T[i][3]
    Tm=(Tm[0]+Tmean1,Tm[1]+Tmean2,Tm[2]+Tmean3,Tm[3]+Tmean4)
Tm=(Tm[0]/len(T),Tm[1]/len(T),Tm[2]/len(T),Tm[3]/len(T))
#==============================================================================
## Calculo la distancia de los puntos al medio y elijo el más centrado
#==============================================================================
for i in range(len(T)):
    dT1[i]=norm([T1m[0]-T1[i][0],T1m[1]-T1[i][1]])
    dT2[i]=norm([T2m[0]-T2[i][0],T2m[1]-T2[i][1],T2m[2]-T2[i][2]])
    dT[i]=norm([Tm[0]-T[i][0],Tm[1]-T[i][1],Tm[2]-T[i][2],Tm[3]-T[i][3]])
i=np.argmin(dT1);T1m=(T1[i][0],T1[i][1])
i=np.argmin(dT2);T2m=(T2[i][0],T2[i][1],T2[i][2])
i=np.argmin(dT);Tm=(T[i][0],T[i][1],T[i][2],T[i][3])
#==============================================================================
## Calculo la distancia de los puntos al centro de la esfera
#==============================================================================
for i in range(len(T)):
    dT1[i]=norm([T1m[0]-T1[i][0],T1m[1]-T1[i][1]])
    dT2[i]=norm([T2m[0]-T2[i][0],T2m[1]-T2[i][1],T2m[2]-T2[i][2]])
    dT[i]=norm([Tm[0]-T[i][0],Tm[1]-T[i][1],Tm[2]-T[i][2],Tm[3]-T[i][3]])
#==============================================================================
## Calculo la distancia de puntos Random al centro de la esfera
#==============================================================================
N1=1000 ##Cantidad de puntos Random
N2=10000 ##Cantidad de puntos Random
N=100000 ##Cantidad de puntos Random
dR1, dR2, dR, R1, R2, R=np.zeros(N1),np.zeros(N2),np.zeros(N),np.zeros((N1,2)),np.zeros((N2,3)),np.zeros((N,4))
for i in range(N1):
    R1[i,0]=np.random.rand()*(np.max([np.max(T),np.max(T1),np.max(T2)])-np.min([np.min(T),np.min(T1),np.min(T2)]))+np.min([np.min(T),np.min(T1),np.min(T2)])
    R1[i,1]=np.random.rand()*(np.max([np.max(T),np.max(T1),np.max(T2)])-np.min([np.min(T),np.min(T1),np.min(T2)]))+np.min([np.min(T),np.min(T1),np.min(T2)])
    dR1[i]=norm([T1m[0]-R1[i,0],T1m[1]-R1[i,1]])
for i in range(N2):
    R2[i,0]=np.random.rand()*(np.max([np.max(T),np.max(T1),np.max(T2)])-np.min([np.min(T),np.min(T1),np.min(T2)]))+np.min([np.min(T),np.min(T1),np.min(T2)])
    R2[i,1]=np.random.rand()*(np.max([np.max(T),np.max(T1),np.max(T2)])-np.min([np.min(T),np.min(T1),np.min(T2)]))+np.min([np.min(T),np.min(T1),np.min(T2)])
    R2[i,2]=np.random.rand()*(np.max([np.max(T),np.max(T1),np.max(T2)])-np.min([np.min(T),np.min(T1),np.min(T2)]))+np.min([np.min(T),np.min(T1),np.min(T2)])
    dR2[i]=norm([T2m[0]-R2[i,0],T2m[1]-R2[i,1],T2m[2]-R2[i,2]])
for i in range(N):
    R[i,0]=np.random.rand()*(np.max([np.max(T),np.max(T1),np.max(T2)])-np.min([np.min(T),np.min(T1),np.min(T2)]))+np.min([np.min(T),np.min(T1),np.min(T2)])
    R[i,1]=np.random.rand()*(np.max([np.max(T),np.max(T1),np.max(T2)])-np.min([np.min(T),np.min(T1),np.min(T2)]))+np.min([np.min(T),np.min(T1),np.min(T2)])
    R[i,2]=np.random.rand()*(np.max([np.max(T),np.max(T1),np.max(T2)])-np.min([np.min(T),np.min(T1),np.min(T2)]))+np.min([np.min(T),np.min(T1),np.min(T2)])
    R[i,3]=np.random.rand()*(np.max([np.max(T),np.max(T1),np.max(T2)])-np.min([np.min(T),np.min(T1),np.min(T2)]))+np.min([np.min(T),np.min(T1),np.min(T2)])
    dR[i]=norm([Tm[0]-R[i,0],Tm[1]-R[i,1],Tm[2]-R[i,2],Tm[3]-R[i,3]])
    print(i)
#==============================================================================
# Calculo cotas de los radios
#==============================================================================
dr1minx,dr1maxx,dr1miny,dr1maxy=np.min(R1[:,0]),np.max(R1[:,0]),np.min(R1[:,1]),np.max(R1[:,1])
dr1=min(abs(T1m[0]-dr1minx),abs(T1m[0]-dr1maxx),abs(T1m[1]-dr1miny),abs(T1m[1]-dr1maxy))
dr2minx,dr2maxx,dr2miny,dr2maxy,dr2minz,dr2maxz=np.min(R2[:,0]),np.max(R2[:,0]),np.min(R2[:,1]),np.max(R2[:,1]),np.min(R2[:,2]),np.max(R2[:,2])
dr2=min(abs(T2m[0]-dr2minx),abs(T2m[0]-dr2maxx),abs(T2m[1]-dr2miny),abs(T2m[1]-dr2maxy),abs(T2m[2]-dr2minz),abs(T2m[2]-dr2maxz))
drminx,drmaxx,drminy,drmaxy,drminz,drmaxz,drminw,drmaxw=np.min(R[:,0]),np.max(R[:,0]),np.min(R[:,1]),np.max(R[:,1]),np.min(R[:,2]),np.max(R[:,2]),np.min(R[:,3]),np.max(R[:,3])
dr=min(abs(Tm[0]-drminx),abs(Tm[0]-drmaxx),abs(Tm[1]-drminy),abs(Tm[1]-drmaxy),abs(Tm[2]-drminz),abs(Tm[2]-drmaxz),abs(Tm[3]-drminw),abs(Tm[3]-drmaxw))
#==============================================================================
## Cuento puntos dentro de la esfera en función del radio r hasta la cota
#==============================================================================
resR=1000 ## Cantidad de puntos para el radio
r1=np.logspace(np.log10(dr1)-1,np.log10(dr1),resR);C1=np.zeros(len(r1));CR1=np.zeros(len(r1))
for j in range(len(r1)):
    for i in range(len(T)):
        if dT1[i]<r1[j]:
            C1[j]=C1[j]+1
    for i in range(N1):
        if dR1[i]<r1[j]:
            CR1[j]=CR1[j]+1

r2=np.logspace(np.log10(dr2)-1,np.log10(dr2),resR);C2=np.zeros(len(r2));CR2=np.zeros(len(r2))
for j in range(len(r2)):
    for i in range(len(T)):
        if dT2[i]<r2[j]:
            C2[j]=C2[j]+1
    for i in range(N2):
        if dR2[i]<r2[j]:
            CR2[j]=CR2[j]+1
             
r=np.logspace(np.log10(dr)-1,np.log10(dr),resR);C=np.zeros(len(r));CR=np.zeros(len(r))
for j in range(len(r)):
    for i in range(len(T)):
        if dT[i]<r[j]:
            C[j]=C[j]+1
    for i in range(N):
        if dR[i]<r[j]:
            CR[j]=CR[j]+1
#==============================================================================
# Normalizo radio y C
#==============================================================================
rr,rr1,rr2=r/np.max(r),r1/np.max(r1),r2/np.max(r2)
C,C1,C2,CR,CR1,CR2=C/np.max(C),C1/np.max(C1),C2/np.max(C2),CR/np.max(CR),CR1/np.max(CR1),CR2/np.max(CR2)
#==============================================================================
# Ajusto lineal
#==============================================================================
p0=[4,5];xd=np.log10(rr[~(C==0)]);yd=np.log10(C[~(C==0)])
p=sc.optimize.leastsq(fc.res, p0, args=(xd,yd), Dfun=None, full_output=0, col_deriv=0, ftol=1.49012e-08, xtol=1.49012e-08, gtol=0.0, maxfev=0, epsfcn=None, factor=100, diag=None)
p0=[2,5];xd=np.log10(rr1[~(C1==0)]);yd=np.log10(C1[~(C1==0)])
p1=sc.optimize.leastsq(fc.res, p0, args=(xd,yd), Dfun=None, full_output=0, col_deriv=0, ftol=1.49012e-08, xtol=1.49012e-08, gtol=0.0, maxfev=0, epsfcn=None, factor=100, diag=None)
p0=[3,5];xd=np.log10(rr2[~(C2==0)]);yd=np.log10(C2[~(C2==0)])
p2=sc.optimize.leastsq(fc.res, p0, args=(xd,yd), Dfun=None, full_output=0, col_deriv=0, ftol=1.49012e-08, xtol=1.49012e-08, gtol=0.0, maxfev=0, epsfcn=None, factor=100, diag=None)
p0=[4,5];xd=np.log10(rr[~(CR==0)]);yd=np.log10(CR[~(CR==0)])
pR=sc.optimize.leastsq(fc.res, p0, args=(xd,yd), Dfun=None, full_output=0, col_deriv=0, ftol=1.49012e-08, xtol=1.49012e-08, gtol=0.0, maxfev=0, epsfcn=None, factor=100, diag=None)
p0=[2,5];xd=np.log10(rr1[~(CR1==0)]);yd=np.log10(CR1[~(CR1==0)])
pR1=sc.optimize.leastsq(fc.res, p0, args=(xd,yd), Dfun=None, full_output=0, col_deriv=0, ftol=1.49012e-08, xtol=1.49012e-08, gtol=0.0, maxfev=0, epsfcn=None, factor=100, diag=None)
p0=[3,5];xd=np.log10(rr2[~(CR2==0)]);yd=np.log10(CR2[~(CR2==0)])
pR2=sc.optimize.leastsq(fc.res, p0, args=(xd,yd), Dfun=None, full_output=0, col_deriv=0, ftol=1.49012e-08, xtol=1.49012e-08, gtol=0.0, maxfev=0, epsfcn=None, factor=100, diag=None)
#==============================================================================
# Calculo D=A*r^n
#==============================================================================
D=np.zeros(len(rr))
for i in range(len(D)):
    D[i]=(10**p[0][1])*(rr[i]**p[0][0])
D1=np.zeros(len(rr1))
for i in range(len(D1)):
    D1[i]=(10**p1[0][1])*(rr1[i]**p1[0][0])
D2=np.zeros(len(rr2))
for i in range(len(D2)):
    D2[i]=(10**p2[0][1])*(rr2[i]**p2[0][0])
    
DR=np.zeros(len(rr))
for i in range(len(DR)):
    DR[i]=(10**pR[0][1])*(rr[i]**pR[0][0])
DR1=np.zeros(len(rr1))
for i in range(len(DR1)):
    DR1[i]=(10**pR1[0][1])*(rr1[i]**pR1[0][0])
DR2=np.zeros(len(rr2))
for i in range(len(DR2)):
    DR2[i]=(10**pR2[0][1])*(rr2[i]**pR2[0][0])
#==============================================================================
# Ploteo C vs r
#==============================================================================
plt.figure()    
plt.subplot(2,3,1)
fig1, =plt.loglog(rr1,C1,'o',label='2',ms=1);plt.hold(True);fig2, =plt.loglog(rr2,C2,'o',label='3',ms=1);fig3, =plt.loglog(rr,C,'o',label='4',ms=1)
fig4, =plt.loglog(rr1,D1,'r',label='5');fig4, =plt.loglog(rr2,D2,'r',label='5');fig4, =plt.loglog(rr,D,'r',label='Ajustes')
plt.legend(handles=[fig1, fig2, fig3, fig4]);plt.xlabel('r');plt.ylabel('C');plt.axis([0.1, 1, 0.01, 1.5])
plt.subplot(2,3,2)
fig1, =plt.loglog(rr1,CR1,'o',label='2',ms=1);plt.hold(True);fig2, =plt.loglog(rr2,CR2,'o',label='3',ms=1);fig3, =plt.loglog(rr,CR,'o',label='4',ms=1)
fig4, =plt.loglog(rr1,DR1,'r',label='5');fig4, =plt.loglog(rr2,DR2,'r',label='5');fig4, =plt.loglog(rr,DR,'r',label='Ajustes')
plt.legend(handles=[fig1, fig2, fig3, fig4])
plt.xlabel('r');plt.ylabel('C Rand');plt.axis([0.1, 1, 0.01, 1.5])
plt.subplot(2,3,4)
fig1, =plt.semilogx(rr1,C1,'o',label='2',ms=1);plt.hold(True);fig2, =plt.plot(rr2,C2,'o',label='3',ms=1);fig3, =plt.plot(rr,C,'o',label='4',ms=1)
fig4, =plt.plot(rr1,D1,'r',label='5');fig4, =plt.plot(rr2,D2,'r',label='5');fig4, =plt.plot(rr,D,'r',label='Ajustes')
plt.legend(handles=[fig1, fig2, fig3, fig4]);plt.xlabel('r');plt.ylabel('C');plt.axis([0.1, 1, 0.01, 1])
plt.subplot(2,3,5)
fig1, =plt.plot(rr1,CR1,'o',label='2',ms=1);plt.hold(True);fig2, =plt.plot(rr2,CR2,'o',label='3',ms=1);fig3, =plt.plot(rr,CR,'o',label='4',ms=1)
fig4, =plt.plot(rr1,DR1,'r',label='5');fig4, =plt.plot(rr2,DR2,'r',label='5');fig4, =plt.plot(rr,DR,'r',label='Ajustes')
plt.legend(handles=[fig1, fig2, fig3, fig4]);plt.xlabel('r');plt.ylabel('C Rand');plt.axis([0.1, 1, 0.01, 1])
plt.subplot(2,3,3)
fig1, =plt.plot([2,3,4],[p1[0][0],p2[0][0],p[0][0]],label='Datos')
plt.hold(True)
fig2, =plt.plot([2,3,4],[pR1[0][0],pR2[0][0],pR[0][0]],label='Random')
plt.legend(handles=[fig1, fig2]);plt.xlabel('De');plt.ylabel('d');
##==============================================================================
## Ploteo M-H
##==============================================================================
#plt.figure()
#plt.subplot(1,2,1);jet=plt.get_cmap('jet');cNorm  = colors.Normalize(vmin=1, vmax=5);scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
#for i in range(len(T1)):
#    colorVal = scalarMap.to_rgba(G1[i])
#    fig=plt.plot(T1[i][0],T1[i][1], 'o',color=colorVal,ms=5)
#    plt.hold(True)
#plt.plot(T1m[0],T1m[1],'b*',ms=10);plt.xlabel('CC Madre');plt.ylabel('CC Hija');plt.title('Dif W3')
#plt.subplot(1,2,2)
#for i in range(N1):
#    fig=plt.plot(R1[i][0],R1[i][1], 'ro',ms=5)
#    plt.hold(True)
#plt.xlabel('CC Madre');plt.ylabel('CC Hija');plt.title('Random')
##==============================================================================
## Ploteo A-M-H
##==============================================================================
#fig=plt.figure()
#jet=plt.get_cmap('jet');cNorm  = colors.Normalize(vmin=1, vmax=5);scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
#for i in range(len(T2)):
#    colorVal = scalarMap.to_rgba(G2[i])
#    ax = fig.gca(projection='3d')
#    ax.scatter(T2[i][0],T2[i][1],T2[i][2], 'o',color=colorVal,s=10)
#    plt.hold(True)
#ax.scatter(T2m[0],T2m[1],T2m[2], 'b*',s=20)
#plt.title('Dif W3');ax.set_xlabel('CC Abuela');ax.set_ylabel('CC Madre');ax.set_zlabel('CC Hija')
#fig=plt.figure()
#for i in range(N2):
#    ax = fig.gca(projection='3d')
#    ax.scatter(R2[i,0],R2[i,1],R2[i,2], 'o',c='r',s=10)
#    plt.hold(True)
#plt.title('Random');ax.set_xlabel('CC Abuela');ax.set_ylabel('CC Madre');ax.set_zlabel('CC Hija')

##==============================================================================
## Ploteo VisA-A-M-H
##==============================================================================
#fig=plt.figure()
#jet=plt.get_cmap('jet');cNorm  = colors.Normalize(vmin=1, vmax=5);scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
#for i in range(len(T)):
#    colorVal = scalarMap.to_rgba(G2[i])
#    ax = fig.gca(projection='3d')
#    ax.scatter(T[i][0],T[i][1],T[i][2], 'o',color=colorVal,s=10)
#    plt.hold(True)
#ax.scatter(Tm[0],Tm[1],Tm[2], 'b*',s=20)
#plt.title('Dif W3');ax.set_xlabel('CC Abuela');ax.set_ylabel('CC Madre');ax.set_zlabel('CC Hija')
#fig=plt.figure()
#for i in range(N):
#    ax = fig.gca(projection='3d')
#    ax.scatter(R[i,0],R[i,1],R[i,2], 'o',c='r',s=10)
#    plt.hold(True)
#plt.title('Random');ax.set_xlabel('CC Abuela');ax.set_ylabel('CC Madre');ax.set_zlabel('CC Hija')
