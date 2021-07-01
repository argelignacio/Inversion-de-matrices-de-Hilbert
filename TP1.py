# -*- coding: utf-8 -*-
"""
Created on Wed May 19 12:23:58 2021

@author: Usuario
"""

from typing import ClassVar
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import imshow
from sympy import *
import math 
import pandas as pd
#punto a
'constantes a usar en el problema'
N = 4

def calculoerror(N,Hinversa):
    'Defino la matriz de hilbert'
    arreglos =[]
    
    for i in range(N):
        subarreglo = [0 for i in range(N)]
    
        for j in range(N):
            constanteHilbert = (i+1)+(j+1)-1
            subarreglo[j] = (1/(constanteHilbert))
        arreglos.append(subarreglo)
    
    hilbertMatrix = np.array(arreglos) 
    
    matrizU = hilbertMatrix.copy()
    
    #defino las matrices L 
    matrizL= np.identity(N)
    'defino la matriz L'
    columna=0
    while columna < N:
        for fila in range (N):
            if fila>columna:
                m_ia=(matrizU[fila][columna])/(matrizU[columna][columna])
                for j in range(1,N-columna):
                    matrizU[fila][columna+j] = matrizU[fila][columna+j]-m_ia*matrizU[columna][columna+j]
                matrizL[fila][columna]=m_ia
        columna+=1
    'defino la matriz U'
    for fil in range(N):
            for col in range(N):
                if (col < fil) :
                    matrizU[fil,col] = 0  
    
    x= np.zeros((N,N))
    def sumatoria1(i,y):
        acum=0
        for j in range(i):
            acum+=matrizL[i][j]*y[j]
        return(acum)
    def sumatoria2(i,xk):
        acum=0
        for j in range(i+1,N):
            acum+=matrizU[i][j]*xk[j]
        return(acum)      
    canonico=0
    while canonico < N:
        #creo el termino independiente
        b=np.zeros(N) #obs deberia quedar transpuesto
        for i in range(N):
            if i==canonico:
                b[i]=1        
        #calculo Ly=b con este b
        y=np.zeros(N) #transpuesto deberia
        for i in range(N):
            for j in range(N):
                if i>=j:
                    y[i]=(b[i]/(matrizL[i][i]))-(1/(matrizL[i][i])*sumatoria1(i,y))
        #calculo Ux=y
        xk=np.zeros(N)
        
        for i in range(N-1,-1,-1):
            for j in range(N-1,-1,-1):      
                if i<=j:
                    xk[i]=(y[i]/(matrizU[i][i]))-(1/(matrizU[i][i])*sumatoria2(i,xk))
        
        x[:,canonico]=xk
        canonico+=1
        #calculo norma infinito

    difmatrices = x-Hinversa
    normainf1 = normaInfinito(difmatrices, N)
    normainf2 = normaInfinito(Hinversa, N)
    Error = normainf1/normainf2
    logerror = math.log10(Error)
    numeroCondicion = normaInfinito(hilbertMatrix,N)*normaInfinito(Hinversa,N)
    logNumeroCondicion = math.log10(numeroCondicion)
    plt.stem(N,logNumeroCondicion, markerfmt='ro',linefmt="ro") 
    plt.stem(N,logerror)
    dfinv = pd.DataFrame(x)
    csv = 'matrizinv'+str(N)+'.csv'
    dfinv.to_csv(csv)
    return(logerror,numeroCondicion)


def calculoerrorConPerturbaciones(N,Hinversa):
    'Defino la matriz de hilbert'
    arreglos =[]
    
    for i in range(N):
        subarreglo = [0 for i in range(N)]
    
        for j in range(N):
            constanteHilbert = (i+1)+(j+1)-1
            subarreglo[j] = (1/(constanteHilbert))
        arreglos.append(subarreglo)
    
    hilbertMatrix = np.array(arreglos)  
    matrizU = hilbertMatrix.copy()
    
    #defino las matrices L 
    matrizL= np.identity(N)
    'defino la matriz L'
    columna=0
    while columna < N:
        for fila in range (N):
            if fila>columna:
                m_ia=(matrizU[fila][columna])/(matrizU[columna][columna])
                for j in range(1,N-columna):
                    matrizU[fila][columna+j] = matrizU[fila][columna+j]-m_ia*matrizU[columna][columna+j]
                matrizL[fila][columna]=m_ia
        columna+=1
    'defino la matriz U'
    for fil in range(N):
            for col in range(N):
                if (col < fil) :
                    matrizU[fil,col] = 0
    x= np.zeros((N,N))
    def sumatoria1(i,y):
        acum=0
        for j in range(i):
            acum+=matrizL[i][j]*y[j]
        return(acum)
    def sumatoria2(i,xk):
        acum=0
        for j in range(i+1,N):
            acum+=matrizU[i][j]*xk[j]
        return(acum)
    inversas=[]
    logaritmoError=[]
    for p in range(1,8,2):      
        canonico=0
        while canonico < N:
            #creo el termino independiente
            b=np.zeros(N) #obs deberia quedar transpuesto
            for i in range(N):
                if i==canonico:
                    b[i]=1+((-1)**(i+1))*10**(-p)
                else:
                    b[i]=((-1)**(i+1))*10**(-p)
            #calculo Ly=b con este b
            y=np.zeros(N) #transpuesto deberia
            for i in range(N):
                for j in range(N):
                    if i>=j:
                        y[i]=(b[i]/(matrizL[i][i]))-(1/(matrizL[i][i])*sumatoria1(i,y))
            #calculo Ux=y
            xk=np.zeros(N)
            
            for i in range(N-1,-1,-1):
                for j in range(N-1,-1,-1):      
                    if i<=j:
                        xk[i]=(y[i]/(matrizU[i][i]))-(1/(matrizU[i][i])*sumatoria2(i,xk))
            
            x[:,canonico]=xk
            canonico+=1
            #calculo norma infinito
        
        difmatrices = x-Hinversa
        normainf1 = normaInfinito(difmatrices, N)
        normainf2 = normaInfinito(Hinversa, N)
        Error = normainf1/normainf2
        logerror = math.log10(Error)
        numeroCondicion = normaInfinito(hilbertMatrix,N)*normaInfinito(Hinversa,N)
        logNumeroCondicion = math.log10(numeroCondicion)
        logaritmoError.append(logerror)
    return(logaritmoError)   

'calculo norma infinito'
def normaInfinito(matriz,N):
    max=0
    for fila in range (N):
        sumafila= 0
        for columna in range (N):
            sumafila += abs(matriz[fila][columna])
        if max < sumafila:
            max = sumafila
    return(max)

'ingreso de las matrices exactas'
Hinversa4=np.array([[16,-120,240,-140],[-120,1200,-2700,1680],[240,-2700,6480,-4200],[-140,1680,-4200,2800]])
Hinversa5=np.array([[25,-300,1050,-1400,630],[-300,4800,-18900,26880,-12600],[1050,-18900,79380,-117600,56700],[-1400,26880,-117600,179200,-88200],[630,-12600,56700,-88200,44100]])
Hinversa6=np.array([[36,-630,3360,-7560,7560,-2772],[-630,14700,-88200,211680,-220500,83160],[3360,-88200,564480,-1411200,1512000,-582120],[-7560,211680,-1411200,3628800,-3969000,1552320],[7560,-220500,1512000,-3969000,4410000,-1746360],[-2772,83160,-582120,1552320,-1746360,698544]])
Hinversa7=np.array([[49,-1176,8820,-29400,48510,-38808,12012],[-1176,37632,-317520,1128960,-1940400,1596672,-504504],[8820,-317520,2857680,-10584000,18711000,-15717240,5045040],[-29400,1128960,-10584000,40320000,-72765000,62092800,-20180160],[48510,-1940400,18711000,-72765000,133402500,-115259760,37837800],[-38808,1596672,-15717240,62092800,-115259760,100590336,-33297264],[12012,-504504,5045040,-20180160,37837800,-33297264,11099088]])
Hinversa8=np.array([[64 ,-2016, 20160, -92400, 221760, -288288, 192192, -51480],[-2016, 84672, -952560, 4656960, -11642400, 15567552, -10594584, 2882880],[20160, -952560, 11430720, -58212000, 149688000, -204324120, 141261120, -38918880],[-92400, 4656960, -58212000, 304920000, -800415000, 1109908800, -776936160, 216216000],[221760, -11642400, 149688000, -800415000, 2134440000,-2996753760, 2118916800, -594594000],[-288288, 15567552, -204324120, 1109908800, -2996753760, 4249941696, -3030051024,856215360],[192192, -10594584, 141261120, -776936160, 2118916800,-3030051024, 2175421248, -618377760],[-51480, 2882880, -38918880, 216216000, -594594000, 856215360, -618377760, 176679360]])

logaritmoErrorExacto=[0,0,0,0,0]
numeroCondicionMatrices=[0,0,0,0,0]
logaritmoErrorExacto[0],numeroCondicionMatrices[0]=calculoerror(4,Hinversa4)
logaritmoErrorExacto[1],numeroCondicionMatrices[1]=calculoerror(5,Hinversa5)
logaritmoErrorExacto[2],numeroCondicionMatrices[2]=calculoerror(6,Hinversa6)
logaritmoErrorExacto[3],numeroCondicionMatrices[3]=calculoerror(7,Hinversa7)
logaritmoErrorExacto[4],numeroCondicionMatrices[4]=calculoerror(8,Hinversa8)
plt.xticks(range(4,9))
plt.title('log10(error) en funcion de N')
plt.xlabel('N')
plt.ylabel('log(error)')
plt.axhline(y=0)
plt.show()


logaritmoError4=calculoerrorConPerturbaciones(4,Hinversa4)
logaritmoError5=calculoerrorConPerturbaciones(5,Hinversa5)
logaritmoError6=calculoerrorConPerturbaciones(6,Hinversa6)
logaritmoError7=calculoerrorConPerturbaciones(7,Hinversa7)
logaritmoError8=calculoerrorConPerturbaciones(8,Hinversa8)
logaritmoError=[logaritmoError4,logaritmoError5,logaritmoError6,logaritmoError7,logaritmoError8]

for i in range(4):
    p=(2*i)+1
    plt.ylim(-14)
    plt.xlim(3.9,8.1)
    plt.xticks(range(4,9))
    x= (4,5,6,7,8)
    y= (logaritmoError4[i],logaritmoError5[i],logaritmoError6[i],logaritmoError7[i],logaritmoError8[i])
    plt.plot(x,y,'purple')
    plt.plot(x, y, 'o')
    z=(logaritmoErrorExacto[0],logaritmoErrorExacto[1],logaritmoErrorExacto[2],logaritmoErrorExacto[3],logaritmoErrorExacto[4])
    plt.plot(x,z,'red')
    plt.plot(x, z, 'o')
    plt.title('log10(error) en funcion de N con p='+str(p))
    plt.xlabel('N')
    plt.ylabel('log(error)')
    plt.axhline(y=0)
    plt.show()

for i in range(4):
    print("Numero de condicion de: "+str(i+4)+" "+str(numeroCondicionMatrices[i]))
'''
for j in range(5):
    archivo = 'matrizinv'+str(j+4)+'.csv'
    with open(archivo) as file:
        file.readline()

        for line in file:
            line.replace('\n','')
            linea = line.split(',')
            linea.pop(0)
            for i in range(len(linea)):
                linea[i]=(round(float(linea[i]),5))
                linea[i] = str('{:.1f}'.format(linea[i]))
            print('&'.join(linea))'''
