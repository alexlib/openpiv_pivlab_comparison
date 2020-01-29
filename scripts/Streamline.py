# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 16:41:07 2020

@author: arkajyoti.ghoshal
"""

from pivpy import io, graphics, pivpy
import pkg_resources as pkg
import matplotlib.pyplot as plt
import numpy as np
import timeit

import pandas as pd
import matplotlib.gridspec as gridspec
start = timeit.default_timer()

with open ('C:\PhD\PIV\Anupam_hak\Original frames\Ext_srch\pass120.txt') as f:
    r=pd.read_csv(f, header= None, delimiter='\t')
    
    
x=list(r.iloc[:,0])
y=list(r.iloc[:,1])
u=list(r.iloc[:,2])
v=list(r.iloc[:,3])


a,b,c,d=(max(x),max(y),min(x),min(y))

xx=np.linspace(c,a,48)
#good print(xx)
yy=(np.flipud(np.linspace(d,b,42)))
yz=((np.linspace(d,b,42)))
#print(yz-yy)
X,Y=np.meshgrid(xx,yy)
X,Y1=np.meshgrid(xx,yz)
#print(888888888888,X)
#meshgrid made with reversed y
#tot number of points = 48x42 as both x and y repeats itself at 48,42 intervals
#print(V)
#plt.plot(X,Y)
#plt.show
#gs = gridspec.GridSpec(nrows=3, ncols=2, height_ratios=[1, 1, 2])
'''fig,ax0 = plt.subplots(figsize=(12,12))
ax0 = fig.add_subplot(gs[0, 0])
ax0.streamplot(X, Y, U, V, density=[1.4, 1.4])'''
#plt.plot(X,Y, marker='.', color='k', linestyle='none')
#M=X+2
#print(X,M)
U1= np.array(u)
#print(U1)
V1=np.array(v)
U=U1.reshape((42,48)) #Y,X
V=V1.reshape((42,48))
speed=np.sqrt(U**2+V**2)
lw = 5*speed / speed.max()

fig,ax = plt.subplots(figsize=(48,30))
#ax0 = fig.add_subplot(gs[0, 0])
#ax[0,0].streamplot(X, Y, U, V, density=[1, 1],linewidth=2*lw)
#ax[1,1].streamplot(X, Y1, U, V, density=[1, 1])
#ax0.quiver(X,Y,U,V,color='r')

#strm = ax[0,1].streamplot(X, Y, U, V, color=U, linewidth=2*lw, cmap='autumn')
#fig.colorbar(strm.lines)
#ax1.set_title('Varying Color')


ax.streamplot(X, Y, U, V, density=1, color='k', linewidth=2*lw)