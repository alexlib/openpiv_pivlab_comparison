# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 16:14:40 2020

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

with open ('C:\PhD\PIV\Anupam_hak\Original frames\Ext_srch\\pass120.txt') as f:
    r=pd.read_csv(f, header= None, names=['x','y','u','v','m'],delimiter='\t')
    
    

    
#x= df[['x']]
#y=df[['y']]
#u=df[['u']]
#v=df[['v']]

tmp=r.values
#plt.quiver(x,y,u,v)
    #print(tmp)
x,y,u,v = tmp[:,0],tmp[:,1],tmp[:,2],tmp[:,3] #(each col as a list)
#print(x)
   
rows = np.unique(y).shape[0]
cols = np.unique(x).shape[0]
x1 = x.reshape(rows,cols)
y1 = y.reshape(rows,cols)
u1 = u.reshape(rows,cols)
v1 = v.reshape(rows,cols)

d = io.from_arrays(x1,y1,u1,v1,np.ones_like(u1))
fig,axes=plt.subplots(figsize=(28,20))
graphics.quiver(d.isel(t=0),nthArr=3, arrScale=10) #arrScale scales arrows

d.piv.vec2scal(property='curl')
fig, ax = plt.subplots(figsize=(56,40))
graphics.contour_plot(d)

from scipy.ndimage.filters import gaussian_filter

d.piv.vorticity()
tmp2 =  gaussian_filter(d.isel(t=0)['w'],0.7)

fig, ax = plt.subplots(figsize=(8,6))

levels = np.linspace(np.min(tmp),np.max(tmp), 10)
#c = ax.contourf(r.x,r.y,tmp, levels=levels,
                 #cmap = plt.get_cmap('RdYlBu'))