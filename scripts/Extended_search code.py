# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 16:39:46 2020

@author: arkajyoti.ghoshal
"""

import numpy as np
import multiprocessing
import pylab
#print(multiprocessing.cpu_count())
#basic code
from openpiv import tools, process, validation, filters, scaling 

import matplotlib.pyplot as plt
from openpiv import  pyprocess
folder='C:\PhD\PIV\Anupam_hak\Original frames\\'
#name=str(folder[:-1]+'\\Ext_srch\\pass1.txt')
for i in range(0,9):
    j=i+1
    
    xx=str(0)+str(0)+str(0)+str(i)
    xy=str(0)+str(0)+str(0)+str(j)
    frame_a  = tools.imread( str(folder+'PIV-Sample-Ghoshal'+str(xx)+'.tif' ))
    frame_b  = tools.imread( str(folder+'PIV-Sample-Ghoshal'+str(xy)+'.tif' ))
    name=str(folder[:-1]+'\Ext_srch\\pass1'+str(i)+'.txt')

    #fig,ax = plt.subplots(1,2)
#ax[0].imshow(frame_a,cmap=plt.cm.gray)
#ax[1].imshow(frame_b,cmap=plt.cm.gray)
    winsize = 40 # pixels
    searchsize = 80  # pixels, search in image B
    overlap = 30 # pixels
    dt = .06 # sec
    u0, v0, sig2noise = process.extended_search_area_piv( frame_a.astype(np.int32), frame_b.astype(np.int32), window_size=winsize, overlap=overlap, dt=dt, search_area_size=searchsize, sig2noise_method='peak2peak' )
    x, y = process.get_coordinates( image_size=frame_a.shape, window_size=winsize, overlap=overlap )
    u1, v1, mask = validation.sig2noise_val( u0, v0, sig2noise, threshold = 2 )
    u2, v2 = filters.replace_outliers( u1, v1, method='localmean', max_iter=10, kernel_size=2)
    x, y, u3, v3 = scaling.uniform(x, y, u2, v2, scaling_factor = 100 )

    tools.save(x, y, u3, v3, mask, name,delimiter='\t' )
    #tools.display_vector_field(name, scale=100, width=0.0025)
    fig,ax = plt.subplots(figsize=(12,12))
#ax.imshow(plt.imread('C:\PhD\Images\FD2_1\\'+str(1)+' '+'('+str(px)+')'+'.tif')[::-1,:],cmap=plt.cm.gray,origin='lower')
    ax.quiver(x*100,y*100,u3*100,v3*100,color='green',units='xy',angles='xy',scale=.5)
    plt.savefig(folder+'\Ext_srch\\pass1img'+str(xx)+'.png')


