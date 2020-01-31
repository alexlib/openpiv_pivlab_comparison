# -*- coding: utf-8 -*-
"""
Created on Sun Jul 21 22:45:53 2019

@author: Theo
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import scipy.misc
import scipy
import glob
from scipy.optimize import curve_fit, minimize
from scipy.stats import norm
import matplotlib.mlab as mlab
from scipy import signal
from scipy import asarray as ar, exp


def Trans_PIVlab(Data):

    Data = np.delete(Data[:], np.where(
        np.equal(np.isnan(Data[:, :, 2]), True)), axis=1)
    Data[:, :, 1] = Data[:, :, 1] * (-1)
    Data[:, :, 3] = Data[:, :, 3] * (-1)
    Data[:,
         :,
         1] = np.flipud(Data[:,
                             :,
                             1]) + np.absolute(np.min(Data[:,
                                                           :,
                                                           1])) + np.absolute(np.max(Data[:,
                                                                                          :,
                                                                                          1]))
    Data[:, :, 3] = np.flipud(Data[:, :, 3])
    Data[:, :, [2, 3]] = np.absolute(Data[:, :, [2, 3]])
    return(Data)


def Trans_PIVlab_2(Data):

    Data = np.delete(Data[:], np.where(
        np.equal(np.isnan(Data[:, :, 2]), True)), axis=1)
    # ata[:,:,1]=Data[:,:,1]*(-1)
    Data[:, :, 3] = Data[:, :, 3] * (-1)
    #Data[:,:,1]=np.flipud(Data[:,:,1])+np.absolute(np.min(Data[:,:,1]))+ np.absolute(np.max(Data[:,:,1]))
    # Data[:,:,3]=np.flipud(Data[:,:,3])
    # Data[:,:,[2,3]]=np.absolute(Data[:,:,[2,3]])
    return(Data)


def Trans_OpenPIV(Data):
    Data = np.delete(Data[:], np.where(
        np.equal(np.isnan(Data[:, :, 2]), True)), axis=1)
    Data = Data[:, :, [0, 1, 2, 3]]
    Data[:, :, 3] = Data[:, :, 3]  # *(-1)#correction of old results
    Data[:, :, [2, 3]] = np.absolute(Data[:, :, [2, 3]])
    return(Data)


def averaged_profil_gaus_fit(Data_Set, DATA_AVG, profil_position=0):
    x = profil_position
    "setting position x position for profil"

    "calculating average"

    #Data=np.delete(Data_Set, np.where(np.not_equal(Data_Set[:,:,0],x)), axis=1)
    Data_avg_p = np.delete(DATA_AVG, np.where(
        np.not_equal(DATA_AVG[:, 0], x)), axis=0)
    "deleting data that is not needed"

    vel_u = Data_avg_p[:, 2]
    cor_x = Data_avg_p[:, 1]
    " stuff for gaus fit"

    mean = np.mean(cor_x)
    sigma = 100  # *np.std(Data[:,:,2])
    "calculating mean an standard deviation"

    def gaus(x, a, x0, sigma):
        return a * exp(-(x - x0)**2 / (2 * sigma**2))
    "function whith gaus equation"

    popt, pcov = curve_fit(
        gaus, cor_x, vel_u, p0=[
            1, mean, sigma], maxfev=1000)
    "optimation of gaus fit"

    cor_max = cor_x[np.argmax(gaus(cor_x, *popt))]
    "get index of extrempoint of gausfit"
    return(gaus(cor_x, *popt), cor_max, popt)
    # return(cor_y,vel_u,gaus(cor_y,*popt),cor_max,Data_avg)


def instant_vel(Data, Data_num, position=0):

    Data = Data[Data_num, :, :]
    Data = np.delete(Data, np.where(
        np.not_equal(Data[:, 0], position)), axis=0)

    return(Data)


# Load PIVlab data
PIVlab_DATA = np.array([np.array(np.loadtxt(f, skiprows=1, dtype=float, delimiter="\t", unpack=False)) for f in glob.glob('vectorfields_folder_pivlab')])
# arrange PIVlab data for comparison
PIVlab_DATA_BACK = Trans_PIVlab_2(PIVlab_DATA)
PIVlab_DATA = Trans_PIVlab(PIVlab_DATA)
# claculate averaged data
PIVlab_DATA_AVG = np.mean(PIVlab_DATA, axis=0)

# load OpenPIV data
OpenPIV_DATA = np.array(
    [
        np.array(
            np.loadtxt(
                g,
                skiprows=1,
                dtype=float,
                delimiter="\t",
                unpack=False)) for g in glob.glob('vectorfields_folder_openpiv')])
# arrange OpenPIV data for comparison
OpenPIV_DATA_BACK = OpenPIV_DATA
OpenPIV_DATA = Trans_OpenPIV(OpenPIV_DATA)
# claculate averaged data
OpenPIV_DATA_AVG = np.mean(OpenPIV_DATA, axis=0)

# select the vectorfiled for the non-time averaged comparison
VECTORFIELD_NUMBER = 40
# select the X_CORRDINATES of the displacement profiles (in pixel)

PIVlab_POSITION = (201, 401, 601, 801)
OpenPIV_POSITION = (200, 400, 600, 800)


#  Gaussian fit for Position 1
U_VEL_PIVlab_P1_GAUS, CORR_Y_MAXIMUM_PIVlab_P1, POPT_PIVlab_1 = averaged_profil_gaus_fit(
    PIVlab_DATA, PIVlab_DATA_AVG, profil_position=PIVlab_POSITION[0])
U_VEL_OpenPIV_P1_GAUS, CORR_Y_MAXIMUM_OpenPIV_P1, POPT_OpenPIV_1 = averaged_profil_gaus_fit(
    OpenPIV_DATA, OpenPIV_DATA_AVG, profil_position=OpenPIV_POSITION[0])
#  Gaussian fit for Position 2
U_VEL_PIVlab_P2_GAUS, CORR_Y_MAXIMUM_PIVlab_P2, POPT_PIVlab_2 = averaged_profil_gaus_fit(
    PIVlab_DATA, PIVlab_DATA_AVG, profil_position=PIVlab_POSITION[1])
U_VEL_OpenPIV_P2_GAUS, CORR_Y_MAXIMUM_OpenPIV_P2, POPT_OpenPIV_2 = averaged_profil_gaus_fit(
    OpenPIV_DATA, OpenPIV_DATA_AVG, profil_position=OpenPIV_POSITION[1])
#  Gaussian fit for Position 3
U_VEL_PIVlab_P3_GAUS, CORR_Y_MAXIMUM_PIVlab_P3, POPT_PIVlab_3 = averaged_profil_gaus_fit(
    PIVlab_DATA, PIVlab_DATA_AVG, profil_position=PIVlab_POSITION[2])
U_VEL_OpenPIV_P3_GAUS, CORR_Y_MAXIMUM_OpenPIV_P3, POPT_OpenPIV_3 = averaged_profil_gaus_fit(
    OpenPIV_DATA, OpenPIV_DATA_AVG, profil_position=OpenPIV_POSITION[2])
#  Gaussian fit for Position 4
U_VEL_PIVlab_P4_GAUS, CORR_Y_MAXIMUM_PIVlab_P4, POPT_PIVlab_4 = averaged_profil_gaus_fit(
    PIVlab_DATA, PIVlab_DATA_AVG, profil_position=PIVlab_POSITION[3])
U_VEL_OpenPIV_P4_GAUS, CORR_Y_MAXIMUM_OpenPIV_P4, POPT_OpenPIV_4 = averaged_profil_gaus_fit(
    OpenPIV_DATA, OpenPIV_DATA_AVG, profil_position=OpenPIV_POSITION[3])


SMALL_SIZE = 20
MEDIUM_SIZE = 30
BIGGER_SIZE = 35

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
plt.close('all')
fig_1 = plt.figure(1, figsize=(20, 10))
# Position 1
ax1_4 = plt.plot((OpenPIV_DATA_AVG[np.squeeze(np.where(OpenPIV_DATA_AVG[:, 0] == OpenPIV_POSITION[0])), 1] -
                  CORR_Y_MAXIMUM_OpenPIV_P1) /
                 np.power(np.log(2) *
                          np.power(POPT_OpenPIV_1[2], 2), 0.5), OpenPIV_DATA_AVG[np.squeeze(np.where(OpenPIV_DATA_AVG[:, 0] == OpenPIV_POSITION[0])), 2] /
                 OpenPIV_DATA_AVG[np.squeeze(np.where(OpenPIV_DATA_AVG[:, 0] == OpenPIV_POSITION[0])), 2].max(), 'r.:', label='$Profil$ $positon: ' +
                 str(OpenPIV_POSITION[0]) +
                 r'\,px$')  # averaged velocity
# Position 2
ax2_4 = plt.plot((OpenPIV_DATA_AVG[np.squeeze(np.where(OpenPIV_DATA_AVG[:, 0] == OpenPIV_POSITION[1])), 1] -
                  CORR_Y_MAXIMUM_OpenPIV_P2) /
                 np.power(np.log(2) *
                          np.power(POPT_OpenPIV_2[2], 2), 0.5), OpenPIV_DATA_AVG[np.squeeze(np.where(OpenPIV_DATA_AVG[:, 0] == OpenPIV_POSITION[1])), 2] /
                 OpenPIV_DATA_AVG[np.squeeze(np.where(OpenPIV_DATA_AVG[:, 0] == OpenPIV_POSITION[1])), 2].max(), 'y.:', label='$Profil$ $positon: ' +
                 str(OpenPIV_POSITION[1]) +
                 r'\,px$')  # averaged velocity
# Position 3
ax3_4 = plt.plot((OpenPIV_DATA_AVG[np.squeeze(np.where(OpenPIV_DATA_AVG[:, 0] == OpenPIV_POSITION[2])), 1] -
                  CORR_Y_MAXIMUM_OpenPIV_P3) /
                 np.power(np.log(2) *
                          np.power(POPT_OpenPIV_3[2], 2), 0.5), OpenPIV_DATA_AVG[np.squeeze(np.where(OpenPIV_DATA_AVG[:, 0] == OpenPIV_POSITION[2])), 2] /
                 OpenPIV_DATA_AVG[np.squeeze(np.where(OpenPIV_DATA_AVG[:, 0] == OpenPIV_POSITION[2])), 2].max(), 'm.:', label='$Profil$ $positon: ' +
                 str(OpenPIV_POSITION[2]) +
                 r'\,px$')  # averaged velocity
# Position 4
ax4_4 = plt.plot((OpenPIV_DATA_AVG[np.squeeze(np.where(OpenPIV_DATA_AVG[:, 0] == OpenPIV_POSITION[3])), 1] -
                  CORR_Y_MAXIMUM_OpenPIV_P4) /
                 np.power(np.log(2) *
                          np.power(POPT_OpenPIV_4[2], 2), 0.5), OpenPIV_DATA_AVG[np.squeeze(np.where(OpenPIV_DATA_AVG[:, 0] == OpenPIV_POSITION[3])), 2] /
                 OpenPIV_DATA_AVG[np.squeeze(np.where(OpenPIV_DATA_AVG[:, 0] == OpenPIV_POSITION[3])), 2].max(), 'k.:', label='$Profil$ $positon: ' +
                 str(OpenPIV_POSITION[3]) +
                 r'\,px$')  # averaged velocity
plt.xlabel(r'$\hat{Y}$')
plt.ylabel(r'$\hat{U}$')
plt.legend()
plt.show()
plt.savefig('figure_1.pdf')

# plt.close('all')
PIVlab_DATA_V = np.squeeze(PIVlab_DATA[VECTORFIELD_NUMBER, :, :])
OpenPIV_DATA_V = np.squeeze(OpenPIV_DATA[VECTORFIELD_NUMBER, :, :])
# plt.close('all')
fig_2 = plt.figure(2, figsize=(10, 10))
# Position 1
Ax1_1 = plt.plot((PIVlab_DATA_V[np.squeeze(np.where(PIVlab_DATA_V[:,
                                                                  0] == PIVlab_POSITION[0])),
                                1] - CORR_Y_MAXIMUM_PIVlab_P1),
                 PIVlab_DATA_V[np.squeeze(np.where(PIVlab_DATA_V[:,
                                                                 0] == PIVlab_POSITION[0])),
                               2],
                 'b.:',
                 label='$PIVlab$')  # instantaneous velocity
Ax2_2 = plt.plot((OpenPIV_DATA_V[np.squeeze(np.where(OpenPIV_DATA_V[:,
                                                                    0] == OpenPIV_POSITION[0])),
                                 1] - CORR_Y_MAXIMUM_OpenPIV_P1),
                 OpenPIV_DATA_V[np.squeeze(np.where(OpenPIV_DATA_V[:,
                                                                   0] == OpenPIV_POSITION[0])),
                                2],
                 'r.:',
                 label='$OpenPIV$')  # instantaneous velocity
plt.ylim(0, 10)
plt.xlabel(r'$ y\ position \ in \ pixel$')
plt.ylabel(r'$\mid u \mid\ in \ pixel\,/\,frame$')
plt.legend()
plt.savefig('figure_2.pdf')
fig_3 = plt.figure(3, figsize=(10, 10))
# Position 2
Ax3_1 = plt.plot((PIVlab_DATA_V[np.squeeze(np.where(PIVlab_DATA_V[:,
                                                                  0] == PIVlab_POSITION[1])),
                                1] - CORR_Y_MAXIMUM_PIVlab_P2),
                 PIVlab_DATA_V[np.squeeze(np.where(PIVlab_DATA_V[:,
                                                                 0] == PIVlab_POSITION[1])),
                               2],
                 'b.:',
                 label='$PIVlab$')  # instantaneous velocity
Ax4_2 = plt.plot((OpenPIV_DATA_V[np.squeeze(np.where(OpenPIV_DATA_V[:,
                                                                    0] == OpenPIV_POSITION[1])),
                                 1] - CORR_Y_MAXIMUM_OpenPIV_P2),
                 OpenPIV_DATA_V[np.squeeze(np.where(OpenPIV_DATA_V[:,
                                                                   0] == OpenPIV_POSITION[1])),
                                2],
                 'r.:',
                 label='$OpenPIV$')  # instantaneous velocity
plt.ylim(0, 10)
plt.xlabel(r'$ y\ position \ in \ pixel$')
plt.ylabel(r'$\mid u \mid\ in \ pixel\,/\,frame$')
plt.legend()
plt.savefig('figure_3.pdf')
fig_4 = plt.figure(4, figsize=(10, 10))
# Position 2
Ax5_1 = plt.plot((PIVlab_DATA_V[np.squeeze(np.where(PIVlab_DATA_V[:,
                                                                  0] == PIVlab_POSITION[2])),
                                1] - CORR_Y_MAXIMUM_PIVlab_P3),
                 PIVlab_DATA_V[np.squeeze(np.where(PIVlab_DATA_V[:,
                                                                 0] == PIVlab_POSITION[2])),
                               2],
                 'b.:',
                 label='$PIVlab$')  # instantaneous velocity
Ax6_2 = plt.plot((OpenPIV_DATA_V[np.squeeze(np.where(OpenPIV_DATA_V[:,
                                                                    0] == OpenPIV_POSITION[2])),
                                 1] - CORR_Y_MAXIMUM_OpenPIV_P3),
                 OpenPIV_DATA_V[np.squeeze(np.where(OpenPIV_DATA_V[:,
                                                                   0] == OpenPIV_POSITION[2])),
                                2],
                 'r.:',
                 label='$OpenPIV$')  # instantaneous velocity
plt.ylim(0, 10)
plt.xlabel(r'$ y\ position \ in \ pixel$')
plt.ylabel(r'$\mid u \mid\ in \ pixel\,/\,frame$')
plt.legend()
plt.savefig('figure_4.pdf')
fig_5 = plt.figure(5, figsize=(10, 10))
# Position 2
Ax7_1 = plt.plot((PIVlab_DATA_V[np.squeeze(np.where(PIVlab_DATA_V[:,
                                                                  0] == PIVlab_POSITION[3])),
                                1] - CORR_Y_MAXIMUM_PIVlab_P4),
                 PIVlab_DATA_V[np.squeeze(np.where(PIVlab_DATA_V[:,
                                                                 0] == PIVlab_POSITION[3])),
                               2],
                 'b.:',
                 label='$PIVlab$')  # instantaneous velocity
Ax8_2 = plt.plot((OpenPIV_DATA_V[np.squeeze(np.where(OpenPIV_DATA_V[:,
                                                                    0] == OpenPIV_POSITION[3])),
                                 1] - CORR_Y_MAXIMUM_OpenPIV_P4),
                 OpenPIV_DATA_V[np.squeeze(np.where(OpenPIV_DATA_V[:,
                                                                   0] == OpenPIV_POSITION[3])),
                                2],
                 'r.:',
                 label='$OpenPIV$')  # instantaneous velocity
plt.ylim(0, 10)
plt.xlabel(r'$ y\ position \ in \ pixel$')
plt.ylabel(r'$\mid u \mid\ in \ pixel\,/\,frame$')
plt.legend()
plt.savefig('figure_5.pdf')
plt.show()


PIVlab_DATA_V_BACK = np.squeeze(PIVlab_DATA_BACK[VECTORFIELD_NUMBER, :, :])
OpenPIV_DATA_V_BACK = np.squeeze(OpenPIV_DATA_BACK[VECTORFIELD_NUMBER, :, :])

# plt.close('all')
fig_6 = plt.figure(6, figsize=(10, 10))
# Position 1
AX1_1 = plt.plot((PIVlab_DATA_V[np.squeeze(np.where(PIVlab_DATA_V[:,
                                                                  0] == PIVlab_POSITION[0])),
                                1] - CORR_Y_MAXIMUM_PIVlab_P1),
                 PIVlab_DATA_V_BACK[np.squeeze(np.where(PIVlab_DATA_V_BACK[:,
                                                                           0] == PIVlab_POSITION[0])),
                                    3],
                 'b.:',
                 label='$PIVlab$')  # instantaneous velocity
AX2_2 = plt.plot((OpenPIV_DATA_V[np.squeeze(np.where(OpenPIV_DATA_V[:,
                                                                    0] == OpenPIV_POSITION[0])),
                                 1] - CORR_Y_MAXIMUM_OpenPIV_P1),
                 OpenPIV_DATA_V_BACK[np.squeeze(np.where(OpenPIV_DATA_V_BACK[:,
                                                                             0] == OpenPIV_POSITION[0])),
                                     3],
                 'r.:',
                 label='$OpenPIV$')  # instantaneous velocity
plt.ylim(-5, 5)
plt.xlabel(r'$ y\ position \ in \ pixel$')
plt.ylabel(r'$v\ in \ pixel\,/\,frame$')
plt.legend()
plt.savefig('figure_6.pdf')
fig_7 = plt.figure(7, figsize=(10, 10))
# Position 2
AX3_1 = plt.plot((PIVlab_DATA_V[np.squeeze(np.where(PIVlab_DATA_V[:,
                                                                  0] == PIVlab_POSITION[1])),
                                1] - CORR_Y_MAXIMUM_PIVlab_P1),
                 PIVlab_DATA_V_BACK[np.squeeze(np.where(PIVlab_DATA_V_BACK[:,
                                                                           0] == PIVlab_POSITION[1])),
                                    3],
                 'b.:',
                 label='$PIVlab$')  # instantaneous velocity
AX4_2 = plt.plot((OpenPIV_DATA_V[np.squeeze(np.where(OpenPIV_DATA_V[:,
                                                                    0] == OpenPIV_POSITION[1])),
                                 1] - CORR_Y_MAXIMUM_OpenPIV_P1),
                 OpenPIV_DATA_V_BACK[np.squeeze(np.where(OpenPIV_DATA_V_BACK[:,
                                                                             0] == OpenPIV_POSITION[1])),
                                     3],
                 'r.:',
                 label='$OpenPIV$')  # instantaneous velocity
plt.ylim(-5, 5)
plt.xlabel(r'$ y\ position \ in \ pixel$')
plt.ylabel(r'$v\ in \ pixel\,/\,frame$')
plt.legend()
plt.savefig('figure_7.pdf')
fig_8 = plt.figure(8, figsize=(10, 10))
# Position 2
AX5_1 = plt.plot((PIVlab_DATA_V[np.squeeze(np.where(PIVlab_DATA_V[:,
                                                                  0] == PIVlab_POSITION[2])),
                                1] - CORR_Y_MAXIMUM_PIVlab_P1),
                 PIVlab_DATA_V_BACK[np.squeeze(np.where(PIVlab_DATA_V_BACK[:,
                                                                           0] == PIVlab_POSITION[2])),
                                    3],
                 'b.:',
                 label='$PIVlab$')  # instantaneous velocity
AX6_2 = plt.plot((OpenPIV_DATA_V[np.squeeze(np.where(OpenPIV_DATA_V[:,
                                                                    0] == OpenPIV_POSITION[2])),
                                 1] - CORR_Y_MAXIMUM_OpenPIV_P1),
                 OpenPIV_DATA_V_BACK[np.squeeze(np.where(OpenPIV_DATA_V_BACK[:,
                                                                             0] == OpenPIV_POSITION[2])),
                                     3],
                 'r.:',
                 label='$OpenPIV$')  # instantaneous velocity
plt.ylim(-5, 5)
plt.xlabel(r'$ y\ position \ in \ pixel$')
plt.ylabel(r'$v\ in \ pixel\,/\,frame$')
plt.legend()
plt.savefig('figure_8.pdf')
fig_9 = plt.figure(9, figsize=(10, 10))
# Position 2
AX7_1 = plt.plot((PIVlab_DATA_V[np.squeeze(np.where(PIVlab_DATA_V[:,
                                                                  0] == PIVlab_POSITION[3])),
                                1] - CORR_Y_MAXIMUM_PIVlab_P1),
                 PIVlab_DATA_V_BACK[np.squeeze(np.where(PIVlab_DATA_V_BACK[:,
                                                                           0] == PIVlab_POSITION[3])),
                                    3],
                 'b.:',
                 label='$PIVlab$')  # instantaneous velocity
AX8_2 = plt.plot((OpenPIV_DATA_V[np.squeeze(np.where(OpenPIV_DATA_V[:,
                                                                    0] == OpenPIV_POSITION[3])),
                                 1] - CORR_Y_MAXIMUM_OpenPIV_P1),
                 OpenPIV_DATA_V_BACK[np.squeeze(np.where(OpenPIV_DATA_V_BACK[:,
                                                                             0] == OpenPIV_POSITION[3])),
                                     3],
                 'r.:',
                 label='$OpenPIV$')  # instantaneous velocity
plt.ylim(-5, 5)
plt.xlabel(r'$ y\ position \ in \ pixel$')
plt.ylabel(r'$v\ in \ pixel\,/\,frame$')
plt.legend()
plt.savefig('figure_9.pdf')
plt.show()
