# 6/7/2019
# find peaks for the average data (and all BPMs)
# modify plot (plot all types)

import time
import operator
import h5py
import heapq
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
#from scipy.signal import find_peaks

from cothread.catools import *
from cothread import Timedout

t0 = time.time()
# ----------------------------------------------
### get BPM data (from .h5, need to modify to get live data from PV)
##
##fid = h5py.File('SR_AllIDBPMs_FA_20190430_1405_44_220mA_10s_ops1.h5', 'r')
###fid = h5py.File('SR_AllIDBPMs_FA_20190222_1623_38_FOFBoffNoPinger1.h5', 'r')
##x_all = fid.get('faX')
##y_all = fid.get('faY')
##fid.close
##
##x_all = np.array(x_all)
##y_all = np.array(y_all)

# ----------------------------------------------
# create list of BPM PVs
print "start getting BPM FA data..."
fa_recordLen = 100000

p_index = ['1','2','3','4','5','6']
Cell_index = ['30','01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29']

prefix = []

for i in Cell_index:
	for j in p_index:
		prefix.append('SR:C'+i+'-BI{BPM:'+j+'}')

p_index=['7','8','9','10']
Cell_index=['04','07','12','19']
for i in Cell_index:
	for j in p_index:
		prefix.append('SR:C'+i+'-BI{BPM:'+j+'}')

p_index=['7','8']
Cell_index=['02','03','08','10','11','16','18','21','28']
for i in Cell_index:
	for j in p_index:
		prefix.append('SR:C'+i+'-BI{BPM:'+j+'}')

p_index=['7','8','9']
Cell_index=['05','17','23']
for i in Cell_index:
	for j in p_index:
		prefix.append('SR:C'+i+'-BI{BPM:'+j+'}')

# get BPM data from PV
FA_X=[]
FA_Y=[]
for i in range(len(prefix)):
	FA_X.append(prefix[i]+'FA-X')
	FA_Y.append(prefix[i]+'FA-Y')
    
x_all = caget(FA_X, count=fa_recordLen)
y_all = caget(FA_Y, count=fa_recordLen)

x_all = np.asarray(x_all)
y_all = np.asarray(y_all)

print "faX data:"
print x_all
print "faX size:"
print x_all.shape
print "\nfaY data:"
print y_all
print "faY size:"
print y_all.shape
print "\n"

# ----------------------------------------------
t = time.time() - t0
print "finish getting BPM FA data at t = %.f"  % (t)
print "start calculation..."

# generate disp bpm index (start from 0 to 179)
disp = []
non_disp = []
for i in range(1,181):
    if np.mod(i,6) == 3 or np.mod(i,6) == 4:
        disp.append(i-1)
    else:
        non_disp.append(i-1)
    
print "dispersive bpm index is (start from 0 to 179)"
print disp
print "\n"
print "non dispersive bpm index is"
print non_disp
print "\n"

# initialize
rf = 499.68e6
fs = rf/1320/38
n = x_all.shape
r = n[0] # number of BPM (rows 233)
c = n[1] # number of samples (columns 100000)

n_max = 5 # max n peaks
prom_x = 0.5e4 # prominence for finding x peaks
prom_y = 0.5e4 # prominence for finding y peaks
dist = 5 # dist between peaks

# pre-allowcate output
Pxx_all = []
Pyy_all = []
pks_freq_x = [0]*r
pks_freq_y = [0]*r
pks_hight_x = [0]*r
pks_hight_y = [0]*r
pks_n_freq_x = [0]*r
pks_n_freq_y = [0]*r
pks_n_hight_x = [0]*r
pks_n_hight_y = [0]*r
for i in range(0,r):
    x = x_all[i]
    y = y_all[i] 

    # ---------------- horizontal ------------------
    f, P = signal.welch(x, fs, window='hann', nperseg = c)
    Pxx_all = np.append(Pxx_all, P) # in 1D
    
    # find peaks of psd
    loc, _ = find_peaks(P, prominence = prom_x, distance = dist)
    pks = P[loc]

    # pick max N peaks
    loc_n = (-pks).argsort()
    loc_n = np.sort(loc_n[0:n_max])
    loc_n = loc_n.astype(int)

    # collect peaks data for all horizontal bpm
    pks_freq_x[i] = f[loc] # all pks frequency
    pks_hight_x[i] = pks # all pks hight
    pks_n_freq_x[i] = f[loc[loc_n]] # n max pks frequency
    pks_n_hight_x[i] = pks[loc_n] # n max pks hight
                     
    # ------------------ vertical -------------------
    _, P = signal.welch(y, fs, window='hann', nperseg = c)
    Pyy_all = np.append(Pyy_all, P) # in 1D
                     
    # find peaks of psd
    loc, _ = find_peaks(P, prominence = prom_y, distance = dist)
    pks = P[loc]

    # pick max N peaks
    loc_n = (-pks).argsort()
    loc_n = np.sort(loc_n[0:n_max])
    loc_n = loc_n.astype(int)

    # collect peaks data for all vertical bpm
    pks_freq_y[i] = f[loc] # all pks frequency
    pks_hight_y[i] = pks # all pks hight
    pks_n_freq_y[i] = f[loc[loc_n]] # n max pks frequency
    pks_n_hight_y[i] = pks[loc_n] # n max pks hight

t = time.time() - t0
print "checkpoint after loop t = %.f" % (t)

# reshape psd (1D -> 2D)
Pxx_all = np.transpose(Pxx_all.reshape(r,len(f)))
Pyy_all = np.transpose(Pyy_all.reshape(r,len(f)))

# group disp, non disp, insertion device (id)
Pxx_non_disp = Pxx_all[:,non_disp]
Pxx_disp = Pxx_all[:,disp]
Pxx_id = Pxx_all[:,180:r+1]
Pyy = Pyy_all[:,0:180]
Pyy_id = Pyy_all[:,180:r+1]

# find integrated PSD (2D)
df = f[1]-f[0]
int_Pxx_all = np.sqrt(np.cumsum(Pxx_all, axis = 0)*df)
int_Pyy_all = np.sqrt(np.cumsum(Pyy_all, axis = 0)*df)
int_Pxx_non_disp = int_Pxx_all[:,non_disp]
int_Pxx_disp = int_Pxx_all[:,disp]
int_Pxx_id = int_Pxx_all[:,180:r+1]
int_Pyy = int_Pyy_all[:,0:180]
int_Pyy_id = int_Pyy_all[:,180:r+1]

# PSD mean of all bpm for different types
Pxx_disp_mean = np.mean(Pxx_disp, axis=1)
Pxx_non_disp_mean = np.mean(Pxx_non_disp, axis=1)
Pxx_id_mean = np.mean(Pxx_id, axis=1)
Pyy_mean = np.mean(Pyy, axis=1)
Pyy_id_mean = np.mean(Pyy_id, axis=1)

# find int. PSD mean
int_Pxx_disp_mean = np.sqrt(np.cumsum(Pxx_disp_mean)*df)
int_Pxx_non_disp_mean = np.sqrt(np.cumsum(Pxx_non_disp_mean)*df)
int_Pxx_id_mean = np.sqrt(np.cumsum(Pxx_id_mean)*df)
int_Pyy_mean = np.sqrt(np.cumsum(Pyy_mean)*df)
int_Pyy_id_mean = np.sqrt(np.cumsum(Pyy_id_mean)*df)

# find int. PSD of specific ranges of frequency
# 0 <= f0 < f1 < f2 < f3 <= 5000 for FA data
f0 = 0.1
f1 = 1
f2 = 500
f3 = 5000
i0 = np.where(f>=f0)
i1 = np.where(f>=f1)
i2 = np.where(f>=f2)
i3 = np.where(f>=f3)
i0 = i0[0][0]
i1 = i1[0][0]
i2 = i2[0][0]
i3 = i3[0]
if not i3.size:
    i3 = f.size -1

# diff int. PSD. from f0 to f1
int_Pxx_disp_mean_diff01 = int_Pxx_disp_mean[i1-1] - int_Pxx_disp_mean[i0-1]
int_Pxx_non_disp_mean_diff01 = int_Pxx_non_disp_mean[i1-1] - int_Pxx_non_disp_mean[i0-1]
int_Pxx_id_mean_diff01 = int_Pxx_id_mean[i1-1] - int_Pxx_id_mean[i0-1]
int_Pyy_mean_diff01 = int_Pyy_mean[i1-1] - int_Pyy_mean[i0-1]
int_Pyy_id_mean_diff01 = int_Pyy_id_mean[i1-1] - int_Pyy_id_mean[i0-1]

# diff int. PSD. from f1 to f2
int_Pxx_disp_mean_diff12 = int_Pxx_disp_mean[i2-1] - int_Pxx_disp_mean[i1-1]
int_Pxx_non_disp_mean_diff12 = int_Pxx_non_disp_mean[i2-1] - int_Pxx_non_disp_mean[i1-1]
int_Pxx_id_mean_diff12 = int_Pxx_id_mean[i2-1] - int_Pxx_id_mean[i1-1]
int_Pyy_mean_diff12 = int_Pyy_mean[i2-1] - int_Pyy_mean[i1-1]
int_Pyy_id_mean_diff121 = int_Pyy_id_mean[i2-1] - int_Pyy_id_mean[i1-1]

# diff int. PSD. from f2 to f3
int_Pxx_disp_mean_diff23 = int_Pxx_disp_mean[i3-1] - int_Pxx_disp_mean[i2-1]
int_Pxx_non_disp_mean_diff23 = int_Pxx_non_disp_mean[i3-1] - int_Pxx_non_disp_mean[i2-1]
int_Pxx_id_mean_diff23 = int_Pxx_id_mean[i3-1] - int_Pxx_id_mean[i2-1]
int_Pyy_mean_diff23 = int_Pyy_mean[i3-1] - int_Pyy_mean[i2-1]
int_Pyy_id_mean_diff23 = int_Pyy_id_mean[i3-1] - int_Pyy_id_mean[i2-1]

# int PSD from f0 to f1
Pxx_disp_mean_sum01 = np.sqrt(np.sum(Pxx_disp_mean[i0:i1])*df)
Pxx_non_disp_mean_sum01 = np.sqrt(np.sum(Pxx_non_disp_mean[i0:i1])*df)
Pxx_id_mean_sum01 = np.sqrt(np.sum(Pxx_id_mean[i0:i1])*df)
Pyy_mean_sum01 = np.sqrt(np.sum(Pyy_mean[i0:i1])*df)
Pyy_id_mean_sum01 = np.sqrt(np.sum(Pyy_id_mean[i0:i1])*df)

# int PSD from f1 to f2
Pxx_disp_mean_sum12 = np.sqrt(np.sum(Pxx_disp_mean[i1:i2])*df)
Pxx_non_disp_mean_sum12 = np.sqrt(np.sum(Pxx_non_disp_mean[i1:i2])*df)
Pxx_id_mean_sum12 = np.sqrt(np.sum(Pxx_id_mean[i1:i2])*df)
Pyy_mean_sum12 = np.sqrt(np.sum(Pyy_mean[i1:i2])*df)
Pyy_id_mean_sum12 = np.sqrt(np.sum(Pyy_id_mean[i1:i2])*df)

# int PSD from f2 to f3
Pxx_disp_mean_sum23 = np.sqrt(np.sum(Pxx_disp_mean[i2:i3])*df)
Pxx_non_disp_mean_sum23 = np.sqrt(np.sum(Pxx_non_disp_mean[i2:i3])*df)
Pxx_id_mean_sum23 = np.sqrt(np.sum(Pxx_id_mean[i2:i3])*df)
Pyy_mean_sum23 = np.sqrt(np.sum(Pyy_mean[i2:i3])*df)
Pyy_id_mean_sum23 = np.sqrt(np.sum(Pyy_id_mean[i2:i3])*df)

### test
##int_Pxx_disp_mean_sum01_2 = np.sqrt((int_Pxx_disp_mean[i1-1])**2 - (int_Pxx_disp_mean[i0-1])**2)
##int_Pxx_non_disp_mean_sum01_2 = np.sqrt(int_Pxx_non_disp_mean[i1-1]**2 - int_Pxx_non_disp_mean[i0-1]**2)

# find peaks for average PSD
for i in range(0,5):
    if i == 0:
        P = Pxx_disp_mean
        prom = prom_x
    elif i == 1:
        P = Pxx_non_disp_mean
        prom = prom_x
    elif i == 2:
        P = Pxx_id_mean
        prom = prom_x
    elif i == 3:
        P = Pyy_mean
        prom = prom_y
    elif i == 4:
        P = Pyy_id_mean
        prom = prom_y
        
    loc, _ = find_peaks(P, prominence = prom, distance = dist)
    pks = P[loc]

    # pick max N peaks
    loc_n = (-pks).argsort()
    loc_n = np.sort(loc_n[0:n_max])
    loc_n = loc_n.astype(int)

    # collect peaks data for all horizontal bpm
    if i == 0:
        Pxx_disp_mean_pks_freq = f[loc] # all pks locations
        Pxx_disp_mean_pks_hight = pks # all pks hight
        Pxx_disp_mean_pks_n_freq = f[loc[loc_n]] # n max pks locations
        Pxx_disp_mean_pks_n_hight = pks[loc_n] # n max pks hight
    elif i == 1:
        Pxx_non_disp_mean_pks_freq = f[loc] # all pks locations
        Pxx_non_disp_mean_pks_hight = pks # all pks hight
        Pxx_non_disp_mean_pks_n_freq = f[loc[loc_n]] # n max pks locations
        Pxx_non_disp_mean_pks_n_hight = pks[loc_n] # n max pks hight
    elif i == 2:
        Pxx_id_mean_pks_freq = f[loc] # all pks locations
        Pxx_id_mean_pks_hight = pks # all pks hight
        Pxx_id_mean_pks_n_freq = f[loc[loc_n]] # n max pks locations
        Pxx_id_mean_pks_n_hight = pks[loc_n] # n max pks hight
    elif i == 3:
        Pyy_mean_pks_freq = f[loc] # all pks locations
        Pyy_mean_pks_hight = pks # all pks hight
        Pyy_mean_pks_n_freq = f[loc[loc_n]] # n max pks locations
        Pyy_mean_pks_n_hight = pks[loc_n] # n max pks hight
    elif i == 4:
        Pyy_id_mean_pks_freq = f[loc] # all pks locations
        Pyy_id_mean_pks_hight = pks # all pks hight
        Pyy_id_mean_pks_n_freq = f[loc[loc_n]] # n max pks locations
        Pyy_id_mean_pks_n_hight = pks[loc_n] # n max pks hight


t = time.time() - t0
print "finish calculation at t = %.f" % (t)

#----------------------- plot -------------------------
print "start plotting"

# plot PSD of all non-disp bpm (Horizontal)
plt.figure()
plt.subplot(2,1,1) # PSD
plt.semilogx(f, Pxx_non_disp)
plt.grid(True)
plt.title('All Non-Dispersive BPMs (Horizontal)')
plt.ylabel('PSD [nm^2/Hz]')
plt.subplot(2,1,2) # Int. PSD 
plt.semilogx(f, int_Pxx_non_disp)
plt.grid(True)
plt.xlabel('Frequency [Hz]')
plt.ylabel('Int PSD [nm]')

# plot PSD of all disp bpm (Horizontal)
plt.figure()
plt.subplot(2,1,1) # PSD
plt.semilogx(f, Pxx_disp)
plt.grid(True)
plt.title('All Dispersive BPMs (Horizontal)')
plt.ylabel('PSD [nm^2/Hz]')
plt.subplot(2,1,2) # Int. PSD 
plt.semilogx(f, int_Pxx_disp)
plt.grid(True)
plt.xlabel('Frequency [Hz]')
plt.ylabel('Int PSD [nm]')

# plot PSD of all id bpm (Horizontal)
plt.figure()
plt.subplot(2,1,1) # PSD
plt.semilogx(f, Pxx_id)
plt.grid(True)
plt.title('All ID BPMs (Horizontal)')
plt.ylabel('PSD [nm^2/Hz]')
plt.subplot(2,1,2) # Int. PSD 
plt.semilogx(f, int_Pxx_id)
plt.grid(True)
plt.xlabel('Frequency [Hz]')
plt.ylabel('Int PSD [nm]')

# plot PSD of all normal Vertical bpm
plt.figure()
plt.subplot(2,1,1) # PSD
plt.semilogx(f, Pyy)
plt.grid(True)
plt.title('All Normal BPMs (Vertical)')
plt.ylabel('PSD [nm^2/Hz]')
plt.subplot(2,1,2) # Int. PSD 
plt.semilogx(f, int_Pyy)
plt.grid(True)
plt.xlabel('Frequency [Hz]')
plt.ylabel('Int PSD [nm]')

# plot PSD of id Vertical bpm
plt.figure()
plt.subplot(2,1,1) # PSD
plt.semilogx(f, Pyy_id)
plt.grid(True)
plt.title('All ID BPMs (Vertical)')
plt.ylabel('PSD [nm^2/Hz]')
plt.subplot(2,1,2) # Int. PSD 
plt.semilogx(f, int_Pyy_id)
plt.grid(True)
plt.xlabel('Frequency [Hz]')
plt.ylabel('Int PSD [nm]')

# plot PSD of mean disp bpm (Horizontal)
plt.figure()
plt.subplot(2,1,1) # PSD
plt.semilogx(f, Pxx_disp_mean)
plt.semilogx(Pxx_disp_mean_pks_freq, Pxx_disp_mean_pks_hight, 'x', label = 'all peaks')
plt.semilogx(Pxx_disp_mean_pks_n_freq, Pxx_disp_mean_pks_n_hight, 'o', label = 'max '+str(n_max)+' peaks')
plt.grid(True)
plt.legend(loc='upper left')
plt.title('Average Dispersive BPMs (Horizontal)')
plt.ylabel('PSD [nm^2/Hz]')
plt.subplot(2,1,2) # Int. PSD 
plt.semilogx(f, int_Pxx_disp_mean)
plt.grid(True)
plt.xlabel('Frequency [Hz]')
plt.ylabel('Int PSD [nm]')

# plot PSD of mean non-disp bpm (Horizontal)
plt.figure()
plt.subplot(2,1,1) # PSD
plt.semilogx(f, Pxx_non_disp_mean)
plt.semilogx(Pxx_non_disp_mean_pks_freq, Pxx_non_disp_mean_pks_hight, 'x', label = 'all peaks')
plt.semilogx(Pxx_non_disp_mean_pks_n_freq, Pxx_non_disp_mean_pks_n_hight, 'o', label = 'max '+str(n_max)+' peaks')
plt.grid(True)
plt.legend(loc='upper left')
plt.title('Average Non-Dispersive BPMs (Horizontal)')
plt.ylabel('PSD [nm^2/Hz]')
plt.subplot(2,1,2) # Int. PSD 
plt.semilogx(f, int_Pxx_non_disp_mean)
plt.grid(True)
plt.xlabel('Frequency [Hz]')
plt.ylabel('Int PSD [nm]')

# plot PSD of mean id bpm (Horizontal)
plt.figure()
plt.subplot(2,1,1) # PSD
plt.semilogx(f, Pxx_id_mean)
plt.semilogx(Pxx_id_mean_pks_freq, Pxx_id_mean_pks_hight, 'x', label = 'all peaks')
plt.semilogx(Pxx_id_mean_pks_n_freq, Pxx_id_mean_pks_n_hight, 'o', label = 'max '+str(n_max)+' peaks')
plt.grid(True)
plt.legend(loc='upper left')
plt.title('Average ID BPMs(Horizontal)')
plt.ylabel('PSD [nm^2/Hz]')
plt.subplot(2,1,2) # Int. PSD 
plt.semilogx(f, int_Pxx_id_mean)
plt.grid(True)
plt.xlabel('Frequency [Hz]')
plt.ylabel('Int PSD [nm]')

# plot PSD of mean normal Vertical bpm
plt.figure()
plt.subplot(2,1,1) # PSD
plt.semilogx(f, Pyy_mean)
plt.semilogx(Pyy_mean_pks_freq, Pyy_mean_pks_hight, 'x', label = 'all peaks')
plt.semilogx(Pyy_mean_pks_n_freq, Pyy_mean_pks_n_hight, 'o', label = 'max '+str(n_max)+' peaks')
plt.grid(True)
plt.legend(loc='upper left')
plt.title('Average Normal BPMs (Vertical)')
plt.ylabel('PSD [nm^2/Hz]')
plt.subplot(2,1,2) # Int. PSD 
plt.semilogx(f, int_Pyy_mean)
plt.grid(True)
plt.xlabel('Frequency [Hz]')
plt.ylabel('Int PSD [nm]')

# plot PSD of mean normal Vertical bpm
plt.figure()
plt.subplot(2,1,1) # PSD
plt.semilogx(f, Pyy_id_mean)
plt.semilogx(Pyy_id_mean_pks_freq, Pyy_id_mean_pks_hight, 'x', label = 'all peaks')
plt.semilogx(Pyy_id_mean_pks_n_freq, Pyy_id_mean_pks_n_hight, 'o', label = 'max '+str(n_max)+' peaks')
plt.grid(True)
plt.legend(loc='upper left')
plt.title('Average ID BPMs (Vertical)')
plt.ylabel('PSD [nm^2/Hz]')
plt.subplot(2,1,2) # Int. PSD 
plt.semilogx(f, int_Pyy_id_mean)
plt.grid(True)
plt.xlabel('Frequency [Hz]')
plt.ylabel('Int PSD [nm]')

plt.show()
plt.close()

t = time.time() - t0
print "compleated at time = %.f" % (t)

