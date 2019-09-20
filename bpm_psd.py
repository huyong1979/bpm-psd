'''
This script gets FA data from all Storage Ring BPMs and then process them to get 
all kinds of PSD(Power Spectral Density)-related results:

1) PSD(based on welch) and PSD's peaks for each BPM in X and Y planes; 

2) integral PSD for for each BPM in X and Y planes; 

3) averaged(mean) PSDs of 5 different types of BPMs: horizontal dispersive BPMs 
(only in X plane, no Y plane), horizontal non-dispersive BPMs, horizontal ID BPMs, 
vertical BPMs, and vertical ID BPMs;  

4) integrated PSDs of #3 (based on numpy's cumulative sum 'np.cumsum');

5) summed PSDs of #3 for specific frequency ranges (based on numpy's 'np.sum');   

'''

#08/14/2019,yhu: based on Sukho's scripts: test04_getBPM.py, test03.py
##this is set in st.cmd: PATH="/home/skongtawong/anaconda2/bin:$PATH"

import time
import datetime
t0 = time.time()
print "\n%s: start..."%datetime.datetime.now()
import sys
sys.path.append('/usr/lib/python2.7/dist-packages')
import cothread
from cothread.catools import caget, caput, DBR_CHAR_STR
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import traceback
import h5py
import math
#import operator
#import heapq

#ver = caget('SR:C22-FOFB{CC}FpgaFirmVer-I')
#print ver

#create lists for BPM and PSD PVs
prefix = []
psd_prefix = []
invalid_hPV = [] #yhu: invalid horizontal/vertical BPMs
invalid_vPV = []
#regular BPMs: 6*30
p_index = ['1','2','3','4','5','6']
cell_index =['30','01','02','03','04','05','06','07','08','09','10','11','12','13',
'14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29']
for i in cell_index:
    for j in p_index:
        prefix.append('SR:C'+i+'-BI{BPM:'+j+'}')
        psd_prefix.append('SR:C'+i+'-APHLA{BPM:'+j+'}')
        invalid_hPV.append('SR:C'+i+'-APHLA{BPM:'+j+'}'+'PSD:BadX-Cmd')
        invalid_vPV.append('SR:C'+i+'-APHLA{BPM:'+j+'}'+'PSD:BadY-Cmd')
#ID BPMs
p_index  = ['7','8','9','10']
cell_index=['04','07','12','19']
for i in cell_index:
    for j in p_index:
        prefix.append('SR:C'+i+'-BI{BPM:'+j+'}')
        psd_prefix.append('SR:C'+i+'-APHLA{BPM:'+j+'}')
        invalid_hPV.append('SR:C'+i+'-APHLA{BPM:'+j+'}'+'PSD:BadX-Cmd')
        invalid_vPV.append('SR:C'+i+'-APHLA{BPM:'+j+'}'+'PSD:BadY-Cmd')
p_index=['7','8']
cell_index=['02','03','08','10','11','16','18','21','28']
for i in cell_index:
    for j in p_index:
        prefix.append('SR:C'+i+'-BI{BPM:'+j+'}')
        psd_prefix.append('SR:C'+i+'-APHLA{BPM:'+j+'}')
        invalid_hPV.append('SR:C'+i+'-APHLA{BPM:'+j+'}'+'PSD:BadX-Cmd')
        invalid_vPV.append('SR:C'+i+'-APHLA{BPM:'+j+'}'+'PSD:BadY-Cmd')
p_index=['7','8','9']
cell_index=['05','17','23']
for i in cell_index:
    for j in p_index:
        prefix.append('SR:C'+i+'-BI{BPM:'+j+'}')
        psd_prefix.append('SR:C'+i+'-APHLA{BPM:'+j+'}')
        invalid_hPV.append('SR:C'+i+'-APHLA{BPM:'+j+'}'+'PSD:BadX-Cmd')
        invalid_vPV.append('SR:C'+i+'-APHLA{BPM:'+j+'}'+'PSD:BadY-Cmd')

#a few small functions
def print_status(message):
    print "%s: %s"%(datetime.datetime.now(), str(message))
    #put array of characters as string
    caput('SR-APHLA{BPM}PSD:Status-Wf',str(message),datatype=DBR_CHAR_STR)

def print_status_and_exit(message):
    print_status(message)
    sys.exit(1)

def check_settings(pvList, value):
    values = caget(pvList)
    wrong_settings = []
    for i in range(len(pvList)):
        if values[i] != value:
            wrong_settings.append(pvList[i])
    if len(wrong_settings) != 0:
        message = "No calculation because of wrong settings for %d BPMs: "\
                  %len(wrong_settings) + str(wrong_settings[0])+ " ..."
        print_status_and_exit(message)

# pv lists
FA_X=[] #FA data in X  plane; a list of waveform PVs
FA_Y=[]
PSD_X = [] # PSD in X plane; a list of waveform PVs
PSD_Y = []
PSD_INTX = [] #integral PSD in X plane
PSD_INTY = []
machine_type = [] # BPM settings
trig_src = []
wfm_sel = []
fa_len = []
tx_status = []
for i in range(len(prefix)):
    FA_X.append(prefix[i]+'FA-X')
    FA_Y.append(prefix[i]+'FA-Y') 
    machine_type.append(prefix[i]+'Loc:Machine-SP') #SR:5
    trig_src.append(prefix[i]+'Trig:TrigSrc-SP') #Ext:1
    wfm_sel.append(prefix[i]+'DDR:WfmSel-SP') #FA Wfm:2
    fa_len.append(prefix[i]+'Burst:FaEnableLen-SP') #100000 
    tx_status.append(prefix[i]+'DDR:TxStatus-I') #Wfm Ready:0 
    PSD_X.append(psd_prefix[i]+'PSD:X-Wf') 
    PSD_Y.append(psd_prefix[i]+'PSD:Y-Wf') 
    PSD_INTX.append(psd_prefix[i]+'PSD:IntX-Wf') 
    PSD_INTY.append(psd_prefix[i]+'PSD:IntY-Wf') 

# get BPM FA data from PVs or file
fa_recordLen = caget('SR-APHLA{BPM}PSD:Len-SP')
if caget('SR-APHLA{BPM}PSD:LiveData-Cmd') == 1: # Data Source: Live Data
    if caget('SR:C03-BI{DCCT:1}I:Real-I') < 0.01: # if too low beam current 
        print_status_and_exit("No beam, no calculation")
    if caget('INJ{TOC}OpControl-Sel') == 2 and caget('INJ{TOC-SM}Cnt:Next-I') < 20: # if inject beam in 20 sec
        caput('SR-APHLA{BPM}PSD:Counter-Calc_.PROC', 1)
        print_status_and_exit("Top-off inj. is coming in 20-sec so no calculation")

    check_settings(machine_type, 5)
    #check_settings(trig_src, 1)
    check_settings(wfm_sel, 2)
    #check_settings(fa_len, 100000)
    if caget("SR-APHLA{BPM}PSD:IgnoreMatlabTrigger-Cmd") == 1: # Ignore Matlab trigger
        caput('SR:C21-PS{Pinger}Ping-Cmd', 1) # send event 35 (Pinger) trigger
        print_status("Ignore matlab triggering. Sending Pinger trigger, wait...")
        #time.sleep(11.0); # DO NOT use time.sleep()
        cothread.Sleep(12)
    else: # Use Matlab trigger
        print_status("Use Matlab trigger")
        cothread.Sleep(3)
    for i in range(0,10):
        if sum(caget(tx_status)) > 0:
            cothread.Sleep(1);
        else:
            break
        if i==9:
            print_status_and_exit("Timeout: FA data are not ready")

    print_status("Finally read and process live data...")
    x_all = caget(FA_X, count=fa_recordLen)
    y_all = caget(FA_Y, count=fa_recordLen)
else: # Data Source: Data from file
    path = '/home/skongtawong/Documents/Guimei/20190805_beamstudy_new_bpm_psd_and_gain/'
    filename = 'SR_AllIDBPMs_FA_20190805_0701_06_fofb_on_newPI999992_350Hz_pinger_off_02.h5'
    message = "read data from " + path+filename + " and process the data ..."
    print_status(str(message))
    fid = h5py.File(path+filename, 'r')
    x_all = fid.get('faX')
    y_all = fid.get('faY')
    fid.close

x_all = np.asarray(x_all)[:,0:fa_recordLen]/1000.0 #unit: um
y_all = np.asarray(y_all)[:,0:fa_recordLen]/1000.0
size = x_all.shape
# generate disp, non-disp, id bpm index (start from 0)
disp = []
non_disp = []
for i in range(1,181):
    if np.mod(i,6) == 3 or np.mod(i,6) == 4:
        disp.append(i-1)
    else:
        non_disp.append(i-1)  
id_bpm = range(180,size[0]) 
norm_bpm = range(0,180)

# initialize 
#rf = 499.68e6 # would be good to read from PV
rf = caget('RF{Osc:1}Freq:I')
fs = rf/1320/38
n = x_all.shape
r = n[0] # number of BPMs (rows 223)
c = n[1] # number of samples (columns 100000)
# for finding peaks
n_max = 5 # get 5 peaks for now
# prominence for finding peaks (um^2/Hz)
(prom_x,prom_y)=caget(['SR-APHLA{BPM}PSD:PromX-SP','SR-APHLA{BPM}PSD:PromY-SP']) 
prom_dist = caget('SR-APHLA{BPM}PSD:PromDist-SP') #distance between peaks (Hz)

# pre-allowcate output
Pxx_all = [] #individual BPM PSD in X plane
Pyy_all = []
pks_freq_x = [0]*r #peak's frequency in X plane; no PV for this 
pks_freq_y = [0]*r
pks_hight_x = [0]*r #peak's amplitude/height
pks_hight_y = [0]*r
pks_n_freq_x = [0]*r
pks_n_freq_y = [0]*r
pks_n_hight_x = [0]*r
pks_n_hight_y = [0]*r
#start and end frequencies for finding peaks
f0_pks, f1_pks = caget(['SR-APHLA{BPM}PSD:Freq0-SP','SR-APHLA{BPM}PSD:Freq1-SP'])
n_perseg = c
df = fs/n_perseg #~0.1Hz
i_f0 = int(math.ceil(f0_pks/df)) #index of f0
i_f1 = int(math.ceil(f1_pks/df))
dist = int(math.ceil(prom_dist/df))

#main loop: PSD calculation
def get_PSD_and_peaks(fa_data, prom):
    P_all = []
    pks_freq, pks_hight, pks_n_freq, pks_n_hight = [0]*r, [0]*r, [0]*r, [0]*r
    for i in range(0,r): #for each BPM
        x = fa_data[i]
        #if nperseg = 1024, len(f)=513
        f, P = signal.welch(x, fs, window='hann', nperseg=n_perseg) 
        P_all = np.append(P_all, P) # in 1D
        f2 = f[i_f0:i_f1]
        P2 = P[i_f0:i_f1]
        # find peaks/locations of individual BPM's psd; no PVs for these peaks
        loc, _ = find_peaks(P2[i_f0:i_f1], prominence = prom, distance = dist)
        pks = P2[loc]
        # pick max N peaks
        loc_n = (-pks).argsort()
        loc_n = loc_n[0:n_max]
        #loc_n = np.sort(loc_n[0:n_max]) # sort frequency from minimum
        loc_n = loc_n.astype(int)
        # collect peaks data for all horizontal bpm
        pks_freq[i] = f2[loc] # all pks frequency
        pks_hight[i] = pks # all pks hight
        pks_n_freq[i] = f2[loc[loc_n]] # n max pks frequency
        pks_n_hight[i] = pks[loc_n] # n max pks hight  
    return P_all, pks_freq, pks_hight, pks_n_freq, pks_n_hight, f

Pxx_all, pks_freq_x, pks_hight_x, pks_n_freq_x, pks_n_hight_x, f = get_PSD_and_peaks(x_all, prom_x)
Pyy_all, pks_freq_y, pks_hight_y, pks_n_freq_y, pks_n_hight_y, f = get_PSD_and_peaks(y_all, prom_y)

#caput results
# reshape psd (1D -> 2D)
Pxx_all = np.transpose(Pxx_all.reshape(r,len(f)))
Pyy_all = np.transpose(Pyy_all.reshape(r,len(f)))
Pxx_all_4ca = np.transpose(Pxx_all)
Pyy_all_4ca = np.transpose(Pyy_all)
#print Pxx_all_4ca.shape
#print Pxx_all_4ca
caput(PSD_X, np.sqrt(Pxx_all_4ca[0:r,1:])) #PSD of individual BPM
caput(PSD_Y, np.sqrt(Pyy_all_4ca[0:r,1:]))
caput('SR-APHLA{BPM}Freq-Wf', f[1:]) #skip f[0] which is 0 (DC)
# get invalid BPMs
badx = [] #bad BPM x
bady = [] #[101,125] # start from 0
invalid_x = caget(invalid_hPV)#len(invalid_x) == len(invalid_y): 223
invalid_y = caget(invalid_vPV)
for i in range(len(invalid_x)):
    if invalid_x[i] == 1:
        badx.append(i) 
    if invalid_y[i] == 1:
        bady.append(i) 
print badx
print bady

# function to remove invalid BPM
def remove_BPM(allBPM, badBPM):
    goodBPM = list(set(allBPM)-set(badBPM))
    goodBPM.sort()
    return goodBPM

good_non_disp_x = remove_BPM(non_disp, badx)
good_disp_x = remove_BPM(disp, badx)
good_id_x = remove_BPM(id_bpm, badx)
good_norm_y = remove_BPM(norm_bpm, bady)
good_id_y = remove_BPM(id_bpm, bady)

# pick valid BPMs for disp, non disp, insertion device (id)
Pxx_non_disp = Pxx_all[:,good_non_disp_x] # no PVs for this
Pxx_disp = Pxx_all[:,good_disp_x]
Pxx_id = Pxx_all[:,good_id_x]
Pyy = Pyy_all[:,good_norm_y]
Pyy_id = Pyy_all[:,good_id_y]

# start frequency for integrating PSD
sf = caget("SR-APHLA{BPM}PSD:StartFreq-SP")#start freq for integral PSD
i_sf = np.where(f>=sf)
i_sf = i_sf[0][0]
f_int = f[i_sf:]
# find integrated PSD (2D)
caput('SR-APHLA{BPM}Freq-Wf_', f_int)# frequency for integral PSD
#df = f[1]-f[0] #~0.1Hz
int_Pxx_all = np.sqrt(np.cumsum(Pxx_all[i_sf:], axis = 0)*df)
int_Pyy_all = np.sqrt(np.cumsum(Pyy_all[i_sf:], axis = 0)*df)
int_Pxx_all_4ca = np.transpose(int_Pxx_all)
int_Pyy_all_4ca = np.transpose(int_Pyy_all)
caput(PSD_INTX, int_Pxx_all_4ca[0:r,1:]) #integral PSD of individual BPM
caput(PSD_INTY, int_Pyy_all_4ca[0:r,1:])

int_Pxx_non_disp = np.sqrt(np.cumsum(Pxx_non_disp, axis = 0)*df) # no PV for this
int_Pxx_disp = np.sqrt(np.cumsum(Pxx_disp, axis = 0)*df)
int_Pxx_id = np.sqrt(np.cumsum(Pxx_id, axis = 0)*df)
int_Pyy = np.sqrt(np.cumsum(Pyy, axis = 0)*df)
int_Pyy_id = np.sqrt(np.cumsum(Pyy_id, axis = 0)*df)

# PSD mean of all valid BPMs for different types
Pxx_disp_mean = np.mean(Pxx_disp, axis=1)
Pxx_non_disp_mean = np.mean(Pxx_non_disp, axis=1)
Pxx_id_mean = np.mean(Pxx_id, axis=1)
Pyy_mean = np.mean(Pyy, axis=1)
Pyy_id_mean = np.mean(Pyy_id, axis=1)
pvs = [
'SR-APHLA{BPM}PSD:X_DISP_MEAN-Wf', 'SR-APHLA{BPM}PSD:X_NON_DISP_MEAN-Wf',
'SR-APHLA{BPM}PSD:X_ID_MEAN-Wf',   'SR-APHLA{BPM}PSD:Y_MEAN-Wf',
'SR-APHLA{BPM}PSD:Y_ID_MEAN-Wf']
values = [
np.sqrt(Pxx_disp_mean)[1:], np.sqrt(Pxx_non_disp_mean)[1:], np.sqrt(Pxx_id_mean)[1:],
np.sqrt(Pyy_mean)[1:],      np.sqrt(Pyy_id_mean)[1:] ]
caput(pvs, values)

# find int. PSD mean
int_Pxx_disp_mean = np.sqrt(np.cumsum(Pxx_disp_mean[i_sf:])*df)
int_Pxx_non_disp_mean = np.sqrt(np.cumsum(Pxx_non_disp_mean[i_sf:])*df)
int_Pxx_id_mean = np.sqrt(np.cumsum(Pxx_id_mean[i_sf:])*df)
int_Pyy_mean = np.sqrt(np.cumsum(Pyy_mean[i_sf:])*df)
int_Pyy_id_mean = np.sqrt(np.cumsum(Pyy_id_mean[i_sf:])*df)
pvs = [
'SR-APHLA{BPM}PSD:IntX_DISP_MEAN-Wf', 'SR-APHLA{BPM}PSD:IntX_NON_DISP_MEAN-Wf',
'SR-APHLA{BPM}PSD:IntX_ID_MEAN-Wf',   'SR-APHLA{BPM}PSD:IntY_MEAN-Wf',
'SR-APHLA{BPM}PSD:IntY_ID_MEAN-Wf']
values = [
int_Pxx_disp_mean, int_Pxx_non_disp_mean, int_Pxx_id_mean,
int_Pyy_mean,      int_Pyy_id_mean]
caput(pvs, values)

# find int. PSD of specific ranges of frequency
# 0 <= f0 < f1 < f2 < f3 <= 5000 for FA data
#f0 = 0.1,f1 = 1,f2 = 500,f3 = 5000
f0,f1,f2,f3 = caget(['SR-APHLA{BPM}PSD:F0-SP','SR-APHLA{BPM}PSD:F1-SP',
'SR-APHLA{BPM}PSD:F2-SP','SR-APHLA{BPM}PSD:F3-SP'])
#print f0, f1, f2, f3
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
# int PSD from f0 to f1 (sum before sqrt)
int_Pxx_disp_mean_f0f1 = np.sqrt(np.sum(Pxx_disp_mean[i0:i1])*df)
int_Pxx_non_disp_mean_f0f1 = np.sqrt(np.sum(Pxx_non_disp_mean[i0:i1])*df)
int_Pxx_id_mean_f0f1 = np.sqrt(np.sum(Pxx_id_mean[i0:i1])*df)
int_Pyy_mean_f0f1 = np.sqrt(np.sum(Pyy_mean[i0:i1])*df)
int_Pyy_id_mean_f0f1 = np.sqrt(np.sum(Pyy_id_mean[i0:i1])*df)
# int PSD from f1 to f2
int_Pxx_disp_mean_f1f2 = np.sqrt(np.sum(Pxx_disp_mean[i1:i2])*df)
int_Pxx_non_disp_mean_f1f2 = np.sqrt(np.sum(Pxx_non_disp_mean[i1:i2])*df)
int_Pxx_id_mean_f1f2 = np.sqrt(np.sum(Pxx_id_mean[i1:i2])*df)
int_Pyy_mean_f1f2 = np.sqrt(np.sum(Pyy_mean[i1:i2])*df)
int_Pyy_id_mean_f1f2 = np.sqrt(np.sum(Pyy_id_mean[i1:i2])*df)
# int PSD from f2 to f3
int_Pxx_disp_mean_f2f3 = np.sqrt(np.sum(Pxx_disp_mean[i2:i3])*df)
int_Pxx_non_disp_mean_f2f3 = np.sqrt(np.sum(Pxx_non_disp_mean[i2:i3])*df)
int_Pxx_id_mean_f2f3 = np.sqrt(np.sum(Pxx_id_mean[i2:i3])*df)
int_Pyy_mean_f2f3 = np.sqrt(np.sum(Pyy_mean[i2:i3])*df)
int_Pyy_id_mean_f2f3 = np.sqrt(np.sum(Pyy_id_mean[i2:i3])*df)
pvs = [
'SR-APHLA{BPM}PSD:IntX_DISP_MEAN_F0F1-I',    'SR-APHLA{BPM}PSD:IntX_NON_DISP_MEAN_F0F1-I',
'SR-APHLA{BPM}PSD:IntX_ID_MEAN_F0F1-I',      'SR-APHLA{BPM}PSD:IntY_MEAN_F0F1-I',
'SR-APHLA{BPM}PSD:IntY_ID_MEAN_F0F1-I',      'SR-APHLA{BPM}PSD:IntX_DISP_MEAN_F1F2-I',
'SR-APHLA{BPM}PSD:IntX_NON_DISP_MEAN_F1F2-I','SR-APHLA{BPM}PSD:IntX_ID_MEAN_F1F2-I',
'SR-APHLA{BPM}PSD:IntY_MEAN_F1F2-I',         'SR-APHLA{BPM}PSD:IntY_ID_MEAN_F1F2-I',
'SR-APHLA{BPM}PSD:IntX_DISP_MEAN_F2F3-I',    'SR-APHLA{BPM}PSD:IntX_NON_DISP_MEAN_F2F3-I',
'SR-APHLA{BPM}PSD:IntX_ID_MEAN_F2F3-I',      'SR-APHLA{BPM}PSD:IntY_MEAN_F2F3-I',
'SR-APHLA{BPM}PSD:IntY_ID_MEAN_F2F3-I']
values = [
int_Pxx_disp_mean_f0f1,     int_Pxx_non_disp_mean_f0f1, int_Pxx_id_mean_f0f1,
int_Pyy_mean_f0f1,          int_Pyy_id_mean_f0f1,       int_Pxx_disp_mean_f1f2,
int_Pxx_non_disp_mean_f1f2, int_Pxx_id_mean_f1f2,       int_Pyy_mean_f1f2,
int_Pyy_id_mean_f1f2,       int_Pxx_disp_mean_f2f3,     int_Pxx_non_disp_mean_f2f3,
int_Pxx_id_mean_f2f3,       int_Pyy_mean_f2f3,           int_Pyy_id_mean_f2f3]
caput(pvs, values)

# find peaks for averaged PSD
def get_peaks(f, P, prom):
    #if i_f0 and i_f1 are used, P*_mean (Pxx_disp_mean, etc.) needs to be sliced too 
    f2 = f[i_f0:i_f1]
    P2 = P[i_f0:i_f1]
    #loc, _ = find_peaks(P[i_f0:i_f1], prominence = prom, distance = dist)
    #loc, _ = find_peaks(P, prominence = prom, distance = dist)
    #pks = P[loc] # no PVs for f[loc] and pks

    loc, _ = find_peaks(P2, prominence = prom, distance = dist)
    pks = P2[loc] # no PVs for f[loc] and pks
    # pick max N peaks
    loc_n = (-pks).argsort()
    loc_n = loc_n[0:n_max]
    loc_n = loc_n.astype(int)   
    if len(loc_n) < n_max:
        f_n = [df]*n_max # avoid 0 frequency for logplot
        pks_n = [0]*n_max
    else:
        f_n = f2[loc[loc_n]] # n max pks locations 
        pks_n = pks[loc_n]  # n max pks hight
    return f_n, pks_n 

Pxx_disp_mean_pks_n_freq,     Pxx_disp_mean_pks_n_hight     = get_peaks(f, Pxx_disp_mean, prom_x)
Pxx_non_disp_mean_pks_n_freq, Pxx_non_disp_mean_pks_n_hight = get_peaks(f, Pxx_non_disp_mean, prom_x)
Pxx_id_mean_pks_n_freq,       Pxx_id_mean_pks_n_hight       = get_peaks(f, Pxx_id_mean, prom_x,)
Pyy_mean_pks_n_freq,          Pyy_mean_pks_n_hight          = get_peaks(f, Pyy_mean, prom_y)
Pyy_id_mean_pks_n_freq,       Pyy_id_mean_pks_n_hight       = get_peaks(f, Pyy_id_mean, prom_y)
pvs = [
'SR-APHLA{BPM}PSD:X_DISP_MEAN_PKS_N_F-Wf',     'SR-APHLA{BPM}PSD:X_DISP_MEAN_PKS_N_H-Wf',
'SR-APHLA{BPM}PSD:X_NON_DISP_MEAN_PKS_N_F-Wf', 'SR-APHLA{BPM}PSD:X_NON_DISP_MEAN_PKS_N_H-Wf',
'SR-APHLA{BPM}PSD:X_ID_MEAN_PKS_N_F-Wf',       'SR-APHLA{BPM}PSD:X_ID_MEAN_PKS_N_H-Wf',
'SR-APHLA{BPM}PSD:Y_MEAN_PKS_N_F-Wf',          'SR-APHLA{BPM}PSD:Y_MEAN_PKS_N_H-Wf',
'SR-APHLA{BPM}PSD:Y_ID_MEAN_PKS_N_F-Wf',       'SR-APHLA{BPM}PSD:Y_ID_MEAN_PKS_N_H-Wf']
#print len(Pxx_disp_mean_pks_n_freq)
#print len(Pyy_id_mean_pks_n_hight)
values = [
Pxx_disp_mean_pks_n_freq,     np.sqrt(Pxx_disp_mean_pks_n_hight), 
Pxx_non_disp_mean_pks_n_freq, np.sqrt(Pxx_non_disp_mean_pks_n_hight),
Pxx_id_mean_pks_n_freq,       np.sqrt(Pxx_id_mean_pks_n_hight),
Pyy_mean_pks_n_freq,          np.sqrt(Pyy_mean_pks_n_hight),
Pyy_id_mean_pks_n_freq,       np.sqrt(Pyy_id_mean_pks_n_hight)]
caput(pvs, values)
#print "Avg Horizontal n peaks: "
#print Pxx_non_disp_mean_pks_n_freq
#print np.sqrt(Pxx_non_disp_mean_pks_n_hight)
t = time.time() - t0
message = "Done! t = %.3f seconds" %(t)
print_status(str(message))
caput('SR-APHLA{BPM}PSD:LoopTime-I', t)


# plot if enabled and manually stop IOC and type ./st.cmd
if caget('SR-APHLA{BPM}PSD:Plot-Cmd') == 1:
    print "start plotting"
    
    # plot PSD of all non-disp bpm (Horizontal)
    plt.figure()
    plt.subplot(2,1,1) # PSD
    plt.semilogx(f, np.sqrt(Pxx_non_disp))
    plt.grid(True)
    plt.title('All Non-Dispersive BPMs (Horizontal)')
    plt.ylabel('PSD [nm/sqrt(Hz)]')
    plt.subplot(2,1,2) # Int. PSD 
    plt.semilogx(f, int_Pxx_non_disp)
    plt.grid(True)
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Int PSD [nm]')
    
    ## plot PSD of all disp bpm (Horizontal)
    #plt.figure()
    #plt.subplot(2,1,1) # PSD
    #plt.semilogx(f, np.sqrt(Pxx_disp))
    #plt.grid(True)
    #plt.title('All Dispersive BPMs (Horizontal)')
    #plt.ylabel('PSD [nm/sqrt(Hz)]')
    #plt.subplot(2,1,2) # Int. PSD 
    #plt.semilogx(f, int_Pxx_disp)
    #plt.grid(True)
    #plt.xlabel('Frequency [Hz]')
    #plt.ylabel('Int PSD [nm]')
    
    # plot PSD of all id bpm (Horizontal)
    #plt.figure()
    #plt.subplot(2,1,1) # PSD
    #plt.semilogx(f, np.sqrt(Pxx_id))
    #plt.grid(True)
    #plt.title('All ID BPMs (Horizontal)')
    #plt.ylabel('PSD [nm/sqrt(Hz)]')
    #plt.subplot(2,1,2) # Int. PSD 
    #plt.semilogx(f, int_Pxx_id)
    #plt.grid(True)
    #plt.xlabel('Frequency [Hz]')
    #plt.ylabel('Int PSD [nm]')
    
    # plot PSD of all normal Vertical bpm
    plt.figure()
    plt.subplot(2,1,1) # PSD
    plt.semilogx(f, np.sqrt(Pyy))
    plt.grid(True)
    plt.title('All Normal BPMs (Vertical)')
    plt.ylabel('PSD [nm/sqrt(Hz)]')
    plt.subplot(2,1,2) # Int. PSD 
    plt.semilogx(f, int_Pyy)
    plt.grid(True)
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Int PSD [nm]')
    
    # plot PSD of id Vertical bpm
    #plt.figure()
    #plt.subplot(2,1,1) # PSD
    #plt.semilogx(f, np.sqrt(Pyy_id))
    #plt.grid(True)
    #plt.title('All ID BPMs (Vertical)')
    #plt.ylabel('PSD [nm/sqrt(Hz)]')
    #plt.subplot(2,1,2) # Int. PSD 
    #plt.semilogx(f, int_Pyy_id)
    #plt.grid(True)
    #plt.xlabel('Frequency [Hz]')
    #plt.ylabel('Int PSD [nm]')
    
    ## plot PSD of mean disp bpm (Horizontal)
    #plt.figure()
    #plt.subplot(2,1,1) # PSD
    #plt.semilogx(f, np.sqrt(Pxx_disp_mean))
    #plt.semilogx(Pxx_disp_mean_pks_freq, np.sqrt(Pxx_disp_mean_pks_hight), 'x', label = 'all peaks')
    #plt.semilogx(Pxx_disp_mean_pks_n_freq, np.sqrt(Pxx_disp_mean_pks_n_hight), 'o', label = 'max '+str(n_max)+' peaks')
    #plt.grid(True)
    #plt.legend(loc='upper left')
    #plt.title('Average Dispersive BPMs (Horizontal)')
    #plt.ylabel('PSD [nm/sqrt(Hz)]')
    #plt.subplot(2,1,2) # Int. PSD 
    #plt.semilogx(f, int_Pxx_disp_mean)
    #plt.grid(True)
    #plt.xlabel('Frequency [Hz]')
    #plt.ylabel('Int PSD [nm]')
    
    # plot PSD of mean non-disp bpm (Horizontal)
    plt.figure()
    plt.subplot(2,1,1) # PSD
    plt.semilogx(f, np.sqrt(Pxx_non_disp_mean))
    plt.semilogx(Pxx_non_disp_mean_pks_freq, np.sqrt(Pxx_non_disp_mean_pks_hight), 'x', label = 'all peaks')
    plt.semilogx(Pxx_non_disp_mean_pks_n_freq, np.sqrt(Pxx_non_disp_mean_pks_n_hight), 'o', label = 'max '+str(n_max)+' peaks')
    plt.grid(True)
    plt.legend(loc='upper left')
    plt.title('Average Non-Dispersive BPMs (Horizontal)')
    plt.ylabel('PSD [nm/sqrt(Hz)]')
    plt.subplot(2,1,2) # Int. PSD 
    plt.semilogx(f, int_Pxx_non_disp_mean)
    plt.grid(True)
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Int PSD [nm]')
    
    ## plot PSD of mean id bpm (Horizontal)
    #plt.figure()
    #plt.subplot(2,1,1) # PSD
    #plt.semilogx(f, np.sqrt(Pxx_id_mean))
    #plt.semilogx(Pxx_id_mean_pks_freq, np.sqrt(Pxx_id_mean_pks_hight), 'x', label = 'all peaks')
    #plt.semilogx(Pxx_id_mean_pks_n_freq, np.sqrt(Pxx_id_mean_pks_n_hight), 'o', label = 'max '+str(n_max)+' peaks')
    #plt.grid(True)
    #plt.legend(loc='upper left')
    #plt.title('Average ID BPMs(Horizontal)')
    #plt.ylabel('PSD [nm^2/Hz]')
    #plt.subplot(2,1,2) # Int. PSD 
    #plt.semilogx(f, int_Pxx_id_mean)
    #plt.grid(True)
    #plt.xlabel('Frequency [Hz]')
    #plt.ylabel('Int PSD [nm]')
    
    # plot PSD of mean normal Vertical bpm
    plt.figure()
    plt.subplot(2,1,1) # PSD
    plt.semilogx(f, np.sqrt(Pyy_mean))
    plt.semilogx(Pyy_mean_pks_freq, np.sqrt(Pyy_mean_pks_hight), 'x', label = 'all peaks')
    plt.semilogx(Pyy_mean_pks_n_freq, np.sqrt(Pyy_mean_pks_n_hight), 'o', label = 'max '+str(n_max)+' peaks')
    plt.grid(True)
    plt.legend(loc='upper left')
    plt.title('Average Normal BPMs (Vertical)')
    plt.ylabel('PSD [nm/sqrt(Hz)]')
    plt.subplot(2,1,2) # Int. PSD 
    plt.semilogx(f, int_Pyy_mean)
    plt.grid(True)
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Int PSD [nm]')
    
    ## plot PSD of mean id Vertical bpm
    #plt.figure()
    #plt.subplot(2,1,1) # PSD
    #plt.semilogx(f, np.sqrt(Pyy_id_mean))
    #plt.semilogx(Pyy_id_mean_pks_freq, np.sqrt(Pyy_id_mean_pks_hight), 'x', label = 'all peaks')
    #plt.semilogx(Pyy_id_mean_pks_n_freq, np.sqrt(Pyy_id_mean_pks_n_hight), 'o', label = 'max '+str(n_max)+' peaks')
    #plt.grid(True)
    #plt.legend(loc='upper left')
    #plt.title('Average ID BPMs (Vertical)')
    #plt.ylabel('PSD [nm/sqrt(Hz)]')
    #plt.subplot(2,1,2) # Int. PSD 
    #plt.semilogx(f, int_Pyy_id_mean)
    #plt.grid(True)
    #plt.xlabel('Frequency [Hz]')
    #plt.ylabel('Int PSD [nm]')
    
    plt.show()
    plt.close()
    
    #t = time.time() - t0
    #print "completed at time = %.f" % (t)

