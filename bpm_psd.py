'''
This script gets FA data from all Storage Ring BPMs and then process them to get 
all kinds of PSD (Power Spectral Density)-related results:

1) PSD (based on welch) and PSD's peaks for each BPM in X and Y planes; 

2) integral PSD for for each BPM in X and Y planes; 

3) averaged (mean) PSDs of 5 different types of BPM data; 
   in X plane: dispersive BPMs(60-bad), non-dispersive(120-bad), ID BPMs(43-bad);
   in Y plane: non-ID/regular BPMs(180-bad), ID BPMs;  

4) integrated PSDs of #3 (based on numpy's cumulative sum 'np.cumsum');

5) summed PSDs of #3 for specific frequency ranges (based on numpy's 'np.sum');   

'''

#11/18/2019: use conda enviroment; epicsEnvSet("EPICS_BASE"...) is required.
  #epicsEnvSet("PATH", "/opt/conda_envs/ap-2019-2.0/bin:$PATH")
  #epicsEnvSet("EPICS_BASE", "/usr/lib/epics")

#08/14/2019,yhu: based on Sukho's scripts: test04_getBPM.py, test03.py
##this is set in st.cmd: PATH="/home/skongtawong/anaconda2/bin:$PATH"

import os
import time
import datetime
t0 = time.time()
print("\n%s: start a new cycle..."%datetime.datetime.now())
import sys
#sys.path.append('/usr/lib/python2.7/dist-packages') # for cothread
import cothread
from cothread.catools import caget, caput, DBR_CHAR_STR
import numpy as np
from scipy import signal
from scipy.signal import find_peaks
import h5py
import math
import matplotlib.pyplot as plt
#import operator
#import heapq
from noise_locator import noise_locator, locate_n_peaks, plot_mesh, plot_n_pks
from save_data import save_data

#ver = caget('SR:C22-FOFB{CC}FpgaFirmVer-I')
#print(ver)

#create BPM and PSD PV lists
bpm_index  = [ ['1','2','3','4','5','6'], #group1: 6 regular BPMs: 6*30
               ['7','8','9','10'], #group2: 4 ID BPMs    
               ['7','8'], #group3: 2 ID BPMs    
               ['7','8','9'] #group4: 3 ID BPMs      
             ]
cell_index = [ ['30','01','02','03','04','05','06','07','08','09','10','11','12','13',
'14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29'],
               ['04','07','12','19'], #these cells have 4 ID BPMs
               ['02','03','08','10','11','16','18','21','28'], # 2 ID BPMs/cell
               ['05','17','23'] # 3 ID BPMs per cell
             ]
#List Comprehension works faster?
def create_pvlist(str1, str2):
    return (['SR:C'+i+'-'+str(str1)+'{BPM:'+j+'}'+str(str2) 
  for n in range(len(bpm_index)) for i in cell_index[n] for j in bpm_index[n]])
#fa_x: FA data in X  plane; a list of waveform PVs
fa_x =         create_pvlist('BI', 'FA-X')
fa_y =         create_pvlist('BI', 'FA-Y')     
fa_a =         create_pvlist('BI', 'FA-A') 
fa_s =         create_pvlist('BI', 'FA-S')   
prefix =       create_pvlist('BI', '')
machine_type = create_pvlist('BI', 'Loc:Machine-SP') #should be "5"
trig_src =     create_pvlist('BI', 'Trig:TrigSrc-SP')     
wfm_sel =      create_pvlist('BI', 'DDR:WfmSel-SP') #should be "2"
fa_len =       create_pvlist('BI', 'Burst:FaEnableLen-SP')      
tx_status =    create_pvlist('BI', 'DDR:TxStatus-I') #"0" means ready
#psd_x: PSD in X plane; a list of waveform PVs
psd_x =    create_pvlist('APHLA', 'PSD:X-Wf')   
psd_y =    create_pvlist('APHLA', 'PSD:Y-Wf')    
psd_intx = create_pvlist('APHLA', 'PSD:IntX-Wf') #integral PSD in X plane 
psd_inty = create_pvlist('APHLA', 'PSD:IntY-Wf') 
pv_x =     create_pvlist('APHLA', 'PSD:BadX-Cmd') 
pv_y =     create_pvlist('APHLA', 'PSD:BadY-Cmd') 

#a few small functions
def update_status(message):
    print("%s: %s"%(datetime.datetime.now(), str(message)))
    #put array of characters as string
    caput('SR-APHLA{BPM}PSD:Status-Wf',str(message),datatype=DBR_CHAR_STR)

def update_status_and_exit(message):
    update_status(message)
    sys.exit(1)

def check_settings(pvs, value):
    values = caget(pvs)
    wrong_settings = []
    for i in range(len(pvs)):
        if values[i] != value:
            wrong_settings.append(pvs[i])
    if len(wrong_settings) != 0:
        message = "No calculation because of wrong settings for %d BPMs: "\
                  %len(wrong_settings) + str(wrong_settings[0])+ " ..."
        update_status_and_exit(message)

# get BPM FA data from PVs or a file
fa_recordLen = caget('SR-APHLA{BPM}PSD:Len-SP')
if caget('SR-APHLA{BPM}PSD:LiveData-Cmd') == 1: # Data Source: Live Data
    if caget('SR:C03-BI{DCCT:1}I:Real-I') < 0.01: # if too low beam current 
        update_status_and_exit("No beam, no calculation, waiting for a new cycle...")
    # if beam is injected in 20 sec
    if caget('INJ{TOC}OpControl-Sel') == 2 and caget('INJ{TOC-SM}Cnt:Next-I') < 20: 
        caput('SR-APHLA{BPM}PSD:Counter-Calc_.PROC', 1)
        update_status_and_exit("Top-off inj. is coming in 20-sec so no calculation"+
", waiting for a new cycle")
    check_settings(machine_type, 5)
    #check_settings(trig_src, 1)
    check_settings(wfm_sel, 2)
    #check_settings(fa_len, 100000)
    if caget("SR-APHLA{BPM}PSD:IgnoreMatlabTrigger-Cmd") == 1: # Ignore Matlab trigger
        caput('SR:C21-PS{Pinger}Ping-Cmd', 1) # send event 35 (Pinger) trigger
        update_status("Ignore matlab triggering. Sending Pinger trigger, wait...")
        #time.sleep(11.0); # DO NOT use time.sleep()
        cothread.Sleep(12)
    else: # Use Matlab trigger
        update_status("Use Matlab trigger ...")
        cothread.Sleep(3)
    for i in range(0,10):
        if sum(caget(tx_status)) > 0:
            cothread.Sleep(1)
        else:
            break
        if i==9:
            update_status_and_exit("Timeout: FA data are not ready")

    update_status("Finally read and process live data...")
    x_all = caget(fa_x, count=fa_recordLen)
    y_all = caget(fa_y, count=fa_recordLen)
    a_all = caget(fa_a, count=fa_recordLen)
    s_all = caget(fa_s, count=fa_recordLen)
    fa_xyas = [x_all, y_all, a_all, s_all]
else: # Data Source: Data from file
    #path = '/home/skongtawong/Documents/Guimei/FAData/20191010_fofb_onoff_changepump_p1/' 
    #file_name = path + 'SR_AllIDBPMs_FA_20191010_0311_13_on01.h5'
    file_name = caget('SR-APHLA{BPM}PSD:File4DataSource-SP', datatype=DBR_CHAR_STR)
    update_status("read data from " + file_name + " and process the data ...")
    '''with...as does not work here because x_all will be '<Closed HDF5 dataset>'
    try:
        with h5py.File(file_name, 'r') as fid:
            x_all = fid.get('faX')
            y_all = fid.get('faY')
            print(type(x_all))
    except Exception as e:
        print("{}: {}".format(type(e),e))
        update_status_and_exit("Error: can not read " + file_name)   
    '''  
    if not os.path.isfile(file_name):
        update_status_and_exit("Error: can not read " + file_name)
    fid = h5py.File(file_name, 'r')
    x_all = fid.get('faX')
    y_all = fid.get('faY')
    fid.close
    
x_all = np.asarray(x_all)[:,0:fa_recordLen]/1000.0 #unit: from nm to um
y_all = np.asarray(y_all)[:,0:fa_recordLen]/1000.0

# generate disp, non-disp, id bpm index (start from 0)
size = x_all.shape
disp = []
non_disp = []
for i in range(1,181):
    if np.mod(i,6) == 3 or np.mod(i,6) == 4:
        disp.append(i-1)
    else:
        non_disp.append(i-1)  
id_bpm = range(180,size[0]) 
norm_bpm = range(0,180)
#print(len(disp), len(non_disp), len(id_bpm)) #60, 120, 43

# initialize 
#rf = 499.68e6 # would be good to read from PV
rf = caget('RF{Osc:1}Freq:I')
fs = rf/1320/38
n = x_all.shape
r = n[0] # number of BPMs (rows 223)
c = n[1] # number of samples (columns 100000)
n_perseg = c #samples per segment in signal.welch
df = fs/n_perseg #resolution: ~0.1Hz
#start and end frequencies for finding peaks
f0_pks, f1_pks = caget(['SR-APHLA{BPM}PSD:Freq0-SP','SR-APHLA{BPM}PSD:Freq1-SP'])
i_f0 = int(math.ceil(f0_pks/df)) #index of f0
i_f1 = int(math.ceil(f1_pks/df))
# prominence for finding peaks (um^2/Hz)
(prom_x,prom_y)=caget(['SR-APHLA{BPM}PSD:PromX-SP','SR-APHLA{BPM}PSD:PromY-SP']) 
prom_dist = caget('SR-APHLA{BPM}PSD:PromDist-SP') #distance between peaks (Hz)
dist = int(math.ceil(prom_dist/df))
n_max = 5 # get 5 peaks for now

# pre-allowcate output
Pxx_all = [] #individual BPM PSD in X plane
Pyy_all = []
pks_freq_x = [0]*r #peak's frequency in X plane; no PV for individual PSD peak 
pks_freq_y = [0]*r
pks_hight_x = [0]*r #peak's amplitude/height
pks_hight_y = [0]*r
pks_n_freq_x = [0]*r
pks_n_freq_y = [0]*r
pks_n_hight_x = [0]*r
pks_n_hight_y = [0]*r

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

#caput results (only Pxx_all, Pyy_all)
# reshape psd (1D -> 2D)
Pxx_all = np.transpose(Pxx_all.reshape(r,len(f)))
Pyy_all = np.transpose(Pyy_all.reshape(r,len(f)))
Pxx_all_4ca = np.transpose(Pxx_all)
Pyy_all_4ca = np.transpose(Pyy_all)
#print(Pxx_all_4ca.shape)
caput(psd_x, np.sqrt(Pxx_all_4ca[0:r,1:])) #PSD of individual BPM
caput(psd_y, np.sqrt(Pyy_all_4ca[0:r,1:]))
caput('SR-APHLA{BPM}Freq-Wf', f[1:]) #skip f[0] which is 0 (DC)

# get invalid BPM index
value_x = caget(pv_x)#len(pv_x) == len(pv_y): 223
value_y = caget(pv_y)
badx = [i for i in range(len(value_x)) if value_x[i] == 1] 
bady = [i for i in range(len(value_y)) if value_y[i] == 1] 
badx_pvname = [pv_x[i] for i in badx]
bady_pvname = [pv_y[i] for i in bady]
#print(badx_pvname)
bad_xy = [badx_pvname, bady_pvname]

#by skongtaw: function to remove invalid BPM
def remove_BPM(allBPM, badBPM): # can be all badBPM
    goodBPM = list(set(allBPM)-set(badBPM))
    goodBPM.sort()
    return goodBPM

good_non_disp_x = remove_BPM(non_disp, badx)
good_disp_x = remove_BPM(disp, badx)
good_id_x = remove_BPM(id_bpm, badx)
good_norm_y = remove_BPM(norm_bpm, bady)
good_id_y = remove_BPM(id_bpm, bady)
pvs = ['SR-APHLA{BPM}PSD:X_DISP_GoodBPM-I', 'SR-APHLA{BPM}PSD:X_NON_DISP_GoodBPM-I',
       'SR-APHLA{BPM}PSD:X_ID_GoodBPM-I',   'SR-APHLA{BPM}PSD:Y_GoodBPM-I',
       'SR-APHLA{BPM}PSD:Y_ID_GoodBPM-I']
values = [len(good_disp_x), len(good_non_disp_x), len(good_id_x), 
          len(good_norm_y), len(good_id_y)]
caput(pvs, values)

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
caput(psd_intx, int_Pxx_all_4ca[0:r,1:]) #integral PSD of individual BPM
caput(psd_inty, int_Pyy_all_4ca[0:r,1:])

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
pvs = ['SR-APHLA{BPM}PSD:X_DISP_MEAN-Wf', 'SR-APHLA{BPM}PSD:X_NON_DISP_MEAN-Wf',
       'SR-APHLA{BPM}PSD:X_ID_MEAN-Wf',   'SR-APHLA{BPM}PSD:Y_MEAN-Wf',
       'SR-APHLA{BPM}PSD:Y_ID_MEAN-Wf']
mean_PSDs=[np.sqrt(Pxx_disp_mean)[1:],np.sqrt(Pxx_non_disp_mean)[1:],np.sqrt(Pxx_id_mean)[1:],
        np.sqrt(Pyy_mean)[1:],     np.sqrt(Pyy_id_mean)[1:] ]
caput(pvs, mean_PSDs)

# find int. PSD mean
int_Pxx_disp_mean = np.sqrt(np.cumsum(Pxx_disp_mean[i_sf:])*df)
int_Pxx_non_disp_mean = np.sqrt(np.cumsum(Pxx_non_disp_mean[i_sf:])*df)
int_Pxx_id_mean = np.sqrt(np.cumsum(Pxx_id_mean[i_sf:])*df)
int_Pyy_mean = np.sqrt(np.cumsum(Pyy_mean[i_sf:])*df)
int_Pyy_id_mean = np.sqrt(np.cumsum(Pyy_id_mean[i_sf:])*df)
pvs = ['SR-APHLA{BPM}PSD:IntX_DISP_MEAN-Wf', 'SR-APHLA{BPM}PSD:IntX_NON_DISP_MEAN-Wf',
       'SR-APHLA{BPM}PSD:IntX_ID_MEAN-Wf',   'SR-APHLA{BPM}PSD:IntY_MEAN-Wf',
       'SR-APHLA{BPM}PSD:IntY_ID_MEAN-Wf']
int_mean_PSDs = [int_Pxx_disp_mean, int_Pxx_non_disp_mean, int_Pxx_id_mean,
          int_Pyy_mean,      int_Pyy_id_mean]
caput(pvs, int_mean_PSDs)

# find int. PSD of specific ranges of frequency
# 0 <= f0 < f1 < f2 < f3 <= 5000 for FA data
#f0 = 0.1,f1 = 1,f2 = 500,f3 = 5000
f0,f1,f2,f3 = caget(['SR-APHLA{BPM}PSD:F0-SP','SR-APHLA{BPM}PSD:F1-SP',
                     'SR-APHLA{BPM}PSD:F2-SP','SR-APHLA{BPM}PSD:F3-SP'])
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

# find integral PSD from f_i to f_f (sum before sqrt)
def get_intPSD_fi_ff(i_i, i_f):
    return (np.sqrt(np.sum(Pxx_disp_mean[i_i:i_f])*df), 
            np.sqrt(np.sum(Pxx_non_disp_mean[i_i:i_f])*df),
            np.sqrt(np.sum(Pxx_id_mean[i_i:i_f])*df),
            np.sqrt(np.sum(Pyy_mean[i_i:i_f])*df),
            np.sqrt(np.sum(Pyy_id_mean[i_i:i_f])*df))

(int_Pxx_disp_mean_f0f1, int_Pxx_non_disp_mean_f0f1, int_Pxx_id_mean_f0f1, 
int_Pyy_mean_f0f1, int_Pyy_id_mean_f0f1) = get_intPSD_fi_ff(i0, i1)
(int_Pxx_disp_mean_f0f2, int_Pxx_non_disp_mean_f0f2, int_Pxx_id_mean_f0f2,
int_Pyy_mean_f0f2, int_Pyy_id_mean_f0f2) = get_intPSD_fi_ff(i0, i2)
(int_Pxx_disp_mean_f0f3, int_Pxx_non_disp_mean_f0f3, int_Pxx_id_mean_f0f3,
int_Pyy_mean_f0f3, int_Pyy_id_mean_f0f3) = get_intPSD_fi_ff(i0, i3)
(int_Pxx_disp_mean_f1f2, int_Pxx_non_disp_mean_f1f2, int_Pxx_id_mean_f1f2,
int_Pyy_mean_f1f2, int_Pyy_id_mean_f1f2) = get_intPSD_fi_ff(i1, i2)
(int_Pxx_disp_mean_f1f3, int_Pxx_non_disp_mean_f1f3, int_Pxx_id_mean_f1f3,
int_Pyy_mean_f1f3, int_Pyy_id_mean_f1f3) = get_intPSD_fi_ff(i1, i3)
(int_Pxx_disp_mean_f2f3, int_Pxx_non_disp_mean_f2f3, int_Pxx_id_mean_f2f3,
int_Pyy_mean_f2f3, int_Pyy_id_mean_f2f3) = get_intPSD_fi_ff(i2, i3)

pvs = [ # need to add PVs for _f0f2, _f0f3, _f1f3
'SR-APHLA{BPM}PSD:IntX_DISP_MEAN_F0F1-I', 'SR-APHLA{BPM}PSD:IntX_NON_DISP_MEAN_F0F1-I',
'SR-APHLA{BPM}PSD:IntX_ID_MEAN_F0F1-I',   'SR-APHLA{BPM}PSD:IntY_MEAN_F0F1-I',
'SR-APHLA{BPM}PSD:IntY_ID_MEAN_F0F1-I',  
'SR-APHLA{BPM}PSD:IntX_DISP_MEAN_F0F2-I', 'SR-APHLA{BPM}PSD:IntX_NON_DISP_MEAN_F0F2-I',
'SR-APHLA{BPM}PSD:IntX_ID_MEAN_F0F2-I',   'SR-APHLA{BPM}PSD:IntY_MEAN_F0F2-I',
'SR-APHLA{BPM}PSD:IntY_ID_MEAN_F0F2-I', 
'SR-APHLA{BPM}PSD:IntX_DISP_MEAN_F0F3-I', 'SR-APHLA{BPM}PSD:IntX_NON_DISP_MEAN_F0F3-I',
'SR-APHLA{BPM}PSD:IntX_ID_MEAN_F0F3-I',   'SR-APHLA{BPM}PSD:IntY_MEAN_F0F3-I',
'SR-APHLA{BPM}PSD:IntY_ID_MEAN_F0F3-I', 
'SR-APHLA{BPM}PSD:IntX_DISP_MEAN_F1F2-I', 'SR-APHLA{BPM}PSD:IntX_NON_DISP_MEAN_F1F2-I',
'SR-APHLA{BPM}PSD:IntX_ID_MEAN_F1F2-I',   'SR-APHLA{BPM}PSD:IntY_MEAN_F1F2-I',         
'SR-APHLA{BPM}PSD:IntY_ID_MEAN_F1F2-I',
'SR-APHLA{BPM}PSD:IntX_DISP_MEAN_F1F3-I', 'SR-APHLA{BPM}PSD:IntX_NON_DISP_MEAN_F1F3-I',
'SR-APHLA{BPM}PSD:IntX_ID_MEAN_F1F3-I',   'SR-APHLA{BPM}PSD:IntY_MEAN_F1F3-I',         
'SR-APHLA{BPM}PSD:IntY_ID_MEAN_F1F3-I',
'SR-APHLA{BPM}PSD:IntX_DISP_MEAN_F2F3-I', 'SR-APHLA{BPM}PSD:IntX_NON_DISP_MEAN_F2F3-I',
'SR-APHLA{BPM}PSD:IntX_ID_MEAN_F2F3-I',   'SR-APHLA{BPM}PSD:IntY_MEAN_F2F3-I',
'SR-APHLA{BPM}PSD:IntY_ID_MEAN_F2F3-I']

values = [
int_Pxx_disp_mean_f0f1, int_Pxx_non_disp_mean_f0f1, int_Pxx_id_mean_f0f1,
int_Pyy_mean_f0f1,      int_Pyy_id_mean_f0f1,
int_Pxx_disp_mean_f0f2, int_Pxx_non_disp_mean_f0f2, int_Pxx_id_mean_f0f2,
int_Pyy_mean_f0f2,      int_Pyy_id_mean_f0f2, 
int_Pxx_disp_mean_f0f3, int_Pxx_non_disp_mean_f0f3, int_Pxx_id_mean_f0f3,
int_Pyy_mean_f0f3,      int_Pyy_id_mean_f0f3,     
int_Pxx_disp_mean_f1f2, int_Pxx_non_disp_mean_f1f2, int_Pxx_id_mean_f1f2,       
int_Pyy_mean_f1f2,      int_Pyy_id_mean_f1f2,    
int_Pxx_disp_mean_f1f3, int_Pxx_non_disp_mean_f1f3, int_Pxx_id_mean_f1f3,       
int_Pyy_mean_f1f3,      int_Pyy_id_mean_f1f3,    
int_Pxx_disp_mean_f2f3, int_Pxx_non_disp_mean_f2f3,  int_Pxx_id_mean_f2f3,       
int_Pyy_mean_f2f3,      int_Pyy_id_mean_f2f3]
caput(pvs, values)

# find (5) peaks for averaged PSD
[i_f0_x_disp, i_f1_x_disp, i_f0_x_non_disp, i_f1_x_non_disp, i_f0_x_id, i_f1_x_id, 
i_f0_y, i_f1_y, i_f0_y_id, i_f1_y_id] = [int(math.ceil(val/df)) for val in 
  caget(['SR-APHLA{BPM}PSD:X_DISP_MEAN_Freq0-SP',    'SR-APHLA{BPM}PSD:X_DISP_MEAN_Freq1-SP',
         'SR-APHLA{BPM}PSD:X_NON_DISP_MEAN_Freq0-SP','SR-APHLA{BPM}PSD:X_NON_DISP_MEAN_Freq1-SP',
         'SR-APHLA{BPM}PSD:X_ID_MEAN_Freq0-SP',      'SR-APHLA{BPM}PSD:X_ID_MEAN_Freq1-SP',
         'SR-APHLA{BPM}PSD:Y_MEAN_Freq0-SP',         'SR-APHLA{BPM}PSD:Y_MEAN_Freq1-SP',
         'SR-APHLA{BPM}PSD:Y_ID_MEAN_Freq0-SP',      'SR-APHLA{BPM}PSD:Y_ID_MEAN_Freq1-SP'])]

[prom_x_disp, prom_x_non_disp, prom_x_id, prom_y, prom_y_id] = caget(
['SR-APHLA{BPM}PSD:X_DISP_MEAN_Prom-SP', 'SR-APHLA{BPM}PSD:X_NON_DISP_MEAN_Prom-SP',
 'SR-APHLA{BPM}PSD:X_ID_MEAN_Prom-SP',   'SR-APHLA{BPM}PSD:Y_MEAN_Prom-SP',
 'SR-APHLA{BPM}PSD:Y_ID_MEAN_Prom-SP'])

[dist_x_disp, dist_x_non_disp, dist_x_id, dist_y, dist_y_id] = [int(math.ceil(val/df)) for val in 
  caget(['SR-APHLA{BPM}PSD:X_DISP_MEAN_Dist-SP', 'SR-APHLA{BPM}PSD:X_NON_DISP_MEAN_Dist-SP',
       'SR-APHLA{BPM}PSD:X_ID_MEAN_Dist-SP',   'SR-APHLA{BPM}PSD:Y_MEAN_Dist-SP',
       'SR-APHLA{BPM}PSD:Y_ID_MEAN_Dist-SP'])]

def get_peaks(f, P, prom, dist, _i_f0, _i_f1):
    #if _i_f0 and _i_f1 are used, P*_mean (Pxx_disp_mean, etc.) needs to be sliced too 
    f2 = f[_i_f0:_i_f1]
    P2 = P[_i_f0:_i_f1]
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

Pxx_disp_mean_pks_n_freq,     Pxx_disp_mean_pks_n_hight     = get_peaks(f, 
  Pxx_disp_mean,     prom_x_disp,     dist_x_disp,     i_f0_x_disp,     i_f1_x_disp)
Pxx_non_disp_mean_pks_n_freq, Pxx_non_disp_mean_pks_n_hight = get_peaks(f, 
  Pxx_non_disp_mean, prom_x_non_disp, dist_x_non_disp, i_f0_x_non_disp, i_f1_x_non_disp)
Pxx_id_mean_pks_n_freq,       Pxx_id_mean_pks_n_hight       = get_peaks(f, 
  Pxx_id_mean,       prom_x_id,       dist_x_id,       i_f0_x_id,       i_f1_x_id)
Pyy_mean_pks_n_freq,          Pyy_mean_pks_n_hight          = get_peaks(f, 
  Pyy_mean,          prom_y,          dist_y,          i_f0_y,          i_f1_y)
Pyy_id_mean_pks_n_freq,       Pyy_id_mean_pks_n_hight       = get_peaks(f, 
  Pyy_id_mean,       prom_y_id,       dist_y_id,       i_f0_y_id,       i_f1_y_id)
pvs = [
'SR-APHLA{BPM}PSD:X_DISP_MEAN_PKS_N_F-Wf',     'SR-APHLA{BPM}PSD:X_DISP_MEAN_PKS_N_H-Wf',
'SR-APHLA{BPM}PSD:X_NON_DISP_MEAN_PKS_N_F-Wf', 'SR-APHLA{BPM}PSD:X_NON_DISP_MEAN_PKS_N_H-Wf',
'SR-APHLA{BPM}PSD:X_ID_MEAN_PKS_N_F-Wf',       'SR-APHLA{BPM}PSD:X_ID_MEAN_PKS_N_H-Wf',
'SR-APHLA{BPM}PSD:Y_MEAN_PKS_N_F-Wf',          'SR-APHLA{BPM}PSD:Y_MEAN_PKS_N_H-Wf',
'SR-APHLA{BPM}PSD:Y_ID_MEAN_PKS_N_F-Wf',       'SR-APHLA{BPM}PSD:Y_ID_MEAN_PKS_N_H-Wf']
values = [
Pxx_disp_mean_pks_n_freq,     np.sqrt(Pxx_disp_mean_pks_n_hight), 
Pxx_non_disp_mean_pks_n_freq, np.sqrt(Pxx_non_disp_mean_pks_n_hight),
Pxx_id_mean_pks_n_freq,       np.sqrt(Pxx_id_mean_pks_n_hight),
Pyy_mean_pks_n_freq,          np.sqrt(Pyy_mean_pks_n_hight),
Pyy_id_mean_pks_n_freq,       np.sqrt(Pyy_id_mean_pks_n_hight)]
caput(pvs, values)

mean_peaks_f = [Pxx_disp_mean_pks_n_freq, Pxx_non_disp_mean_pks_n_freq, 
Pxx_id_mean_pks_n_freq, Pyy_mean_pks_n_freq, Pyy_id_mean_pks_n_freq]


# noise locator
# find corrector strenght of all frequencies
update_status("Noise locator...")
corr_all_x, corr_all_y = noise_locator(x_all, y_all, good_non_disp_x, good_norm_y) 
print(corr_all_x.shape)#(90, 100000)
print(corr_all_x[0,:])

# locate noise for disp hor, non-disp hor, id hor, vert, id vert (5 peaks)
corr_disp_mean_pks_x,     f_out_disp_mean_pks_x     = locate_n_peaks(corr_all_x, f, Pxx_disp_mean_pks_n_freq) 
corr_non_disp_mean_pks_x, f_out_non_disp_mean_pks_x = locate_n_peaks(corr_all_x, f, Pxx_non_disp_mean_pks_n_freq) 
corr_id_mean_pks__x,      f_out_id_mean_pks__x      = locate_n_peaks(corr_all_x, f, Pxx_id_mean_pks_n_freq) 
corr_mean_pks_y,          f_out_mean_pks_y          = locate_n_peaks(corr_all_y, f, Pyy_mean_pks_n_freq) 
corr_id_mean_pks_y,       f_out_id_mean_pks_y       = locate_n_peaks(corr_all_y, f, Pyy_id_mean_pks_n_freq) 

pvs = [] # a list of a list of waveforms (a list of 2D array)
for t in ['X_DISP','X_NON_DISP','X_ID','Y','Y_ID']: # five types
    pvs.append(['SR-APHLA{CORR}Noise:'+t+'_Freq'+str(n)+'-Wf'for n in [0,1,2,3,4]])
#values: a list of 2D array (5*(5*90))
values = [np.transpose(corr_disp_mean_pks_x), np.transpose(corr_non_disp_mean_pks_x), 
          np.transpose(corr_id_mean_pks__x),  np.transpose(corr_mean_pks_y), 
          np.transpose(corr_id_mean_pks_y)]
#caput(pvs, values) does not work because caput only works on 1D/2D array   
for (pv, value) in zip(pvs, values):
    caput(pv, value)

# can put any frequencies that we are interested to find the locations
f_arb = np.array([52.2, 53.2, 54.5, 59.9, 38, 273, 58]) 
corr_arb_x,               f_out_arb_x               = locate_n_peaks(corr_all_x, f, f_arb)
# test output
#print(corr_non_disp_mean_pks_x.shape) #shape: (90, 5)
#print(f_out_non_disp_mean_pks_x)
#plt.figure(1)
#plot_n_pks(corr_non_disp_mean_pks_x, f_out_non_disp_mean_pks_x, 'Horizontal non-disp noise location (max 5 peaks)')


# save all types of live data to .h5 file:
if caget('SR-APHLA{BPM}PSD:LiveData-Cmd') == 1: # Data Source: Live Data
    update_status("Saving data to .h5 file...")
    save_data(fa_xyas, prefix, bad_xy, mean_PSDs, int_mean_PSDs, mean_peaks_f)

t = time.time() - t0
update_status("Done! Waiting for a new cycle..." )
caput('SR-APHLA{BPM}PSD:LoopTime-I', t)
if caget('SR-APHLA{BPM}PSD:LiveData-Cmd') == 1: #update TS only for live data
    caput('SR-APHLA{BPM}PSD:Timestamp-Sts', str(datetime.datetime.now()))
