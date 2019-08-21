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
#import operator
#import heapq

#ver = caget('SR:C22-FOFB{CC}FpgaFirmVer-I')
#print ver

#create prefix lists for BPM and PSD PVs
prefix = []
psdPrefix = []#yhu
#regular BPMs: 6*30
p_index = ['1','2','3','4','5','6']
Cell_index =['30','01','02','03','04','05','06','07','08','09','10','11','12','13','14','15',
'16','17','18','19','20','21','22','23','24','25','26','27','28','29']
for i in Cell_index:
    for j in p_index:
        prefix.append('SR:C'+i+'-BI{BPM:'+j+'}')
        psdPrefix.append('SR:C'+i+'-APHLA{BPM:'+j+'}')
#ID BPMs
p_index=['7','8','9','10']
Cell_index=['04','07','12','19']
for i in Cell_index:
    for j in p_index:
        prefix.append('SR:C'+i+'-BI{BPM:'+j+'}')
        psdPrefix.append('SR:C'+i+'-APHLA{BPM:'+j+'}')
p_index=['7','8']
Cell_index=['02','03','08','10','11','16','18','21','28']
for i in Cell_index:
    for j in p_index:
        prefix.append('SR:C'+i+'-BI{BPM:'+j+'}')
        psdPrefix.append('SR:C'+i+'-APHLA{BPM:'+j+'}')
p_index=['7','8','9']
Cell_index=['05','17','23']
for i in Cell_index:
    for j in p_index:
        prefix.append('SR:C'+i+'-BI{BPM:'+j+'}')
        psdPrefix.append('SR:C'+i+'-APHLA{BPM:'+j+'}')

#a few small functions
def just_print(message):
    print "%s: %s"%(datetime.datetime.now(), str(message))
    #put array of characters as string
    caput('SR-APHLA{BPM}PSD:Status-Wf',str(message),datatype=DBR_CHAR_STR)

def print_and_exit(message):
    just_print(message)
    sys.exit()

def check_settings(pvList, value):
    values = caget(pvList)
    wrongSettings = []
    for i in range(len(pvList)):
        if values[i] != value:
            wrongSettings.append(pvList[i])
    if len(wrongSettings) != 0:
            print_and_exit("Wrong settings: "+str(wrongSettings[:2])+" ...")

# pv lists
FA_X=[]
FA_Y=[]
PSD_X = []
PSD_Y = []
PSD_INTX = [] #integral PSD for X plane
PSD_INTY = []
MachineType = []
TrigSrc = []
WfmSel = []
FaLen = []
TxStatus = []
for i in range(len(prefix)):
    FA_X.append(prefix[i]+'FA-X')
    FA_Y.append(prefix[i]+'FA-Y') 
    MachineType.append(prefix[i]+'Loc:Machine-SP') #SR:5
    TrigSrc.append(prefix[i]+'Trig:TrigSrc-SP') #Ext:1
    WfmSel.append(prefix[i]+'DDR:WfmSel-SP') #FA Wfm:2
    FaLen.append(prefix[i]+'Burst:FaEnableLen-SP') #100000 
    TxStatus.append(prefix[i]+'DDR:TxStatus-I') #Wfm Ready:0 
    PSD_X.append(psdPrefix[i]+'PSD:X-Wf') 
    PSD_Y.append(psdPrefix[i]+'PSD:Y-Wf') 
    PSD_INTX.append(psdPrefix[i]+'PSD:IntX-Wf') 
    PSD_INTY.append(psdPrefix[i]+'PSD:IntY-Wf') 

# get BPM FA data from PVs or file
fa_recordLen = caget('SR-APHLA{BPM}PSD:Len-SP')
#print fa_recordLen
if caget('SR-APHLA{BPM}PSD:LiveData-Cmd') == 1:  
    if caget('SR:C03-BI{DCCT:1}I:Real-I') < 0.01:
        print_and_exit("No beam, no calculation")
    if caget('INJ{TOC}OpControl-Sel')==2 and caget('INJ{TOC-SM}Cnt:Next-I')<20:    
        caput('SR-APHLA{BPM}PSD:Counter-Calc_.PROC', 1)
        print_and_exit("Top-off injection is coming with 20-sec so no calculation")

    check_settings(MachineType, 5)
    check_settings(TrigSrc, 1)
    check_settings(WfmSel, 2)
    #check_settings(FaLen, 100000)
    caput('SR:C21-PS{Pinger}Ping-Cmd', 1)
    just_print("Pinger trigger is sent, wait...")
    #time.sleep(11.0); # DO NOT use time.sleep()
    cothread.Sleep(12)
    for i in range(0,10):
        if sum(caget(TxStatus)) > 0:
            cothread.Sleep(1);
        else:
            break
        if i==9:
            print_and_exit("Timeout: FA data are not ready")

    just_print("Finally read and process live data...")
    x_all = caget(FA_X, count=fa_recordLen)
    y_all = caget(FA_Y, count=fa_recordLen)
else:
    path = '/home/skongtawong/Documents/Guimei/20190805_beamstudy_new_bpm_psd_and_gain/'
    filename = 'SR_AllIDBPMs_FA_20190805_0701_06_fofb_on_newPI999992_350Hz_pinger_off_02.h5'
    message = "read and process data from " + path+filename
    just_print(str(message))
    fid = h5py.File(path+filename, 'r')
    x_all = fid.get('faX')
    y_all = fid.get('faY')
    fid.close

x_all = np.asarray(x_all)[:,0:fa_recordLen]/1000.0 #unit: um
y_all = np.asarray(y_all)[:,0:fa_recordLen]/1000.0
#print "faX size:"
#print x_all.shape #(223, 100000)

# generate disp bpm index (start from 0 to 179)
disp = []
non_disp = []
for i in range(1,181):
    if np.mod(i,6) == 3 or np.mod(i,6) == 4:
        disp.append(i-1)
    else:
        non_disp.append(i-1)   
#print "dispersive bpm index is (start from 0 to 179)"
#print disp

# initialize 
#rf = 499.68e6 # would be good to read from PV
rf = caget('RF{Osc:1}Freq:I')
#print rf
fs = rf/1320/38
n = x_all.shape
r = n[0] # number of BPMs (rows 223)
c = n[1] # number of samples (columns 100000)
# for finding peaks
n_max = 5 # max n peaks
# prominence for finding x peaks (um^2/Hz)
(prom_x, prom_y) = caget(['SR-APHLA{BPM}PSD:PromX-SP','SR-APHLA{BPM}PSD:PromY-SP']) 
#print prom_x, prom_y
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

#main loop: PSD calculation
for i in range(0,r):
    x = x_all[i]
    y = y_all[i] 
    # ---------------- horizontal ------------------
    f, P = signal.welch(x, fs, window='hann', nperseg = c)
    Pxx_all = np.append(Pxx_all, P) # in 1D
    # find peaks of individual BPM's psd; no PVs for these peaks
    loc, _ = find_peaks(P, prominence = prom_x, distance = dist)
    pks = P[loc]
    # pick max N peaks
    loc_n = (-pks).argsort()
    loc_n = loc_n[0:n_max]
    #loc_n = np.sort(loc_n[0:n_max]) # sort frequency from minimum
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
    loc_n = loc_n[0:n_max]
    #loc_n = np.sort(loc_n[0:n_max]) # sort frequency from minimum
    loc_n = loc_n.astype(int)
    # collect peaks data for all vertical bpm
    pks_freq_y[i] = f[loc] # all pks frequency
    pks_hight_y[i] = pks # all pks hight
    pks_n_freq_y[i] = f[loc[loc_n]] # n max pks frequency
    pks_n_hight_y[i] = pks[loc_n] # n max pks hight
#t = time.time() - t0
#print "checkpoint after loop t = %.f" % (t)

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
# disp, non disp, insertion device (id)
badx = []
bady = [101,125] # start from 0
Pxx_non_disp = Pxx_all[:,non_disp] # no PVs for this
Pxx_disp = Pxx_all[:,disp]
Pxx_id = Pxx_all[:,180:r+1]
Pyy = Pyy_all[:,0:180]
Pyy = np.delete(Pyy,bady,1)
Pyy_id = Pyy_all[:,180:r+1]

# find integrated PSD (2D)
df = f[1]-f[0] #~0.1Hz
int_Pxx_all = np.sqrt(np.cumsum(Pxx_all, axis = 0)*df)
int_Pyy_all = np.sqrt(np.cumsum(Pyy_all, axis = 0)*df)
int_Pxx_all_4ca = np.transpose(int_Pxx_all)
int_Pyy_all_4ca = np.transpose(int_Pyy_all)
caput(PSD_INTX, int_Pxx_all_4ca[0:r,1:])
caput(PSD_INTY, int_Pyy_all_4ca[0:r,1:])
int_Pxx_non_disp = np.sqrt(np.cumsum(Pxx_non_disp, axis = 0)*df) # no PVs for this
int_Pxx_disp = np.sqrt(np.cumsum(Pxx_disp, axis = 0)*df)
int_Pxx_id = np.sqrt(np.cumsum(Pxx_id, axis = 0)*df)
int_Pyy = np.sqrt(np.cumsum(Pyy, axis = 0)*df)
int_Pyy_id = np.sqrt(np.cumsum(Pyy_id, axis = 0)*df)

# PSD mean of all bpm for different types
Pxx_disp_mean = np.mean(Pxx_disp, axis=1)
Pxx_non_disp_mean = np.mean(Pxx_non_disp, axis=1)
Pxx_id_mean = np.mean(Pxx_id, axis=1)
Pyy_mean = np.mean(Pyy, axis=1)
Pyy_id_mean = np.mean(Pyy_id, axis=1)
pvs = ['SR-APHLA{BPM}PSD:X_DISP_MEAN-Wf','SR-APHLA{BPM}PSD:X_NON_DISP_MEAN-Wf',
'SR-APHLA{BPM}PSD:X_ID_MEAN-Wf','SR-APHLA{BPM}PSD:Y_MEAN-Wf','SR-APHLA{BPM}PSD:Y_ID_MEAN-Wf']
values = [np.sqrt(Pxx_disp_mean)[1:],np.sqrt(Pxx_non_disp_mean)[1:],np.sqrt(Pxx_id_mean)[1:],np.sqrt(Pyy_mean)[1:],np.sqrt(Pyy_id_mean)[1:]]
caput(pvs, values)

# find int. PSD mean
sf = caget("SR-APHLA{BPM}PSD:StarFreq-SP")#start freq for integral PSD
Isf = np.where(f>=sf)
Isf = Isf[0][0]
f_int = f[Isf:]
caput('SR-APHLA{BPM}Freq-Wf_', f_int)# frequency for integral PSD
int_Pxx_disp_mean = np.sqrt(np.cumsum(Pxx_disp_mean[Isf:])*df)
int_Pxx_non_disp_mean = np.sqrt(np.cumsum(Pxx_non_disp_mean[Isf:])*df)
int_Pxx_id_mean = np.sqrt(np.cumsum(Pxx_id_mean[Isf:])*df)
int_Pyy_mean = np.sqrt(np.cumsum(Pyy_mean[Isf:])*df)
int_Pyy_id_mean = np.sqrt(np.cumsum(Pyy_id_mean[Isf:])*df)
pvs = ['SR-APHLA{BPM}PSD:IntX_DISP_MEAN-Wf','SR-APHLA{BPM}PSD:IntX_NON_DISP_MEAN-Wf',
'SR-APHLA{BPM}PSD:IntX_ID_MEAN-Wf','SR-APHLA{BPM}PSD:IntY_MEAN-Wf','SR-APHLA{BPM}PSD:IntY_ID_MEAN-Wf']
values = [int_Pxx_disp_mean,int_Pxx_non_disp_mean,int_Pxx_id_mean,int_Pyy_mean,int_Pyy_id_mean]
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
pvs = ['SR-APHLA{BPM}PSD:IntX_DISP_MEAN_F0F1-I','SR-APHLA{BPM}PSD:IntX_NON_DISP_MEAN_F0F1-I',
'SR-APHLA{BPM}PSD:IntX_ID_MEAN_F0F1-I','SR-APHLA{BPM}PSD:IntY_MEAN_F0F1-I','SR-APHLA{BPM}PSD:IntY_ID_MEAN_F0F1-I',
'SR-APHLA{BPM}PSD:IntX_DISP_MEAN_F1F2-I','SR-APHLA{BPM}PSD:IntX_NON_DISP_MEAN_F1F2-I',
'SR-APHLA{BPM}PSD:IntX_ID_MEAN_F1F2-I','SR-APHLA{BPM}PSD:IntY_MEAN_F1F2-I','SR-APHLA{BPM}PSD:IntY_ID_MEAN_F1F2-I',
'SR-APHLA{BPM}PSD:IntX_DISP_MEAN_F2F3-I','SR-APHLA{BPM}PSD:IntX_NON_DISP_MEAN_F2F3-I',
'SR-APHLA{BPM}PSD:IntX_ID_MEAN_F2F3-I','SR-APHLA{BPM}PSD:IntY_MEAN_F2F3-I','SR-APHLA{BPM}PSD:IntY_ID_MEAN_F2F3-I']
values = [int_Pxx_disp_mean_f0f1,int_Pxx_non_disp_mean_f0f1,int_Pxx_id_mean_f0f1,int_Pyy_mean_f0f1,int_Pyy_id_mean_f0f1,
int_Pxx_disp_mean_f1f2,int_Pxx_non_disp_mean_f1f2,int_Pxx_id_mean_f1f2,int_Pyy_mean_f1f2,int_Pyy_id_mean_f1f2,
int_Pxx_disp_mean_f2f3,int_Pxx_non_disp_mean_f2f3,int_Pxx_id_mean_f2f3,int_Pyy_mean_f2f3,int_Pyy_id_mean_f2f3]
caput(pvs, values)

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
    loc_n = loc_n[0:n_max]
    loc_n = loc_n.astype(int)   
    if len(loc_n) < n_max:
        f_n = [df]*n_max # avoid 0 frequency for logplot
        pks_n = [0]*n_max
    else:
        pks_n = pks[loc_n]
        f_n = f[loc[loc_n]]
    #print pks_n
    #print f_n
    # collect peaks data for all horizontal bpm
    if i == 0:
        Pxx_disp_mean_pks_freq = f[loc] # all pks locations
        Pxx_disp_mean_pks_hight = pks # all pks hight
        Pxx_disp_mean_pks_n_freq = f_n # n max pks locations
        Pxx_disp_mean_pks_n_hight = pks_n # n max pks hight
    elif i == 1:
        Pxx_non_disp_mean_pks_freq = f[loc] # all pks locations
        Pxx_non_disp_mean_pks_hight = pks # all pks hight
        Pxx_non_disp_mean_pks_n_freq = f_n # n max pks locations
        Pxx_non_disp_mean_pks_n_hight = pks_n # n max pks hight
    elif i == 2:
        Pxx_id_mean_pks_freq = f[loc] # all pks locations
        Pxx_id_mean_pks_hight = pks # all pks hight
        Pxx_id_mean_pks_n_freq = f_n # n max pks locations
        Pxx_id_mean_pks_n_hight = pks_n # n max pks hight
    elif i == 3:
        Pyy_mean_pks_freq = f[loc] # all pks locations
        Pyy_mean_pks_hight = pks # all pks hight
        Pyy_mean_pks_n_freq = f_n # n max pks locations
        Pyy_mean_pks_n_hight = pks_n # n max pks hight
    elif i == 4:
        Pyy_id_mean_pks_freq = f[loc] # all pks locations
        Pyy_id_mean_pks_hight = pks # all pks hight
        Pyy_id_mean_pks_n_freq = f_n # n max pks locations
        Pyy_id_mean_pks_n_hight = pks_n # n max pks hight
pvs = ['SR-APHLA{BPM}PSD:X_DISP_MEAN_PKS_N_F-Wf','SR-APHLA{BPM}PSD:X_DISP_MEAN_PKS_N_H-Wf',
'SR-APHLA{BPM}PSD:X_NON_DISP_MEAN_PKS_N_F-Wf','SR-APHLA{BPM}PSD:X_NON_DISP_MEAN_PKS_N_H-Wf',
'SR-APHLA{BPM}PSD:X_ID_MEAN_PKS_N_F-Wf','SR-APHLA{BPM}PSD:X_ID_MEAN_PKS_N_H-Wf',
'SR-APHLA{BPM}PSD:Y_MEAN_PKS_N_F-Wf','SR-APHLA{BPM}PSD:Y_MEAN_PKS_N_H-Wf',
'SR-APHLA{BPM}PSD:Y_ID_MEAN_PKS_N_F-Wf','SR-APHLA{BPM}PSD:Y_ID_MEAN_PKS_N_H-Wf']
#print len(Pxx_disp_mean_pks_n_freq)
#print len(Pyy_id_mean_pks_n_hight)
values = [Pxx_disp_mean_pks_n_freq,Pxx_disp_mean_pks_n_hight,Pxx_non_disp_mean_pks_n_freq,Pxx_non_disp_mean_pks_n_hight,
Pxx_id_mean_pks_n_freq,Pxx_id_mean_pks_n_hight,Pyy_mean_pks_n_freq,Pyy_mean_pks_n_hight,Pyy_id_mean_pks_n_freq,Pyy_id_mean_pks_n_hight]
caput(pvs, values)
#print "Avg Horizontal n peaks: "
#print Pxx_non_disp_mean_pks_n_freq
#print np.sqrt(Pxx_non_disp_mean_pks_n_hight)
t = time.time() - t0
message = "Done! t = %.f seconds" %(t)
just_print(str(message))
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

