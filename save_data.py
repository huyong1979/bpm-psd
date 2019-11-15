'''
Save these types of live data to .h5 file:
  1) faX (x_all): all BPM FA data in X plane;
  2) faY (y_all): all BPM FA data in Y plane;
  3) faA (a_all): all BPM FA data for button A;
  4) faS (s_all): all BPM FA data for sum of 4 buttons;
  5) beamCur (beam_cur): SR beam current;
  6) nBunch (n_bunch): number of bunches;
  7) prefix: prefix for BPM PVs;
  8) s (s_BPM): BPM location in Z plane;
  9) badX (badx_pvname): invalid BPM PV names in X plane;
  10)badY (bady_pvname): invalid BPM PV names in Y plane;
  11)xDispMeanPSD (Pxx_disp_mean): averaged PSD of dispersive BPMs in X plane;
  12)xNonDispMeanPSD (*non_disp*): averaged PSD of non-dispersive BPMs in X plane;
  13)xIDMeanPSD (Pxx_id_mean): averaged PSD of ID BPMs in X plane;
  14)yMeanPSD (Pyy_mean): averaged PSD of non-ID BPMs in Y plane;
  15)yIDMeanPSD (Pyy_id_mean): averaged PSD of ID BPMs in Y plane;
  16)intXDispMeanPSD (int_Pxx_disp_mean): integral PSD of xDispMeanPSD;
  17)intXNonDispMeanPSD (int_Pxx_non_disp_mean): ...;
  18)intXIDMeanPSD (int_Pxx_id_mean): ...;
  19)intYMeanPSD (int_Pyy_mean): ...;
  20)intYIDMeanPSD (int_Pyy_id_mean): integral PSD of yIDMeanPSD;;
  21)xDispMeanPeaksFreq (Pxx_disp_mean_pks_n_freq): 5 peak frequencies of xDispMeanPSD;
  22)xNonDispMeanPeaksFreq (Pxx_non_disp_mean_pks_n_freq): ...;
  23)xIDMeanPeaksFreq (Pxx_id_mean_pks_n_freq): ...;
  24)yMeanPeaksFreq (Pyy_mean_pks_n_freq): ...;
  25)yIDMeanPeaksFreq (Pyy_id_mean_pks_n_freq): 5 peak frequencies of yIDMeanPSD;
'''

import time
import datetime
t0 = time.time()
import sys
sys.path.append('/usr/lib/python2.7/dist-packages')
from cothread.catools import caget
import numpy as np
import h5py

path = '/epics/data/bpm_psd_data/' #the directory where .h5 file is saved
s_BPM=np.array([  
         4.935   ,   7.46002 ,   13.1446 ,   15.3773 ,   20.2472 ,
         22.8109 ,   29.9886 ,   32.5523 ,   38.3018 ,   40.5345 ,
         45.3368 ,   47.8618 ,   57.7322 ,   60.2572 ,   65.9418 ,
         68.1745 ,   73.0444 ,   75.6081 ,   82.7858 ,   85.3495 ,
         91.099  ,   93.3317 ,   98.134  ,  100.659  ,  110.529  ,
        113.054  ,  118.739  ,  120.972  ,  125.842  ,  128.405  ,
        135.583  ,  138.147  ,  143.896  ,  146.129  ,  150.931  ,
        153.456  ,  163.327  ,  165.852  ,  171.536  ,  173.769  ,
        178.639  ,  181.202  ,  188.38   ,  190.944  ,  196.693  ,
        198.926  ,  203.728  ,  206.253  ,  216.124  ,  218.649  ,
        224.333  ,  226.566  ,  231.436  ,  234.     ,  241.177  ,
        243.741  ,  249.491  ,  251.723  ,  256.526  ,  259.051  ,
        268.921  ,  271.446  ,  277.131  ,  279.363  ,  284.233  ,
        286.797  ,  293.975  ,  296.538  ,  302.288  ,  304.52   ,
        309.323  ,  311.848  ,  321.718  ,  324.243  ,  329.928  ,
        332.16   ,  337.03   ,  339.594  ,  346.772  ,  349.335  ,
        355.085  ,  357.318  ,  362.12   ,  364.645  ,  374.515  ,
        377.04   ,  382.725  ,  384.958  ,  389.828  ,  392.391  ,
        399.569  ,  402.133  ,  407.882  ,  410.115  ,  414.917  ,
        417.442  ,  427.313  ,  429.838  ,  435.522  ,  437.755  ,
        442.625  ,  445.188  ,  452.366  ,  454.93   ,  460.679  ,
        462.912  ,  467.714  ,  470.239  ,  480.11   ,  482.635  ,
        488.319  ,  490.552  ,  495.422  ,  497.986  ,  505.163  ,
        507.727  ,  513.477  ,  515.709  ,  520.512  ,  523.037  ,
        532.907  ,  535.432  ,  541.117  ,  543.349  ,  548.219  ,

        550.783  ,  557.961  ,  560.524  ,  566.274  ,  568.506  ,
        573.309  ,  575.834  ,  585.704  ,  588.229  ,  593.914  ,
        596.146  ,  601.016  ,  603.58   ,  610.758  ,  613.321  ,
        619.071  ,  621.304  ,  626.106  ,  628.631  ,  638.501  ,
        641.026  ,  646.711  ,  648.944  ,  653.814  ,  656.377  ,
        663.555  ,  666.119  ,  671.868  ,  674.101  ,  678.903  ,
        681.428  ,  691.299  ,  693.824  ,  699.508  ,  701.741  ,
        706.611  ,  709.175  ,  716.352  ,  718.916  ,  724.665  ,
        726.898  ,  731.7    ,  734.225  ,  744.096  ,  746.621  ,
        752.305  ,  754.538  ,  759.408  ,  761.972  ,  769.149  ,
        771.713  ,  777.463  ,  779.695  ,  784.498  ,  787.023  ])

def save_data(fa_xyas, prefix, bad_xy, mean_PSDs, int_mean_PSDs, mean_peaks_f):
    '''Save all kinds of live data to .h5 file'''
    #x_all, y_all, a_all, s_all = fa_xyas[0], fa_xyas[1], fa_xyas[2], fa_xyas[3];
    [x_all, y_all, a_all, s_all] = [fa for fa in fa_xyas]
    [badx_pvname, bady_pvname] = [bad for bad in bad_xy]
    #print(bady_pvname)
    [Pxx_disp_mean, Pxx_non_disp_mean, Pxx_id_mean, Pyy_mean, Pyy_id_mean] = \
        [psd for psd in mean_PSDs]  

    [int_Pxx_disp_mean, int_Pxx_non_disp_mean, int_Pxx_id_mean,
     int_Pyy_mean, int_Pyy_id_mean] = [psd for psd in int_mean_PSDs]  

    [Pxx_disp_mean_pks_n_freq, Pxx_non_disp_mean_pks_n_freq, Pxx_id_mean_pks_n_freq,
     Pyy_mean_pks_n_freq, Pyy_id_mean_pks_n_freq] = [f for f in mean_peaks_f] 

    beam_cur=caget('SR:C03-BI{DCCT:1}I:Total-I')
    n_bunch=caget('SR:C16-BI{FPM:1}NbrBunches-I') 

    values = [x_all, y_all, a_all, s_all] + [beam_cur, n_bunch, prefix, s_BPM] \
    + [badx_pvname, bady_pvname] \
    + [Pxx_disp_mean, Pxx_non_disp_mean, Pxx_id_mean, Pyy_mean, Pyy_id_mean] \
    + [int_Pxx_disp_mean, int_Pxx_non_disp_mean, int_Pxx_id_mean,
       int_Pyy_mean, int_Pyy_id_mean]\
    + [Pxx_disp_mean_pks_n_freq, Pxx_non_disp_mean_pks_n_freq, Pxx_id_mean_pks_n_freq, 
       Pyy_mean_pks_n_freq, Pyy_id_mean_pks_n_freq]
    
    fields = ['faX', 'faY', 'faA', 'faS', 'beamCur', 'nBunch', 'prefix', 's',
    'badX', 'badY', 'xDispMeanPSD', 'xNonDispMeanPSD', 'xIDMeanPSD', 
    'yMeanPSD', 'yIDMeanPSD', 'intXDispMeanPSD', 'intXNonDispMeanPSD',
    'intXIDMeanPSD', 'intYMeanPSD', 'intYIDMeanPSD', 'xDispMeanPeaksFreq',
    'xNonDispMeanPeaksFreq', 'xIDMeanPeaksFreq', 'yMeanPeaksFreq', 'yIDMeanPeaksFreq']
    
    #file_name = path + "bpm-fa-psd_" + time.strftime("%Y%b%d-%H%M%S") + ".h5"
    file_name = path + "bpm-fa-psd_" + time.strftime("%Y%b%d") + ".h5"
    try:
        with h5py.File(file_name, 'w') as f:
            for (field, value) in zip(fields, values):
                f[field] = value
    except IOError as e:
        print e.message

    print("%s: data saved in %s"%(datetime.datetime.now(), file_name))
