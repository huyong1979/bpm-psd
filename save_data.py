'''
Save 29 types of data to .h5 file as 6 groups:
  1. environment:
     1)beam_current: Storage Ring (SR) beam current;
     2)beam_bunch: number of bunches in the SR; 
     3)bpm_prefix: prefix for SR BPM PV names;
     4)bpm_location: BPM location in Z plane;
     5)bpm_badx: invalid BPM PV names in X plane;
     6)bpm_bady: invalid BPM PV names in Y plane;
  2. FA:
     1)X: all SR BPM FA data (223*100000) in X plane;
     2)Y: all SR BPM FA data (223*100000) in Y plane;
     3)S: all BPM FA data for sum of 4 buttons;
  3. averaged_PSDs:
     1)x_dispersive: averaged PSD of dispersive BPMs in X plane;
     2)x_non-dispersive: averaged PSD of non-dispersive BPMs in X plane;  
     3)x_ID: averaged PSD of ID BPMs in X plane;  
     4)y: averaged PSD of non-ID BPMs in Y plane;  
     5)y_ID: averaged PSD of ID BPMs in Y plane;  
  4. integral_PSDs:
     1)x_dispersive: integral PSD of averaged_PSDs/x_dispersive;
     2)x_non-dispersive: ...
     3)x_ID: ...  
     4)y: ...
     5)y_ID: ...
  5. peak_frequencies:
     1)x_dispersive: 5 peak frequencies of averaged_PSDs/x_dispersive;
     2)x_non-dispersive: ...
     3)x_ID: ...  
     4)y: ...
     5)y_ID: ...
  6. corrector_locations:
     1)x_dispersive: corrector locations for 5 peak frequencies;
     2)x_non-dispersive: ...
     3)x_ID: ...  
     4)y: ...
     5)y_ID: ...

Some data are passed from bpm_psd.py to the function save_data()  
'''

import numpy as np
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

import time
import datetime
#t0 = time.time()
#sys.path.append('/usr/lib/python2.7/dist-packages')
from cothread.catools import caget, caput, DBR_CHAR_STR
import h5py

def save_data(prefix, bad_xy, fa_xys, PSDs, int_PSDs, peaks_f, corr_locs):
    '''Save all kinds of data to .h5 file as groups'''
    keys = [
'environment/beam_current', 'environment/beam_bunch', 'environment/bpm_prefix', 
'environment/bpm_location', 'environment/bpm_badx',   'environment/bpm_bady',
'FA/X', 'FA/Y', 'FA/S',  
'averaged_PSDs/x_dispersive', 'averaged_PSDs/x_non-dispersive',
'averaged_PSDs/x_ID',         'averaged_PSDs/y', 'averaged_PSDs/y_ID',
'integral_PSDs/x_dispersive', 'integral_PSDs/x_non-dispersive',
'integral_PSDs/x_ID',         'integral_PSDs/y', 'integral_PSDs/y_ID',
'peak_frequencies/x_dispersive', 'peak_frequencies/x_non-dispersive',
'peak_frequencies/x_ID',         'peak_frequencies/y', 'peak_frequencies/y_ID',
'corrector_locations/x_dispersive', 'corrector_locations/x_non-dispersive',
'corrector_locations/x_ID', 'corrector_locations/y', 'corrector_locations/y_ID']

    beam_cur=caget('SR:C03-BI{DCCT:1}I:Total-I')
    n_bunch=caget('SR:C16-BI{FPM:1}NbrBunches-I') 

    #How to handle a list of strings in Python 3:
    #see https://github.com/h5py/h5py/issues/892
    prefix = np.array(prefix, dtype='S')
    [badx, bady] = [bad for bad in bad_xy]
    #print(badx) #['SR:C16-APHLA{BPM:6}PSD:BadX-Cmd', ...]
    badx = np.array(badx, dtype='S')
    #print(badx) #[b'SR:C16-APHLA{BPM:6}PSD:BadX-Cmd' b...]
    bady = np.array(bady, dtype='S')

    values =  [beam_cur, n_bunch, prefix, s_BPM,badx, bady] \
            + [fa for fa in fa_xys] \
            + [psd for psd in PSDs] \
            + [psd for psd in int_PSDs] \
            + [f for f in peaks_f] \
            + [loc for loc in corr_locs]

    #path = '/epics/data/bpm_psd_data/' #the directory where .h5 file is saved
    path = caget('SR-APHLA{BPM}PSD:Path-SP')
    file_name = str(path) + "bpm-fa-psd_" + time.strftime("%Y%b%d-%Hh%Mm") + ".h5"
    #file_name = str(path) + "test" + time.strftime("%Y%b%d-%Ham") + ".h5"
    hf = h5py.File(file_name, 'w')
    for (key, value) in zip (keys, values):
        hf[key] = value
        hf[key].attrs['beam_current'] = beam_cur
        hf[key].attrs['beam_bunches'] = n_bunch
 
    hf.close()
    print("%s: %d types of data are saved in %s"%(datetime.datetime.now(),
                                                len(keys), file_name))
    caput('SR-APHLA{BPM}PSD:h5Name-Wf', file_name, datatype=DBR_CHAR_STR)
