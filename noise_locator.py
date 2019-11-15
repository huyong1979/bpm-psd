# Noise locator from FA data
# 10/18/2019
import datetime
import h5py
import numpy as np
import matplotlib.pyplot as plt

def noise_locator(x_all, y_all, goodx, goody): # for both planes
    # get ORM form .h5
    # note: ORM start from cell 1, not cell 30. each cell has 8 bpm, bpm 1-6 are for norm BPM, bpm 7-8 are for ID, (begin from 1)
    #print 'start noise locator...'
    fid1 = h5py.File('ORM2017.h5', 'r')
    ORM = fid1.get('ORM')
    fid1.close
    ORM = np.array(ORM) #shape: 480*180
    
    n_corr = 90 # fast corrector
    n_bpm = 240 # all bpm in one plane from CC
    n_bpm_norm = 180 # non-ID bpm
    bpm_one_cell = 6 # non-ID bpm per cell
    
    # generate bpm index for ORM (start from 0)
    orm_id = []
    orm_norm = []
    for i in range(1,(n_bpm+1)):
        if np.mod(i,8) == 7 or np.mod(i,8) == 0:
            i
            orm_id.append(i-1)
        else:
            orm_norm.append(i-1)
    Rx = ORM[orm_norm, 0:n_corr] #180*90
    Ry = ORM[n_bpm:, n_corr:] #240*90
    Ry = Ry[orm_norm, 0:n_corr] #180*90
    
    # shift cell 30 to the top (6 rows)
    def rotate_c30(Rin):
        Rout = np.zeros(Rin.shape)
        Rout[0:bpm_one_cell, :] = Rin[(n_bpm_norm-bpm_one_cell):n_bpm_norm, :]
        Rout[bpm_one_cell:, :] = Rin[0:(n_bpm_norm-bpm_one_cell), :]
        return Rout
    
    # switch BPM domain of ORM to 30,1,2,...,29
    Rx = rotate_c30(Rx)
    Ry = rotate_c30(Ry)

    # reconstruct data: f, x, y, ORM
    x = x_all[goodx, :]
    y = y_all[goody, :]
    Rx = Rx[goodx, :]
    Ry = Ry[goody, :]
    
    # find inverse with Tikhonov
    def Inv_Tikhonov(M, beta):
        n_bpm, n_corr = np.shape(M)
        U, s, VT = np.linalg.svd(M, full_matrices=False)
        UT = np.transpose(U)
        Sinv = np.diag(s/(s*s+beta))
        V = np.transpose(VT)
        Rinv = np.matmul(V, np.matmul(Sinv, UT))
        return Rinv
    
    # calculate corrector strength
    Rx_inv = Inv_Tikhonov(Rx, 0.01)
    Ry_inv = Inv_Tikhonov(Ry, 0.01)
    ffx = np.fft.fft(x)
    ffy = np.fft.fft(y)
    Cx = np.abs(np.matmul(Rx_inv, ffx))
    Cy = np.abs(np.matmul(Ry_inv, ffy))
    print('%s: noise locator complete...'%datetime.datetime.now())
    return Cx, Cy

def locate_n_peaks(corr_all, f_in, f_n_pks): #for one plane
    # Find index of pks
    n_pks = len(f_n_pks)
    i_n_pks = np.zeros([n_pks]) # allowcate for index of f_n_pks
    for i in range(0,len(f_n_pks)):
        tmp = np.where(f_in >= f_n_pks[i])
        i_n_pks[i] = tmp[0][0]    
        
    i_n_pks = i_n_pks.astype(int)    
    f_out = f_in[i_n_pks]
    corr_n_pks = corr_all[:, i_n_pks]
    # norm_corr_n_pks = corr_n_pks/np.max(corr_n_pks, axis = 0)
    return corr_n_pks, f_out

def plot_n_pks(corr_n_pks, f_pks, name):
    n_corr, n_subplot = corr_n_pks.shape
    corr_no = range(1, n_corr+1)
    for i in range(0, n_subplot):
        plt.subplot(n_subplot, 1, i+1)
        plt.plot(corr_no, corr_n_pks[:,i])
        lg = np.round(f_pks[i], 2)
        lg = str(lg) + ' Hz'
        if i == 0:
            plt.title(name)
        plt.legend([lg])
        plt.xlabel('Corr No.')
        plt.ylabel('Corr strength [au]')
        plt.xlim(1, n_corr)
        plt.show()
        
def plot_mesh(corr_all, f_in, f_i, f_f): #plot topview of noise locator from f_i to f_f
    tmp = np.where(f_in >= f_i)
    i_i = tmp[0][0]    
    tmp = np.where(f_in >= f_f)
    i_f = tmp[0][0]
    # find the real frequency
    f_i = f_in[i_i] 
    f_f = f_in[i_f]
    plt.imshow(corr_all[:, i_i:i_f], aspect = 'auto', extent=[f_i, f_f, 0, 90])
    plt.ylabel('Corr No.')
    plt.xlabel('Frequency [Hz]')
    plt.show()
