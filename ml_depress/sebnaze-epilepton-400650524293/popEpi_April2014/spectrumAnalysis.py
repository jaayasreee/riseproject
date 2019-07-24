'''
Created on 1 feb. 2013

@author: Squirel
'''

import pylab as pb
import h5py
from scipy.io import loadmat
from numpy import fft, angle, abs, load
import numpy as np

"""
y: timeserie; Fs: sampling freq; fmin, fmax: minimum & maximum freqs; dt: step interval
"""
def getPowSpectrum(y,Fs, fmin, fmax, dt):
    n=len(y)
    Y=fft.fft(y)
    
    return abs(Y)**2, fft.fftfreq(n, dt)



def highPassFilter(signal):
    #path="C:/Users/Squirel/workspace/HC-Septum-Network/popEpi_"
    #path="./traces_22feb/"
    """x1plot = load("./traces_22feb/x1plot_40_HMR_40_ML_CpES_50_gx1x2_15_gx2x1_15_x0_300_10s.npy")
    x2plot = load("./traces_22feb/x2plot_40_HMR_40_ML_CpES_50_gx1x2_15_gx2x1_15_x0_300_10s.npy")
    tplot  = load("./traces_22feb/tplot_40_HMR_40_ML_CpES_50_gx1x2_15_gx2x1_15_x0_300_10s.npy")"""
    
    #epi = x1plot.mean(axis=0)+x2plot.mean(axis=0)
    spectrum = fft.fftn(signal)
    freqs = fft.fftfreq(signal.size, d=100)
    
    bp=spectrum[:]
    for i in range(len(bp)):
        if i <= 100: 
            bp[i]=0
            bp[-i]=0
    ibp=fft.ifftn(bp)
    return ibp
    """fig=pb.figure()
    ax1=pb.subplot(211)
    ax1.plot(freqs, abs(bp),'.b')
    ax2=pb.subplot(212)
    ax2.plot(tplot, ibp)
    pb.show()"""
    
def bandPassFilter(signal):
    spectrum = fft.fft(signal)
    freqs = fft.fftfreq(signal.size, d=1000)
    
    bp=spectrum[:]
    for i in range(len(bp)/2):
        if i < 100 or i > 10000 : 
            bp[i]=0
            bp[-i]=0
    ibp=fft.ifft(bp)
    return ibp
    
    
def SE_spectrum_analysis(sim_exp):
    #load matlab file where timeserie is stored
    #f = h5py.File('C:/Users/Squirel/Desktop/data_antoine/status_epilepticus/timeserie_SE_4j_rat8l_sr250Hz.mat', 'r')
    if sim_exp=='exp':
        f = h5py.File('C:/Users/Squirel/Desktop/data_antoine/status_epilepticus/EEG_rat8l.mat', 'r')
        data = f.get('EEG/data')
        data = pb.array(data) #transform to numpy array --> now each element of data is an array of one float
        sf = 250. # sampling freq (Hz)
        beg=120000; lgth = 20.; #stating time (s) and length (s) of the selected period
        spectrum_ymax = 0.2 #maximum value displayed in plot (to have all the same scale)
        ts_ylim = (-0.0020, 0.0010)
    else: #sim_exp=='sim'
        #radical='_40_HMR_40_ML_CpES1_0_CpES2_0_x0_35_noise3_0_noise2_30_noise1_60_gx1x2_20_gx2x1_20_gx1x1_20_gx2x2_20_r40_20s' #phase I
        #radical='_40_HMR_40_ML_CpES1_20_CpES2_20_x0_20_noise3_0_noise2_30_noise1_60_gx1x2_20_gx2x1_20_gx1x1_20_gx2x2_20_r40_20s' # phase II
        #radical='_40_HMR_40_ML_CpES1_80_CpES2_80_x0_20_noise3_0_noise2_30_noise1_60_gx1x2_20_gx2x1_20_gx1x1_20_gx2x2_20_r40_20s' #phase III
        radical='_40_HMR_40_ML_CpES1_40_CpES2_40_x0_35_noise3_0_noise2_30_noise1_60_gx1x2_20_gx2x1_20_gx1x1_20_gx2x2_20_r40_100s' #phase IV     !!! set beg=20 instead of 10s
        x1_plot = np.load('./traces_dec13/x1plot'+radical+'.npy')
        x2_plot = np.load('./traces_dec13/x2plot'+radical+'.npy')
        #downsampling to 1kHz
        x1_plot =  x1_plot[:,0:-1:10]
        x2_plot = x2_plot[:,0:-1:10]
        
        data = -(0.75*x1_plot.mean(axis=0) + 0.25*x2_plot.mean(axis=0)) #transform to numpy array --> now each element of data is an array of one float
        sf = 1000. # sampling freq (Hz)
        beg=20.; lgth = 80.; #stating time (s) and length (s) of the selected period
        spectrum_ymax = 300 #maximum value displayed in plot (to have all the same scale)
        ts_ylim = (-0.5, 1.75)
    
    minfreq=1.; maxfreq=50.; #Hz
    #pb.plot(range(lgth*sf), data[beg*sf:beg*sf+lgth*sf]); pb.show() #if splot needed 
    signal = data[beg*sf:beg*sf+lgth*sf]
    #signal = pb.sin(2*pb.pi*pb.arange(0,15,0.004))
    #pb.plot(signal); pb.show()
    spectrum = fft.fftn(signal)
    delta_f = 1.0 / ((lgth*sf)/2) * lgth
    #pb.plot(pb.arange(minfreq, maxfreq, delta_f),abs(spectrum[minfreq*delta_f:maxfreq*delta_f]), "."); pb.show()
    #pb.plot(pb.arange(minfreq, maxfreq, 1/lgth), abs(spectrum[minfreq*lgth:maxfreq*lgth]), "."); pb.show()
    
    spec = pb.arange(minfreq,maxfreq,1.0)
    for i in range(len(spec)):
        spec[i] = pb.mean(abs(spectrum[(minfreq+i)*lgth:(minfreq+i)*lgth + lgth]))
    #pb.plot(spec)
    #pb.show()
    
    #plots:
    fig = pb.figure()
    ax1 = pb.subplot2grid((2,2),(0,0),colspan=2)
    ax2 = pb.subplot2grid((2,2),(1,0))
    ax3 = pb.subplot2grid((2,2),(1,1))
    
    ax1.plot(signal)
    ax1.set_ylim(ts_ylim)
    ax2.plot(pb.arange(minfreq, maxfreq, 1/lgth), abs(spectrum[minfreq*lgth:maxfreq*lgth]), ".")
    ax2.set_ylim(0,spectrum_ymax)
    ax3.plot(spec)
    pb.show()
    
    
def okeedokee_spectrum():
    f = h5py.File('C:/Users/Squirel/Desktop/data_abdel/Data/okeedokee_dCA3_period02.mat', 'r') #load the data
    data = f.get('ts_dCA3')
    data = pb.array(data) #transform to numpy array --> now each element of data is an array of one float
    sf = 1250. # sampling freq (Hz)
    beg=40; lgth = 10.; #stating time (s) and length (s) of the selected period
    minfreq=1.; maxfreq=50.; #Hz 
    signal = data[beg*sf:beg*sf+lgth*sf] #troncate
    spectrum = fft.fftn(signal) #extract spectra
    spec = pb.arange(minfreq, maxfreq,1.0) #
    for i in range(len(spec)):
        spec[i] = pb.mean(abs(spectrum[(minfreq+i)*lgth:(minfreq+i)*lgth]))
    #display
    fig = pb.figure()
    ax1 = pb.subplot2grid((1,3),(0,1),colspan=2)
    ax2 = pb.subplot2grid((1,3),(0,0))
    
    
    ax1.plot(signal)
    ax2.plot(pb.arange(minfreq, maxfreq, 1/lgth), abs(spectrum[minfreq*lgth:maxfreq*lgth]), ".")
    pb.show()
    
if __name__ == "__main__":
    #highPassFilter()
    #SE_spectrum_analysis('exp')
    okeedokee_spectrum()