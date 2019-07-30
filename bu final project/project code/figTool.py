import numpy as np
import pylab as pb
import matplotlib as mpl
from spectrumAnalysis import highPassFilter, bandPassFilter


def zoom(beg, end, x1_plot, y1_plot, z_plot, x2_plot, y2_plot, t_plot, KOP_plot, radical):
    #resize sample according zoom interval
    x1_plot = x1_plot[:,beg/0.05:end/0.05]
    x2_plot = x2_plot[:,beg/0.05:end/0.05]
    y1_plot = y1_plot[:,beg/0.05:end/0.05]
    y2_plot = y2_plot[:,beg/0.05:end/0.05]
    z_plot = z_plot[:,beg/0.05:end/0.05]
    t_plot = t_plot[beg/0.05:end/0.05]
    KOP_plot = KOP_plot[0, beg/0.05:end/0.05] #0 because k only needed, no psi
    
    nbn1=x1_plot.shape[0]
    nbn2=x2_plot.shape[0]
    x1bar_plot = pb.zeros(x1_plot.shape[1])
    x2bar_plot = pb.zeros(x2_plot.shape[1]) 
    zbar_plot = pb.zeros(z_plot.shape[1])
    
    for i in range(x1bar_plot.size):
        x1bar_plot[i]=pb.mean(x1_plot[:,i])
        x2bar_plot[i]=pb.mean(x2_plot[:,i])
        zbar_plot[i]=pb.mean(z_plot[:,i])
    
    #plotting    
    fig = pb.figure(figsize=(20,10))
    
    ax1=pb.subplot(5,1,1)
    ax1.set_title("x1 (thick black=x1bar)")
    ax2=pb.subplot(5,1,2)
    ax2.set_title("x2 (thick black=x2bar)")
    ax3=pb.subplot(5,1,3)
    ax3.set_title("x1bar - x2bar")
    ax4=pb.subplot(5,1,4)
    ax4.set_title("Z (thick black=zbar)")
    ax5=pb.subplot(5,1,5)
    ax5.set_title("Amplitude of the Kuramoto Order parameter")
    
    for i in range(nbn1):
        #time series pop1
        ax1.plot(t_plot, x1_plot[i,:]) #i -> all neurons, 0 -> only neuron 0 ... 
        #time series z
        ax4.plot(t_plot, z_plot[i,:], label=None)
    for j in range(nbn2):
        #time series pop2
        ax2.plot(t_plot, x2_plot[j,:])
        
    #draw time series
    ax1.plot(t_plot, x1bar_plot, 'black', linewidth=1.5)
    ax2.plot(t_plot, x2bar_plot, 'black', linewidth=1.5)
    ax3.plot(t_plot, x2bar_plot - x1bar_plot, label='x2bar - x1bar')
    ax3.legend(prop={'size':10})
    ax4.plot(t_plot, zbar_plot, 'black', linewidth=2., label="zbar")
    ax4.legend(prop={'size':10})
    ax5.plot(t_plot, KOP_plot[:])
    ax5.legend(prop={'size':10})
    
    fig.savefig("epilepton"+radical+"_zoom.png", dpi=200)
    pb.show()
    
    
def synchPlotPop2():
    #load data
    #radical='_40_HMR_40_ML_CpES1_100_CpES2_100_x0_30_noise3_0_noise2_40_noise1_80_gx1x2_20_gx2x1_20_gx1x1_10_gx2x2_10_r1_100s'
    radical='_40_HMR_40_ML_CpES1_60_CpES2_60_x0_30_noise3_0_noise2_40_noise1_80_gx1x2_30_gx2x1_30_gx1x1_20_gx2x2_20_r40_80s'
    x1_plot = np.load('./traces_march25/x1plot'+radical+'.npy')
    x2_plot = np.load('./traces_march25/x2plot'+radical+'.npy')
    z_plot =  np.load('./traces_march25/zplot'+radical+'.npy')
    #KOP2_plot = np.load('./traces_07feb/KOP2noRest_40_HMR_40_ML_var_CpES_x0_2500to12000ms.npy')
    #x2_plot = x2_plot[:,2500:12000]
    #KOP2_plot = KOP2_plot[:,2500:12000]
    
    #get minimums and maximums before resizing
    x1max=x1_plot.max()
    x1min=x1_plot.min()
    x2max=x2_plot.max()
    x2min=x2_plot.min()
    
    #downsampling
    x1_plot = x1_plot[:,0:-1:10]
    x2_plot = x2_plot[:,0:-1:10]
    z_plot = z_plot[:,0:-1:10]
    
    x1_plot_nf = x1_plot.copy() #non-filtered
    x2_plot_nf = x2_plot.copy() #non-filtered
    
    #filter signal before show
    for i in range(x1_plot.shape[0]):
        x1_plot[i,:]=bandPassFilter(x1_plot[i,:])
        x2_plot[i,:]=bandPassFilter(x2_plot[i,:])
    
    beg=000
    end=80000
    # fs=1kHz: sampling freq !!! downsampling above
    #resize sample according above interval
    x1_plot = x1_plot[:,beg:end]
    x2_plot = x2_plot[:,beg:end]
    #y1_plot = y1_plot[:,beg:end]
    #y2_plot = y2_plot[:,beg:end]
    z_plot = z_plot[:,beg:end]
    #t_plot = t_plot[beg:end]
    x1_plot_nf = x1_plot_nf[:,beg:end]
    x2_plot_nf = x2_plot_nf[:,beg:end]

    #create colormaps
    cdict1 = {'red':[(0,1,1),
                    (0.52,1,1),
                    (0.52,0,0),
                    (1,0,0)],
             'green':[(0,1,1),
                     (0.52,1,1),
                     (0.52,0,0),
                     (1,0,0)],
             'blue':[(0,1,1),
                    (0.52,1,1),
                    (0.52,0,0),
                    (1,0,0)]}
    my_cmap1 = mpl.colors.LinearSegmentedColormap('my_colormap1',cdict1,256)
    cdict2 = {'red':[(0,1,1),
                    (0.8,1,1),
                    (0.8,0,0),
                    (1,0,0)],
             'green':[(0,1,1),
                    (0.8,1,1),
                    (0.8,0,0),
                    (1,0,0)],
             'blue':[(0,1,1),
                    (0.8,1,1),
                    (0.8,0,0),
                    (1,0,0)]}
    my_cmap2 = mpl.colors.LinearSegmentedColormap('my_colormap2',cdict2,256)
    
    fig = pb.figure(figsize=(30,10))
    ax1 = fig.add_subplot(611)
    ax2 = fig.add_subplot(612)
    ax3 = fig.add_subplot(613)
    ax4 = fig.add_subplot(614)
    ax5 = fig.add_subplot(615)
    ax6 = fig.add_subplot(616)

    img1 = ax1.imshow(x1_plot, origin='lower', aspect='auto', interpolation="nearest", cmap=my_cmap1, vmin=x1min, vmax=x1max)
    ax2.plot(range(x1_plot.shape[1]), x1_plot[20,:], label='x1')
    ax2.plot(range(x1_plot.shape[1]), x1_plot_nf.mean(axis=0), label='x1')
    ax2.legend(prop={'size':10})
    ax2.set_xlim(0,x1_plot.shape[1])

    img2 = ax3.imshow(x2_plot, origin='lower', aspect='auto', interpolation="nearest", cmap=my_cmap2, vmin=x2min, vmax=x2max)
    ax4.plot(range(x2_plot.shape[1]), x2_plot[20,:], label='x2')
    ax4.plot(range(x2_plot.shape[1]), x2_plot_nf.mean(axis=0), label='x2')
    ax4.legend(prop={'size':10})
    ax4.set_xlim(0,x2_plot.shape[1])
    
    ax5.plot(range(x2_plot.shape[1]), -(0.8*x1_plot.mean(axis=0) + 0.2*x2_plot.mean(axis=0)), label='x1bar + x2bar')
    ax5.legend(prop={'size':10})
    ax5.set_xlim(0,x2_plot.shape[1])
    
    ax6.plot(range(z_plot.shape[1]), z_plot.mean(axis=0), 'k', linewidth=2)
    
    pb.show()
    fig.savefig("spontaneous_seizure_sept10_CpES_100_noise2_40_x0_30_gxixj_20_gxixi_10_r1_cut60s_threshold_HD.svg", dpi=200)
    


def plotNbSpikes_Cp():
    gxixis = [0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2] #values of intrapopulation couplings
    nb_pts = len(gxixis) #nb of samples in the x axis (intrapop cp)
    values = pb.zeros(nb_pts)
    for idx in range(nb_pts):
        #load data
        radical='_40_HMR_40_ML_CpES1_80_CpES2_80_x0_30_noise3_0_noise2_20_noise1_40_gx1x2_10_gx2x1_10_gx1x1_%i_gx2x2_%i_r200_20s' %(gxixis[idx]*100,gxixis[idx]*100)
        x1_plot = np.load('./traces_4oct/x1plot'+radical+'.npy')
        x2_plot = np.load('./traces_4oct/x2plot'+radical+'.npy')
        z_plot =  np.load('./traces_4oct/zplot'+radical+'.npy')
        #downsampling
        x1_plot = x1_plot[:,0:-1:10]
        x2_plot = x2_plot[:,0:-1:10]
        z_plot = z_plot[:,0:-1:10]
        #keep non filtered
        #x1_plot_nf = x1_plot.copy() 
        #x2_plot_nf = x2_plot.copy()
        
        beg=0000
        end=20000
        fs=0.1 # fs: sampling freq
        #resize sample according above interval
        x1_plot = x1_plot[:,beg:end]
        x2_plot = x2_plot[:,beg:end]
        y1_plot = y1_plot[:,beg:end]
        y2_plot = y2_plot[:,beg:end]
        z_plot = z_plot[:,beg:end]
        t_plot = t_plot[beg:end]
        
        x1_mean = x1_plot.mean(axis=0)
        x2_mean = x2_plot.mean(axis=0)
        IS = pb.zeros(len(x1_mean)) #binary timeserie to store interictal spikes
        #algo to keep IS only (pop1 silent (<-1), pop2 active (>-1))
        for i in range(len(x1_mean)-1):
            if (x1_mean[i] < -1 and x2_mean[i] > -1):
                IS[i]=1
            else:
                IS[i]=0
        #average over time to have a metric value for each timeserie
        avg = IS.mean()
        values[idx]=avg
        #display
        fig = pb.figure()
        ax1 = fig.add_subplot(311)
        ax2 = fig.add_subplot(312)
        ax3 = fig.add_subplot(313)
        
        ax1.plot(range(len(x1_mean)), x1_mean)
        ax2.plot(range(len(x2_mean)), x2_mean)
        ax3.plot(range(len(x2_mean)), IS)
        pb.show()
    print(values)
    pb.plot(gxixis,values)
    pb.show()
            
def plotLengthSeizure_Cp():
    gxixjs = [0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2] #values of intrapopulation couplings
    nb_pts = len(gxixjs) #nb of samples in the x axis (intrapop cp)
    values = pb.zeros(nb_pts)
    for idx in range(nb_pts):
        #load data
        radical='_40_HMR_40_ML_CpES1_80_CpES2_80_x0_30_noise3_0_noise2_20_noise1_40_gx1x2_%i_gx2x1_%i_gx1x1_10_gx2x2_10_r200_20s' %(gxixjs[idx]*100,gxixjs[idx]*100)
        x1_plot = np.load('./traces_4oct/x1plot'+radical+'.npy')
        x2_plot = np.load('./traces_4oct/x2plot'+radical+'.npy')
        z_plot =  np.load('./traces_4oct/zplot'+radical+'.npy')
        #downsampling
        x1_plot =  x1_plot[:,0:-1:10]
        x2_plot = x2_plot[:,0:-1:10]
        z_plot = z_plot[:,0:-1:10]
        #keep non filtered
        x1_plot_nf = x1_plot.copy() 
        x2_plot_nf = x2_plot.copy()
        
        beg=0000
        end=20000
        fs=0.1 # fs: sampling freq
        #resize sample according above interval
        x1_plot = x1_plot[:,beg:end]
        x2_plot = x2_plot[:,beg:end]
        y1_plot = y1_plot[:,beg:end]
        y2_plot = y2_plot[:,beg:end]
        z_plot = z_plot[:,beg:end]
        t_plot = t_plot[beg:end]
        
        x1_mean = x1_plot.mean(axis=0)
        x2_mean = x2_plot.mean(axis=0)
        S = pb.zeros(len(x1_mean)) #binary timeserie to store seizure
        #algo to keep IS only (pop1 silent (<-1), pop2 active (>-1))
        for i in range(len(x1_mean)-1):
            if (x1_mean[i] > -1):
                S[i]=1
            else:
                S[i]=0
        #average over time to have a metric value for each timeserie
        avg = S.mean()
        values[idx]=avg
        #display
        fig = pb.figure()
        ax1 = fig.add_subplot(311)
        ax2 = fig.add_subplot(312)
        ax3 = fig.add_subplot(313)
        
        ax1.plot(range(len(x1_mean)), x1_mean)
        ax2.plot(range(len(x2_mean)), x2_mean)
        #ax3.plot(range(len(x2_mean)), IS)
        pb.show()

    print(values) 
    pb.plot(gxixjs,values)
    pb.show()

def plotLengthSeizure_noise():
    noise2s = [1, 1.5, 2, 2.5, 3, 3.5, 4] #values of intrapopulation couplings
    nb_pts = len(noise2s) #nb of samples in the x axis (intrapop cp)
    values = pb.zeros(nb_pts)
    for idx in range(nb_pts):
        #load data
        radical='_40_HMR_40_ML_CpES1_60_CpES2_60_x0_30_noise3_0_noise2_%i_noise1_%i_gx1x2_20_gx2x1_20_gx1x1_20_gx2x2_20_r200_20s' %(noise2s[idx]*10,noise2s[idx]*20)
        x1_plot = np.load('./traces_nov04/x1plot'+radical+'.npy')
        x2_plot = np.load('./traces_nov04/x2plot'+radical+'.npy')
        z_plot =  np.load('./traces_nov04/zplot'+radical+'.npy')
        #downsampling
        x1_plot =  x1_plot[:,0:-1:10]
        x2_plot = x2_plot[:,0:-1:10]
        z_plot = z_plot[:,0:-1:10]
        #keep non filtered
        x1_plot_nf = x1_plot.copy() 
        x2_plot_nf = x2_plot.copy()
        
        beg=0000
        end=20000
        fs=0.1 # fs: sampling freq
        #resize sample according above interval
        x1_plot = x1_plot[:,beg:end]
        x2_plot = x2_plot[:,beg:end]
        #y1_plot = y1_plot[:,beg:end]
        #y2_plot = y2_plot[:,beg:end]
        z_plot = z_plot[:,beg:end]
        #t_plot = t_plot[beg:end]
        
        x1_mean = x1_plot.mean(axis=0)
        x2_mean = x2_plot.mean(axis=0)
        S = pb.zeros(len(x1_mean)) #binary timeserie to store seizure
        #algo to keep IS only (pop1 silent (<-1), pop2 active (>-1))
        for i in range(len(x1_mean)-1):
            if (x1_mean[i] > -1):
                S[i]=1
            else:
                S[i]=0
        #also to count the number of seizures occurring in the timeserie
        nb_seizures=0
        for i in range(len(S)-1):
            if (S[i]==0 and S[i+1]==1):
                nb_seizures = nb_seizures + 1
        #average over time to have a metric value for each timeserie
        result = sum(S)/nb_seizures
        values[idx]=result
        """
        #display
        fig = pb.figure()
        ax1 = fig.add_subplot(311)
        ax2 = fig.add_subplot(312)
        ax3 = fig.add_subplot(313)
        
        ax1.plot(range(len(x1_mean)), x1_mean)
        ax2.plot(range(len(x2_mean)), x2_mean)
        ax3.plot(range(len(x2_mean)), IS)
        pb.show()
        """
    print(values)    
    pb.plot(noise2s,values)
    pb.show()


if __name__=="__main__":
    synchPlotPop2()
    #plotNbSpikes_Cp()
    #plotLengthSeizure_Cp()
    #plotLengthSeizure_noise()