'''
Created on 21 aout 2012

@author: Squirel
'''
import numpy as np
import pylab as pb
import os
import cmath
#t_plot=pb.arange(0,3000,0.05)
import pop2n
import gc
#access memory usage
def memory():
    import os
    from wmi import WMI
    w = WMI('.')
    result = w.query("SELECT WorkingSet FROM Win32_PerfRawData_PerfProc_Process WHERE IDProcess=%d" % os.getpid())
    return int(result[0].WorkingSet)

################
## FOR POP1 : ##
################
def get_mean_nullclines(x1span,z):
    a=1.; b=3.; c=1.; d=5.; m=0.6; I1=3.1; 
    x1nullcline=[]
    y1nullcline=[]
    for x in x1span:
        if x<0: x1nullcline.append(a*x**3-b*x**2+z-I1)
        else: x1nullcline.append(-x*(m+0.6*(z-4)**2)+z-I1)
        y1nullcline.append(c-d*x**2)
    return x1nullcline, y1nullcline


def get_rightFixedPt(z):
    a=1.; b=3.; c=1.; d=5.; m=0.6; I1=3.1;
    # dx^2-(m+0.6(z-4)^2)x+z-I1-c=0  --> AX^2+BX+C=0 
    A=d; B=-(m+0.6*(z-4)**2); C=z-I1-c
    delta = B**2-4*A*C
    if(delta==0): 
        x1=-B/(2*A)
    else:
        if delta>0 : 
            x1=(-B+pb.sqrt(delta))/(2*A)
            if x1<0: x1=(-B-pb.sqrt(delta))/(2*A)
        else:
            return 0,0 #if delta<0 --> no root 
    y1=z-I1-x1*(m+0.6*(z-4)**2)
    return x1,y1

#Kuramoto Order Parameter  r*exp(i*psi) = 1/N*sum(exp(i*theta))    CAUTION: x1plot & y1plot must be re-referenced (i.e. rotating around (0,0))  before call
def get_KOP(x1plot,y1plot,index, nbn):
    sum=0
    for n in range(nbn):
        x=x1plot[n]; y=y1plot[n];
        #get polar coordinate from (x,y) euclidean coordinates
        r,theta=cmath.polar(complex(x,y))
        sum+=pb.exp(complex(0,theta))
        
    return cmath.polar(sum/nbn)

###############
## FOR POP2: ##
############### 


###############
## make a phase plan movie from existing traces, in euclidean space, with nullclines and fixed points that serve as origin for the synchronization measure.
###############
def makeMoviePop1and2():
    
    beg=520000; end=570000; fs=10 # fs: sampling freq in kHz !! not downsampled
    #load data an resize sample according interval above
    radical='_40_HMR_40_ML_CpES1_100_CpES2_100_x0_30_noise3_0_noise2_40_noise1_80_gx1x2_20_gx2x1_20_gx1x1_10_gx2x2_10_r1_100s'
    print memory()
    x1_p = np.load('./traces_sept10/x1plot'+radical+'.npy', mmap_mode='r+')
    print memory()
    x1_plot = x1_p[:,beg:end].copy()
    del(x1_p)
    gc.collect()
    print memory()
    x2_p = np.load('./traces_sept10/x2plot'+radical+'.npy', mmap_mode='r+')
    x2_plot = x2_p[:,beg:end] #x2is stored reduced
    del(x2_p)
    gc.collect()
    print memory()
    y1_p = np.load('./traces_sept10/y1plot'+radical+'.npy', mmap_mode='r+')
    y1_plot = y1_p[:,beg:end]
    del(y1_p)
    gc.collect()
    print memory()
    y2_p = np.load('./traces_sept10/y2plot'+radical+'.npy', mmap_mode='r+')
    y2_plot = y2_p[:,beg:end]
    del(y2_p)
    gc.collect()
    print memory()
    z_p = np.load('./traces_sept10/zplot'+radical+'.npy', mmap_mode='r+')
    z_plot = z_p[:,beg:end]
    del(z_p)
    gc.collect()
    print memory()
    t_p = np.load('./traces_sept10/tplot'+radical+'.npy', mmap_mode='r+')
    t_plot = t_p[beg:end]
    del(t_p)
    gc.collect()
    print memory()
    
    #get minimums and maximums before resizing
    x1max=x1_plot.max()
    x1min=x1_plot.min()
    y1max=y1_plot.max()
    y1min=y1_plot.min()
    x2max=x2_plot.max()
    x2min=x2_plot.min()
    y2max=y2_plot.max()
    y2min=y2_plot.min()
    
    #compute mean time series over population
    nbn=x1_plot.shape[0]
    x1bar_plot = pb.zeros(x1_plot.shape[1])
    y1bar_plot = pb.zeros(y1_plot.shape[1]) #normally x1_plot y1_plot have same shape
    x2bar_plot = pb.zeros(x2_plot.shape[1])
    y2bar_plot = pb.zeros(y2_plot.shape[1])
    zbar_plot = pb.zeros(z_plot.shape[1])
    #compute the means:
    for i in range(x1bar_plot.size):
        x1bar_plot[i]=pb.mean(x1_plot[:,i])
        y1bar_plot[i]=pb.mean(y1_plot[:,i])
        x2bar_plot[i]=pb.mean(x2_plot[:,i])
        y2bar_plot[i]=pb.mean(y2_plot[:,i])
        zbar_plot[i]=pb.mean(z_plot[:,i])
    
    x1span=pb.arange(x1min-0.1, x1max+0.1, 0.01)
    x2span=pb.arange(x2min-0.1, x2max+0.1, 0.01)
    #y2span=pb.arange(y2min-0.1, y2max+0.1, 0.01)
    #compute nullclines
    x1nc, y1nc = get_mean_nullclines(x1span, zbar_plot[i])
    x2nc, y2nc = pop2n.pop2n().get_nullclines(x2span)
    #compute pop2 fixed pt (because static, pop1 fxed pt is dynamic so is computed at each time step 
    x2_fixedPt = x2span[np.abs(x2span-0).argmin()]
    y2_fixedPt = y2nc[np.where(x2span==x2span[np.abs(x2span-0).argmin()])]
    
    #scaling factor to compute the normalized phase. -1<x<1.  (*100 in order that most points out of the unit circle)
    x1_scaling = 2/(x1max-x1min)*100
    y1_scaling = 2/(y1max-y1min)*100
    x2_scaling = 2/(x2max-x2min)*100
    y2_scaling = 2/(y2max-y2min)*100
        
    movie_plot_indexes = pb.arange(0,t_plot.size,5)
    
    fig = pb.figure(figsize=(15,12))
    ax1=pb.subplot2grid((3,2),(0,0)); ax1.hold(True)
    ax2=pb.subplot2grid((3,2),(0,1)); ax2.hold(True)
    ax3=pb.subplot2grid((3,2),(1,0), polar=True); ax3.hold(True)
    ax4=pb.subplot2grid((3,2),(0,1), polar=True); ax4.hold(True)
    ax5=pb.subplot2grid((3,2),(2,0), colspan=2); ax5.hold(True)
    ax5.plot(t_plot,  -(0.8*x1_plot.mean(axis=0) + 0.2*x2_plot.mean(axis=0)), 'k')
    pt = ax5.plot(t_plot[0], -(0.8*x1_plot[:,0].mean(axis=0) + 0.2*x2_plot[:,0].mean(axis=0)), 'ob')
    count=1;
    for i in movie_plot_indexes:
        #fig.suptitle("t=%i" %(beg+i))
        ax1.cla()
        ax2.cla()
        ax3.cla()
        ax4.cla()
        pt.pop(0).remove() #remove dot from timeserie at each frame of the video
        #plot nullclines
        ax1.plot(x1span,x1nc,x1span,y1nc)
        ax2.plot(x2span,x2nc, x2span,y2nc)
        
        #get "right" fixed point
        x1_fixedPt,y1_fixedPt = get_rightFixedPt(zbar_plot[i])
        #x1_fixedPt,y1_fixedPt = 0 , 0 # if we consider the center of the limit circle on the phase plane to be the origin
        #plot "right" fixed point
        ax1.plot(x1_fixedPt,y1_fixedPt, 'or')
        ax2.plot(x2_fixedPt,y2_fixedPt, 'or')
        #        #get kuramoto order parameter
        r,psi = get_KOP(x1_scaling*(x1_plot[:,i]-x1_fixedPt), y1_scaling*(y1_plot[:,i]-y1_fixedPt), i, nbn)
        r2,psi2 = get_KOP(x2_scaling*(x2_plot[:,i]-x2_fixedPt), y2_scaling*(y2_plot[:,i]-y2_fixedPt), i, nbn)
        #plot KOP vector from fixed point
        #ax1.arrow(x1_fixedPt, y1_fixedPt, r*pb.cos(psi)/x1_scaling, r*pb.sin(psi)/y1_scaling)
        #ax2.arrow(x2_fixedPt, y2_fixedPt, r2*pb.cos(psi2)/x2_scaling, r2*pb.sin(psi2)/y2_scaling)
        ax3.arrow(0, 0, psi, r/2) #/ such that the arrow (in fact line) does not reach the points on the circle
        ax4.arrow(0, 0, psi2, r2/2)
                  
        #plot individual pts 
        ax1.plot(x1_plot[:,i], y1_plot[:,i], 'ob')
        ax1.plot(x1bar_plot[i], y1bar_plot[i], 'ok' )
        ax2.plot(x2_plot[:,i], y2_plot[:,i], 'ob')
        ax2.plot(x2bar_plot[i], y2bar_plot[i], 'ok' )
        for n in range(nbn):
            r1,theta1=cmath.polar(complex(x1_scaling*(x1_plot[n,i]-x1_fixedPt),y1_scaling*(y1_plot[n,i]-y1_fixedPt)))
            r2,theta2=cmath.polar(complex(x2_scaling*(x2_plot[n,i]-x2_fixedPt),y2_scaling*(y2_plot[n,i]-y2_fixedPt)))
            ax3.plot(theta1,1,'o', markersize=10)
            ax4.plot(theta2,1,'o', markersize=10)
        pt = ax5.plot(t_plot[i], -(0.8*x1_plot[:,i].mean(axis=0) + 0.2*x2_plot[:,i].mean(axis=0)), 'ob', markersize=10)
        #set graph boundaries (and labels if needed)
        ax1.set_ylim(y1min-0.1,y1max+0.1)
        ax1.set_xlim(x1min-0.1,x1max+0.1)
        ax2.set_ylim(y2min-0.1,y2max+0.1)
        ax2.set_xlim(x2min-0.1,x2max+0.1)
        ax3.set_yticks([])
        ax3.set_ylim(0,1.05) #1.05 is the r value, it is to set the diameter of the circle
        ax4.set_yticks([])
        ax4.set_ylim(0,1.05) #1.05 is the r value, it is to set the diameter of the circle
        
        ax1.set_title("Population 1")
        ax2.set_title("Population 2")
        ax3.set_title("Phase Population 1 neurons.")
        ax4.set_title("Phase Population 2 neurons.")
        ax5.set_title("Time serie whole system")
        
        #create img from figure frame and save
        fname = '_tmp%06d.png' %count
        print 'Saving frame ', fname
        fig.savefig(fname)
        count += 1
        
    print 'Making movie - may take a while'
    os.system("C:/MPlayer/mencoder.exe mf://_tmp*.png -mf fps=30:type=png -ovc lavc -lavcopts vcodec=mpeg4 -oac copy -o PhaseSpaceCircle_pop1_2_x0_30_noise2_40_sept10.avi")
    os.system("del .\_tmp*.png")


def drawPolarFrame():
    #load data
    radical='_40_HMR_40_ML_CpES1_0_CpES2_100_I2_90_x0_40_noise_100_s8_15s'
    #x1_plot = np.load('./traces_11april/x1plot'+radical+'.npy')
    x2_plot = np.load('./traces_11april/x2plot'+radical+'.npy')
    #y1_plot = np.load('./traces_11april/y1plot'+radical+'.npy')
    y2_plot = np.load('./traces_11april/y2plot'+radical+'.npy')
    #z_plot = np.load('./traces_11april/zplot'+radical+'.npy')
    t_plot = np.load('./traces_11april/tplot'+radical+'.npy')
    
    #x1max=x1_plot.max()
    #x1min=x1_plot.min()
    #y1max=y1_plot.max()
    #y1min=y1_plot.min()
    
    ts=42000 #timestep wanted
    beg=38500; end=42000; dt=0.05
    
    #get the values at the given timestep ts
    #x1_plot = x1_plot[:,ts]
    x2_plot = x2_plot[:,beg:end]*20 #x2is stored reduced
    #y1_plot = y1_plot[:,ts]
    y2_plot = y2_plot[:,beg:end]
    #z_plot = z_plot[:,ts]
    t_plot = t_plot[beg:end]
    count=1
    for ti in pb.arange(0,t_plot.size,2):
        fig=pb.figure()
        ax1=fig.add_subplot(111, polar=True)
        ax1.hold(True)
        ax1.set_yticks([])
        ax1.set_ylim(0,1.05) #1.05 is the r value, it is to set the diameter of the circle
        ax1.set_title("Phase Population 2 neurons. t=%.1fms" %(beg+ti))    
        for i in range(x2_plot.shape[0]):
            r,theta=cmath.polar(complex(x2_plot[i,ti],y2_plot[i,ti]))
            ax1.plot(theta,1, 'o')
        
        #create img from figure frame and save
        fname = '_tmp%06d.png' %count
        print 'Saving frame ', fname
        fig.savefig(fname)
        count += 1
        
    print 'Making movie - may take a while'
    os.system("C:/MPlayer/mencoder.exe mf://_tmp*.png -mf fps=30:type=png -ovc lavc -lavcopts vcodec=mpeg4 -oac copy -o Phase_Circle_pop2.avi")
    os.system("del .\_tmp*.png")
    

#############
## SCRIPT: ##
#############

if __name__=="__main__":
    makeMoviePop1and2()
