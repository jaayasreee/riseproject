import numpy as np
import pylab as pb
#import os
import cmath
#t_plot=pb.arange(0,3000,0.05)
import paramSpace_4D_v5 as paramSpace4D
import pop2n
import pop1n
import sys
from spectrumAnalysis import bandPassFilter

def get_mean_nullclines(x1span,z, m, I1):
    a=1.
    b=3.
    c=1.
    d=5.
    m=m
    I1=I1
    x1nullcline=[]
    y1nullcline=[]
    for x in x1span:
        if x<0: x1nullcline.append(a*x**3-b*x**2+z-I1)
        else: x1nullcline.append(-x*(m+0.6*(z-4)**2)+z-I1)
        y1nullcline.append(c-d*x**2)
    return x1nullcline, y1nullcline


def get_rightFixedPt(z):
    a=1.
    b=3.
    c=1.
    d=5.
    m=0.6
    I1=3.1
    # dx^2-(m+0.6(z-4)^2)x+z-I1-c=0  --> AX^2+BX+C=0 
    A=d
    B=-(m+0.6*(z-4)**2)
    C=z-I1-c
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
def get_KOP(x1s,y1s,nbn):
    sum=0
    for n in range(nbn):
        x=x1s[n]
        y=y1s[n]
        #get polar coordinate from (x,y) euclidean coordinates
        r,theta=cmath.polar(complex(x,y))
        sum+=pb.exp(complex(0,theta))
        
    return cmath.polar(sum/nbn)

#Kuramoto Order Parameter  r*exp(i*psi) = 1/N*sum(exp(i*theta))    CAUTION: x1plot & y1plot must be re-referenced (i.e. rotating around (0,0))  before call
def get_KOP_noRest(x1s,y1s,nbn,baseVal):
    sum=0
    for n in range(nbn):
        #x=x1s[n]; y=y1s[n];
        #get polar coordinate from (x,y) euclidean coordinates
        r,theta=cmath.polar(complex(x1s[n],y1s[n]))   #returns   -pi < theta < pi
        if theta > -pb.pi/2 :
            sum += pb.exp(complex(0,theta))
        else:
            sum += baseVal
        
    return cmath.polar(sum/nbn)


def displayKOP():
    KOP1=np.load("./traces_01mai/KOP1_40_HMR_40_ML_CpES1_100_CpES2_100_x0_25_noiseRatio_160_noise2_5_15s.npy")
    KOP2=np.load("./traces_01mai/KOP2_40_HMR_40_ML_CpES1_100_CpES2_100_x0_25_noiseRatio_160_noise2_5_15s.npy")
    
    print(KOP1[0,:].mean())
    print(KOP2[0,:].mean())
    print(KOP1[0,:].mean() + KOP2[0,:].mean())
    
    fig=pb.figure()
    ax1=fig.add_subplot(311)
    ax1.plot(KOP1[0,:])
    ax2=fig.add_subplot(312)
    ax2.plot(KOP2[0,:])
    ax3=fig.add_subplot(313)
    ax3.plot(KOP1[0,:] + KOP2[0,:])
    
    pb.show()
    
def computeKOP2map():
    out_img_grid_KOP2=[]
    for out_row_index in range(len(paramSpace4D.out_row_val)):
        row_grid_KOP2=[]
        for out_col_index in range(len(paramSpace4D.out_col_val)):
            in_row_length=len(paramSpace4D.in_row_val)
            in_col_length=len(paramSpace4D.in_col_val)
            in_img_grid_KOP2=pb.zeros((in_row_length, in_col_length))
            for in_row_index in range(in_row_length):
                for in_col_index in range(in_col_length):
                    radical="_%i_HMR_%i_ML_CpES_%i_gx1x2_%i_gx2x1_%i_x0_%i_%is" \
                        %(paramSpace4D.nbn1, paramSpace4D.nbn2, int(paramSpace4D.out_col_val[out_col_index]*100), int(paramSpace4D.in_row_val[in_row_index]*100), int(paramSpace4D.in_col_val[in_col_index]*100), int(paramSpace4D.out_row_val[out_row_index]*100), int(paramSpace4D.t_stop/1000))
                    #load data
                    KOP2_plot = np.load('traces_22feb/KOP2'+radical+".npy")
                    in_img_grid_KOP2[in_row_index, in_col_index] = KOP2_plot[0,:].mean()
            row_grid_KOP2.append(in_img_grid_KOP2)
        out_img_grid_KOP2.append(row_grid_KOP2)
    np.save('outImgGridKOP2_v3', out_img_grid_KOP2)
    
    
def computeKOPs(noise2, noiseRatio, x0):
    radical="_40_HMR_40_ML_CpES1_100_CpES2_100_x0_%i_noiseRatio_%i_noise2_%i_15s" %(int(float(x0)*10), int(float(noiseRatio)*10), int(float(noise2)*10))
    #load data
    x1_plot = np.load('traces_01mai/x1plot'+radical+".npy")
    y1_plot = np.load('traces_01mai/y1plot'+radical+".npy")
    z_plot =  np.load('traces_01mai/zplot'+radical+".npy")
    x2_plot = np.load('traces_01mai/x2plot'+radical+".npy")
    y2_plot = np.load('traces_01mai/y2plot'+radical+".npy")
    t_plot = np.load('traces_01mai/tplot'+radical+".npy")
    
    nbn=x1_plot.shape[0] #extract from data the number of neurons used in the simulation
    #get minimums and maximums before resizing
    x1max=x1_plot.max()
    x1min=x1_plot.min()
    y1max=y1_plot.max()
    y1min=y1_plot.min()
    x2max=x2_plot.max()*20 #x2is stored reduced
    x2min=x2_plot.min()*20
    y2max=y2_plot.max()
    y2min=y2_plot.min()
    
    x2span=pb.arange(x2min-0.1, x2max+0.1, 0.01)
    y2span=pb.arange(y2min-0.1, y2max+0.1, 0.01)

    #scaling factor to compute the normalized phase. -1<x<1.  (*100 in order that most points out of the unit circle)
    x1_scaling = 2/int((x1max-x1min))
    y1_scaling = 2/int((y1max-y1min))
    x2_scaling = 2/int((x2max-x2min))
    y2_scaling = 2/int((y2max-y2min))
    
    KOP1_plot=pb.zeros((2,t_plot.size)) #storage of KOP1
    KOP2_plot=pb.zeros((2,t_plot.size)) #storage of KOP2
    #get nullclines and fixed point of pop2  ---!!--- pop1 fixed points needs to be computed each time step
    x2nc, y2nc = pop2n.pop2n().get_nullclines(x2span)
    #get fixed point coordinates
    x2_fixedPt = x2span[np.abs(x2span-0).argmin()]
    y2_fixedPt = y2nc[np.where(x2span==x2span[np.abs(x2span-0).argmin()])]
    x1_fixedPt = -0.5
    y1_fixedPt = -4
        
    #compute KOPs for each timestep
    for i in range(t_plot.size):
        #take z mean and use it to get fxed pt from pop1
        print(i)
        z = z_plot[:,i].mean()
        #x1_fixedPt, y1_fixedPt = get_rightFixedPt(z)
        #get kuramoto order parameters
        r1,psi1 = get_KOP_noRest(x1_scaling*(x1_plot[:,i]-x1_fixedPt), y1_scaling*(y1_plot[:,i]-y1_fixedPt), nbn, 0.35) #base value of pop1 kop with lots of noise is 0.35
        r2,psi2 = get_KOP_noRest(x2_scaling*(x2_plot[:,i]-x2_fixedPt), y2_scaling*(y2_plot[:,i]-y2_fixedPt), nbn, 0.15) #base value of pop2 kop with lots of noise is 0.15
        KOP1_plot[:,i] = r1,psi1 #store KOP1
        KOP2_plot[:,i] = r2,psi2 #store KOP2
    #save KOP2    
    np.save('./traces_01mai/KOP2'+radical+'.npy', KOP2_plot)
    np.save('./traces_01mai/KOP1'+radical+'.npy', KOP1_plot)
    
def computeCoastline(noise2, noiseRatio, x0):
    radical="_40_HMR_40_ML_CpES1_100_CpES2_100_x0_%i_noiseRatio_%i_noise2_%i_15s" %(int(float(x0)*10), int(float(noiseRatio)*10), int(float(noise2)*10))
    #load data
    x1_plot = np.load('traces_01mai/x1plot'+radical+".npy")
    y1_plot = np.load('traces_01mai/y1plot'+radical+".npy")
    z_plot =  np.load('traces_01mai/zplot'+radical+".npy")
    x2_plot = np.load('traces_01mai/x2plot'+radical+".npy")
    y2_plot = np.load('traces_01mai/y2plot'+radical+".npy")
    t_plot = np.load('traces_01mai/tplot'+radical+".npy")
    
    nbn=x1_plot.shape[0] #extract from data the number of neurons used in the simulation
    
    ts=bandPassFilter(x1_plot.mean(axis=0) + x2_plot.mean(axis=0))
    cli = paramSpace4D.get_coastline(ts)
    
    np.save('./traces_01mai/CL'+radical+'.npy',cli)
    
    
if __name__ == "__main__":
    #displayKOP()
    print(sys.argv[1], sys.argv[2], sys.argv[3], computeKOPs(sys.argv[1], sys.argv[2], sys.argv[3]))
    #computeCoastline(sys.argv[1], sys.argv[2], sys.argv[3])
    #computeCoastline(0.5,16,3)