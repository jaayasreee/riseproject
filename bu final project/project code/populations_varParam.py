#import matplotlib
#matplotlib.use('Agg')
import numpy as np
from random import uniform
import pylab as pb

import pop1n
import pop2n 

import KOP_v5 as KOP
#import drawGraph2
from spectrumAnalysis import bandPassFilter

import argparse 

parser = argparse.ArgumentParser(description='Launch Epileptor derived population equations - ex: python populations_args --t_stop 5000 --CpES 0.1 --CpCS 0.9')
parser.add_argument('--t_stop', action='store', dest='t_stop', default=3000, help='time of simulation')
parser.add_argument('--CpES1', action='store', dest='CpES1', default=0.8, help='Coupling % of electrical synapses (i.e. gap junctions) [0-1] within population 1 neurons')
parser.add_argument('--CpES2', action='store', dest='CpES2', default=0.8, help='Coupling % of electrical synapses (i.e. gap junctions) [0-1] within population 2 neurons')
parser.add_argument('--CpCS', action='store', dest='CpCS', default=1., help='Coupling % of chemical synapses [0-1]')
parser.add_argument('--m', action='store', dest='m', default=0.8, help='m parameter in the epileptor equations')
parser.add_argument('--x0', action='store', dest='x0', default=-2., help='x0 in the epileptor equations')
parser.add_argument('--r', action='store', dest='r', default=0.0001, help='r parameter in population 1 neurons equations')
parser.add_argument('--s', action='store', dest='s', default=8., help='s parameter in population 1 neurons equations, reflecting electrical impact rate on slow variable')
parser.add_argument('--nbn1', action='store', dest='nbn1', default=10, help='number of neurons in population 1')
parser.add_argument('--nbn2', action='store', dest='nbn2', default=10, help='number of neurons in population 2')
parser.add_argument('--g_x1x1', action='store', dest='g_x1x1', default=0.2, help='collateral synaptic (maximum) conductance between neurons from population 1')
parser.add_argument('--g_x2x2', action='store', dest='g_x2x2', default=0.2, help='collateral synaptic (maximum) conductance between neurons from population 2')
parser.add_argument('--g_x1x2', action='store', dest='g_x1x2', default=0.2, help='synaptic (maximum) conductance between neurons from population 1 to population 2')
parser.add_argument('--g_x2x1', action='store', dest='g_x2x1', default=0.2, help='fast synaptic (maximum) conductance between neurons from population 2 to population')
parser.add_argument('--g_x2x1_slow', action='store', dest='g_x2x1_slow', default=0., help='slow synaptic (maximum) conductance between neurons from population 2 to population 1')
parser.add_argument('--g_x2x2_slow', action='store', dest='g_x2x2_slow', default=0., help='slow synaptic (maximum) conductance between neurons from population 2 to population 1')
parser.add_argument('--I2', action='store', dest='I2', default=0.8, help='baseline input current in population 2 neurons')
parser.add_argument('--I1', action='store', dest='I1', default=3.1, help='baseline input current in population 1 neurons')
parser.add_argument('--c2', action='store', dest='c2', default=0.3, help='scaling factor of z injection in pop2 neurons')
parser.add_argument('--noise1', action='store', dest='noise1', default=0.5, help='noise amplitude that is introduced in population 1 neurons at each time step')
parser.add_argument('--noiseRatio', action='store', dest='noiseRatio', default=5, help='noise factor for population 1 as compared to population 2')
parser.add_argument('--noise2', action='store', dest='noise2', default=0.3, help='noise amplitude that is introduced in population 2 neurons at each time step')
parser.add_argument('--noise3', action='store', dest='noise3', default=0.0, help='noise amplitude that is introduced in slow variable at each time step')

args = parser.parse_args()

#extract parameters
t_stop = int(args.t_stop)
CpES1 = float(args.CpES1)
CpES2 = float(args.CpES2)
CpCS = float(args.CpCS)
m = float(args.m)
r = float(args.r)
s = float(args.s)
x0 = float(args.x0)
g_x1x1 = float(args.g_x1x1)
g_x2x2 = float(args.g_x2x2)
g_x1x2 = float(args.g_x1x2)
g_x2x1 = float(args.g_x2x1)
g_x2x1_slow = float(args.g_x2x1_slow)
g_x2x2_slow = float(args.g_x2x2_slow)
c2 = float(args.c2)
I2 = float(args.I2)
I1 = float(args.I1)
nbn1= int(args.nbn1) #nb of neurons in pop 1
nbn2= int(args.nbn2) #nb of neurons in pop 2
#noise=float(args.noise1)
noiseRatio=float(args.noiseRatio)
noise2=float(args.noise2)
noise3=float(args.noise3)
noise1=noise2*20
dt=0.05 # time step for simulation
fs=1000 # sampling frequency for storage = 1 kHz
n=2 # number of time steps needed to have a sampling freq of 1kHz (for storage and plotting)

"""#set up parameters that change over time
CpES_variable=pb.zeros(t_stop/dt)
CpES_variable[0:int((t_stop/3)/dt)] = 0.1
CpES_variable[int((t_stop/3+1)/dt):int((2*t_stop/3)/dt)] = pb.linspace(0.1, 1, int(t_stop/3/dt))
CpES_variable[int((2*t_stop/3+1)/dt):t_stop/dt] = pb.linspace(1, 0.8, int((t_stop/3)/dt))
x0_variable=pb.zeros(t_stop/dt)
x0_variable[0:(t_stop/4)/dt] = pb.linspace(-4.5, -2, t_stop/4/dt)
x0_variable[(t_stop/4)/dt:(3*t_stop/4)/dt] = -2
x0_variable[(3*t_stop/4)/dt:t_stop/dt] = pb.linspace(-2, -4.5, t_stop/4/dt)
"""

###########################
# Allocate memory
###########################

#arrays for plotting (store values at each ms in order to have 1kHz "sampling freq.")    
x1_plot=pb.zeros((nbn1,t_stop/dt/n))
x2_plot=pb.zeros((nbn2,t_stop/dt/n))
y1_plot=pb.zeros((nbn1,t_stop/dt/n))
y2_plot=pb.zeros((nbn2,t_stop/dt/n))
x1bar_plot=pb.zeros(t_stop/dt/n)
x2bar_plot=pb.zeros(t_stop/dt/n)
y1bar_plot=pb.zeros(t_stop/dt/n)
y2bar_plot=pb.zeros(t_stop/dt/n)
z_plot=pb.zeros((nbn1,t_stop/dt/n))
zbar_plot=pb.zeros(t_stop/dt/n)
KOP_plot = pb.zeros((2,t_stop/dt/n)) #Kuramoto Order Parameter /// line0: k  line1:psi

#arrays to store means between samples
x1_nsamples=pb.zeros((nbn1,n))
x2_nsamples=pb.zeros((nbn2,n))
y1_nsamples=pb.zeros((nbn1,n))
y2_nsamples=pb.zeros((nbn2,n))
z_nsamples=pb.zeros((nbn1,n))
x1bar_nsamples=pb.zeros(n)
x2bar_nsamples=pb.zeros(n)
y1bar_nsamples=pb.zeros(n)
y2bar_nsamples=pb.zeros(n)
zbar_nsamples=pb.zeros(n)
KOP_nsamples=pb.zeros((2,n))

# arrays to store values of each neurons at one time step
#x1_samples=pb.zeros(nbn1)
#x2_samples=pb.zeros(nbn2)
#y1_samples=pb.zeros(nbn1)
#y2_samples=pb.zeros(nbn2)
#z_samples=pb.zeros(nbn1)

##################################################
# Set Initial Conditions and Instanciate objects
##################################################

# set random initial conditions
x1_init = []; y1_init= []; x2_init = []; y2_init= []; z_init=[];
pop1=[]; pop2=[]; r1=[]; r2=[];  
for i in range(nbn1):
    #initial conditions pop1 neurons
    x1_init.append(uniform(-1.,1.5))
    y1_init.append(uniform(-5.,0.))
    z_init.append(uniform(3.,3.))
    
    #instanciation of pop1 neurons
    pop1.append(pop1n.pop1n(m=m, x0=x0, CpES=CpES1, CpCS=CpCS, g_x1x1=g_x1x1, g_x2x1=g_x2x1, g_x2x1_slow=g_x2x1_slow, I1=I1, r=r, s=s, noise=noise1, noise3=noise3))
    pop1[i].x1=x1_init[i]
    pop1[i].y1=y1_init[i]
    pop1[i].z=z_init[i]

for j in range(nbn2):
    #initial conditions pop2 neurons
    x2_init.append(uniform(-1.25,1.))
    y2_init.append(uniform(0.,1.))
    
    #instanciation of pop2 neurons   
    pop2.append(pop2n.pop2n(CpES=CpES2, CpCS=CpCS, g_x2x2=g_x2x2, g_x1x2=g_x1x2, I2=I2, g_x2x2_slow=g_x2x2_slow, c2=c2, noise=noise2))
    pop2[j].x2=x2_init[j]
    pop2[j].y2=y2_init[j]
    
#connections between neurons
for i in range(nbn1):
    pop1[i].connect_syn_pop2n(pop2[:])
    pop1[i].connect_gap(pop1[:])
for j in range(nbn2):
    pop2[j].connect_syn_pop1n(pop1[:])
    pop2[j].connect_syn_pop2n(pop2[:])
    pop2[j].connect_gap(pop2[:])

x1bar=pb.average(x1_init)
x2bar=pb.average(x2_init)
zbar=pb.average(z_init)

#get pop2 fixed point coordinates. (pop1 need to be computed at each timestep)


############################
# Integration loop
############################
count_samples=0  
for ti in pb.arange(t_stop/dt):
    for i in range(nbn1):
        #pop1[i].x0=x0_variable[ti]
        #pop1[i].CpES=CpES_variable[ti]
        x1_nsamples[i,count_samples],y1_nsamples[i,count_samples],z_nsamples[i,count_samples]=pop1[i].euler(dt, 0, x1bar,x2bar,zbar, ti) # 4th parameter: x2bar 
    for j in range(nbn2):
        #pop2[j].CpES=CpES_variable[ti]
        x2_nsamples[j,count_samples],y2_nsamples[j,count_samples]=pop2[j].euler(dt, 0, x1bar,x2bar,zbar, ti) # 5th parameter: zbar
    #update means
    x1bar=pb.average(x1_nsamples[:,count_samples])
    x2bar=pb.average(x2_nsamples[:,count_samples])
    zbar=pb.average(z_nsamples[:,count_samples])
    #zbar=zbar + r*(4*(x1bar - x0) - zbar)*dt    

    #store means for plotting
    x1bar_nsamples[count_samples]=x1bar
    x2bar_nsamples[count_samples]=x2bar
    y1bar_nsamples[count_samples]=pb.average(y1_nsamples[:,count_samples])
    y2bar_nsamples[count_samples]=pb.average(y2_nsamples[:,count_samples])
    zbar_nsamples[count_samples]=zbar

    # calculate KURAMOTO ORDER PARAMETER
    #x_fixedPt,y_fixedPt = KOP.get_rightFixedPt(zbar) #get "right" fixed point coordinates
    #k,psi = KOP.get_KOP(x1_nsamples[:,count_samples]-x_fixedPt, y1_nsamples[:,count_samples]-y_fixedPt, nbn1)  #get kuramoto order parameter
    #KOP_nsamples[:,count_samples] = k,psi  #store KOP
        
    count_samples+=1  
    if count_samples==n: #each n samples, store values (correspond to sampling freq.)
        x1_plot[:,ti/n]=x1_nsamples.mean(axis=1)
        x2_plot[:,ti/n]=x2_nsamples.mean(axis=1)
        y1_plot[:,ti/n]=y1_nsamples.mean(axis=1)
        y2_plot[:,ti/n]=y2_nsamples.mean(axis=1)
        z_plot[:,ti/n]=z_nsamples.mean(axis=1)
        x1bar_plot[ti/n]=x1bar_nsamples.mean()
        x2bar_plot[ti/n]=x2bar_nsamples.mean()
        y1bar_plot[ti/n]=y1bar_nsamples.mean()
        y2bar_plot[ti/n]=y2bar_nsamples.mean()
        zbar_plot[ti/n]=zbar_nsamples.mean()
        KOP_plot[:,ti/n]=np.mean(KOP_nsamples, axis=1)
        count_samples=0 # reset counter 

    
t_plot=pb.arange(0,t_stop/dt/n)
t_plot_dt=pb.arange(0,t_stop,dt)

########################################
# save & plot
########################################

#save data
radical="_%i_HMR_%i_ML_CpES1_%i_CpES2_%i_x0_%i_noise1_%i_noise2_%i_noise3_%i_%is" %(nbn1, nbn2, int(CpES1*100), int(CpES2*100), int(-x0*10), int(noise1*10), int(noise2*100), int(noise3*100),int(t_stop/1000))
np.save('x1plot'+radical+'.npy', x1_plot)
np.save('x2plot'+radical+'.npy', x2_plot)
np.save('y1plot'+radical+'.npy', y1_plot)
np.save('y2plot'+radical+'.npy', y2_plot)
np.save('zplot'+radical+'.npy', z_plot)
np.save('tplot'+radical+'.npy', t_plot)
#np.save('KOP'+radical+'.npy', KOP_plot)


fig = pb.figure(figsize=(30,25)) #pb.figure(figsize=(20,10))
#fig.hold(True)
fig.suptitle('%i pop1 neurons; %i pop2 neurons; I1=%.2f; I2=%.2f; x0=%.2f; r=%.5f, s=%i; CpES1=%.2f; CpES2=%.2f; g_x1x1=%.2f; g_x1x2=%.2f; g_x2x2=%.2f; g_x2x1=%.2f; c2=%.3f; noise1=%.2f; noise2=%.2f; noise3=%.2f' %(nbn1, nbn2,I1,I2, x0, r, s, CpES1, CpES2, g_x1x1, g_x1x2, g_x2x2, g_x2x1, c2, noise1, noise2, noise3)) 
    
ax1=fig.add_subplot(611); ax1.hold(True)
ax2=fig.add_subplot(613); ax2.hold(True)
ax3=fig.add_subplot(615); ax3.hold(True)
ax4=fig.add_subplot(616); ax4.hold(True)
ax5=fig.add_subplot(612); ax5.hold(True)
ax6=fig.add_subplot(614); ax5.hold(True)

for i in range(nbn1):
        #time series pop1
        ax5.plot(t_plot, x1_plot[i,:], label=None)
        #drawGraph2.scatter_plot(t_plot, x1_plot[i], 0.0, i, ax1_1)
        #time series z
        ax4.plot(t_plot, z_plot[i,:], label=None)
for j in range(nbn2): 
        #time series pop2
        ax6.plot(t_plot, x2_plot[j,:], label=None)
        
#draw mean values in time series
ax1.imshow(x1_plot, interpolation='nearest', aspect='auto')
ax5.plot(t_plot, x1bar_plot, 'black', linewidth=1.5, label="x1bar")
ax5.set_xlim(right=t_plot.size)
#ax1.legend(prop={'size':10})
ax1.set_title("Population 1 neurons")
#ax1_1.set_xlim(left=0)
#ax1_1.set_ylim((-0.5, nbn1-0.5))
ax2.imshow(x2_plot, interpolation='nearest', aspect='auto')
ax6.plot(t_plot, x2bar_plot, 'black', linewidth=1.5, label="x2bar")
ax6.set_xlim(right=t_plot.size)
#ax2.legend(prop={'size':10})
ax2.set_title("Population 2 neurons")
#plot energy at last
#ax3.plot(t_plot, highPassFilter(x2bar_plot + x1bar_plot), 'purple', label='x2bar + x1bar filtered')
ax3.plot(t_plot, 0.2*x2bar_plot + 0.8*x1bar_plot, 'purple', label='x2bar + x1bar')
ax3.set_xlim(0,t_plot.size)
ax3.legend(prop={'size':10})
ax3.set_title("Whole system")
#draw z mean 
ax4.plot(t_plot, zbar_plot, 'black', linewidth=1.5, label="zbar")
ax4.set_xlim(0,t_plot.size)
#ax4.legend(prop={'size':10})
ax4.set_title("Slow variable")
#draw x0 and CpES
#ax4.plot(t_plot_dt, CpES_variable, 'blue' ,linewidth=2.5, label="CpES")
#ax4.plot(t_plot_dt, x0_variable, 'red',linewidth=2.5, label="x0")
#ax5.legend(prop={'size':10})
#set figure name & save & show
figname="epilepton"+radical+".png"
fig.savefig(figname, dpi=200)
#figTool.zoom(10000, 12000, x1_plot, y1_plot, z_plot, x2_plot, y2_plot, t_plot, KOP_plot, radical)
pb.show()