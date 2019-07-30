import matplotlib
import numpy as np
from random import uniform
import pylab as pb
from numpy import arange, mean, exp, zeros, tanh, cosh
from numpy.random import rand

import argparse 

import spectrumAnalysis as sA
import fast_synapse_T 
import slow_synapse_T 
import KOP_v5 as KOP
import pop1n
import pop2n
import paramSpace3D_v6 
import paramSpace_4D_v5
import figTool


parser = argparse.ArgumentParser(description='Launch Epileptor derived population equations - ex: python populations_args --t_stop 5000 --CpES 0.1 --CpCS 0.9')
parser.add_argument('--t_stop', action='store', dest='t_stop', default=1000, help='time of simulation')
parser.add_argument('--CpES1', action='store', dest='CpES1', default=1, help='Coupling % of electrical synapses (i.e. gap junctions) [0-1] within population 1 neurons')
parser.add_argument('--CpES2', action='store', dest='CpES2', default=1, help='Coupling % of electrical synapses (i.e. gap junctions) [0-1] within population 2 neurons')
parser.add_argument('--CpCS', action='store', dest='CpCS', default=1., help='Coupling % of chemical synapses [0-1]')
parser.add_argument('--m', action='store', dest='m', default=0.8, help='m parameter in the epileptor equations')
parser.add_argument('--x0', action='store', dest='x0', default=-3., help='x0 in the epileptor equations')
parser.add_argument('--r', action='store', dest='r', default=0.0002, help='r parameter in population 1 neurons equations')
parser.add_argument('--s', action='store', dest='s', default=8., help='s parameter in population 1 neurons equations, reflecting electrical impact rate on slow variable')
parser.add_argument('--nbn1', action='store', dest='nbn1', default=80, help='number of neurons in population 1')
parser.add_argument('--nbn2', action='store', dest='nbn2', default=20, help='number of neurons in population 2')
parser.add_argument('--g_x1x1', action='store', dest='g_x1x1', default=0.1, help='collateral synaptic (maximum) conductance between neurons from population 1')
parser.add_argument('--g_x2x2', action='store', dest='g_x2x2', default=0.1, help='collateral synaptic (maximum) conductance between neurons from population 2')
parser.add_argument('--g_x1x2', action='store', dest='g_x1x2', default=0.2, help='synaptic (maximum) conductance between neurons from population 1 to population 2')
parser.add_argument('--g_x2x1', action='store', dest='g_x2x1', default=0.2, help='fast synaptic (maximum) conductance between neurons from population 2 to population')
parser.add_argument('--g_x2x1_slow', action='store', dest='g_x2x1_slow', default=0., help='slow synaptic (maximum) conductance between neurons from population 2 to population 1')
parser.add_argument('--g_x2x2_slow', action='store', dest='g_x2x2_slow', default=0., help='slow synaptic (maximum) conductance between neurons from population 2 to population 1')
parser.add_argument('--I2', action='store', dest='I2', default=0.6, help='baseline input current in population 2 neurons')
parser.add_argument('--I1', action='store', dest='I1', default=3.1, help='baseline input current in population 1 neurons')
parser.add_argument('--c2', action='store', dest='c2', default=0.3, help='scaling factor of z injection in pop2 neurons')
parser.add_argument('--noise1', action='store', dest='noise1', default=0.5, help='noise amplitude that is introduced in population 1 neurons at each time step')
parser.add_argument('--noiseRatio', action='store', dest='noiseRatio', default=5, help='noise factor for population 1 as compared to population 2')
parser.add_argument('--noise2', action='store', dest='noise2', default=0.4, help='noise amplitude that is introduced in population 2 neurons at each time step')
parser.add_argument('--noise3', action='store', dest='noise3', default=0.0, help='noise amplitude that is introduced in slow variable at each time step')

args, unknown = parser.parse_known_args()

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
noise=float(args.noise1)
noiseRatio=float(args.noiseRatio)
noise2=float(args.noise2)
noise3=float(args.noise3)
noise1=noise2*20.
dt=0.05 # time step for simulation
fs=1000 # sampling frequency for storage = 1 kHz
n=2 # number of time steps needed to have a sampling freq of 1kHz (for storage and plotting)

n_e=nbn1; n_i=nbn2
#set other parameters for both populations
#---pop1
a=1.
b=3.
c=1.
d=5.
#---pop2
aa=8.
V1=-1.2
V2=18.
V3=12.
V4=17.4
phi=0.067
gCa=4.
gK=8.
gL=2.
E_K=-84.
E_L=-60.
E_Ca=120.
Cm=20.
#---synapses
Vt=2.
Kd=5.
Tmax=1. 
Ee=0.
alpha_e=1.1
beta_e=0.19
Ei=-80.
alpha_i=5.
beta_i=0.18
GsEE=g_x1x1
GsEI=g_x1x2
GsIE=g_x2x1
GsII=g_x2x2 #synaptic coupling
Ce1 = CpES1
Ce2=CpES2

##################################################
# Set Initial Conditions and Instanciate objects
##################################################

# set random initial conditions
x1 = 2.5*pb.rand(nbn1,1).ravel()-1.
y1 = 5.*pb.rand(nbn1,1).ravel()-5.
x2 = 2.*pb.rand(nbn2,1).ravel()-1.
y2 = pb.rand(nbn2,1).ravel()
z = 0.1*pb.rand(nbn1,1).ravel()+3.45
uE=zeros((n_e,1)).ravel()
uI=zeros((n_i,1)).ravel()
#define integration parameter if manual integration (i.e euler) 
t_end=t_stop
dt=0.05
t = arange(0,t_end,dt)
#create storage arrays
x_1s = zeros((nbn1,t_end))
x_2s = zeros((nbn2,t_end))
y_1s = zeros((nbn1,t_end))
n_s = zeros((nbn2,t_end))
z_s = zeros((nbn1,t_end))

#integrate
for t_step in t:
    #equations
    #transmitter release
    Te = Tmax / (1.+exp(-(mean(x1)*50.-Vt)/Kd)) #exc
    Ti = Tmax / (1.+exp(-(mean(x2)*50.-Vt)/Kd)) #inh
    #opening and closing of gates
    duE = (alpha_e*Te*(1.-uE)-beta_e*uE) #exc gates
    duI = (alpha_i*Ti*(1.-uI)-beta_i*uI) #inh gates
    #---synapses exc-exc
    IsynEE = -GsEE*mean(uE)*(x1*50.-Ee) #synaptic current
    #---synapses exc-inh
    IsynEI = -GsEI*mean(uE)*(x2*50.-Ee) #synaptic current
    #---synapses inh-exc
    IsynIE = -GsIE*mean(uI)*(x1*50.-Ei) #synaptic current
    #---synapses inh-inh
    IsynII = -GsII*mean(uI)*(x2*50.-Ei) #synaptic current
    #---pop1 exc
    dx1 = y1 - x1**3 + 3.*x1**2 - z + I1 + Ce1*(mean(x1)-x1) + IsynEE/50. + IsynIE/50. + noise1*(2.*(rand(n_e,1)-0.5)).ravel()
    dy1 = c - 5.*x1**2 - y1
    #---pop2 inh
    V=x2*20.
    m_inf=0.5*(1+tanh((V-V1)/V2))
    n_inf=0.5*(1+tanh((V-V3)/V4))
    tau_n=1./cosh((V-V3)/(2*V4))
    dn = phi*(n_inf-n)/tau_n; #slow 
    I_in= I2*50. + IsynII + IsynEI + Ce2*(mean(x2)-x2)*50. - 0.3*(mean(z)-3.)*50.
    dV = ((I_in - gCa*m_inf*(V-E_Ca) - gK*n*(V-E_K) - gL*(V-E_L))/Cm + 50.*noise2*(2.*(rand(n_i,1)-0.5)).ravel())
    #---slow variable exogenous
    dz = r*(s*(x1+mean(x2)-x0)-mean(z))
    
    #update variables
    x1=x1+dx1*dt
    x2=(V+dV*dt)/20.
    y1=y1+dy1*dt
    n=n+dn*dt
    z=z+dz*dt
    uE=uE+duE*dt
    uI=uI+duI*dt
    #store
    if(pb.mod(t_step,1)==0):
        x_1s[:,int(t_step)]=x1.ravel()
        x_2s[:,int(t_step)]=x2.ravel()
        z_s[:,int(t_step)]=z.ravel()
        y_1s[:,int(t_step)]=y1.ravel()
        n_s[:,int(t_step)]=n.ravel()
        print(t_step)  

   
t_plot=pb.arange(0,t_stop)
t_plot_dt=pb.arange(0,t_stop,dt)

#save data
radical="_%i_HMR_%i_ML_CpES1_%i_CpES2_%i_x0_%i_noise1_%i_noise2_%i_noise3_%i_gx1x2_%i_gx2x1_%i_gx1x1_%i_gx2x2_%i_r%i_%is" %(nbn1, nbn2, int(CpES1*100), int(CpES2*100), int(-x0*10), int(noise1*10), int(noise2*100), int(noise3*100), int(g_x1x2*100), int(g_x2x1*100), int(g_x1x1*100), int(g_x2x2*100), int(r*10000000),int(t_stop/1000))
np.save('x1s'+radical+'.npy', x_1s)
np.save('x2s'+radical+'.npy', x_2s)
np.save('y1s'+radical+'.npy', y_1s)
np.save('ns'+radical+'.npy', n_s)
np.save('zs'+radical+'.npy', z_s)
np.save('t'+radical+'.npy', t_plot)
#np.save('KOP'+radical+'.npy', KOP_plot)'''


fig = pb.figure(figsize=(30,25)) #pb.figure(figsize=(20,10))
#fig.hold(True)
fig.suptitle('%i pop1 neurons; %i pop2 neurons; I1=%.2f; I2=%.2f; x0=%.2f; r=%.5f, s=%i; CpES1=%.2f; CpES2=%.2f; g_x1x1=%.2f; g_x1x2=%.2f; g_x2x2=%.2f; g_x2x1=%.2f; c2=%.3f; noise1=%.2f; noise2=%.2f; noise3=%.2f' %(nbn1, nbn2,I1,I2, x0, r, s, CpES1, CpES2, g_x1x1, g_x1x2, g_x2x2, g_x2x1, c2, noise1, noise2, noise3)) 
    
ax1=fig.add_subplot(611)
ax2=fig.add_subplot(613)
ax3=fig.add_subplot(615)
ax4=fig.add_subplot(616)
ax5=fig.add_subplot(612)
ax6=fig.add_subplot(614)
#useless loops since it is vectorized now but whatever...
for i in range(nbn1):
        #time series pop1
        ax5.plot(t_plot, x_1s[i,:], label=None)
        drawGraph2.scatter_plot(t_plot, x1_plot[i], 0.0, i, ax1_1)
        #time series z
        ax4.plot(t_plot, z_s[i,:], label=None)
for j in range(nbn2): 
        #time series pop2
        ax6.plot(t_plot, x_2s[j,:], label=None)
        
#draw mean values in time series
ax1.imshow(x_1s, interpolation='nearest', aspect='auto')
ax5.plot(t_plot, mean(x_1s,axis=0), 'black', linewidth=1.5, label="x1bar")
ax5.set_xlim(right=t_plot.size)
ax1.legend(prop={'size':10})
ax1.set_title("Population 1 neurons")
ax1.set_xlim(left=0)
ax1.set_ylim((-0.5, nbn1-0.5))
ax2.imshow(x_2s, interpolation='nearest', aspect='auto')
ax6.plot(t_plot, mean(x_2s,axis=0), 'black', linewidth=1.5, label="x2bar")
ax6.set_xlim(right=t_plot.size)
ax2.legend(prop={'size':10})
ax2.set_title("Population 2 neurons")
#plot energy at last
ax3.plot(t_plot, sA.highPassFilter(x2bar_plot + x1bar_plot), 'purple', label='x2bar + x1bar filtered')
ax3.plot(t_plot, 0.2*mean(x_2s,axis=0) + 0.8*mean(x_1s,axis=0), 'cyan', label='x2bar + x1bar')
ax3.set_xlim(0,t_plot.size)
ax3.legend(prop={'size':10})
ax3.set_title("Whole system")
#draw z mean 
ax4.plot(t_plot, mean(z_s,axis=0), 'black', linewidth=1.5, label="zbar")
ax4.set_xlim(0,t_plot.size)
ax4.legend(prop={'size':10})
ax4.set_title("Slow variable")
#draw x0 and CpES
ax4.plot(t_plot_dt, CpES_variable, 'blue' ,linewidth=2.5, label="CpES")
ax4.plot(t_plot_dt, x0_variable, 'red',linewidth=2.5, label="x0")
ax5.legend(prop={'size':10})
#set figure name & save & show
figname="epilepton"+radical+".png"
fig.savefig(figname, dpi=200)
figTool.zoom(10000, 12000, x1_plot, y1_plot, z_plot, x2_plot, y2_plot, t_plot, KOP_plot, radical)
pb.show()