import matplotlib
#matplotlib.use('Agg')
import pylab as pb
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
 
nbn1=40
nbn2=40
t_stop=15000

CpESs=pb.array([1])
noise3=pb.array([0,2,4])
x0s=pb.array([2.5,3,3.5,4,4.5])
noise2=pb.array([0,0.2,0.4,0.6,0.8,1])

def computeKOPgrid():
    #compute KOP grid
    KOPgrid = pb.zeros((len(noise2),len(noise3),len(x0s)))
    i,j,k=0,0,0
    for i in range(len(noise2)):
        for j in range(len(noise3)):
            for k in range(len(x0s)):
                radical="_%i_HMR_%i_ML_CpES1_100_CpES2_100_x0_%i_noise3_%i_noise2_%i_noise1_%i_15s" %(nbn1, nbn2, int(x0s[k]*10), int(noise3[j]), int(noise2[i]*100),  int(noise2[i]*200),)
                KOP1 = np.load('traces_may22/KOP1'+radical+".npy")
                KOP2 = np.load('traces_may22/KOP2'+radical+".npy")
                KOPgrid[i,j,k] = KOP1[0,:].mean() + KOP2[0,:].mean()
    
    kopmin = KOPgrid.min()
    kopmax = KOPgrid.max()
    print KOPgrid
    #plot
    fig = pb.figure()
    ax = fig.add_subplot(311, projection='3d')
    #k=0 
    X,Y = pb.meshgrid(noise2,noise3)
    for k in range(len(x0s)):
        #img = ax.contourf(KOPgrid[:,:,k], vmin=kopmin, vmax=kopmax), offset=x0s[k])
        img = ax.contourf(X, Y, KOPgrid[:,:,k].T, zdir='z', levels=pb.linspace(kopmin, kopmax, 50), vmin=kopmin, vmax=kopmax, offset=x0s[k])
        #cbar = pb.colorbar(img)
        #fig.savefig("param3Dspace_CpES1_%i_CpES2_%i_x0_%i.png" %(int(CpES1s[i]*100), int(CpES2s[j]*100), int(x0s[k]*10)))
    ax.set_xlabel("noise2")
    ax.set_xticks(pb.linspace(min(noise2),max(noise2),5))
    ax.set_ylabel("noise3")
    ax.set_yticks(pb.linspace(min(noise3),max(noise3),5))
    ax.set_zlabel("x0")
    ax.set_zticks(pb.linspace(min(x0s),max(x0s),5))
    ax.set_zlim(min(x0s), max(x0s))
    
    ax2=fig.add_subplot(312)
    index_noise3=0;
    X2,Y2=pb.meshgrid(x0s, noise2)
    img2 = ax2.contourf(X2, Y2,KOPgrid[:,index_noise3,:],levels=pb.linspace(kopmin, kopmax, 50), vmin=kopmin, vmax=kopmax)
#    img = ax2.imshow(KOPgrid[index_noise2,:,:].T, interpolation="nearest", origin="lower")
    ax2.set_xlabel("|x0|")
    ax2.set_ylabel("noise_1_2")
    cbar=pb.colorbar(img2)
    
    noiseGrid=pb.zeros((len(noise2),len(x0s)))
    for n in range(len(noise2)):
        noiseGrid[n,:] = KOPgrid[n,0,:]
    ngmin=noiseGrid.min()
    ngmax=noiseGrid.max()
    
    ax3=fig.add_subplot(313)
    X3,Y3=pb.meshgrid(x0s, range(len(noise2)))
    img3 = ax3.contourf(X3, Y3, noiseGrid, levels=pb.linspace(ngmin, ngmax, 50), vmin=ngmin, vmax=ngmax)
    #img2 = ax3.imshow(noiseGrid.T, interpolation="nearest", origin="lower")
    cbar2=pb.colorbar(img3)
    ax3.set_xlabel("|x0|")
    ax3.set_ylabel("noise_1_2_3")
    
    print noiseGrid.T
    
    pb.show()
    
"""  DEPRECATED    
def computeCSgrid():
    #compute KOP grid
    CSgrid = pb.zeros((len(noise2),len(noiseRatios),len(x0s)))
    i,j,k=0,0,0
    for i in xrange(len(noise2)):
        for j in xrange(len(noiseRatios)):
            for k in xrange(len(x0s)):
                radical="_%i_HMR_%i_ML_CpES1_100_CpES2_100_x0_%i_noiseRatio_%i_noise2_%i_15s" %(nbn1, nbn2, int(x0s[k]*10), int(noiseRatios[j]*10), int(noise2[i]*10))
                CS = np.load('traces_01mai/CL'+radical+".npy")
                CSgrid[i,j,k] = pb.real(CS)
    
    CSmin = CSgrid.min()
    CSmax = CSgrid.max()
    print CSgrid
    #plot
    fig = pb.figure()
    ax = fig.add_subplot(211, projection='3d')
    k=0 
    for k in range(len(x0s)):
        #img = ax.contourf(KOPgrid[:,:,k], vmin=kopmin, vmax=kopmax), offset=x0s[k])
        img = ax.contourf(noise2, noiseRatios, CSgrid[:,:,k], zdir='z', levels=pb.linspace(CSmin, CSmax, 50), vmin=CSmin, vmax=CSmax, offset=x0s[k])
        #cbar = pb.colorbar(img)
        #fig.savefig("param3Dspace_CpES1_%i_CpES2_%i_x0_%i.png" %(int(CpES1s[i]*100), int(CpES2s[j]*100), int(x0s[k]*10)))
    ax.set_xlabel("noise2")
    ax.set_xticks(pb.linspace(min(noise2),max(noise2),5))
    ax.set_ylabel("noiseRatio")
    ax.set_yticks(pb.linspace(min(noiseRatios),max(noiseRatios),5))
    ax.set_zlabel("x0")
    ax.set_zticks(pb.linspace(min(x0s), max(x0s),5))
    ax.set_zlim(min(x0s),max(x0s))
    
    ax2=fig.add_subplot(212)
    index_noise2=0;
    X,Y=pb.meshgrid(x0s, noiseRatios)
    img = ax2.contourf(X, Y,CSgrid[index_noise2,:,:],levels=pb.linspace(CSmin, CSmax, 50), vmin=CSmin, vmax=CSmax)
    #img = ax2.imshow(CSgrid[index_noise2,:,:], interpolation='nearest', origin='lower')
    cbar=pb.colorbar(img)
    
    
    pb.show()
    
def computeCLKOPgrid():
    #compute KOP grid
    CLKOPgrid = pb.zeros((len(noise2),len(noiseRatios),len(x0s)))
    i,j,k=0,0,0
    for i in xrange(len(noise2)):
        for j in xrange(len(noiseRatios)):
            for k in xrange(len(x0s)):
                radical="_%i_HMR_%i_ML_CpES1_100_CpES2_100_x0_%i_noiseRatio_%i_noise2_%i_15s" %(nbn1, nbn2, int(x0s[k]*10), int(noiseRatios[j]*10), int(noise2[i]*10))
                KOP1 = np.load('traces_01mai/KOP1'+radical+".npy")
                KOP2 = np.load('traces_01mai/KOP2'+radical+".npy")
                CL = np.load('traces_01mai/CL'+radical+".npy")
                CLKOPgrid[i,j,k] = (KOP1[0,:].mean() + KOP2[0,:].mean())* pb.real(CL)
    
    clkopmin = CLKOPgrid.min()
    clkopmax = CLKOPgrid.max()
    #plot
    fig = pb.figure()
    ax = fig.add_subplot(211, projection='3d')
    #k=0 
    for k in xrange(len(x0s)):
        #img = ax.contourf(KOPgrid[:,:,k], vmin=kopmin, vmax=kopmax), offset=x0s[k])
        img = ax.contourf(noise2, noiseRatios, CLKOPgrid[:,:,k], zdir='z', levels=pb.linspace(clkopmin, clkopmax, 50), vmin=clkopmin, vmax=clkopmax, offset=x0s[k])
        #cbar = pb.colorbar(img)
        #fig.savefig("param3Dspace_CpES1_%i_CpES2_%i_x0_%i.png" %(int(CpES1s[i]*100), int(CpES2s[j]*100), int(x0s[k]*10)))
    ax.set_xlabel("noise2")
    ax.set_xticks(pb.linspace(min(noise2),max(noise2),5))
    ax.set_ylabel("noiseRatio")
    ax.set_yticks(pb.linspace(min(noiseRatios),max(noiseRatios),5))
    ax.set_zlabel("x0")
    ax.set_zticks(pb.linspace(min(x0s),max(x0s),5))
    ax.set_zlim(min(x0s), max(x0s))
    
    ax2=fig.add_subplot(212)
    index_noise2=0;
    X,Y=pb.meshgrid(x0s, noiseRatios)
    img = ax2.contourf(X, Y,CLKOPgrid[index_noise2,:,:],levels=pb.linspace(clkopmin, clkopmax, 50), vmin=clkopmin, vmax=clkopmax)
    #img = ax2.imshow(KOPgrid[index_noise2,:,:])
    ax2.set_xlabel("|x0|")
    ax2.set_ylabel("noise_1")
    cbar=pb.colorbar(img)
    
    
    pb.show()
"""    
    
if __name__ == "__main__":
    computeKOPgrid()