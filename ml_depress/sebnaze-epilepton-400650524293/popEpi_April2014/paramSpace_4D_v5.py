'''
Created on 10 janv. 2013

@author: Squirel

CAREFUL: this code is hell of a shit to maintain or understand due to row/col -> x/y relationships that differ from outer frame to inner fame
'''

import matplotlib
#matplotlib.use('Agg')
import numpy as np
import pylab as pb
from numpy import load
#import sys
#matplotlib.RcParams.update({'font.sans-serif':'Tahoma'})


nbn1=40
nbn2=40
#CpES=0.6
t_stop=20000

out_row="x0" # = out_y axis
out_col="CpES" # = out_x
in_row="gxixi"  # = in_y
in_col="gxixj"  # = in_x

out_row_val=pb.array([2, 2.5, 3, 3.5, 4])
out_row_val=out_row_val[::-1] #reverse array to have down -> axes up rather than down  
out_col_val=pb.array([0, 0.2, 0.4, 0.6, 0.8, 1]) 
in_row_val=pb.array([0, 0.2, 0.4, 0.6])
#in_row_val=in_row_val[::-1] #reverse array to have down -> axes up rather than down  
in_col_val=pb.array([0, 0.1, 0.2, 0.3])


def drawSubplot(out_row_index, out_col_index):
    in_row_length=len(in_row_val)
    in_col_length=len(in_col_val)
    in_img_grid=pb.zeros((in_row_length, in_col_length))
    in_img_grid_mean=pb.zeros((in_row_length, in_col_length))
    in_img_grid_KOP=pb.zeros((in_row_length, in_col_length))
    for in_row_index in range(in_row_length):
        for in_col_index in range(in_col_length):
            radical="_%i_HMR_%i_ML_CpES1_%i_CpES2_%i_x0_%i_noise3_0_noise2_20_noise1_40_gx1x2_%i_gx2x1_%i_gx1x1_%i_gx2x2_%i_%is" %(nbn1, nbn2, int(out_col_val[out_col_index]*100), int(out_col_val[out_col_index]*100), int(out_row_val[out_row_index]*100), int(in_col_val[in_col_index]*100), int(in_col_val[in_col_index]*100), int(in_row_val[in_row_index]*100), int(in_row_val[in_row_index]*100), int(t_stop/1000))
            #load data
            #path='C://Users/Squirel/workspace/HC-Septum-Network/popEpi_HMR_ML/'
            #x1_plot = np.load('traces_22feb/x1plot'+radical+".npy")
            #x2_plot = np.load('traces_22feb/x2plot'+radical+".npy")
            KOP1_plot = np.load('traces_july05/KOP1'+radical+".npy")
            KOP2_plot = np.load('traces_july05/KOP2'+radical+".npy")
            #epi=x1_plot+x2_plot
            #in_img_grid_mean[in_row_index, in_col_index]=epi.mean() #for the mean metric
            #in_img_grid[in_row_index, in_col_index] = get_coastline(epi.mean(axis=0)) #give the mean neural activity as a timeserie to compute the coastline
            in_img_grid_KOP[in_row_index, in_col_index] = KOP1_plot[0,:].mean() + KOP2_plot[0,:].mean()
    #img=ax.imshow(in_img_grid, interpolation=None, norm=None, vmin=-2.4, vmax=-1.9)
    #pb.colorbar(img)
    return in_img_grid, in_img_grid_mean, in_img_grid_KOP 


def get_coastline(timeserie):
    """careful: consider all timeseries have same sampling rate"""
    x=timeserie; y=range(len(timeserie));
    coastline_index=0
    for i in range(len(timeserie)-1):
        coastline_index += pb.sqrt( (x[i]-y[i])**2 + (x[i+1]-y[i+1])**2 )
    return coastline_index

def computeGrid():
    """compute grid with "mean" time series"""
    out_img_grid=[]
    out_img_grid_mean=[]
    out_img_grid_KOP=[]
    #loop to load data create grid
    for out_row_index in range(len(out_row_val)):
        row_grid=[]  #each cell containing 1 column
        row_grid_mean=[]
        row_grid_KOP=[]
        for out_col_index in range(len(out_col_val)):
            coastline_grid, mean_grid, KOP_grid = drawSubplot(out_row_index, out_col_index)
            row_grid.append(coastline_grid)
            row_grid_mean.append(mean_grid)
            row_grid_KOP.append(KOP_grid)
        out_img_grid.append(row_grid)
        out_img_grid_mean.append(row_grid_mean)
        out_img_grid_KOP.append(row_grid_KOP)
    #np.save('outImgGrid_v3.npy', out_img_grid)
    np.save('outImgGridMean_v5.npy', out_img_grid_mean)
    np.save('outImgGridKOP_v5.npy', out_img_grid_KOP)
    return out_img_grid

def drawFig(out_img_grid_KOP):
    #figure properties
    fig=pb.figure(figsize=(10,6))
    row_length=len(out_row_val)
    col_length=len(out_col_val)
    gs=matplotlib.gridspec.GridSpec(row_length, col_length, wspace=.05, hspace=.05) # +1 for the colormap area
    """figure processing"""
    #out_img_grid = out_img_grid.tolist()
    out_img_grid_KOP = out_img_grid_KOP.tolist()
    #get min and max to set up colormap scale, and loop all again to plot images
    #VMIN=pb.array(out_img_grid).min()
    #VMAX=pb.array(out_img_grid).max()
    KOPMIN=pb.array(out_img_grid_KOP).min()
    KOPMAX=pb.array(out_img_grid_KOP).max()
    for out_row_index in range(len(out_row_val)):
        for out_col_index in range(len(out_col_val)):
            ax=fig.add_subplot(gs[out_row_index, out_col_index])
            ax.hold(True)
            """out_img_grid[out_row_index][out_col_index] = ax.imshow(out_img_grid[out_row_index][out_col_index], \
	    								interpolation="nearest", norm=None, vmin=vmin, vmax=vmax, \
	    								extent=[in_col_val.min(), in_col_val.max(), in_row_val.min(), in_row_val.max()], \
									aspect='auto', origin='lower')"""
            
            """out_img_grid[out_row_index][out_col_index] = ax.contourf(out_img_grid[out_row_index][out_col_index], \
                                                                     levels=pb.linspace(VMIN,VMAX,100),  #100 steps for contours
                                                                     #linewidths=pb.linspace(0.2,2,100), \
                                                                     #colors='k', \
                                                                     #linestyles='solid', \
                                                                     extent=[in_col_val.min(), in_col_val.max(), in_row_val.min(), in_row_val.max()], \
                                                                     vmin=VMIN, vmax=VMAX)"""
            
            out_img_grid_KOP[out_row_index][out_col_index] = ax.contourf(out_img_grid_KOP[out_row_index][out_col_index], \
                                                                     levels=pb.linspace(KOPMIN,KOPMAX,50),  #100 steps for contours
                                                                     #alpha=1,
                                                                     #linewidths=pb.linspace(0.1,1,50), \
                                                                     #colors='k', \
                                                                     #cmap=pb.get_cmap('gray'),
                                                                     extent=[in_col_val.min(), in_col_val.max(), in_row_val.min(), in_row_val.max()], \
                                                                     vmin=KOPMIN, vmax=KOPMAX)
            
            #out_img_grid[out_row_index][out_col_index].set_clim(vmin=VMIN, vmax=VMAX)
            
            #set ticks and labels of axis (no number in legend for graphs not on borders):
            #for the x axis:
            xticks=[0,0.1,0.2] #from in_col_val
            ax.set_xticks(xticks)
            if(out_row_index < row_length-1):
                #ax.xaxis.set_major_locator(pb.MaxNLocator(nbins=len(in_col_val)))  # set number of ticks on x axis 
                ax.set_xticklabels([])
            else:
                ax.set_xlabel(in_col+"\n\n"+out_col+"="+str(out_col_val[out_col_index]))#, fontsize='xx-large')
                #ax.xaxis.set_major_locator(pb.MaxNLocator(nbins=len(in_col_val)))  # set number of ticks on x axis
                ax.set_xticklabels(map(str, xticks))#, fontsize='xx-large') #in_col_val
                
            #for the y axis:
            yticks=[0,0.2,0.4] #from in_row_val
            ax.set_yticks(yticks)    
            if(out_col_index > 0):
                #ax.yaxis.set_major_locator(pb.MaxNLocator(nbins=len(in_row_val)))  # set number of ticks on y axis
                ax.set_yticklabels([])
            else:
                ax.set_ylabel("|"+out_row+"|="+str(out_row_val[out_row_index])+'      '+in_row, rotation='horizontal')#, fontsize='xx-large')
                #ax.yaxis.set_major_locator(pb.MaxNLocator(nbins=len(in_row_val)))  # set number of ticks on y axis
                ax.set_yticklabels(map(str,yticks))#, fontsize='xx-large') #in_row_val

    #setup colorbar
    ax=pb.axes([0.92, 0.1, 0.01, 0.8])  #guess [left, bottom, width, heigth]. in percents
    #ticksValues=[1.3,1.5,1.7, 1.9]
    cbar=pb.colorbar( out_img_grid_KOP[0][0], ax)#, ticks=ticksValues)
    #cbar.ax.set_yticklabels(map(str,ticksValues), fontsize='xx-large')
    #fig.savefig("paramExplore_"+out_row+"_"+out_col+"_"+in_row+"_"+in_col+"_invAxes.png", dpi=200)
    pb.show()
    return fig

if __name__ == "__main__":
    #out_img_grid = computeGrid()
    out_img_grid_KOP = load('C://Users/Squirel/workspace/HC-Septum-Network/popEpi_variableParamaters/outImgGridKOP_v5.npy') #v3: initial, KOP but inacurate.  v4: KOPnoRest.   v5:KOP good.
    drawFig(out_img_grid_KOP)
            
