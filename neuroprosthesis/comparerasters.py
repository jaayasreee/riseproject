"""
This script generates an output figure for the NEST paper -- comparing the 
rasters in the healthy, damaged, and restored cases. It's mostly copied-and-
pasted out of plotspklfp.py, and bears an even stronger resemblance to
ckintf/ck/1103/comparerasters.py.

It also plots LFP time series, so the name is not very informative. Sorry.

Usage: just run!, but choose cortical or thalamic with the switch on line 22.

Version: 2012dec10
"""


# Import useful tools -- not sure how many of these will actually be used
import time; tic=time.clock()
from scipy import loadtxt, size, shape, zeros, mod, floor, array, ceil
from pylab import figure, plot, scatter, xlabel, ylabel, xlim, ylim, hold, figtext, r_, show

whichlfp=4 # Which LFP to use, default 4 = from all layers
calculate=1 # Whether or not to recalculate

if calculate:
    dr='/home/cliffk/bill/ckintf/ck/1109' # Directory
    inputstems=['./neuropros1','./neuropros2','./neuropros3'];
    nfiles=size(inputstems) # Number of files
    layernames=['L2/3','L4','L5','L6']
    killtime=0 # Number of time points to remove from beginning and end
    fs=200 # Sampling rate
    results=list() # Initialize list to store results
    changedindices=1; # WARNING 2011mar21 -- this is needed since Bill/Sam seem to have changed which indices correspond to which cells 
    popdata=list() # For storing spikes per population
    lfpdata=list()
    for i in range(nfiles):
        popdata.append([]) # Add an element to popdata to store this sim's entire raster
        print "Loading spike data for file %i..." % i
        allspkdata=loadtxt("%s-spk.txt" % inputstems[i]) # Set spike filename (if used)) # Load the file
        cellscale=ceil(max(allspkdata[:,1])/470) # Calculate how many more cells than normal there are -- 470 is standard, so any more is unusual
        popinfo=[["I2L", "I2", "E2", "I4L", "I4", "E4", "I5L", "I5", "E5R", "E5B", "I6L", "I6", "E6"],array([23, 20, 18, 17, 16, 15, 14, 13, 12, 11, 10, 8, 7])+changedindices,array([13, 25, 150, 14, 20, 30, 13, 25, 65, 17, 13, 25, 60])*cellscale] # List cell population names, numbers, and number of cells in each
        ncells=int(sum(popinfo[2])) # Calculate total number of cells
        npops=size(popinfo,1) # Get the number of cell populations
        maxtime=3 # Maximum time is the time of the latest spike in seconds, rounded up to the nearest integer
        numcols=int(allspkdata[-1,3]+1) # Number of columns is the last number listed in the file, plus one
        nspikes=size(allspkdata,0) # Find out how many spikes across all populations and columns
        tmp2=list() # Re-initialize the spikes-per-population for this file
        for j in range(npops):
            tmp2=allspkdata[allspkdata[:,2]==popinfo[1][j],:] # Pull out the spikes from a particular population
            popdata[i].append(tmp2) # Add one population to popinfo
            
        # And now for LFPs...
        print "Loading LFP data..."
        lfporig=loadtxt("%s-lfp.txt" % inputstems[i]) # Set spike filename (if used)) # Load the file
        [rows,cols]=shape(lfporig) # Find out how big the array is
        nlayers=int(cols/numcols) # Since layers and columns are interlaced, dividing the total columns by the number of model columns gives layers
        tmp=zeros((rows-2*killtime,nlayers,cols/nlayers)) # Dimensions are: | time | layer | column |
        for j in range(cols):# Reshape the array from 2D to 3D (time, layer, column)
            L=mod(j,nlayers) # Which layer is it? Loop with period nlayers.
            C=int(floor(j/nlayers)) # Which column is it? Increment every nlayers.
            tmp[:,L,C]=lfporig[:,j] # Reshape array, eliminating the first <killtime> points
        lfpdata.append(tmp)
    
    # Normalize LFPs
    lfpmax=0;
    for i in range(nfiles): lfpmax=max(lfpmax,lfpdata[i][:,whichlfp,0].max())
    for i in range(nfiles): lfpdata[i]/=lfpmax
    
print "Plotting..."
figh=figure(figsize=(8,11)) # Open figure and store its handle
figh.subplots_adjust(left=0.16) # Less space on left
figh.subplots_adjust(right=0.98) # Less space on right
figh.subplots_adjust(top=0.95) # Less space on top
figh.subplots_adjust(bottom=0.06) # Less space on bottom
figh.subplots_adjust(wspace=0.25) # More space between
figh.subplots_adjust(hspace=0.2) # More space between
figtext(0.03,0.85,'A',size='xx-large')
figtext(0.03,0.62,'B',size='xx-large')
figtext(0.03,0.39,'C',size='xx-large')
figtext(0.03,0.16,'D',size='xx-large')


# Plot rasters
colors=array([[0.42,0.67,0.84],[0.50,0.80,1.00],[0.90,0.32,0.00],[0.34,0.67,0.67],[0.42,0.82,0.83],[0.90,0.59,0.00],[0.33,0.67,0.47],[0.42,0.83,0.59],[0.90,0.76,0.00],[1.00,0.85,0.00],[0.71,0.82,0.41],[0.57,0.67,0.33],[1.00,0.38,0.60]]) # Colors for each cell population
spc=0 # Only one column, so pick it
for i in range(nfiles):
    thisaxis=figh.add_subplot(4,1,i+1)
    hold(True)
    for j in range(npops):
        if size(popdata[i][j][:,3]==spc):
            thisdata=popdata[i][j][popdata[i][j][:,3]==spc,:] # Pull out spikes from column spc
            # Reorder spikes so the E2 is at the top of the plot and I6 is at the bottom
            scatter(thisdata[:,0]/1000.,thisdata[:,1],s=5,c=colors[j],edgecolors='none')
            xlim(0,maxtime) # To convert to seconds
            ylim(0,ncells+10) # Just larger than the number of cells in the model
            ylabel("Neuron number")
    thisaxis.set_xticks(r_[0:maxtime+1]) # Don't show half seconds

# Plot LFPs
lfpcolors=[[0,0,1],[1,0,0],[0,0.5,0]] # Use colors from comparecausality
for i in range(nfiles):
    offset=20 # Don't show the first few points
    thisaxis=figh.add_subplot(4,1,4)
    hold(True)
    timeaxis=r_[0:maxtime*fs]/float(fs) # Set time axis
    npts=size(timeaxis,0)
    tmp=plot(timeaxis,lfpdata[i][0+offset:npts+offset,whichlfp,0],linewidth=1,c=lfpcolors[i]) # The middle index -- default 4 -- indicates which LFP to use
    xlim(0,maxtime)
    ylim(0.5,1.0) # Cortex: big voltage change
    if i==nfiles-1: xlabel('Time (s)')
    ylabel('Normalized voltage')
    thisaxis.set_xticks(r_[0:maxtime+1]) # Don't show half seconds


toc=time.clock()
print 'Done; elapsed time was %0.1f seconds.' % (toc-tic)
show()
