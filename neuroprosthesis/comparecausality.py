"""
This script generates an output figure for the NEST paper -- comparing the 
causality in the healthy, damaged, and restored cases. It's been updated.

Usage: just run!

Version: 2012dec10
"""


# Import useful tools -- not sure how many of these will actually be used
import time; tic=time.clock()
from scipy import loadtxt, size, shape, zeros, mod, floor, mean
from pylab import figure, plot, xlabel, ylabel, legend, xlim, ylim, show, hold, squeeze, sqrt
from bsmart import granger # Load the Granger calculation tool


nscenarios=3; # How many scenarios are being evaluated?
calc=1; # Should results be recalculated?

if calc: # Save some time if just plotting adjustments
    inputstems=['./neuropros1','./neuropros2','./neuropros3'];
    nfiles=size(inputstems) # Number of files
    layerA=1 # Which layers to compare? First layer default: 1=L4
    layerB=2 # Which layers to cmopare? Second layer default: 2=L5
    killtime=200 # Number of time points to remove from beginning and end
    results=[] # Initialize list to store results
    for i in range(nfiles):
        print "Loading LFP data for file %i..." % i
        lfpfilename="%s-lfp.txt" % inputstems[i] # Set LFP filename (if used)
        lfporig=loadtxt(lfpfilename) # Load the file
        [rows,cols]=shape(lfporig) # Find out how big the array is
        nlayers=5 # Since layers and columns are interlaced, dividing the total columns by the number of model columns gives layers
        lfpdata=zeros((rows-2*killtime,nlayers,cols/nlayers)) # Dimensions are: | time | layer | column |
        for j in range(cols):# Reshape the array from 2D to 3D (time, layer, column)
            L=mod(j,nlayers) # Which layer is it? Loop with period nlayers.
            C=int(floor(j/nlayers)) # Which column is it? Increment every nlayers.
            lfpdata[:,L,C]=lfporig[killtime:-killtime,j] # Reshape array, eliminating the first <killtime> points
            lfpdata[:,L,C]-=mean(lfpdata[:,L,C])
            
        print 'Calculating Granger...'
        F,pp,cohe,Fx2y,Fy2x,Fxy=granger(lfpdata[:,layerA,0],lfpdata[:,layerB,0],order=20) # Pick out LFPs from layers (default 2/3 and 5)
        results.append(Fx2y)

print "Plotting..."
figh=figure(figsize=(10,6)) # Open figure and store its handle

# Plot Granger spectra
hold(True)
labels=list()
labels=['Baseline','With damage','With prosthesis']
alldata=list()
colors=[[0,0,1],[1,0,0],[0,0.5,0]]
for i in range(nscenarios):
    data=squeeze(results[i])
    plot(F,data,label=labels[i],linewidth=2,c=colors[i])
    alldata.append(data)
    xlim(0,50)
    ylim(0,0.5) # Need to define explicitly since otherwise 
xlabel('Frequency (Hz)')
ylabel('Granger causality')
legend()

toc=time.clock()
print 'Done; elapsed time was %0.1f seconds.' % (toc-tic)


# Do stats
endpt=size(alldata[0],0)/1.2 # Only go out to 50 Hz
startpt=size(alldata[0],0)/12 # Start at 5 Hz
totals=zeros((nscenarios))
for i in range(nscenarios): totals[i]=alldata[i][startpt:endpt].sum()

discrepancy=zeros((nscenarios-1))
for i in range(nscenarios-1): discrepancy[i]=100*(totals[0]-totals[i+1])/totals[0]

rmse=zeros((nscenarios-1))
for i in range(nscenarios-1): rmse[i]=sqrt(sum((1-alldata[i+1][0:endpt]/alldata[0][0:endpt])**2))

# Shorten names
b=totals[0]
d=totals[1]
p=totals[2]

print "File: %s" % inputstems[0]
print "Percentage preserved with damage: %f" % (100*d/b)
print "Percentage preserved with prosthesis: %f" % (100*p/b)

# Do more stats
from scipy import stats
start=size(alldata[0],0)/12 # Start at 5 Hz
corrbd=stats.spearmanr(alldata[0][startpt:endpt],alldata[1][startpt:endpt]) # Correlation between baseline and damage
corrbp=stats.spearmanr(alldata[0][startpt:endpt],alldata[2][startpt:endpt]) # Correlation between baseline and prosthesis

show()
