#!/usr/bin/env python

# plots stimulus paramaters
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm

## Options ####################################################################
oFreq       = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
oFreqLabels = [1, 2.5, 5, 7.5, 10, 12.5, 15, 20, 30, 60]
oTime = 210
oRotPeriod = 21
oFlickPeriod = 30
oResolution = 500

## Generate Stimulus Vectors ##################################################
nFlickCycles = oTime / oFlickPeriod
nFlickUnits  = oFlickPeriod / len(oFreq)

flickCycle = []
for f in oFreq:
	flickSection = np.repeat(f, nFlickUnits*oResolution)
	flickCycle   = np.concatenate([flickCycle, flickSection])

time  = np.arange((oTime)*oResolution)
rot   = np.sin(2*np.pi*time / (oRotPeriod*oResolution))
flick = np.tile(flickCycle, nFlickCycles)

## Plot Stimulus ##############################################################
cMap = ListedColormap(plt.cm.jet(np.arange(len(oFreq)).astype(np.float)
	                             /float(len(oFreq)))[:, 0:3])

pnts = np.array([time, rot]).T.reshape(-1,1,2)
segs = np.concatenate([pnts[:-1], pnts[1:]], axis=1)

fig = plt.figure()
axLine = LineCollection(segs, cmap=cMap)
axLine.set_array(flick)
#axLine.set_linewidth(3)
#plt.gca().add_collection(axLine)

#plt.xlim(time.min(), time.max()+1)
#plt.xticks(np.arange(0, len(time)+1, oRotPeriod*oResolution),
#	       np.arange(0, len(time), oRotPeriod),  size=36)
#plt.xlabel('Time (s)', size=36)

#plt.ylim(-1.1, 1.1)
#plt.yticks([-1, 1], (r'$-\pi$', r'$\pi$'), size=44)
#plt.ylabel('Rotation (rads)', size=36)

axCBar = fig.colorbar(axLine)
axCBar.set_ticks(np.linspace(1.5, 9.5, num=10))
axCBar.set_ticklabels(oFreqLabels)
axCBar.ax.tick_params(labelsize=48)

plt.box('on')
fig.set_frameon('False')
fig.set_facecolor('white')
plt.show()
## JDV Jun 23 2013 ############################################################
