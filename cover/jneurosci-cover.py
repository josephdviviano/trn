#!/usr/bin/env python

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import nibabel as nib
import networkx as nx

"""
run in :~/data/MSC_JNEUROSCI/s01 :D.
"""

def calculate_normalizer(ts):
    minimum = np.min(ts)
    maximum = np.max(ts)

    divisor = np.max(np.array(np.abs(minimum), np.abs(maximum)))

    return divisor

# load timeseries data
mean_run = nib.load('mean_RUN_bpass_smoothHead.nii.gz')
# load anatomical ROIs resampled to match data
mask = nib.load('mask_LGN_TRN_resamp_split.nii.gz')
# load regions that modulate with retinotopy (q[FDR] = 0.05)
stat = nib.load('mask_stats_FFT_single_multitaper_LGNTRN_mean_RUN_allpass_smoothHeadLGNTRN.nii.gz')

# reshape data to nvox x ntrs
dims = mean_run.shape

mean_run = mean_run.get_data()
mask = mask.get_data()
stat = stat.get_data()

nvox = dims[0]*dims[1]*dims[2]
ntrs = dims[3]

mean_run = np.reshape(mean_run, (nvox, ntrs))
mask = np.reshape(mask, (nvox, 1))
stat = np.reshape(stat, (nvox, 1))

# set the ROIs from the mask: 
# 1 = rlgn, 
# 2 = llgn, 
# 3 = rvtrn, 
# 4 = lvtrn, 
# 5 = rdtrn, 
# 6 = ldtrn

rois = [1,2,4]

n_samps = 20

# plot 10 time series from each ROI
for i, val in enumerate(rois):
    # find the values within each ROI that are statistically significant
    idx = np.intersect1d(np.where(mask == val)[0], np.where(stat > 0)[0])

    plt.figure(num=None, figsize=(2, 8))
    # take the first 10, arbitrary
    for n in range(n_samps):
    
        # calculate normalizer (min/max)
        divisor = calculate_normalizer(mean_run[idx[n], :])

        # generate each subplot
        plt.subplot(n_samps, 1, n)
        plt.plot(mean_run[idx[n], :] / divisor, linewidth=2, color='black')
        plt.ylim(-2, 2)
        plt.axis('off')

    plt.suptitle(str(val))
    
    plt.savefig(str(val) + '.svg')
    plt.close()

    # if this is the first one, init the out_ts matrix
    if i == 0:
        out_ts = mean_run[idx[:n_samps], :]
    # otherwise stack them
    else:
        out_ts = np.vstack((out_ts, mean_run[idx[:n_samps], :]))

cmat = np.corrcoef(out_ts)

# plot full graph
plt.imshow(cmat, vmin=-1, vmax=1, cmap=plt.cm.RdBu_r, interpolation='nearest')
plt.axis('off')
plt.savefig('0_lvtrn_mat.svg')

# plot thresholded greaph
cmat[cmat <= 0] = 0
plt.close()
plt.imshow(cmat, vmin=-1, vmax=1, cmap=plt.cm.RdBu_r, interpolation='nearest')
plt.axis('off')
plt.savefig('0_lvtrn_thresholded-mat.svg')

# turn thresholded matrix into a networkx graph object 
g = nx.Graph(cmat)
nx.write_gexf(g, '0_lvtrn_thresholded-mat.gexf')