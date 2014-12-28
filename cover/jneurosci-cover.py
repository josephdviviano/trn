#!/usr/bin/env python

import numpy as np
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

def load_nifti(filename, returns_dims=False):
    """
    Loads a nifti file, returns the data matrix reshaped to 2D and it's original
    dimensions.

    If returns_dims is True, we also return the 3D voxel dimensions, the number
    of voxels, and the number of TRs.
    """
    data = nib.load(filename)

    # retrieve the dimensions
    dims = data.shape
    nvox = dims[0]*dims[1]*dims[2]
    ntrs = dims[3]

    # retrieve the matrix
    data = data.get_data()
    data = np.reshape(data, (nvox, ntrs))

    if returns_dims == True:
        return data, dims, nvox, ntrs
    else:
        return data

# mask ROIs: 1 = rlgn, 2 = llgn, 3 = rvtrn, 4 = lvtrn, 5 = rdtrn, 6 = ldtrn
rois = [1,2,4]
n_samps = 20
outnames = ['01_rlgn-ts', '02_llgn-ts', '03_lvtrn-ts']


# load timeseries data, anatomical ROIs, and statmask from retino (q[FDR]=0.05)
mean_run = load_nifti('mean_RUN_bpass_smoothHead.nii.gz')
mask = load_nifti('mask_LGN_TRN_resamp_split.nii.gz')
stat = load_nifti('mask_stats_FFT_single_multitaper_LGNTRN_mean_RUN_allpass_smoothHeadLGNTRN.nii.gz')

# plot n time series from each ROI
for i, val in enumerate(rois):
    
    # find the values within each ROI that are statistically significant
    idx = np.intersect1d(np.where(mask == val)[0], np.where(stat > 0)[0])

    plt.figure(num=None, figsize=(2, 8))
    # take the first n, arbitrary
    for n in range(n_samps):
    
        # calculate normalizer (min/max)
        divisor = calculate_normalizer(mean_run[idx[n], :])

        # generate each subplot
        plt.subplot(n_samps, 1, n)
        plt.plot(mean_run[idx[n], :] / divisor, linewidth=2, color='black')
        plt.ylim(-2, 2)
        plt.axis('off')

    plt.suptitle(str(val))
    
    plt.savefig(outnames[i] + '.svg')
    plt.close()

    # if this is the first one, init the out_ts matrix
    if i == 0:
        out_ts = mean_run[idx[:n_samps], :]
    # otherwise stack them
    else:
        out_ts = np.vstack((out_ts, mean_run[idx[:n_samps], :]))

# generate correlation matrix
cmat = np.corrcoef(out_ts)

# plot full graph and thresholded graph
plt.imshow(cmat, vmin=-1, vmax=1, cmap=plt.cm.RdBu_r, interpolation='nearest')
plt.axis('off')
plt.savefig('04_full-mat.svg')
plt.close()

cmat[cmat <= 0] = 0
plt.imshow(cmat, vmin=-1, vmax=1, cmap=plt.cm.RdBu_r, interpolation='nearest')
plt.axis('off')
plt.savefig('04_pos-mat.svg')

# turn thresholded matrix into a networkx graph object 
g = nx.Graph(cmat)

# number the nodes by ROI
for n in g.nodes():
    if n < n_samps: 
        g.node[n]['label'] = 1

    elif n >= n_samps and n < n_samps*2:
        g.node[n]['label'] = 2

    else:
        g.node[n]['label'] = 3

nx.write_gexf(g, '05_pos-mat.gexf')
