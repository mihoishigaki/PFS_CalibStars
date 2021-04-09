import numpy as np


import matplotlib.pyplot as plt
import scipy.optimize as so
import sys, logging, glob

import types

from astropy.io import fits

import pfs_calibstars.utility as ut
import pfs_calibstars.simulate as sm


import pandas as pd


def plot_mag_sn_rv_abund(ObjIDs, objtypes, mags, mks, abundfile):

    plt.rcParams["font.size"] = 18
    cmap = plt.get_cmap("tab10")

    fig, ax = plt.subplots(1, 2, figsize = (15, 8))

    # Define a second axis
    ax2 = ax[0].twinx()


    
    for i, ObjID in enumerate(ObjIDs):

        datafile = "../../pfs_calibstars_data/ssp/mag_rv_sn_" + \
            objtypes[i] + ".csv"
        df = pd.read_csv(datafile)

        gmags = ()
        rvs = ()
        snrs = ()
        for j, mag in enumerate(mags[i]):

            filt = df['gmag'] == mag

            gmags = np.append(gmags, df['gmag'][filt])
            rvs = np.append(rvs, df['rvsig'][filt])
            snrs = np.append(snrs, df['SN'][filt])

        ax[0].plot(gmags, snrs, linestyle = "-", marker = mks[i], \
                   mec = cmap(i), mfc = 'none', \
                   ms = 12, label = "S/N (" + objtypes[i] + ")", \
                   lw = 2, alpha = 0.8)

        ax2.plot(gmags, rvs, linestyle = "--", marker = mks[i], \
                   color = cmap(i), ms = 12, label = "RVerr (" + objtypes[i] + ")", \
                 lw = 2, alpha = 0.8)



        
    # Get abund data
    hdul = fits.open(abundfile)
    head = hdul[1].header
    data = hdul[1].data
    gmag_abund = data['gmag']
    feherr_G = data['feherr_G']
    feherr_K = data['feherr_K']
    alphafeerr_G = data['alphafeerr_G']
    alphafeerr_K = data['alphafeerr_K']


    for k in range(0, 2):

        filt = (gmag_abund >= np.min(mags[k])) & (gmag_abund <= np.max(mags[k])) 
        
        if k == 0:
            
            feherr = feherr_G[filt]
            alphaerr = alphafeerr_G[filt]
            lab = "(G-dwarf)"
        else:
            feherr = feherr_K[filt]
            alphaerr = alphafeerr_K[filt]
            lab = "(K-giant)"
            
        ax[1].plot(gmag_abund[filt], feherr, ls = "-", lw = 2, \
                   marker = mks[k], ms = 12, mec = cmap(k), mfc = 'none',\
                   label = "[Fe/H] " + lab, alpha = 0.8)

        ax[1].plot(gmag_abund[filt], alphaerr, ls = ":", lw = 2, \
                 marker = mks[k], ms = 12, color = cmap(k), \
                   label = r"[$\alpha$/Fe] " + lab, alpha = 0.8)



        
    ax[0].set_xlabel("gmag")
    ax[0].set_ylabel("S/N (per resolution, @800nm, ETC-ver1.2)")
    ax[0].set_ylim(0, 105)
    
    ax2.set_ylabel("RV error [km/s]")
    ax2.set_ylim(0, 27)
    
    
    ax[0].legend(loc = 'upper left', bbox_to_anchor = (0.0, 1.00), ncol = 2, \
                 prop = {"size":17}, columnspacing = 1.0, frameon=False)
    ax2.legend(loc = 'upper left', bbox_to_anchor = (0.0, 0.95), ncol = 2, \
               prop = {"size":17}, columnspacing = 1.0, frameon=False)

    ax[1].set_xlabel("gmag")
    ax[1].set_ylabel('Abundance precision [dex]')
    ax[1].set_ylim(0, 0.15)
    ax[1].legend(prop = {"size":17}, frameon=False)

    fig.tight_layout(pad=1.5)

    plt.savefig("../../pfs_calibstars_data/ssp/gmag_sn_rv_abund.png")

    
    return()

