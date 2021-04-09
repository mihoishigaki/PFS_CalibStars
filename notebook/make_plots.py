#!/usr/bin/env python
# coding: utf-8

# In[68]:


get_ipython().run_line_magic('matplotlib', 'inline')

get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')

import imp
import os, glob,sys,pprint
import numpy as np
import matplotlib.pyplot as plt

import synphot
from astropy.utils.data import get_pkg_data_filename
from astropy.io import ascii,fits
import pandas as pd


from pfsspecsim import pfsetc
from pfsspecsim import pfsspec


# Please specify a place where the codes are located 
pfs_calibstars_dir = '/Users/ishigakimiho/PFS/Github/pfs_calibstars'
sys.path.append(pfs_calibstars_dir)


database_dir = "../../pfs_calibstars_data/database"
output_dir = "../../pfs_calibstars_data/outputs"

from pfs_calibstars import *
import pfs_calibstars as cs
import pfs_calibstars.utility as ut



# In[2]:


#dir(cs)


# In[50]:


# Get input files
ObjIDs = ["00081", "00092"]
#ObjIDs = ["00092"]
objtypes = ["K-giant", "G-dwarf"]
#objtypes = ["G-dwarf"]
libID = "HDS"
band = "sdss_g"
mags = np.arange(17, 24, 1)
cont_file = '../utility/specregion/contregions_MR.csv'
wcover_file = '../utility/specregion/wcoverage_MR.csv'
outdir = "../../pfs_calibstars_data/ssp"
rv0 = 100.
cs.measure_rv_mags(ObjIDs, objtypes, libID, output_dir, band, mags, outdir, rv0, cont_file, wcover_file)


# In[ ]:





# In[47]:



h5file = '../../pfs_calibstars_data/outputs/HDS/00092/MR/MP0.00_FA0.000_MTangle060_TP1.00/G53-41_Teff5907_logg4.3_MH-1.6_Vrad0.0_sdss_g_22.0_Texp10800.snc.sim.h5'
template = '../../pfs_calibstars_data/database/HDS/00092/Synspec/WL_Flux_Error_G53-41_Teff5907_logg4.3_MH-1.4_Vrad0.0.txt'
rv0 = 100.
cont_file = '../utility/specregion/contregions_MR.csv'
wcover_file = '../utility/specregion/wcoverage_MR.csv'
cs.measure_rv(h5file, template, rv0, cont_file, wcover_file, plot = False)


# In[109]:



ObjIDs = ["00092", "00081"]
objtypes = ["G-dwarf", "K-giant"]
mags = [[18, 19, 20, 21, 22], [18, 19, 20, 21, 22, 23]]
mks = ['o', '^']
abundfile = "../../pfs_calibstars_data/ssp/Ivanna/pfs_fe_alpha_precision_data.fits"
cs.plot_mag_sn_rv_abund(ObjIDs, objtypes, mags, mks, abundfile)



# In[70]:


abundfile = "../../pfs_calibstars_data/ssp/Ivanna/pfs_fe_alpha_precision_data.fits"
hdul = fits.open(abundfile)
head = hdul[1].header
data = hdul[1].data
gmag_abund = data['gmag']
feherr_G = data['feherr_G']
feherr_K = data['feherr_K']
alphafeerr_G = data['alphafeerr_G']
alphafeerr_K = data['alphafeerr_K']


plt.plot(gmag_abund, feherr_G)
plt.show()


# In[30]:


path = output_dir + "/HDS/00081/MR/MP0.00_FA0.000_MTangle060_TP1.00/"
file = "BD+292356_Teff4827_logg1.6_MH-1.7_Vrad0.0_sdss_g_21.0_Texp7200.snc.sim_00081.91.dat"


import matplotlib.pyplot as plt
wv, fl = np.loadtxt(path + file, usecols = (0, 1), unpack = True)


# In[38]:


plt.plot(wv, fl)
plt.ylim(0.0, 100000)
plt.xlim(710., 885.)


# In[ ]:




