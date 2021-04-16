#!/usr/bin/env python
# coding: utf-8

# In[2]:


get_ipython().run_line_magic('matplotlib', 'inline')

get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')

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


import pfs_calibstars as cs



# ## Preparation

# In[20]:


# Set output directory 

output_dir = "../../pfs_calibstars_data/outputs/Synspec"


# Please prepare for a comma-separated table like this. 
#
#    * The header should be the same as below
#    * The column "synfilepath" should specifiy a path to the synthetic spectrum of a given star
#

catalogfile = "../catalogs/Synspec/catalog_Synspec.csv"
df = pd.read_csv(catalog_name, dtype = {'synfilepath': str})
df


# ## Produce many simulated spectra 

# In[ ]:


# Get a list of stars
ObjIDs = df["starname"].values


# Photometric band for magnitudes. At the moment only "sdss_g" is acceptable.
band = "sdss_g"


# Number of realizations of simulated spectra
nreal = 10




for ObjID in ObjIDs:
    
    synfile = ((df["synfilepath"][df["starname"] == ObjID]).values)[0]

    print(synfile)
    for mag in [17., 18., 19., 20.]: 
        
        exptime, nexp = cs.get_exptime("sdss_g", mag, "science")
        nexps = [nexp]

        
        cs.simulate_many_spectra(ObjID, catalogfile, nexps, band, mag, nreal,                                  synfile, output_dir, setting = "Optimistic",                                  write_h5 = False)


# 
