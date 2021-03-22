#!/usr/bin/env python
# coding: utf-8

# In[2]:


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

import pfs_calibstars as cs



# ## Create databases of simulated spectra 

# In[3]:


libID = "HDS"
cs.produce_database(libID, multiproc = True)


# ## Produce many simulated spectra 

# In[9]:


libID = "HDS"
nreal = 100

ObjIDs = ["00092", "00081"]

band = "sdss_g"

for ObjID in ObjIDs:
    
    for mag in [17., 18., 19., 20., 21., 22., 23.]: 
        
        exptime, nexp = cs.get_exptime("sdss_g", mag, "science")
        nexps = [nexp]

        cs.simulate_many_spectra(ObjID, libID, nexps, band, mag, nreal,                                  database_dir, output_dir, setting = "Optimistic",                                  write_h5 = True)


# In[ ]:




