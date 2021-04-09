#!/usr/bin/env python
# coding: utf-8

# In[1]:


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



# In[5]:


synspecpath = database_dir + "/HDS/00081/Synspec/"
synspecfile = synspecpath + "WL_Flux_Error_BD+292356_Teff4827_logg1.6_MH-1.7_Vrad0.0.txt"
w, f = np.loadtxt(synspecfile, usecols = (0, 1), unpack = True)
band = "sdss_g"
mag = 17.

wv, mag = cs.flux2ABmag(w, f, band, mag)

plt.plot(wv, mag)
plt.show()


# In[ ]:




