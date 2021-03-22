# pfs_calibstars


## Abstract

This repository contains python codes to calculate simulated PFS
spectra from theoretical stellar spectra ('pfs_calibstars' directory).  

As inputs, the codes take a theoretical spectrum which mimics an existing
(previously obverved) star ('inputs' directory). Coordinates (RA, DEC) and 
physical parameters (Teff, logg, [Fe/H]) of the stars should be 
prepared in the 'catalogs' directory.  

The two main capabilities of pfs_calibstars are shown in 'notebook/main.ipynb':   

(1) Producing simulated PFS spectra from a theoretical spectrum of 
    an observed magnitude with an observing mode of LR or MR and 
    with either "Optimistic" and "Pessimistic" observing condition.  

(2) Producing simulated PFS spectra from theoretical spectra of 
    various magnitudes with a fixed mode and conditiion.  


To run the code in main.ipynb, please specify paths to the codes and 
the output directories.  


-----------------------  

## Directory structure  


- inputs:              Input theoretical spectra calculated to fit the
                       observed stars.   

- pfs_calibstars :     Python codes to generate simulated PFS spectra with
                       desired settings.    

- notebook:            Example scripts to run the codes.  

- catalogs:            Catalogs of observed stars and their estimated physical
                       parameters.  

- speclib:             Library of synthetic and observed spectra.  

- utility:             Utility data, e.g., photometric filters, isochrones.  





