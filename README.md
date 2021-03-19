# pfs_calibstars

This repository contains python codes to calculate simulated PFS
spectra from theoretical stellar spectra ('pfs_calibstars' directory).
As inputs, the codes take theoretical spectra ('input' directory) calculated
with stellar physical parameters estimated for observed stars
by various spectroscopic instruments (e.g., HDS, SDSS, etc.; 'catalog'
directory).


-----------------------

## Directory structure

|
|-- inputs:             Input theoretical spectra calculated to fit the
|                       observed stars. 
|
|-- outputs:            Outputs of the codes
|
|-- pfs_calibstars :    Python codes to generate simulated PFS spectra with
|                       desired settings. 
|
|-- notebook:           Example scripts to run the codes
|
|-- catalogs:           Catalogs of observed stars and their estimated physical
|                       parameters.
|
|-- speclib:            Library of synthetic and observed spectra
|
|-- utility:            Utility data, e.g., photometric filters, isochrones
|
|-- database:           Database of simulated PFS spectra under the
|                       standard settings
|
