import numpy as np
import pandas as pd


import os, sys
import logging


from multiprocessing import Pool
import multiprocessing as multi



################################################################################
#--- iSpec directory -------------------------------------------------------------
ispec_dir = os.path.dirname(os.path.abspath(__file__)) + "/../../iSpec/iSpec_v20201001"
sys.path.insert(0, os.path.abspath(ispec_dir))
import ispec
#########################################################



pfs_calibstars_dir = os.path.dirname(os.path.abspath(__file__)) + "/../"
sys.path.append(pfs_calibstars_dir)

import pfs_calibstars.utility as ut



def synthesize_spectrum(code="turbospectrum",teff=5771.0,logg=4.44,MH=0.0,alpha=0.0, \
	vsini=0.0,wave_base=380.,wave_top=400.,R=5000., fixed_abundances=None):
    #--- Synthesizing spectrum -----------------------------------------------------
    # Parameters
    #teff = 5771.0
    #logg = 4.44
    #MH = 0.00
    #alpha = ispec.determine_abundance_enchancements(MH)
	microturbulence_vel = ispec.estimate_vmic(teff, logg, MH) # 1.07
	macroturbulence = ispec.estimate_vmac(teff, logg, MH) # 4.21
    #vsini = 1.60 # Sun
	limb_darkening_coeff = 0.6
	resolution = R
	wave_step = 0.1

    # Wavelengths to synthesis
    #regions = ispec.read_segment_regions(ispec_dir + "/input/regions/fe_lines_segments.txt")
	regions = None
    #wave_base = 515.0 # Magnesium triplet region
    #wave_top = 525.0


    # Selected model amtosphere, linelist and solar abundances
    #model = ispec_dir + "/input/atmospheres/MARCS/"
	model = ispec_dir + "/input/atmospheres/MARCS.GES/"
    #model = ispec_dir + "/input/atmospheres/MARCS.APOGEE/"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.APOGEE/"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Castelli/"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Kurucz/"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Kirby/"


	if wave_top <= 1100.:
		atomic_linelist_file = ispec_dir + "/input/linelists/transitions/VALD.300_1100nm/atomic_lines.tsv"
	else:
		atomic_linelist_file = ispec_dir + "/input/linelists/transitions/VALD.1100_2400nm/atomic_lines.tsv"

    #atomic_linelist_file = ispec_dir + "/input/linelists/transitions/GESv5_atom_hfs_iso.420_920nm/atomic_lines.tsv"
    #atomic_linelist_file = ispec_dir + "/input/linelists/transitions/GESv5_atom_nohfs_noiso.420_920nm/atomic_lines.tsv"
	
	isotope_file = ispec_dir + "/input/isotopes/SPECTRUM.lst"

    # Load chemical information and linelist
	atomic_linelist = ispec.read_atomic_linelist(atomic_linelist_file, wave_base=wave_base, wave_top=wave_top)
	atomic_linelist = atomic_linelist[atomic_linelist['theoretical_depth'] >= 0.01] # Select lines that have some minimal contribution in the sun

	isotopes = ispec.read_isotope_data(isotope_file)

	#solar_abundances_file = ispec_dir + "/input/abundances/Grevesse.2007/stdatom.dat"
	#solar_abundances_file = ispec_dir + "/input/abundances/Asplund.2005/stdatom.dat"
	solar_abundances_file = ispec_dir + "/input/abundances/Asplund.2009/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Grevesse.1998/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Anders.1989/stdatom.dat"

    # Load model atmospheres
	modeled_layers_pack = ispec.load_modeled_layers_pack(model)
    # Load SPECTRUM abundances
	solar_abundances = ispec.read_solar_abundances(solar_abundances_file)

    ## Custom fixed abundances
    #fixed_abundances = ispec.create_free_abundances_structure(["C", "N", "O"], chemical_elements, solar_abundances)
    #fixed_abundances['Abund'] = [-3.49, -3.71, -3.54] # Abundances in SPECTRUM scale (i.e., x - 12.0 - 0.036) and in the same order ["C", "N", "O"]
    ## No fixed abundances
    #fixed_abundances = None

	# Validate parameters
	if not ispec.valid_atmosphere_target(modeled_layers_pack, {'teff':teff, 'logg':logg, 'MH':MH, 'alpha':alpha}):
		msg = "The specified effective temperature, gravity (log g) and metallicity [M/H] fall out of theatmospheric models."
		print(msg)

    # Prepare atmosphere model
	atmosphere_layers = ispec.interpolate_atmosphere_layers(modeled_layers_pack, {'teff':teff, 'logg':logg, 'MH':MH, 'alpha':alpha}, code=code)

    # Synthesis
	synth_spectrum = ispec.create_spectrum_structure(np.arange(wave_base, wave_top, wave_step))
	synth_spectrum['flux'] = ispec.generate_spectrum(synth_spectrum['waveobs'], \
            atmosphere_layers, teff, logg, MH, alpha, atomic_linelist, isotopes, solar_abundances, \
            fixed_abundances, microturbulence_vel = microturbulence_vel, \
            macroturbulence=macroturbulence, vsini=vsini, limb_darkening_coeff=limb_darkening_coeff, \
            R=resolution, regions=regions, verbose=1,
            code=code)
    ##--- Save spectrum ------------------------------------------------------------
	logging.info("Saving spectrum...")
	synth_filename = "example_synth_%s.fits" % (code)
	ispec.write_spectrum(synth_spectrum, synth_filename)

	return



def run_synspec(starname, teff, logg, feh, alpha, outdir, mode = "MR", \
	arms = ['blue','red','nir']):
    
	## Custom abundance   e.g., to change C, N, O abundances to the desired values
	solar_abundances_file = ispec_dir + "/input/abundances/Asplund.2009/stdatom.dat"
	solar_abundances = ispec.read_solar_abundances(solar_abundances_file)
	chemical_elements_file = ispec_dir + "/input/abundances/chemical_elements_symbols.dat"
	chemical_elements = ispec.read_chemical_elements(chemical_elements_file)
	#fixed_abundances = ispec.create_free_abundances_structure(["C", "N", "O"], chemical_elements, solar_abundances)
    #fixed_abundances['Abund'] = [-3.49, -3.71, -3.54] # Abundances in SPECTRUM scale (i.e., x - 12.0 - 0.036) and in the same order ["C", "N", "O"]
	fixed_abundances = None

    ## Synthesize a spectrum with the code "turbospectrum" 
	
	
	spectrum_wave = []
	spectrum_flux = []


	## - Wavelength coverage defined at
    ##    https://pfs.ipmu.jp/research/parameters.html
    ## - Mergins of 10nm are applied if possible. 
    ## - Overlaps in wavelength between neighboring arms are removed. 
	
	for arm in arms:
		if arm != 'nir':

			if arm == 'blue':
				wmin = 370.	
				wmax = 650.	
				R = 2300

			elif arm == 'red':

				if mode == 'MR':
					wmin = 700.
					wmax = 895.
					R = 5000.
				else:
					wmin = 650.
					wmax = 970.
					R = 3000.

			synthesize_spectrum(code = "turbospectrum", teff = teff, logg = logg, MH = feh, \
				vsini = 0.0, wave_base = wmin, wave_top = wmax, R=R, fixed_abundances = fixed_abundances)
			spectrum_tmp = ispec.read_spectrum('outflux.txt')
			spectrum_wave = np.hstack((spectrum_wave, spectrum_tmp['waveobs']))
			spectrum_flux = np.hstack((spectrum_flux, spectrum_tmp['flux']))


		elif arm == 'nir':
				
			for i in range(0, 2):
				
				if i==0:
					if mode == "MR":
						wmin_tmp = 930.
					else:
						wmin_tmp = 970.
					wmax_tmp = 1100.
				else:
					wmin_tmp = 1100.
					wmax_tmp = 1270.

				R = 4300.
				synthesize_spectrum(code = "turbospectrum", teff = teff, logg = logg, MH = feh, \
					vsini = 0.0, wave_base = wmin_tmp, wave_top = wmax_tmp, R=R, fixed_abundances = fixed_abundances)
				spectrum_tmp = ispec.read_spectrum('outflux.txt')
				spectrum_wave = np.hstack((spectrum_wave, spectrum_tmp['waveobs']))
				spectrum_flux = np.hstack((spectrum_flux, spectrum_tmp['flux']))


	synth_spectrum = ispec.create_spectrum_structure(spectrum_wave)
	synth_spectrum['flux'] = spectrum_flux
	ispec.write_spectrum(synth_spectrum, outdir + '/' + starname + '.txt')

	
	return()


