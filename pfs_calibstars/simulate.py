#
#    This file contains modules to read/write synthetic spectra,
#    run the PFS exposure time calculator and the simulators. 
#    
#    This file contains following subroutines
#
#          - flux2ABmag: Convert flux (in Flambda, any unit) to desired AB
#            magnitude
#          - read_spectext: Read synthetic spectra 
#
#
#    
#

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os, sys, glob, shutil, logging, pathlib
import pandas as pd
import subprocess as sp
from os import path

from astropy.utils.data import get_pkg_data_filename
from astropy import units as u
from synphot import SpectralElement
from synphot import SourceSpectrum, units
from synphot.models import Empirical1D
from synphot import Observation, specio
from astropy.io import ascii,fits
from astropy.table import Table
from pathlib import Path

from scipy import ndimage
from scipy.ndimage import gaussian_filter1d
from scipy import signal
from scipy import interpolate
import scipy.optimize as so




from pfsspecsim import pfsetc
from pfsspecsim import pfsspec


import pfs_calibstars.utility as ut




from multiprocessing import Pool
import multiprocessing as multi





def produce_database(libID, multiproc = False):

    if libID == 'HDS':

        catalogfile = "../catalogs/"+libID+"/"+"catalog_test_"+libID+".csv"
        
        df = pd.read_csv(catalogfile)

        ids = df["HaKiDaSName"]


        if multiproc == True:
            p = Pool(multi.cpu_count())
            p.map(calc_simspec_HDS, ids)
            p.close()

        else:

            for id in ids:
                calc_simspec_HDS(id)

    return()



def calc_simspec_HDS(id):

    #id, libID = arg[0], arg[1]

    libID = "HDS"
    ObjID = id[7:] 

   
    


    # Initialize subdirectories in the 'database' directory
    initialize_directories(ObjID, libID)

    
    # Get a file for the synthetic spectrum
    synfile = get_synfile(ObjID, libID)
    
  



    # Read catalog of star's basic information
    starname, ra, dec, vmag, kmag, teff, logg, feh, rv = \
        get_stellarinfo(ObjID, libID)
    
    #catalogfile = "../catalogs/"+libID+"/"+"catalog_test_"+libID+".csv"
    #df = pd.read_csv(catalogfile)
    


    
    # Assign magnitudes consistent with the stellar parameters
    gabs, rabs, iabs, zabs, yabs = get_isochrone_info(teff, logg, feh)
    

    
    # A rough estimate of absolute V
    
    Vabs = gabs - 0.03 - 0.42 * (gabs - rabs)
    
    

    #mag = get_mag(ObjID, libID, band)
    
    dmod = vmag - Vabs
    
    gmag = gabs + dmod
    rmag = rabs + dmod
    imag = iabs + dmod
    zmag = zabs + dmod
    ymag = yabs + dmod


    
    # Convirt units of the spectrum from flamda to AB mag
    
    wave, flux = read_spectext(synfile)

    
    band = ["V"]
    wvs, mags = flux2ABmag(wave, flux, band[0], vmag) \
        # The input wave should be in nm
    
    
    # Save magnitudes to a file 
    aa = np.array([wvs, mags])    
    synfilename = (((synfile[:-4]).split('/'))[-1]).replace("WL_Flux_Error_", "")
    magfilename = (synfile[:-4].split('Synspec'))[0]+"ABmag/"+synfilename+"_%s_%.1f.txt"%(band[0],vmag)
    
    np.savetxt(magfilename,aa.T,fmt='%.3f %.3f')

    
    
    # Assign typical number of exposures (1 exposure = 15min)
    objtype = "standard"
    exptime, nexp = get_exptime(band[0], vmag, objtype)
    
    
    # Roop over different settings: 
    for res in ["LR", "MR"]: # spectral resolutionn
        
        
        for condition in ["Optimistic", "Pessimistic"]:
            
            setting_label, MP, field_angle, MTangle, TP = get_setting(condition)
            

            # output directory
            outdir="../database/" + libID + "/" + \
                ObjID + "/ETC_Outputs/" + res + "/" + setting_label + "/"


            
            # Run the ETC code 
            outsncname = run_etc(magfilename, ObjID, res, MP, field_angle, \
                                 MTangle, TP, exptime, nexp, setting_label, outdir)

            
            # Directory for simulated spectra
            simdir="../database/" + libID + "/" + ObjID + "/Simulated/" + res + "/" + setting_label
            
            
            # Run the simulator 
            outsimname, outsimname_fits_list \
            = run_simulator(magfilename, outsncname, TP, res, exptime, nexp, ObjID, \
                            setting_label, simdir, nreal=10)

            

            # Assign photometric data
            for outsimname_fits in outsimname_fits_list:
                
                hdul=fits.open(outsimname_fits, mode='update')

                hdul[3].data['fiberMag'] = [gmag, rmag, imag, zmag, ymag]
                hdul[3].header['RA'] = ra
                hdul[3].header['DEC'] = dec
                
                hdul.flush()
                hdul.close()
            
    return()




def read_spectext(textfile):

    #
    # The input "textfile" should contain wavelength (nm) and
    # flux (arbitrary unit) in the first two columns. 
    # A header line, if present, should be commented out by "#"
    #

    if path.exists(textfile)==False:
        print("The file does not exists!: ",textfile)
        sys.exit()
        
    wv,fl=np.loadtxt(textfile,usecols=(0,1),unpack=True)
    

    return(wv,fl)




def flux2ABmag(wave, flux, band, mag, plot = False):


    
    # "band" should be either "sdss_g" or "V"

    # The input "wave" should be in nm
    

    # Define bandpass:
    if band == "sdss_g":
        bp = SpectralElement.from_file('../utility/photometry/g.dat')
        magvega = -0.08
    elif band == "V":
        bp = SpectralElement.from_filter('johnson_v')
        magvega = 0.03

    # SDSS g-band magnitude of Vega
    #gmagvega=-0.08
    

    # Read the spectrum of Vega
    sp_vega = SourceSpectrum.from_vega()
    wv_vega = sp_vega.waveset
    fl_vega = sp_vega(wv_vega,flux_unit=units.FLAM)  

    ## Convolve with the bandpass
    obs_vega = Observation(sp_vega, bp)
    
    ## Integrated flux
    fluxtot_vega = obs_vega.integrate()
    
    
    # Read the synthetic spectrum
    sp = SourceSpectrum(Empirical1D, points = wave * 10., \
                        lookup_table=flux*units.FLAM)
    wv = sp.waveset
    fl = sp(wv,flux_unit=units.FLAM)

    ## Convolve with the bandpass
    obs = Observation(sp, bp, force='extrap')
    
    ## Integrated g-band flux
    fluxtot = obs.integrate()
    
    
    # Scaling factor to make the flux compatible with the desired magnitude
    dm = mag-magvega
    const = fluxtot_vega*10**(-0.4*dm)/fluxtot
 
    # Scale the original flux by const
    fl_scale = const*fl
 
    # Convert to ABmag
    fl_scale_mag = units.convert_flux(wv, fl_scale, out_flux_unit='abmag')

    
    sp_scaled_mag = SourceSpectrum(Empirical1D, points = wv, lookup_table = fl_scale_mag)

    #aa = np.array([wv[1:-1]/10.,sp_scaled_mag(wv,flux_unit=u.ABmag)[1:-1]])
   
    
    #fluxfilename = ((fluxfile[:-4]).split('/'))[-1]
    #magfilename = (fluxfile[:-4].split('Synspec'))[0]+"ABmag/"+fluxfilename+"_%s_%.1f.txt"%(band,mag)
    
    #np.savetxt(magfilename,aa.T,fmt='%.3f %.3f')

    # Plot
    if plot == True: 
        fig,ax=plt.subplots(2,1,sharex=True)
        ax[0].plot(wave * 10.,flux,linestyle="-",marker="")
        ax[1].plot(wv[1:-1],sp_scaled_mag(wv,flux_unit=u.ABmag)[1:-1],linestyle="-",marker="")
        ax[1].set_xlabel("Angstrom")
        ax[0].set_ylabel("$F_{\lambda}$")
        ax[1].set_ylabel("ABmag")
        ax[1].set_ylim(mag+3.,mag-2.0)
        #plt.show()

    
    return(wv[1:-1]/10., sp_scaled_mag(wv, flux_unit = u.ABmag)[1:-1])


def run_etc(input_name, ObjID, res, MP, field_angle, MTangle, TP, \
            exptime, nexp, setting_label, outdir):
    


    
    if res=='MR':
        mr_mode='Y'
    else:
        mr_mode='N'

    totalexptime=exptime*nexp
    
    #outID=outdir + (input_name[:-4].split("/"))[-1]+"_%2s_MP%.2f_FA%.3f_MTa%.1f_TP%.1f_Texp%.0f"%(res,MP,field_angle,MTangle,TP,totalexptime)


    outID=outdir + (input_name[:-4].split("/"))[-1]+"_Texp%.0f"%(totalexptime)
    
    outnoisename = outID + ".noise.dat"
    outsncname= outID + ".snc.dat"

   
    # If the output filename contains a blank space:
    #outnoisename = outnoisename.replace(' ', '\ ')
    #outsncname = outsncname.replace(' ', '\ ')

    
    etc = pfsetc.Etc()
    etc.set_param('EXP_TIME',exptime)
    etc.set_param('EXP_NUM', nexp)
    etc.set_param('MAG_FILE',input_name)
    etc.set_param('REFF',0.0)  # should be 0 for a point source
    etc.set_param('MOON_PHASE',MP)
    etc.set_param('MOON_TARGET_ANG',MTangle)
    etc.set_param('FIELD_ANG',field_angle)
    etc.set_param('OUTFILE_NOISE',outnoisename)
    etc.set_param('OUTFILE_SNC',outsncname)
    etc.set_param('OUTFILE_SNL','-')
    etc.set_param('OUTFILE_OII','-')
    etc.set_param('MR_MODE',mr_mode)
    etc.set_param('degrade',TP)

    if os.path.exists(outsncname):
        print("The SNC file already exists")
        return(outsncname)
    else:
        etc.run()
    


    # Plot
    #w,snc=np.loadtxt(outsncname,usecols=(2,3),unpack=True)
    #plt.plot(w,snc)
    #plt.show()
    
    
    return(outsncname)


def run_simulator(input_name, sncname, TP, res, exptime, \
                  nexp, ObjID, setting_label, simdir, nreal = 1):
    

    #etcoutname=input_name+'_TP%3.1f'%TP+'_'+res+'_Nexp%02.0f'%num_exp
    #mag='abmag/'+input_name+'_abmag.dat'

      
    sim_name=((sncname[:-4]).split("/"))[-1]+'.sim'

    simlist = glob.glob(simdir + sim_name + "*")
    if len(simlist)!=0:
        for simfile in simlist:
            os.remove(simfile)

    fitslist=glob.glob(simdir+"pfsObject*.fits")
    if len(fitslist)!=0:
        for fits in fitslist:
            os.remove(fits)


    
    sim = pfsspec.Pfsspec()
    sim.set_param('outDir', simdir)
    sim.set_param('etcFile', sncname)
    sim.set_param('MAG_FILE', input_name)
    sim.set_param('EXP_NUM',nexp)
    sim.set_param('asciiTable', sim_name+"_"+ObjID)
    sim.set_param('nrealize',nreal)

    sim.make_sim_spec()

    fitsfiles = \
        glob.glob(simdir+"/pfsObject-000-00000-0,0-0000000000000*.fits")



    #sim_name_fits_list = ()

    #for fitsfile in fitsfiles:


    #    print("fitsfile = ",fitsfile)

        
        #sim_name_short=sim_name.replace("WL_Flux_Error_","")
        #sim_name_short=sim_name_short.replace(".snc.sim","")
        #sim_name_short=sim_name_short.replace("-","M")

        #sim_name_fits="-".join((fitsfile.split("-"))[0:4])+\
        #    "-%016i"%(10000+np.int(ObjID))+\
        #    "-"+"-".join((fitsfile.split("-"))[-2:])
        #shutil.copyfile(fitsfile, sim_name_fits)

        #sim_name_fits_list = np.append(sim_name_fits_list, sim_name_fits)
        
    
    return(simdir+"/"+sim_name, fitsfiles)





def plot_spec(files, w_min, w_max, synthfiles, outdir, \
              synth_norm_filepath = "../speclib/Evan_normalized_synspec/", inset = False, \
              textlabel = "", ymin = 0.3, ymax = 1.25):

    cmap = plt.get_cmap("tab10")

    plt.rcParams["font.size"] = 18
    
    fig = plt.figure(figsize = (11, 4.5))
    ax = fig.add_subplot(111)
    
    if inset == True:
        axins = inset_axes(ax, width="25%", height="45%", loc = 4)
    
    for i,file in enumerate(files):

        #teff = ((((file.split('/'))[10]).split("_"))[4])[4:]
        #logg = ((((file.split('/'))[10]).split("_"))[5])[4:]
        #mh = ((((file.split('/'))[10]).split("_"))[6])[2:]
        
        
        w0, f0, ferr0 = np.loadtxt(file,usecols=(0, 1, 2),unpack=True)

        filt = ((((w0 > 385.) & (w0 < 650.)) | \
            ((w0 > 700.) & (w0 < 900.)) | \
            ((w0 > 950.) & (w0 < 1250.))) & (np.isnan(f0) == False)) & \
            ((w0 >= w_min) & (w0 <= w_max)) & \
            ((w0 < 847.8023) | ((w0 > 851.8023) & (w0 < 852.2091)) | \
             ((w0 > 856.2091) & (w0 < 864.2141)) | (w0 > 868.2141))

        w = w0[filt]
        f = f0[filt]

        f_norm, ferr_norm = normalize_spec(w, f, w0, f0, ferr0)
 
        
        # Put labels for absorption line features
        catw =  849.8023
        if np.min(w0) < catw and np.max(w0) > catw:
            wcens, labs, cols = get_linelist_red()
            
            for kk, wcen in enumerate(wcens):
                #x = [wcen-0.1, wcen+0.1, wcen+0.1, wcen-0.1]
                #y = [0.0, 0.0, 2.0, 2.2]
                x = [wcen, wcen]
                y = [0.0, 2.0]
                ax.plot(x, y, linestyle = "-", linewidth = 2, \
                        color = cmap(cols[kk]), \
                        alpha = 0.5)
                if wcen == 851.5108:
                    continue
                ax.text(wcen, ymax - (ymax - ymin)*0.03, \
                        labs[kk], \
                        rotation = 90., ha = 'center', va = 'top', \
                        fontsize = 15)


        #filt = (w>855.) &  (w<865)
        #std = np.std(f_norm[filt])
        #sn = 1.0/std

        a = np.zeros(4)
        a.fill(0.25)
        f_conv = np.convolve(f_norm, a, 'same')

        #f_conv = gaussian_filter1d(f_norm, 0.5)

        
        ax.plot(w0, f_conv, linestyle = "-", marker = '', \
                color = 'k', alpha = 0.8)
        
        if inset == True:
            #x = [851.4072, 851.4072]
            x = [880.6756, 880.6756]
            y = [0, 2]
            axins.plot(x, y, linestyle = '-', color = cmap(3), \
                       linewidth = 2, alpha = 0.5)
            #x = [851.5108, 851.5108]
            #axins.plot(x, y, linestyle = '-', color = cmap(2), \
            #           linewidth = 2, alpha = 0.5)

            

                    
            axins.plot(w0, f_norm, color = 'k', linestyle = "-", \
                       marker = '', alpha = 0.8)


        
        # Get synthetic spectra
        for j, synthfile in enumerate(synthfiles):


            if j==1:
                ls = '--'
            else:
                ls = '--'


            w0_s, f_s = np.loadtxt(synth_norm_filepath + synthfile, \
                                   delimiter = ',', usecols = (0, 1), \
                                   unpack = True)
            w_s = w0_s/10.

            
            ax.plot(w_s, f_s, linestyle = ls, color = cmap(0), marker = "", linewidth = 1)
            if inset == True:
                axins.plot(w_s, f_s, linestyle = ls, color = cmap(0), marker = "", linewidth = 1)
        
        
        #ax.set_xlim(w_min, w_max)
        #ax.set_ylim(ymin, ymax)

        #ax.set_ylabel("Flux")
    if w_min>840. and w_max < 900.:
        ax.text(w_min + (w_max - w_min)*0.02, ymin + (ymax-ymin)*0.05, \
                textlabel, ha = 'left')



    ax.tick_params(labelleft=True, labelbottom=True)

    ax.set(xlim = (w_min, w_max), ylim = (ymin, ymax))
    ax.set_xlabel("Wavelength [nm]")
    #ax.set_rasterized(True)

    if inset == True:
        axins.set(xlim = (879.6756, 881.6756), ylim = (0.75, 1.1))
        axins.tick_params(labelleft=False, labelbottom=False)
        #axins.set_rasterized(True)
        
    outfig=(file[:-13].split("/"))[-1]


    plt.savefig(outdir + "%s_%.0f_%.0f.sim.png"%(outfig, w_min, w_max))
    #plt.show()



def plot_simulated_spec(ObjIDs, nexps, bands, mags, textlabels, \
                        synthfiles, output_rootdir, insets, libID = "HDS"):

    cmap = plt.get_cmap("tab10")

    plt.rcParams["font.size"] = 18
    

    # Define observational setting
    res = "MR"
    exptime = 900.
    setting = "Optimistic"
    setting_label,MP,field_angle,MTangle,TP = get_setting(setting)


    
    npanels = len(ObjIDs)

    fig, ax = plt.subplots(npanels, 1, figsize = (11, 8.), sharex = True)


    for i, ObjID in enumerate(ObjIDs):
    
        outdir = output_rootdir + "/" + libID + "/" + ObjID + "/" + \
            res + "/" + setting_label + "/"
        if os.path.isdir(outdir) == False:
            Path(outdir).mkdir(parents=True, exist_ok=True)
                    
        
    
        synfile = get_synfile(ObjID, libID, Vrad0flag = True)

    
        wave, flux = np.loadtxt(synfile, usecols = (0, 1), unpack = True)
    
        wvs, ABmags = flux2ABmag(wave, flux, bands[i], mags[i])


        # Save magnitudes to a file 
        aa = np.array([wvs, ABmags])
    
        synfilename = (((synfile[:-4]).split('/'))[-1]).replace("WL_Flux_Error_", "")
        magfilename = outdir + synfilename + "_%s_%.1f.txt"%(bands[i],mags[i])
    
        np.savetxt(magfilename,aa.T,fmt='%.3f %.3f')

        # Ca T region
        w_min = 847.
        w_max = 885.

        
        ymin = 0.2
        ymax = 1.3

    
        for nexp in nexps[i]: 
        
    
            # Run the ETC code
            
            outsncname = run_etc(magfilename, ObjID, res, MP, field_angle,\
                                 MTangle, TP, \
                             exptime, nexp, setting_label, outdir)
    
            # Run the simulator 
            outsimname,outsimname_fits = run_simulator(magfilename,outsncname,\
                                                       TP,res, exptime,nexp,\
                                                       ObjID, setting_label, outdir)




            simspecfile = outsimname + "_" + ObjID+".dat"
            

            w0, f0, ferr0 = np.loadtxt(simspecfile, \
                                       usecols=(0, 1, 2),unpack=True)

            filt = ((((w0 > 385.) & (w0 < 650.)) | \
            ((w0 > 700.) & (w0 < 900.)) | \
            ((w0 > 950.) & (w0 < 1250.))) & (np.isnan(f0) == False)) & \
            ((w0 >= w_min) & (w0 <= w_max)) & \
            ((w0 < 847.8023) | ((w0 > 851.8023) & (w0 < 852.2091)) | \
             ((w0 > 856.2091) & (w0 < 864.2141)) | (w0 > 868.2141))

            w = w0[filt]
            f = f0[filt]

            f_norm, ferr_norm = normalize_spec(w, f, w0, f0, ferr0)
 
        
            # Put labels for absorption line features
            catw =  849.8023
            if np.min(w0) < catw and np.max(w0) > catw:
                wcens, labs, cols = get_linelist_red()
            
                for kk, wcen in enumerate(wcens):
                    #x = [wcen-0.1, wcen+0.1, wcen+0.1, wcen-0.1]
                    #y = [0.0, 0.0, 2.0, 2.2]
                    x = [wcen, wcen]
                    y = [0.0, 2.0]
                    ax[i].plot(x, y, linestyle = "-", linewidth = 2, \
                               color = cmap(cols[kk]), \
                               alpha = 0.5)
                    if wcen == 851.5108:
                        continue
                    if i==0:
                        ax[i].text(wcen, ymax - (ymax - ymin)*0.03, \
                                   labs[kk], \
                                   rotation = 90., ha = 'center', va = 'top', \
                                   fontsize = 15)


            #a = np.zeros(4)
            #a.fill(0.25)
            #f_conv = np.convolve(f_norm, a, 'same')

            #f_conv = gaussian_filter1d(f_norm, 0.5)


      
            # Plot simulated spectrum
            ax[i].plot(w0, f_norm, linestyle = "-", marker = '', \
                       color = 'k', alpha = 0.5)


            

        # Plot synthetic spectra
        filt = (wave > 700.) & (wave < 900.)
        flux_norm, flux_norm_err = normalize_spec(wave[filt], flux[filt], \
                                                      wave, flux, np.zeros_like(flux))
            

        ax[i].plot(wave, flux_norm, linestyle = ":", marker = '', color = cmap(0), alpha = 0.8)
            

      

        if insets[i] == True:
            axins = inset_axes(ax[i], width = "25%", height = "45%", loc = 4)

            axins.plot(w0, f_norm, linestyle = "-", marker = '', color = 'k', alpha = 0.5)

            for synthfile in synthfiles[i]:
                wsyn, fsyn = np.loadtxt(synthfile, usecols = (0, 1), unpack = True)
                filt = (wsyn > 700.) & (wsyn < 900.)
                fsyn_norm, fsyn_norm_err = normalize_spec(wsyn[filt], fsyn[filt], \
                                                          wsyn, fsyn, np.zeros_like(fsyn))
                axins.plot(wsyn, fsyn_norm, color = 'b', linestyle = ":", marker = '', \
                           alpha = 0.8)


            axins.set(xlim = (879.7, 881.7), ylim = (0.7, 1.05))
            axins.tick_params(labelleft = False, labelbottom = False)


            
        ax[i].text(w_min + (w_max - w_min)*0.02, ymin + (ymax-ymin)*0.05, \
                textlabels[i], ha = 'left')



        #ax[i].tick_params(labelleft=True, labelbottom=True)

        ax[i].set_xlim(w_min, w_max)
        ax[i].set_ylim(ymin, ymax)

        
    ax[1].set_xlabel("Wavelength [nm]")

        
    outfig="Simulated_PFSspec.pdf"

    plt.subplots_adjust(hspace=0.02)
    plt.savefig(outfig)
    #plt.show()




    
    return()


def simulate_many_spectra(ObjID, libID, nexps, band, mag, \
                          nreal, output_rootdir,  \
                          setting = "Optimistic", \
                          write_h5 = False, plot = False):


    # Get information about the target

    ## Information from the catalog
    starname, ra, dec, vmag, kmag, teff, logg, feh, rv = \
        get_stellarinfo(ObjID, libID)
    
    ## Assign magnitudes consistent with the stellar parameters
    gabs, rabs, iabs, zabs, yabs = get_isochrone_info(teff, logg, feh)
    
    ## A rough estimate of absolute V
    #Vabs = gabs - 0.03 - 0.42 * (gabs - rabs)
    
    ## Calculate distance modulus
    dmod = mag - gabs

    ## Calculate appearent magnitudes 
    gmag = gabs + dmod
    rmag = rabs + dmod
    imag = iabs + dmod
    zmag = zabs + dmod
    ymag = yabs + dmod


    

    # Get information aboout observational setting 
    
    res = "MR"
    exptime = 900.
       
    setting_label, MP, field_angle, MTangle, TP = get_setting(setting)


    outdir = output_rootdir + "/" + libID + "/" + ObjID + "/" + \
        res + "/" + setting_label + "/"
    
    if os.path.isdir(outdir) == False:
        Path(outdir).mkdir(parents=True, exist_ok=True)
                    
        

    synfile = get_synfile(ObjID, libID, Vrad0flag = True)

    wave, flux = np.loadtxt(synfile, usecols = (0, 1), unpack = True)
    
    wvs, ABmags = flux2ABmag(wave, flux, band, mag)


    # Save magnitudes to a file 
    aa = np.array([wvs, ABmags])

    
    synfilename = (((synfile[:-4]).split('/'))[-1]).replace("WL_Flux_Error_", "")


    magfilename = outdir + synfilename + "_%s_%.1f.txt"%(band,mag)

    
    np.savetxt(magfilename,aa.T,fmt='%.3f %.3f')


    snrs = np.zeros(len(nexps))

    #simh5files = ()
    
    for i, nexp in enumerate(nexps): 
        
    
        # Run the ETC code 
        outsncname = run_etc(magfilename,ObjID,res,MP,field_angle,MTangle, \
                             TP,exptime,nexp, setting_label, outdir)

        snrs[i] = calc_sn(outsncname, res)
        
        # Run the simulator
        
        ## Delete all existing pfsObject*.fits file
        fitsfilelist = glob.glob(outdir + "pfsObject*.fits")
        for file in fitsfilelist:
            os.remove(file)

        
        outsimname, outsimname_fits = run_simulator(magfilename,outsncname,TP, \
                                                   res,exptime, nexp, ObjID, \
                                                   setting_label, outdir, \
                                                   nreal = nreal)
        
        dirname = outsimname.replace(".snc.sim","") + "_pfsObject"
        if os.path.isdir(dirname):
            shutil.rmtree(dirname)
        os.mkdir(dirname)
            
        #fitsfiles = \
            #glob.glob(outdir+"pfsObject-000-00000-0,0-00000000000000*.fits")
        

        #for file in fitsfiles:
        
        for file in outsimname_fits:
            hdul=fits.open(file, mode='update')

            hdul[3].data['fiberMag'] = [gmag, rmag, imag, zmag, ymag]
            hdul[3].header['RA'] = ra
            hdul[3].header['DEC'] = dec
                
            hdul.flush()
            hdul.close()
            
            shutil.move(file, dirname + "/" + (file.split("/"))[-1])
        
        

        if write_h5==True:
            
            textfilelist = glob.glob(outsimname + "*.dat")

            #ids = ()
            #for textfile in textfilelist:
            #    id = (textfile.split(".")[-2])
            #    ids = np.append(ids, id)


            datadictionary = {"snr": snrs[i], "w_unit": "nm", \
                              "f_unit": "10^-17erg/s/cm^2/A", \
                              "e_unit": "10^-17erg/s/cm^2/A"}

            ut.write_to_hdf5(textfilelist, \
                             datadictionary, outsimname + ".h5")

            #print("h5file = ", h5filename)
        
            #simh5files = np.append(simh5files, h5filename)


        if plot==True:
            #fig, ax = plt.subplots(1, 1)
            for k, textfile in enumerate(textfilelist):
                w0, f0, ferr0 = \
                    np.loadtxt(textfile, usecols = (0, 1, 2), unpack = True)

                filt = (w0 > 650.) & (w0 < 900)
                w = w0[filt]
                f = f0[filt]
                ferr = ferr0[filt]

                f_norm, ferr_norm = normalize_spec(w, f, w0, f0, ferr0)
                ax.plot(w0, f_norm)


                if k == 3:
                    break
            
        
    return(outdir, snrs)





def calc_exptime_sn(ObjIDs, types, nexps, band, mags, mode, \
                    output_rootdir, caption = '', xlims = [0., 400]):
    
    libID = "HDS"
    
    
    fig, ax = plt.subplots(1, 1, figsize = (11, 8))


    # Plot S/N constrainits 
    plt.rcParams["font.size"] = 17
    textoffset = 1.
    ax.plot([0.0, xlims[1]], [10.0, 10.0], linestyle = '-', linewidth = 4, \
            color = 'lightgrey')
    ax.text(xlims[1], 10. + textoffset, r"$V_{RV}$", ha = 'right')
    ax.text(xlims[1], 10., r"($\sigma=3-10$km/s, $\sim 3\times 10^5$)", \
            ha = 'right', va = 'top', fontsize = 16)

    ax.plot([0.0, xlims[1]], [20.0, 20.0], linestyle = '-', linewidth = 4, \
            color = 'lightgrey')
    ax.text(xlims[1], 20. + textoffset, r"[Fe/H]", ha = 'right')
    ax.text(xlims[1], 20., r"($\sigma=0.2$dex, $\sim 2\times 10^5$)", \
            ha = 'right', va = 'top', fontsize =16)
    

    ax.plot([0.0, xlims[1]], [60.0, 60.0], linestyle = '-', linewidth = 4, \
            color = 'lightgrey')
    ax.text(xlims[1], 60. + textoffset, r"[$\alpha$/Fe]", \
            ha = 'right')
    ax.text(xlims[1], 60., r"($\sigma=0.2$dex, $\sim 1\times 10^5$)", \
            ha = 'right', va = 'top', fontsize = 16)


    plt.rcParams["font.size"] = 18
    
    for k, ObjID in enumerate(ObjIDs):

        if k==0:
            cmap = plt.get_cmap("Oranges")
            ls = "-"
    
        else:
            cmap = plt.get_cmap("Greens")
            ls = "--"


        outdir = output_rootdir + "/" + libID + "/" + ObjID + "/"
        if os.path.isdir(outdir) == False:
            Path(outdir).mkdir(parents=True, exist_ok=True)
                    
        
        synfile = get_synfile(ObjID, libID, Vrad0flag = False)
 
        snr_meds = np.zeros_like(nexps)



        for j, mag in enumerate(mags[k,:]):

            
            wave, flux = np.loadtxt(synfile, usecols = (0, 1), unpack = True)
            
            wvs, ABmags = flux2ABmag(wave, flux, band, mag)

            # Save magnitudes to a file 
            aa = np.array([wvs, ABmags])
    
            synfilename = (((synfile[:-4]).split('/'))[-1]).replace('WL_Flux_Error_','')
            magfilename = outdir + synfilename+"_%s_%.1f.txt"%(band, mag)
    
            np.savetxt(magfilename,aa.T,fmt='%.3f %.3f')



            
            for i, nexp in enumerate(nexps): 
        
                res = "MR"
                exptime = 900.
       
                setting_label, MP, field_angle, MTangle, TP = \
                    get_setting("Optimistic")
    
                # Run the ETC code
                
               
                outsncname = run_etc(magfilename,ObjID,res,MP,field_angle,MTangle, \
                                     TP,exptime,nexp, setting_label, outdir)

                if types[k] == "K-giant":
                    outsncname_giant = outsncname
                elif types[k] == "G-dwarf":
                    outsncname_dwarf = outsncname
                
            

                wcontmin = 780.
                wcontmax = 820.
                snr_meds[i] = calc_sn(outsncname, res)
            expt_min = nexps * 15.


            ax.plot(expt_min, snr_meds, marker = 'o', linestyle = ls, linewidth = 2, color = cmap(1.0-np.float(j)/4), label = r"%s, $g = %.1f$"%(types[k], mag))



    ax.set_xlim(xlims[0], xlims[1])
    ax.set_xlabel("Total exposure time [min]")
    ax.set_ylabel("S/N (per resolution element)")

    print(expt_min)
    plt.xticks(expt_min, [15, 30, 60, 180, 300, 360])
    plt.legend(prop = {"size": 12})

    plt.savefig(output_rootdir + "/" + libID + "/plot/Exptime_SN.png")

    
    # caption
    #if caption == '':
    #    captiontext = write_caption(outsncname_giant, outsncname_dwarf)


        
    return()


def calc_sn(sncfile, res, wcontmin = 780., wcontmax = 820.):

    # Calc SN
    w, snr = np.loadtxt(sncfile, usecols = (2, 3), unpack = True)

    w_cont_filt = (w >  wcontmin) & (w > wcontmax)

    snr_med = np.median(snr[w_cont_filt])

    if res == "MR":
        angst_per_pix = 0.4
        res = 1.6
    elif res == "LR":
        angst_per_pis = 0.9
        res = 2.7
                    
    snr_per_reselem_factor = np.sqrt(res/angst_per_pix)
                
    snr_med_reselem = snr_med * snr_per_reselem_factor

    return(snr_med_reselem)



def write_caption(outsncname_giant, outsncname_dwarf):


    teffs = np.zeros(2)
    loggs = np.zeros(2)
    mhs = np.zeros(2)

    print(teffs)
    for i, outsncname in enumerate([outsncname_giant, outsncname_dwarf]):
    
        filename = (outsncname.split("/"))[-1]
        print(filename)
        teffs[i] = np.float(((filename.split("Teff"))[1])[0:4])
        loggs[i] = np.float(((filename.split("logg"))[1])[0:3])
        mhs[i] = np.float(((filename.split("MH"))[1])[0:4])

    f = open("../internal/FigureCaptions/caption_proposal.tex")
    texts = f.readlines()
    f.close()
        
    caption = texts %(teffs[0], loggs[0], mhs[0], teffs[1], loggs[1], mhs[1])         

    
    return(caption)




def initialize_directories(ObjID, libID, YorN = "Y"):

    Synspec_Path = "../inputs/" + libID + "/"
    
    #print("Initiallize all relevant directories for Object: "+ObjID+". Is this OK?  (Y or N): \n")
    #YorN=input().strip()

    if YorN == "N":
        sys.exit()


    # If you want to use a setting other than "Optimistic" or "Pessimistic",
    #   edit get_setting() function

    optimistic_setting,MP,field_angle,MTangle,TP = get_setting("Optimistic")
    pessimistic_setting,MP,field_angle,MTangle,TP = get_setting("Pessimistic")
    
    setting_labels=[optimistic_setting,pessimistic_setting]


    
    if os.path.isdir("../database/" + libID + "/" + ObjID):
    
        shutil.rmtree("../database/" + libID + "/" + ObjID)
    
    Path("../database/" + libID + "/" + ObjID).mkdir(parents=True, exist_ok=True)

    

    # Synthetic spectrum
    path_to_copy_of_synspec = "../database/" + libID + "/" + ObjID + "/Synspec/"
    
    Path(path_to_copy_of_synspec).mkdir(parents=True,exist_ok=True) 

    # Copy synthetic spectrum

    ## Get the objname for the ObjID
    ObjName=get_objname(ObjID, libID)

    
    if ObjName !="":
        synfiles=glob.glob(Synspec_Path+"*"+ObjName+"*")
        
        for f in synfiles:
            shutil.copy2(f, path_to_copy_of_synspec)
    else:
        print("Synthetic spectrum not found for ObjId:",ObjID)
        sys.exit()

    # ABmag
    
    Path("../database/" + libID + "/"+ObjID+"/ABmag").mkdir(parents=True,exist_ok=True)
        
    # ETC outputs

    for res in ["MR","LR"]:
    
        for setting_label in setting_labels:
    
            Path("../database/" + libID + "/" + ObjID + \
                 "/ETC_Outputs/"+res+"/"+setting_label).mkdir(parents=True,exist_ok=True)

    
    # Simulated spectra

    for res in ["MR","LR"]:
    
        for setting_label in setting_labels:
    
            Path("../database/" + libID + "/" + ObjID + \
                 "/Simulated/"+res+"/"+setting_label).mkdir(parents=True,exist_ok=True)
    
    return()


def get_objname(ObjID, libID):

    
    df=pd.read_csv("../catalogs/" + libID+"/catalog_HDS.csv")

    filt = df["HaKiDaSName"] == "HaKiDaS%05i"%(np.int(ObjID))

    name=(df["Name_1"][filt]).values

    if np.size(name)==1:
        print("Name of the object "+ObjID+" => "+name)
    else:
        print("Object name for "+ObjID+" not found.\n")
        sys.exit()
        
    return(name[0])
    

def get_mag(ObjID, libID, bands):

    if np.size(bands) == 1 and bands[0] == "V":

        df=pd.read_csv("../catalogs/"+ libID + "/catalog_HDS.csv")

    else:
    
        df=pd.read_csv("../catalogs/" + libID + "/catalog_HDSxPS1.csv")

    filt = df["HaKiDaSName"] == "HaKiDaS%05i"%(np.int(ObjID))


    mags = ()
    
    for band in bands:
    
        cname = band + "mag"
        
        mag=(df[cname][filt]).values

        print(mag)
        if np.size(mag)==1:
            print("Magnitude of the object %s => %.2f at %s"%(ObjID,mag,band))
        else:
            print("Magnitude for %s not found.\n"%(ObjID))
            sys.exit()

        mags = np.hstack((mags, mag[0]))
                
    return(mags)






def get_synfile(ObjID, libID, Vrad0flag = False):


    if Vrad0flag == True:
        synfiles=glob.glob("../database/"+libID+"/"+ObjID+"/Synspec/*Vrad0.0.txt")

    else:
        synfiles=glob.glob("../database/"+libID+"/"+ObjID+"/Synspec/*")

    
    if np.size(synfiles)==0:
        print("No synthetic spectra found for  "+libID+"/"+ObjID+". Please check!")
        sys.exit()
    
    return(synfiles[0])


def get_stellarinfo(ObjID,libID):

    

    if libID == "HDS":

         # Read catalog of star's basic information
        catalogfile = "../catalogs/"+libID+"/"+"catalog_test_"+libID+".csv"
        df = pd.read_csv(catalogfile)
        filt = df['HaKiDaSName'] == "HaKiDaS" + ObjID

        starname = df['Name_1'][filt]
        if len(starname) == 0:
            print("Information on %s not found in the catalog %s\n"%(ObjID, catalogfile))
            sys.exit()
        elif len(starname) > 1:
            print("More than one object found in the catalog %s\n"%(ObjID, catalogfile))
            sys.exit()

            
        ra = np.float(df['RA_1'][filt])
        dec = np.float(df['DEC_1'][filt])
        vmag = np.float(df['Vmag'][filt])
        kmag = np.float(df['Kmag_x'][filt])

        teff = np.float(df['Teff'][filt])
        logg = np.float(df['logg'][filt])
        feh = np.float(df['[Fe/H]'][filt])
        rv = np.float(df['rv_HDS'][filt])

        
    #synfilepath = get_synfile(ObjID, libID)
    #synfilename = (synfilepath.split('/'))[-1]
    #params = synfilename.split('_')
    #teff = np.float64((params[4])[-4:])
    #logg = np.float64((params[5])[-3:])
    #mh = np.float64((params[6])[-4:])
    #vlos = np.float64(((params[7])[:-4])[4:])
    
    
    
    return(starname, ra, dec, vmag, kmag, teff, logg, feh, rv)



def get_exptime(band, mag, objtype):

    # objtype should be either "standard" or "science"

    
    # - One mag brighter => 10**0.4 = 2.5 times more counts
    # 
    # - For a flux-limited case, S/N ~ sqrt(counts)
    #    => S/N is 1.6 times larger for a given exposure time 
    #
    # - Consider an object A gets S/N = 100 with exptime=60min
    #     if another object B is one magnitude brighter, the
    #      exptime to get S/N = 100 is:
    #        t*S/sqrt(t*S)=t0*S0/sqrt(t0*S0)
    #       sqrt(t)*1.6 = sqrt(t0)
    #       t/t0 = (1./1.6)**2
    #       t=0.39t0


    if objtype=="standard":
        
        exptime=60 # sec 

        if band=="V":
            if mag>13.0:
                nexp=20
            elif mag<=13.0 and mag>12.0:
                nexp=8
            elif mag<=12.0 and mag>11.0:
                nexp=3
            elif mag<=11.0:
                nexp=1

    if objtype=="science":
        
        exptime=900 # sec
    
        if band == "V" or band == "sdss_g":

            if mag > 23:
                nexp = 20
            elif mag<=23 and mag>22.:
                nexp=20
            elif mag<=22 and mag>21.:
                nexp=12
            elif mag<=21 and mag>20.:
                nexp = 8
            elif mag<=20.0 and mag>19.0:
                nexp=4
            elif mag<=19.0 and mag>18.0:
                nexp=2
            elif mag<=18:
                nexp=1
                
    return(exptime,nexp)



def get_setting(condition):

    if condition=="Optimistic":

        # Set parameters for the ETC
        MP=0.0 # New moon
        field_angle=0.0  # field center
        MTangle=60.  # Moon-target angluar distance (since MP=0.0, the value has no effect)
        TP=1.0  # throughput change, 1.0 for the original throughput as designed

    elif condition=="Pessimistic":

        # Set parameters for the ETC
        MP=0.25 # Quater moon
        field_angle=0.675  # field center
        MTangle=30.  # Moon-target angluar distance (since MP=0.0, the value has no effect)
        TP=0.8  # throughput change, 1.0 for the original throughput as designed

    label="MP%.2f_FA%.3f_MTangle%03i_TP%.2f"%(MP,field_angle,MTangle,TP)

    return(label,MP,field_angle,MTangle,TP)


def get_isochrone_info(teff, logg, feh):

        
    if feh < -1.75:
        MH = -2.00
    elif feh >= -1.75 and feh < -1.25:
        MH = -1.50
    elif feh >= -1.25 and feh < -0.75:
        MH = -1.00
    elif feh >= -0.75 and feh < -0.25:
        MH = -0.50
    elif feh >= -0.25 and feh < 0.25:
        MH = 0.00
    else:
        MH = 0.50
    
    isochronepath = "../utility/isochrones/PARSEC_Age10Gyrs_MH-2_0.5.dat"

    df = pd.read_table(isochronepath, comment = '#', delimiter = '\s+')


    diffs = np.zeros_like(df['logTe'])
    
    for i, lte in enumerate(df['logTe']):

        if df['MH'][i] != MH:
            diffs[i] = 999.999
            continue

        
        loggiso = np.float64(df['logg'][i])

        diffs[i] = (10**np.float64(lte) - teff)**2 + (loggiso - logg)**2


    filt = (diffs == np.min(diffs))

    gmag = df['gmag'][filt].values 
    rmag = df['rmag'][filt].values 
    imag = df['imag'][filt].values 
    zmag = df['zmag'][filt].values 
    ymag = df['ymag'][filt].values 


    return(gmag[0], rmag[0], imag[0], zmag[0], ymag[0])



def get_linelist_red():

    cmap = plt.get_cmap("tab10")
    
    wcens = ()

    labs = ()
    
    
    labs = np.hstack((labs, "Ca II"))
    wcens = np.hstack((wcens, 8498.023))


    labs = np.hstack((labs, "Fe I"))
    wcens = np.hstack((wcens, 8514.072))


    labs = np.hstack((labs, "Fe I"))
    wcens = np.hstack((wcens, 8515.108)) # blend

    labs = np.hstack((labs, "Ca II"))
    wcens = np.hstack((wcens, 8542.091))

    labs = np.hstack((labs, "Fe I"))
    wcens = np.hstack((wcens, 8582.258))

    labs = np.hstack((labs, "Fe I"))
    wcens = np.hstack((wcens, 8611.804))

    labs = np.hstack((labs, "Fe I"))
    wcens = np.hstack((wcens, 8621.601))

    labs = np.hstack((labs, "Ca II"))
    wcens = np.hstack((wcens, 8662.141))

    labs = np.hstack((labs, "Fe I"))
    wcens = np.hstack((wcens, 8688.626))

    labs = np.hstack((labs, "H I"))
    wcens = np.hstack((wcens, 8750.447))   # (dwarf only)


    labs = np.hstack((labs, "Mg I"))
    wcens = np.hstack((wcens, 8806.756))

    labs = np.hstack((labs, "Fe I"))
    wcens = np.hstack((wcens, 8824.221))


    wcens = wcens/10.

    cols = np.zeros(np.size(labs), dtype = np.int)

    for kk, lab in enumerate(labs):


        if lab == "Ca II":
            cols[kk] = 1
        elif lab == "Fe I":
            cols[kk] = 2
        elif lab == "Mg I":
            cols[kk] = 3
        elif lab == "H I":
            cols[kk] = 9
        else:
            cols[kk] = 7


        
    return(wcens, labs, cols)



def normalize_spec(w, f, w0, f0, ferr0, \
                   medfilt_widh = 10, maxfilt_width = 5, order = 4):

    # w, f : arrays used to define a continuum function
    # w0, f0, ferr0: target flux to be normalized

    
    # Median filter
    f_medfilt = ndimage.median_filter(f, size = 10)
    
    
    # Maximum filter
    f_maxfilt = ndimage.maximum_filter(f_medfilt, size = 5) 
    
    
    
    z = np.polyfit(w, f_maxfilt, order)
    p = np.poly1d(z)
    f_norm = f0 / p(w0)
    ferr_norm = ferr0 / p(w0)
    
    return(f_norm, ferr_norm)

    
