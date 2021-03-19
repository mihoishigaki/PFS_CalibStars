#
#     This file contains modules to define a datamodel for the 
#     standard star database 
#
#     A preliminary list of columns 
#
#     Colums for :
#
#     name                     STR            Name in PFS
#     ra                       DOUBLE         Right Ascension (degree)
#     dec                      DOUBLE         Declination (degree)
#     gmag                     FLOAT          g-band magnitude 
#     rmag                     FLOAT          r-band magnitude
#     imag                     FLOAT          i-band magnitude 
#     SDSS_name                STR            Name in SDSS               
#     PS_name                  STR            Name in PanSTAARS
#     GaiaDR2_name             STR            Name in Gaia DR2
#     spectype                 STR            Spectral type
#     spectype_flag            STR            Flag in spectral type
#     abundance_flag           STR            Flag in chemical abundance
#     teff                     FLOAT          Effective temperature (K)
#     teff_err                 FLOAT          Error in TEFF (K)
#     teff_ref                 STR            Reference of TEFF
#     logg                     FLOAT          Surface gravity (dex)
#     logg_err                 FLOAT          Error in LOGG
#     logg_ref                 STR            Reference of LOGG
#     feh                      FLOAT          [Fe/H] (dex)
#     feh_err                  FLOAT          Error in FEH
#     feh_ref                  STR            Reference of FEH
#     vr_helio                 FLOAT          Heliocentric radial velocity (km/s)
#     vr_helio_err             FLOAT          Error in VR_HELIO (km/s)
#     vr_helio_flag            STR            Flag in VR_HELIO (e.g. variation known)
#     empspec_src              STR            Source of the emperical spectrum
#     empspec_res              INT            Spectral resolution of the emperical spectrum
#     cfe                      FLOAT          [C/Fe]
#     cfe_err                  FLOAT          Error in CFE
#     mgfe                     FLOAT          [Mg/Fe]
#     mgfe_err                 FLOAT          Error in MGFE
#     sife                     FLOAT          [Si/Fe]
#     sife_err                 FLOAT          Error in SIFE
#     cafe                     FLOAT          [Ca/Fe]
#     cafe_err                 FLOAT          Error in CAFE
#     tife                     FLOAT          [Ti/Fe]
#     tife_err                 FLOAT          Error in TIFE
#     srfe                     FLOAT          [Sr/Fe]
#     srfe_err                 FLOAT          Error in SRFE
#     bafe                     FLOAT          [Ba/Fe]
#     bafe_err                 FLOAT          Error in BAFE
#
#


def creat_dataframe():

    header_list=['name','ra','dec','gmag','rmag','imag',\
                 'SDSS_name','PS_name','GaiaDR2_name',\
                 'spectype','spectype_flag','abundance_flag',\
                 'teff','teff_err','teff_ref','logg','logg_err','logg_ref',\
                 'feh','feh_err','feh_ref','vr_helio','vr_helio_err','vr_helio_flag',\
                 'empspec_src','empspec_res',\
                 'cfe','cfe_err',\
                 'mgfe','mgfe_err',\
                 'sife','sife_err',\
                 'cafe','cafe_err',\
                 'tife','tife_err',\
                 'srfe','srfe_err',\
                 'bafe','bafe_err']

    df=df.refindex(columns = header_list)
    

    return(df)


def read_elodie():






    

    return()


