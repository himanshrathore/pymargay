#PYthon based MAnga Resolved GAlaxY analysis

#necessary imports

import numpy as np
import math
from astropy.io import fits
from astropy.table import Table
from astropy.cosmology import Planck13 as cosmo
from astropy import units as u
import matplotlib.pyplot as plt
from photutils import EllipticalAnnulus
from photutils import EllipticalAperture
from photutils import aperture_photometry
import pandas as pd

#######################################################################################################################

#Required catalogs

#global properties table, cross of Pipe3D + Salim
Pipe3D_Salim = Table.read('address1')
#pymorph r-band photometric data
totalrband = Table.read('address2')
#Postion angle correction to pymorph measurements
SPA_corr = Table.read('address3')
#MaNGA target catalog (Wake et al. 2017)
MaNGA_target = Table.read('address4')

########################################################################################################################

#functions for obtaining general properties of a MaNGA galaxy

def get_fits_hdul(MaNGA_ID):
    #retrieves the hdul of the fits file of our object
    
    hdul = fits.open('address5') #address5 should be a function of MaNGA_ID    
    return hdul

def get_coordinates(MaNGA_ID):
    '''
    Gets the equitorial coordinates of the given galaxy from the primary header of the datacube.
    input: Manga_ID of the galaxy
    returns: 2-tuple (ra,dec) of the galaxy, extracted from the primary header of the datacube
    '''
    hdul = get_fits_hdul(MaNGA_ID)
    ra = hdul[0].header['OBJRA']
    dec = hdul[0].header['OBJDEC']
    return ra, dec

def get_redshift(MaNGA_ID):
    #retrieves the redshift mentioned in the Salim catalog for our object
    #Redshifts mentioned in the Salim Catalog are from SDSS, whereas Pipe3D redshifts are those derived by Pipe3D
    
    ID = 'manga' + '-' + MaNGA_ID
    IDs = Pipe3D_Salim['mangaid']
    loc = np.where(IDs == ID)[0][0]
    z = Pipe3D_Salim['col8'][loc]
    return z

def get_luminosity_distance(MaNGA_ID):
    #returns the luminosity distance (in cm) given the MaNGA_ID
    
    z = get_redshift(MaNGA_ID)
    d = cosmo.luminosity_distance(z).to(u.cm)
    d_cm = d.value
    return d_cm

def get_spaxel_area(MaNGA_ID):
    #returns the proper area of each spaxel (in kpc^2) from the angular size of each spaxel mentioned in the datacube header
    
    hdul = get_fits_hdul(MaNGA_ID)
    z = get_redshift(MaNGA_ID)
    spx_scale = hdul[0].header['CD2_2']*3600 #angular scale of each spaxel in arcsec
    spaxel_size = cosmo.kpc_proper_per_arcmin(z)*spx_scale/60 #kiloparsecs corresponding to a spaxel
    spaxel_area = (spaxel_size.value)**2 #spaxel_area in kpc^2
    return spaxel_area

############################################################################################################################

#BPT related functions

def get_bpt_lines_map(MaNGA_ID, snr_cutoff = 3): #Fluxes of the 4 BPT lines
    '''
    inputs: MaNGA_ID of the galaxy, S/N cutoff to be applied to the 4 BPT lines
    returns: The flux maps of the 4 BPT lines, and the corresponding error maps. No dust correction is applied, since these lines have similar wavelengths and we take their ratios. So, dust correction will not play a significant role
    '''
    
    hdul = get_fits_hdul(MaNGA_ID)
    emline_maps = np.copy(hdul[3].data)
    
    #Retrieve the emission line maps

    oiii = np.copy(emline_maps[26,:,:]) #[OIII]5007
    nii = np.copy(emline_maps[46,:,:]) #[NII]6584
    ha = np.copy(emline_maps[45,:,:]) #6562.68
    hb = np.copy(emline_maps[28,:,:]) #4861.32

    #The corresponding noise maps
    oiii_noise = np.copy(emline_maps[26+228,:,:]) #[OIII]5007
    nii_noise = np.copy(emline_maps[46+228,:,:]) #[NII]6584
    ha_noise = np.copy(emline_maps[45+228,:,:]) #6562.68
    hb_noise = np.copy(emline_maps[28+228,:,:]) #4861.32
    
    #To remove negative values and low S/N pixels
    oiii[oiii<=0]= math.nan
    nii[nii<=0]= math.nan
    ha[ha<=0]= math.nan
    hb[hb<=0]= math.nan

    oiii_noise[oiii<=0]= math.nan
    nii_noise[nii<=0]= math.nan
    ha_noise[ha<=0]= math.nan
    hb_noise[hb<=0]= math.nan

    #Set the S/N cutoff
    
    indices = (emline_maps[26,:,:]/emline_maps[26+228,:,:] > snr_cutoff)*(emline_maps[46,:,:]/emline_maps[46+228,:,:] > snr_cutoff)*(emline_maps[45,:,:]/emline_maps[45+228,:,:] > snr_cutoff)*(emline_maps[28,:,:]/emline_maps[28+228,:,:] > snr_cutoff)
    
    oiii[indices == False] = math.nan
    nii[indices == False] = math.nan
    ha[indices == False] = math.nan
    hb[indices == False] = math.nan
    
    oiii_noise[indices == False] = math.nan
    nii_noise[indices == False] = math.nan
    ha_noise[indices == False] = math.nan
    hb_noise[indices == False] = math.nan
    
    return oiii, oiii_noise, nii, nii_noise, ha, ha_noise, hb, hb_noise

def line(x, a, b): #Kewley et al. 2006 line
    return 0.61/(x - a) + b

def liner(x): #Cid Fernandes et al. 2010 line
    return 1.01*x + 0.48

def get_BPT_tags(MaNGA_ID, snr_cutoff = 3):
    #Returns a boolean map for each of the BPT tags
    #Inputs: MaNGA_ID, and S/N cutoff for all the BPT lines (by default taken to be 3)
    
    oiii, oiii_noise, nii, nii_noise, ha, ha_noise, hb, hb_noise = get_bpt_lines_map(MaNGA_ID, snr_cutoff)
    
    #Taking ratios
    ob = np.log10(oiii/hb) #[OIII]/[Hb]
    na = np.log10(nii/ha) #[NII]/[Ha]
    
    SF_con1 = np.log10(oiii/hb) <= line(np.log10(nii/ha), 0.05, 1.3)
    SF_conii = np.log10(nii/ha) < 0.05 
    
    shock_con1 = np.log10(oiii/hb) > line(np.log10(nii/ha), 0.47, 1.19)
    shock_conii = np.log10(nii/ha) <= 0.47
    shock1 = np.logical_and(shock_con1, shock_conii)
    shock2 = np.log10(nii/ha) > 0.47
    shock = np.logical_or(shock1, shock2)
    
    comp1 = np.log10(oiii/hb) > line(np.log10(nii/ha), 0.05, 1.3)
    comp2 = np.log10(oiii/hb) <= line(np.log10(nii/ha), 0.47, 1.19)

    liner_con1 = np.log10(oiii/hb) <= liner(np.log10(nii/ha)) 
    seyfert_con1 = np.log10(oiii/hb) > liner(np.log10(nii/ha))
    
    SF = np.logical_and(SF_con1, SF_conii)
    comp = np.logical_and(comp1, comp2)
    liners = np.logical_and(shock, liner_con1)
    seyfert = np.logical_and(shock, seyfert_con1)
    
    return SF, comp, liners, seyfert

#########################################################################################################################

#functions for obtaining resolved maps

def get_Ha_map(MaNGA_ID, snr_cutoff = 3):
    #Returns the Ha flux map and the Ha noise flux map, with the required S/N ratio 
    #inputs: MaNGA_ID, Ha line S/N cutoff (by default taken to be 3)
    
    #H alpha map
    hdul = get_fits_hdul(MaNGA_ID)
    emline_maps = np.copy(hdul[3].data)
    
    ha = np.copy(emline_maps[45,:,:])
    ha[ha<=0] = math.nan

    #setting s/n cutoff

    ha_noise = np.copy(emline_maps[273,:,:])
    ha_noise[ha_noise<=0] = math.nan
    
    low_snr_indices = np.where(emline_maps[45,:,:]/emline_maps[273,:,:] < snr_cutoff)

    ha[low_snr_indices] = math.nan
    ha_noise[low_snr_indices] = math.nan
    
    return ha , ha_noise

def get_A_Ha_map(MaNGA_ID, snr_cutoff = 3):
    #returns the A_Ha map for the object (attenuation at Ha wavelength) using Balmer decrement
    #inputs: MaNGA_ID, emission line S/N cutoff (by default taken to be 3). This S/N cutoff is for all 4 BPT emission lines
    
    oiii, oiii_noise, nii, nii_noise, ha, ha_noise, hb, hb_noise = get_bpt_lines_map(MaNGA_ID, snr_cutoff)
    EBV = 1.97*np.log10((ha/hb)/2.86) # E(B-V) is the color excess
    
    A_Halpha = EBV * 3.33
    
    return A_Halpha

def apply_dust_correction(A_Ha, Ha):
    #returns the dust corrected Ha flux map
    #inputs: attenuation map of Ha, Ha flux map
    corr_ha = Ha * (10**(0.4*A_Ha)) #corrected halpha flux
    return corr_ha

def get_SFR_map(MaNGA_ID, snr_cutoff = 3):
    #returns the dust corrected SFR map (solar masses/yr) and the sigma SFR map (solar masses/yr/kpc^2) and the respective error maps. No dust correction is applied to the error maps.
    #input: MaNGA_ID, S/N cutoff on the 4 BPT emission lines (by default taken to be 3) 
    
    Ha, Ha_noise = get_Ha_map(MaNGA_ID, snr_cutoff)
    A_Ha = get_A_Ha_map(MaNGA_ID, snr_cutoff)
    corr_ha = apply_dust_correction(A_Ha, Ha)
    
    d_cm = get_luminosity_distance(MaNGA_ID)
    spx_area = get_spaxel_area(MaNGA_ID)
    
    Ha_luminosity = np.copy(corr_ha) * 4 * np.pi * (d_cm**2)*(10**(-16))
    Ha_noise_luminosity = np.copy(Ha_noise) * 4 * np.pi * (d_cm**2)*(10**(-16))
    
    #Star Formation Rate
    SFR = 7.9E-42 * Ha_luminosity #Kennicutt1998, with Saltpeter IMF
    SFR_error = 7.9E-42 * Ha_noise_luminosity
    
    #Star formation rate surface density
    sigma_sfr = SFR/spx_area #Sigma SFR
    sigma_sfr_error = SFR_error/spx_area
    
    return SFR, SFR_error, sigma_sfr, sigma_sfr_error
    
def get_Ha_lum_map(MaNGA_ID, snr_cutoff = 3):
    #input: MaNGA_ID, S/N cutoff on the 4 BPT lines (by default taken to be 3)
    #returns: dust corrected Ha luminosity and its error map. No dust correction is applied to the error map.
    
    Ha, Ha_noise = get_Ha_map(MaNGA_ID, snr_cutoff)
    A_Ha = get_A_Ha_map(MaNGA_ID, snr_cutoff)
    corr_ha = apply_dust_correction(A_Ha, Ha)
    
    d_cm = get_luminosity_distance(MaNGA_ID)
    
    Ha_luminosity = np.copy(corr_ha) * 4 * np.pi * (d_cm**2)*(10**(-16))
    Ha_noise_luminosity = np.copy(Ha_noise) * 4 * np.pi * (d_cm**2)*(10**(-16))
    
    return Ha_luminosity, Ha_noise_luminosity


def get_smd_map(MaNGA_ID):
    #returns the stellar mass surface density map
    
    hdul = get_fits_hdul(MaNGA_ID)
    spx_area = get_spaxel_area(MaNGA_ID)
    
    smd = np.copy(hdul[1].data[19]) #its already logarithmic
    smd = (10**smd)/spx_area #smd in Msun/Kpc^2
    
    return smd


def get_Vband_map(MaNGA_ID):
    #returns the V band image of the galaxy
    
    hdul = get_fits_hdul(MaNGA_ID)
    
    V = np.copy(hdul[1].data[0])
    
    return V


def get_gas_metallicity_map(MaNGA_ID, snr_cutoff = 3):
    #returns the gas metallicity map for the given object and the corresponding error map
    #Inputs: MaNGA_ID, and S/N cutoff for all the BPT lines (by default taken to be 3)
    
    oiii, oiii_noise, nii, nii_noise, ha, ha_noise, hb, hb_noise = get_bpt_lines_map(MaNGA_ID, snr_cutoff)
    
    #Defining ratios
    ob = np.log10(oiii/hb) #[OIII]/Hb
    an = np.log10(ha/nii) # Ha/[NII]
    O3N2 = ob + an
    
    #Calculating gas phase metallicity (Ellison et al. 2018, Marino et al. 2013)
    gm = 8.533 - (0.214*O3N2)
    
    #Calculating error
    gm_error = 0.214 * ((oiii_noise/oiii)**2 + (ha_noise/ha)**2 + (hb_noise/hb)**2 + (nii_noise/nii)**2)**0.5
    
    return gm, gm_error

def get_lum_wt_age_map(MaNGA_ID):
    #returns the luminosity weighted age map of the stellar population
    
    hdul = get_fits_hdul(MaNGA_ID)
    lum_age = np.copy(hdul[1].data[5])
    
    return lum_age

def get_mass_wt_age_map(MaNGA_ID):
    #returns the mass weighted age map of the stellar population
    
    hdul = get_fits_hdul(MaNGA_ID)
    mass_age = np.copy(hdul[1].data[6])
    
    return mass_age

def get_age_error_map(MaNGA_ID):
    #returns the error map of the age of the stellar population

    hdul = get_fits_hdul(MaNGA_ID)
    error_age = np.copy(hdul[1].data[7])
    
    return error_age

def get_lum_wt_sm_map(MaNGA_ID):
    #returns luminosity weighted metallicity of the stellar population. 12 is added to it in order to make the values positive definite 
    
    hdul = get_fits_hdul(MaNGA_ID)
    lum_sm = np.copy(hdul[1].data[8])
    
    return 12 + lum_sm

def get_mass_wt_sm_map(MaNGA_ID):
    #returns mass weighted metallicity of the stellar population. 12 is added to it in order to make the values positive definite 
    
    hdul = get_fits_hdul(MaNGA_ID)
    mass_sm = np.copy(hdul[1].data[9])
    
    return 12 + mass_sm

def get_sm_error_map(MaNGA_ID):
    #returns the metallicity error map of the stellar population. 12 is added to it in order to make the values positive definite     

    hdul = get_fits_hdul(MaNGA_ID)
    error_sm = np.copy(hdul[1].data[10])
    
    return 12 + error_sm

######################################################################################################################

#combined sigma sfr - sigma star arrays

def get_SFR_M_arrays(MaNGA_IDs, snr_cutoff = 3):
    #Returns the combined flattened arrays of SFR and Stellar Mass surface densities for the MaNGA_IDs

    my_sigma_sfr_flattened = np.array([])
    my_sigma_sfr_error_flattened = np.array([])
    my_smd_flattened = np.array([])

    for MaNGA_ID in MaNGA_IDs:

        SFR, SFR_error, sigma_sfr, sigma_sfr_error = get_SFR_map(MaNGA_ID, snr_cutoff)

        SF, comp, liners, seyfert = get_BPT_tags(MaNGA_ID, snr_cutoff)
        SF_comp = np.logical_or(SF, comp) #spaxels that are either SF or comp
        mask1 = np.logical_not(SF_comp) #mask is True when spaxel is neither SF nor comp

        temp_sigma_sfr = np.copy(sigma_sfr)
        temp_sigma_sfr[mask1] = math.nan

        temp_sigma_sfr_error = np.copy(sigma_sfr_error)
        temp_sigma_sfr_error[mask1] = math.nan

        smd = np.log10(get_smd_map(MaNGA_ID))
        temp_smd = np.copy(smd)
        temp_smd[mask1] = math.nan

        temp_sigma_sfr_flattened = temp_sigma_sfr.flatten()
        temp_sigma_sfr_error_flattened = temp_sigma_sfr_error.flatten()
        temp_smd_flattened = temp_smd.flatten()

        mask2 = np.logical_and(np.logical_and(temp_smd_flattened>0, temp_sigma_sfr_flattened>0), temp_sigma_sfr_error_flattened > 0) #removing nan and ambigous values

        my_sigma_sfr_flattened = np.append(my_sigma_sfr_flattened, temp_sigma_sfr_flattened[mask2])
        my_sigma_sfr_error_flattened = np.append(my_sigma_sfr_error_flattened, temp_sigma_sfr_error_flattened[mask2])
        my_smd_flattened = np.append(my_smd_flattened, temp_smd_flattened[mask2])
        #my_sigma_sfr_flattened and my_smd_flattened are of the same size, with no nan or -ve values

    return my_smd_flattened, my_sigma_sfr_flattened, my_sigma_sfr_error_flattened

def get_trimmed_SFR_M_arrays(MaNGA_IDs, snr_cutoff = 3, distance_range = np.array([0,2])):
    #Returns the combined flattened arrays of SFR and Stellar Mass surface densities for the MaNGA_IDs
    
    my_sigma_sfr_flattened = np.array([])
    my_sigma_sfr_error_flattened = np.array([])
    my_smd_flattened = np.array([])

    for MaNGA_ID in MaNGA_IDs:
        
        hdul = get_fits_hdul(MaNGA_ID)
        z = get_redshift(MaNGA_ID)
        spx_scale = hdul[0].header['CD2_2']*3600 #angular scale of each spaxel in arcsec
        kpc_per_spx = ((cosmo.kpc_proper_per_arcmin(z)).value)*(spx_scale/60)
        center = np.array([hdul[0].header['CRPIX1']-hdul[0].header['CRPIX3'], hdul[0].header['CRPIX2']-hdul[0].header['CRPIX3']])
        
        SFR, SFR_error, sigma_sfr, sigma_sfr_error = get_SFR_map(MaNGA_ID, snr_cutoff)

        SF, comp, liners, seyfert = get_BPT_tags(MaNGA_ID, snr_cutoff)
        SF_comp = np.logical_or(SF, comp) #spaxels that are either SF or comp
        mask1 = np.logical_not(SF_comp) #mask is True when spaxel is neither SF nor comp

        temp_sigma_sfr = np.copy(sigma_sfr)
        temp_sigma_sfr[mask1] = math.nan

        temp_sigma_sfr_error = np.copy(sigma_sfr_error)
        temp_sigma_sfr_error[mask1] = math.nan

        smd = np.log10(get_smd_map(MaNGA_ID))
        temp_smd = np.copy(smd)
        temp_smd[mask1] = math.nan
        
        spx_x = np.arange(sigma_sfr.shape[0])
        spx_y = np.arange(sigma_sfr.shape[1])
        
        mask_dist = np.zeros((sigma_sfr.shape[0], sigma_sfr.shape[1]), dtype = "bool")
        
        for i in spx_x:
            
            for j in spx_y:
                
                distance = (np.sqrt((i - center[0])**2 + (j - center[1])**2))*kpc_per_spx
                
                if(distance_range[1] != -1):
                    
                    if((distance >= distance_range[0]) and (distance <= distance_range[1])):
                        mask_dist[i][j] = True
                        
                elif(distance_range[1] == -1):
                    
                    if(distance >= distance_range[0]):
                        mask_dist[i][j] = True
        
        temp_sigma_sfr[mask_dist] = math.nan
        temp_sigma_sfr_error[mask_dist] = math.nan
        temp_smd[mask_dist] = math.nan

        temp_sigma_sfr_flattened = temp_sigma_sfr.flatten()
        temp_sigma_sfr_error_flattened = temp_sigma_sfr_error.flatten()
        temp_smd_flattened = temp_smd.flatten()

        mask2 = np.logical_and(np.logical_and(temp_smd_flattened>0, temp_sigma_sfr_flattened>0), temp_sigma_sfr_error_flattened > 0) #removing nan and ambigous values

        my_sigma_sfr_flattened = np.append(my_sigma_sfr_flattened, temp_sigma_sfr_flattened[mask2])
        my_sigma_sfr_error_flattened = np.append(my_sigma_sfr_error_flattened, temp_sigma_sfr_error_flattened[mask2])
        my_smd_flattened = np.append(my_smd_flattened, temp_smd_flattened[mask2])
        #my_sigma_sfr_flattened and my_smd_flattened are of the same size, with no nan or -ve values

    return my_smd_flattened, my_sigma_sfr_flattened, my_sigma_sfr_error_flattened

#combined gas metallicity - sigma star arrays

def get_Z_M_arrays(MaNGA_IDs, snr_cutoff = 3):
    #Returns the combined flattened arrays of gas phase metallicity (and corresponding errors) and stellar mass surface densities for the given MaNGA_IDs

    my_gm_flattened = np.array([])
    my_smd_flattened = np.array([])
    my_gm_error_flattened = np.array([])

    for MaNGA_ID in MaNGA_IDs:

        gm, gm_error = get_gas_metallicity_map(MaNGA_ID, snr_cutoff)

        SF, comp, liners, seyfert = get_BPT_tags(MaNGA_ID, snr_cutoff)
        SF_comp = np.logical_or(SF, comp) #spaxels that are either SF or comp
        mask1 = np.logical_not(SF_comp) #mask is True when spaxel is neither SF nor comp

        temp_gm = np.copy(gm)
        temp_gm[mask1] = math.nan
        
        temp_gm_error = np.copy(gm_error)
        temp_gm_error[mask1] = math.nan

        smd = np.log10(get_smd_map(MaNGA_ID))
        temp_smd = np.copy(smd)
        temp_smd[mask1] = math.nan

        temp_gm_flattened = temp_gm.flatten()
        temp_smd_flattened = temp_smd.flatten()
        temp_gm_error_flattened = temp_gm_error.flatten()

        mask2 = np.logical_and(np.logical_and(temp_smd_flattened>0, temp_gm_flattened>0), temp_gm_error_flattened > 0) #removing nan and ambigous values

        my_gm_flattened = np.append(my_gm_flattened, temp_gm_flattened[mask2])
        my_smd_flattened = np.append(my_smd_flattened, temp_smd_flattened[mask2])
        my_gm_error_flattened = np.append(my_gm_error_flattened, temp_gm_error_flattened[mask2])
        #my_gm_flattened, my_smd_flattened, my_gm_error_flattened are of the same size, with no nan or -ve values

    return my_smd_flattened, my_gm_flattened, my_gm_error_flattened

def get_trimmed_Z_M_arrays(MaNGA_IDs, snr_cutoff = 3, distance_range = np.array([0,2])):
    
    my_gm_flattened = np.array([])
    my_smd_flattened = np.array([])
    my_gm_error_flattened = np.array([])

    for MaNGA_ID in MaNGA_IDs:
        
        hdul = get_fits_hdul(MaNGA_ID)
        z = get_redshift(MaNGA_ID)
        spx_scale = hdul[0].header['CD2_2']*3600 #angular scale of each spaxel in arcsec
        kpc_per_spx = ((cosmo.kpc_proper_per_arcmin(z)).value)*(spx_scale/60)
        center = np.array([hdul[0].header['CRPIX1']-hdul[0].header['CRPIX3'], hdul[0].header['CRPIX2']-hdul[0].header['CRPIX3']])

        gm, gm_error = get_gas_metallicity_map(MaNGA_ID, snr_cutoff)

        SF, comp, liners, seyfert = get_BPT_tags(MaNGA_ID, snr_cutoff)
        SF_comp = np.logical_or(SF, comp) #spaxels that are either SF or comp
        mask1 = np.logical_not(SF_comp) #mask is True when spaxel is neither SF nor comp
        
        temp_gm = np.copy(gm)
        temp_gm[mask1] = math.nan
        
        temp_gm_error = np.copy(gm_error)
        temp_gm_error[mask1] = math.nan

        smd = np.log10(get_smd_map(MaNGA_ID))
        temp_smd = np.copy(smd)
        temp_smd[mask1] = math.nan
 
        spx_x = np.arange(gm.shape[0])
        spx_y = np.arange(gm.shape[1])
        
        mask_dist = np.zeros((gm.shape[0], gm.shape[1]), dtype = "bool")
        
        for i in spx_x:
            
            for j in spx_y:
                
                distance = (np.sqrt((i - center[0])**2 + (j - center[1])**2))*kpc_per_spx
                
                if(distance_range[1] != -1):
                    
                    if((distance >= distance_range[0]) and (distance <= distance_range[1])):
                        mask_dist[i][j] = True
                        
                elif(distance_range[1] == -1):
                    
                    if(distance >= distance_range[0]):
                        mask_dist[i][j] = True
        
        temp_gm[mask_dist] = math.nan
        temp_gm_error[mask_dist] = math.nan
        temp_smd[mask_dist] = math.nan

        temp_gm_flattened = temp_gm.flatten()
        temp_smd_flattened = temp_smd.flatten()
        temp_gm_error_flattened = temp_gm_error.flatten()

        mask2 = np.logical_and(np.logical_and(temp_smd_flattened>0, temp_gm_flattened>0), temp_gm_error_flattened > 0) #removing nan and ambigous values

        my_gm_flattened = np.append(my_gm_flattened, temp_gm_flattened[mask2])
        my_smd_flattened = np.append(my_smd_flattened, temp_smd_flattened[mask2])
        my_gm_error_flattened = np.append(my_gm_error_flattened, temp_gm_error_flattened[mask2])
        #my_gm_flattened, my_smd_flattened, my_gm_error_flattened are of the same size, with no nan or -ve values

    return my_smd_flattened, my_gm_flattened, my_gm_error_flattened
    

#Getting apertures
        
def get_minor_axis(major_axis, ellip):
    '''
    obtain the minor axis of an ellipse, given the major axis and the ellipticity
    input: major axis and ellipticity
    returns: minor axis length
    '''
    return major_axis - (major_axis*ellip)

def get_inclination(e):
    #Returns the inclination angle of the galaxy in radians, as calculated from ellipticity (e)
    
    if e < 0.8:
        inc = np.arccos(np.sqrt(1 + (((e**2)-(e*2))/0.96)))
    
    else:
        inc = np.pi/2
        
    return inc

def get_photometric_params(MaNGA_ID):   
    
    hdul = get_fits_hdul(MaNGA_ID)
    
    spx_scale = hdul[0].header['CD2_2']*3600 #angular size of each spaxel in arcsec
    FWHM = (hdul[0].header["RFWHM"])/spx_scale  #r band FWHM in spaxels
    center = np.array([hdul[0].header['CRPIX1']-hdul[0].header['CRPIX3'], hdul[0].header['CRPIX2']-hdul[0].header['CRPIX3']])    #IFU Center in spaxels (0 indexed)
    
    #Find the row index of MaNGA_ID in the totalrband data table
    IDs = np.char.strip(totalrband['plateifu_1'])
    IDs = IDs.astype('str')
    loc = np.where(IDs == MaNGA_ID)[0][0]
    
    #Find the row index of MaNGA_ID in the SPA_corr data table
    IDs1 = np.char.strip(SPA_corr['PLATEIFU'])
    IDs1 = IDs1.astype('str')
    loc1 = np.where(IDs1 == MaNGA_ID)[0][0]
    
    
    pos_angle = totalrband['PA_S'][loc] #pos_angle: position angle in degrees from pymorph
    Re_spaxels = totalrband['A_HL_S'][loc]/spx_scale #Re in spaxels from the pymorph
    ellipticity = 1 - (totalrband["BA_S"][loc]) #ellipticity: 1 - (minor_axis/major_axis) from pymorph
    spa = SPA_corr['SPA_R'][loc1]
    inclination = get_inclination(ellipticity)
    
    return center, FWHM, (-1)*(90 - pos_angle - spa)*np.pi/180, Re_spaxels, ellipticity, inclination

def get_apertures_elliptical(MaNGA_ID, method = "Re", binsize = 0.2):
    
    hdul = get_fits_hdul(MaNGA_ID)
    center, FWHM, pos_angle, Re_spaxels, ellipticity, inclination = get_photometric_params(MaNGA_ID)
    spx_scale = hdul[0].header['CD2_2']*3600 #angular scale of each spaxel in arcsec
    z = get_redshift(MaNGA_ID)
    
    First = MaNGA_ID.split('-')[0]
    Last = MaNGA_ID.split('-')[1]

    if(Last.find('127') == 0):
        IFUradius_spaxels = 16/spx_scale

    elif(Last.find('91') == 0):
        IFUradius_spaxels = 13.5/spx_scale

    elif(Last.find('61') == 0):
        IFUradius_spaxels = 11/spx_scale

    elif(Last.find('37') == 0):
        IFUradius_spaxels = 8.5/spx_scale

    elif(Last.find('19') == 0):
        IFUradius_spaxels = 6/spx_scale
        
    
    if(method == "Re"):
        
        '''
        Gets elliptical apertures for a given galaxy
        Inputs: MaNGA ID of the galaxy, radial bin size in units of Re (by default, set to 0.2)
        Outputs: Returns two arrays, first is an array of aperture objects, whose first element is the central elliptical
                aperture, and succesive elements are elliptical annuli. The second array contains the radial distances
                to these apertures, distance is being measured from the center of the IFU to the midpoint of the
                aperture/annuli, along the major axis of the galaxy.
        '''      

        num_bins = np.round((IFUradius_spaxels/(binsize*Re_spaxels)),0)

        start = binsize*Re_spaxels
        limit = num_bins*start
        thickness = binsize*Re_spaxels

        radial_distances = []
        apertures = []
        a_in = start
        a_out = a_in + thickness
        b_out = get_minor_axis(a_out, ellipticity)
        aperture = EllipticalAperture(center, a_in, get_minor_axis(a_in, ellipticity), theta = pos_angle)
        apertures.append(aperture)
        radial_distances.append(a_in/(2*Re_spaxels))
        while(a_out <= limit):
            annulus = EllipticalAnnulus(center, a_in, a_out, get_minor_axis(a_out, ellipticity), theta = pos_angle)
            apertures.append(annulus)
            radial_distances.append((a_in + a_out)/(2*Re_spaxels))
            a_in = a_out
            a_out = a_out + thickness

        radial_distances = [round(r,1) for r in radial_distances]    
        return apertures, radial_distances
    
    elif(method == "Kpc"):
        
        '''
        Gets elliptical apertures for a given galaxy
        Inputs: MaNGA ID of the galaxy, radial bin size in units of Kpc
        Outputs: Returns two arrays, first is an array of aperture objects, whose first element is the central elliptical
                aperture, and succesive elements are elliptical annuli. The second array contains the radial distances
                to these apertures, distance is being measured from the center of the IFU to the midpoint of the
                aperture/annuli, along the major axis of the galaxy.
        '''   
        
        spx_per_kpc =  (cosmo.arcsec_per_kpc_proper(z)).value/spx_scale
        binsize_spx = spx_per_kpc * binsize
        num_bins = np.round((IFUradius_spaxels/binsize_spx),0)
        
        start = binsize_spx
        limit = num_bins*start
        thickness = binsize_spx
        
        radial_distances = []
        apertures = []
        a_in = start
        a_out = a_in + thickness
        b_out = get_minor_axis(a_out, ellipticity)
        aperture = EllipticalAperture(center, a_in, get_minor_axis(a_in, ellipticity), theta = pos_angle)
        apertures.append(aperture)
        radial_distances.append(a_in/(2*spx_per_kpc))
        while(a_out <= limit):
            annulus = EllipticalAnnulus(center, a_in, a_out, get_minor_axis(a_out, ellipticity), theta = pos_angle)
            apertures.append(annulus)
            radial_distances.append((a_in + a_out)/(2*spx_per_kpc))
            a_in = a_out
            a_out = a_out + thickness

        radial_distances = [round(r,1) for r in radial_distances]    
        return apertures, radial_distances
        

#Aperture photometry

def get_aperture_data(image, apertures):
    '''
    Obtain the data contained in each aperture, when they are drawn over an image.
    input: the image and a list of aperture objects.
    returns: an array with elements themselves being arrays of the data within each aperture/annulus.
    '''
    data_list = []
    for my_aperture in apertures:
        aperture_mask = my_aperture.to_mask(method = 'center')
        aperture_data = aperture_mask.multiply(image)
        data_list.append(aperture_data)
    return np.asarray(data_list)

def get_individual_radial_profile(MaNGA_ID, data_list, surface_density = True, statistic = "median"):
    '''
    Obtains the median values over all spaxels within each aperture
    input: a list with elements themselves being array of data within each aperture
    returns: an array of median values in each aperture
    (If surface_density = True, then the final array is multiplied by the deprojection factor (sec(i)))
    '''
    
    center, FWHM, pa, Re_spaxels, ellipticity, inclination = get_photometric_params(MaNGA_ID)
    
    if(statistic == "median"):
    
        values = []
        for aperture_data in data_list:
            data_array = aperture_data.flatten()

            if(len(data_array[data_array>0]) > 0):
                if(surface_density == False):
                    values.append(np.median(data_array[data_array>0]))
                if(surface_density == True):        
                    values.append((np.median(data_array[data_array>0]))/np.cos(inclination))

            else:
                values.append(math.nan)


        return np.asarray(values)
    
    elif(statistic == "mean"):
        
        values = []
        for aperture_data in data_list:
            data_array = aperture_data.flatten()

            if(len(data_array[data_array>0]) > 0):
                if(surface_density == False):
                    values.append(np.mean(data_array[data_array>0]))
                if(surface_density == True):        
                    values.append((np.mean(data_array[data_array>0]))/np.cos(inclination))

            else:
                values.append(math.nan)


        return np.asarray(values)
    
    elif(statistic == "sum"):
        
        values = []
        for aperture_data in data_list:
            data_array = aperture_data.flatten()

            if(len(data_array[data_array>0]) > 0):
                if(surface_density == False):
                    values.append(np.sum(data_array[data_array>0]))
                if(surface_density == True):        
                    values.append((np.sum(data_array[data_array>0]))/np.cos(inclination))

            else:
                values.append(math.nan)


        return np.asarray(values)
    
    elif(statistic == "quad"):
        
        values = []
        for aperture_data in data_list:
            data_array = aperture_data.flatten()

            if(len(data_array[data_array>0]) > 0):
                if(surface_density == False):
                    values.append(np.sqrt(np.sum(data_array[data_array>0]**2))/len(data_array[data_array>0]))
                if(surface_density == True):        
                    values.append(np.sqrt(np.sum(data_array[data_array>0]**2))/(np.cos(inclination)*len(data_array[data_array>0])))

            else:
                values.append(math.nan)


        return np.asarray(values)
    
def get_radial_profile_catalog(MaNGA_IDs, method = "Re", binsize = 0.2, stat = "mean", snr_cutoff = 3, BPT_cut = True, deproject = True):
    '''
    Columns: sno, plateifu, ra, dec, redshift, radial_distances, sigma_star_profile, sigma_sfr_profile
    where,
    radial_distances: array of radial distances of the apertures for that object in units of Re
    sigma_star_profile: array of median sigma_star values for that object in units of log10(Msun/Kpc^2)
    sigma_sfr_profile: array of median sigma_sfr values for that object in units of log10(Msun/yr/Kpc^2)
    '''
    
    mytable = pd.DataFrame({'sno': [], 'plateifu': [], 'ra': [], 'dec': [], 'redshift': [], 'stellar_mass': [], 'radial_distances': [], 'sigma_star_profile': [], 'sigma_star_norm_profile': [], 'sigma_sfr_profile': [], 'sigma_sfr_error_profile': [], 'ssfr_profile': [], 'ssfr_error_profile': [], 'weight': []})
    
    count = 0
    
    for MaNGA_ID in MaNGA_IDs:

        count = count + 1
        
        ra, dec = get_coordinates(MaNGA_ID)
        
        redshift = get_redshift(MaNGA_ID)

        #Find the row index of MaNGA_ID in the Pipe3D_Salim data table
        IDs = np.char.strip(Pipe3D_Salim['plateifu'])
        IDs = IDs.astype('str')
        loc = np.where(IDs == MaNGA_ID)[0][0]
        
        stellar_mass = Pipe3D_Salim["col10"][loc]
        
        #Find the row index of MaNGA_ID in the MaNGA_target data table
        IDs1 = np.char.strip(MaNGA_target['plateifu'])
        IDs1 = IDs1.astype('str')
        loc1 = np.where(IDs1 == MaNGA_ID)[0][0]
        
        weight = MaNGA_target['ESWEIGHT'][loc1]
        
        apertures, radial_distances = get_apertures_elliptical(MaNGA_ID, method, binsize)
        
        SF, comp, liners, seyfert = get_BPT_tags(MaNGA_ID, snr_cutoff)
        SF_comp = np.logical_or(SF, comp) #spaxels that are either SF or comp
        BPT_mask = np.logical_not(SF_comp) #mask is True when spaxel is neither SF nor comp

        smd = np.copy(get_smd_map(MaNGA_ID))
        smd_norm = np.copy(get_smd_map(MaNGA_ID))/(10**(stellar_mass))
        if(BPT_cut == True):
            smd[BPT_mask] = math.nan
            smd_norm[BPT_mask] = math.nan
        smd_datalist = get_aperture_data(smd, apertures)
        smd_norm_datalist = get_aperture_data(smd_norm, apertures)
        if(deproject == True):
            sigma_star_profile = np.log10(get_individual_radial_profile(MaNGA_ID, smd_datalist, surface_density = True, statistic = stat))
            sigma_star_norm_profile = np.log10(get_individual_radial_profile(MaNGA_ID, smd_norm_datalist, surface_density = True, statistic = stat))
        else:
            sigma_star_profile = np.log10(get_individual_radial_profile(MaNGA_ID, smd_datalist, surface_density = False, statistic = stat))
            sigma_star_norm_profile = np.log10(get_individual_radial_profile(MaNGA_ID, smd_norm_datalist, surface_density = False, statistic = stat))

        SFR, SFR_error, sigma_sfr, sigma_sfr_error = get_SFR_map(MaNGA_ID, snr_cutoff)
        
        sigma_sfr[BPT_mask] = math.nan
        sigma_sfr_datalist = get_aperture_data(sigma_sfr, apertures)
        if(deproject == True):
            sigma_sfr_profile = np.log10(get_individual_radial_profile(MaNGA_ID, sigma_sfr_datalist, surface_density = True, statistic = stat)) #sigma_sfr_profile is logarithmic
        else:
            sigma_sfr_profile = np.log10(get_individual_radial_profile(MaNGA_ID, sigma_sfr_datalist, surface_density = False, statistic = stat)) #sigma_sfr_profile is logarithmic
        
        sigma_sfr_error[BPT_mask] = math.nan
        sigma_sfr_error_datalist = get_aperture_data(sigma_sfr_error, apertures)
        if(deproject == True):
            sigma_sfr_error_profile = get_individual_radial_profile(MaNGA_ID, sigma_sfr_error_datalist, surface_density = True, statistic = 'quad')/(10**sigma_sfr_profile)
        else:
            sigma_sfr_error_profile = get_individual_radial_profile(MaNGA_ID, sigma_sfr_error_datalist, surface_density = False, statistic = 'quad')/(10**sigma_sfr_profile)
        
        #ssfr profile
        
        smd1 = np.copy(get_smd_map(MaNGA_ID))
        smd1[BPT_mask] = math.nan
        #smd1_datalist = get_aperture_data(smd1, apertures)
        ssfr_datalist = get_aperture_data(sigma_sfr/smd1, apertures)
        
        ssfr_profile = np.log10(get_individual_radial_profile(MaNGA_ID, ssfr_datalist, surface_density = False, statistic = stat))
        
        #ssfr_profile = np.log10(get_individual_radial_profile(MaNGA_ID, sigma_sfr_datalist, surface_density = False, statistic = "sum")/get_individual_radial_profile(MaNGA_ID, smd1_datalist, surface_density = False, statistic = "sum"))
        
        #ssfr_error_profile
        
        ssfr_error_datalist = get_aperture_data(sigma_sfr_error/smd1, apertures) 
        
        #ssfr_error_profile = get_individual_radial_profile(MaNGA_ID, ssfr_error_datalist, surface_density = False, statistic = 'quad')/(10**ssfr_profile)
        ssfr_error_profile = sigma_sfr_error_profile
        #ssfr_error_profile = get_individual_radial_profile(MaNGA_ID, sigma_sfr_error_datalist, surface_density = False, statistic = "quad")/get_individual_radial_profile(MaNGA_ID, sigma_sfr_datalist, surface_density = False, statistic = "sum")            
        
        newrow = pd.DataFrame({'sno': count, 'plateifu': MaNGA_ID, 'ra': ra, 'dec': dec, 'redshift': redshift, 'stellar_mass': stellar_mass, 'radial_distances' : [np.asarray(radial_distances)], 'sigma_star_profile': [np.asarray(sigma_star_profile)], 'sigma_star_norm_profile': [np.asarray(sigma_star_norm_profile)], 'sigma_sfr_profile': [np.asarray(sigma_sfr_profile)], 'sigma_sfr_error_profile': [np.asarray(sigma_sfr_error_profile)], 'ssfr_profile': [np.asarray(ssfr_profile)], 'ssfr_error_profile': [np.asarray(ssfr_error_profile)], 'weight': weight})
        
        mytable = mytable.append(newrow, ignore_index = True)
       
    return mytable

def get_sigma_map(MaNGA_ID):
    #returns the stellar velocity dispersion map (km/s) and the corresponding error map (km/s)
    
    hdul = get_fits_hdul(MaNGA_ID)
    
    sigma = np.copy(hdul[1].data[15])
    sigma_error = np.copy(hdul[1].data[16])
    
    return sigma, sigma_error

def get_central_sigma(MaNGA_ID, method = "Kpc", binsize = 1, stat = "mean"):
    
    hdul = get_fits_hdul(MaNGA_ID)
    center, FWHM, pos_angle, Re_spaxels, ellipticity, inclination = get_photometric_params(MaNGA_ID)
    spx_scale = hdul[0].header['CD2_2']*3600 #angular scale of each spaxel in arcsec
    z = get_redshift(MaNGA_ID)
    
    spx_per_kpc =  (cosmo.arcsec_per_kpc_proper(z)).value/spx_scale
    binsize_spx = spx_per_kpc * binsize

    a_in = binsize_spx
    
    aperture = [EllipticalAperture(center, a_in, get_minor_axis(a_in, ellipticity), theta = pos_angle)]
    
    sigma, sigma_error = get_sigma_map(MaNGA_ID)
    
    sigma_datalist = get_aperture_data(sigma, aperture)
    sigma_error_datalist = get_aperture_data(sigma_error, aperture)
    
    cen_sigma = get_individual_radial_profile(MaNGA_ID, sigma_datalist, surface_density = False, statistic = "mean")[0]
    cen_sigma_error = get_individual_radial_profile(MaNGA_ID, sigma_error_datalist, surface_density = False, statistic = "mean")[0]
    
    #Find the row index of MaNGA_ID in the MaNGA_target data table
    IDs1 = np.char.strip(MaNGA_target['plateifu'])
    IDs1 = IDs1.astype('str')
    loc1 = np.where(IDs1 == MaNGA_ID)[0][0]
    
    weight = MaNGA_target['ESWEIGHT'][loc1]
    
    return cen_sigma, cen_sigma_error, weight

def get_radial_profile(catalog, qty_colname, radii):
    
    mean_qty = np.array([])
    error_qty = np.array([])

    for i in range(len(radii)): #iterate over each radius

        temp_qty = np.array([]) 

        for j in range(len(catalog)): #iterate over each galaxy

            qty_array = catalog[qty_colname][j]

            if(i <= len(qty_array)-1):
                
                temp_qty = np.append(temp_qty, qty_array[i])

        mean_qty = np.append(mean_qty, np.mean(temp_qty[np.isfinite(temp_qty)]))
        error_qty = np.append(error_qty, np.std(temp_qty[np.isfinite(temp_qty)]))
       
    return mean_qty, error_qty


def get_radially_binned_arrays(catalog, qty_colname, atRadii, remove_nan = True, return_weights = False):
    #obtain the qty_array for all objects in the catalog at a particular radial distance
    
    arr = np.array([])
    weights = np.array([])
    for i in range(len(catalog)):
        r = catalog["radial_distances"][i]
        idx = np.where(r == atRadii)[0]
        if(len(idx)>0):
            index = idx[0]
            qty_array = catalog[qty_colname][i]
            if(index < len(qty_array)):
                arr = np.append(arr, qty_array[index])
                weights = np.append(weights, catalog['weight'][i])
                
    if(remove_nan == True):            
        array = arr[np.isfinite(arr)] #To remove nan's
        weights = weights[np.isfinite(arr)]
    else:
        array = arr
        
    if(return_weights):
        return array, weights
    else:
        return array

#delta_sigma_sfr radial profile

def get_delta_sigma_sfr_profile(MaNGA_ID, main_catalog, control_catalog, cs_sm_bins  = np.array([8.9, 9.4, 9.9, 10.4, 10.9, 11.4])):
    
    IDs = main_catalog['plateifu']
    loc = np.where(IDs == MaNGA_ID)[0][0]

    stellar_mass = main_catalog["stellar_mass"][loc]

    #get the positions of the control galaxies for this main sample galaxy, let that array be control_positions
    for m in range(len(cs_sm_bins)-1):
        if((cs_sm_bins[m] <= stellar_mass) and (cs_sm_bins[m+1] > stellar_mass)):
            idx = m
            break;

    control_positions = np.where(np.logical_and(control_catalog["stellar_mass"] >= cs_sm_bins[idx], control_catalog["stellar_mass"] <= cs_sm_bins[idx+1]))[0]       

    main_radial_distances = main_catalog["radial_distances"][loc]

    delta_sigma_sfr_mean = np.array([])
    delta_sigma_sfr_error = np.array([])

    for r in range(len(main_radial_distances)):
        sigma_sfr_main_r = main_catalog["sigma_sfr_profile"][loc][r]

        if(np.isfinite(sigma_sfr_main_r)):
            temp_array = np.array([])

            for i in control_positions:

                if(r < len(control_catalog["radial_distances"][i])):
                    sigma_sfr_control_r = control_catalog["sigma_sfr_profile"][i][r]
                    if(np.isfinite(sigma_sfr_control_r)):
                        temp_array = np.append(temp_array, sigma_sfr_main_r - sigma_sfr_control_r)
                    else:
                        temp_array = np.append(temp_array, math.nan)


            mean = np.mean(temp_array[np.isfinite(temp_array)])
            error = np.std(temp_array[np.isfinite(temp_array)])

        else:
            mean = math.nan
            error = math.nan

        delta_sigma_sfr_mean = np.append(delta_sigma_sfr_mean, mean)
        delta_sigma_sfr_error = np.append(delta_sigma_sfr_error, error)
        
    n_radial_bins = len(delta_sigma_sfr_mean)
    radial_distances = np.array([])
    
    for i in range(n_radial_bins):
        r = 0.5 + i
        radial_distances = np.append(radial_distances, r)

    return delta_sigma_sfr_mean, delta_sigma_sfr_error, radial_distances

#function to calculate radial profile gradients

def get_avg_gradient(MaNGA_ID, catalog, qty = "sigma_sfr_profile"):
    
    IDs = catalog['plateifu']
    loc = np.where(IDs == MaNGA_ID)[0][0]
    
    radial_distances = catalog["radial_distances"][loc]
    qty_profile = catalog[qty][loc]
    
    #radial_distances and qty_profile are arrays of same length, but qty_profile might contain nan's
    #Goal: find slope of the line joining initial point and final point
    
    if(len(qty_profile[np.isfinite(qty_profile)]) <= 1):
        return 0
    else:
        isfinite_qty_profile = np.where(np.isfinite(qty_profile) == True)[0]
        
        r1 = radial_distances[min(isfinite_qty_profile)]
        q1 = qty_profile[min(isfinite_qty_profile)]
        r2 = radial_distances[max(isfinite_qty_profile)]
        q2 = qty_profile[max(isfinite_qty_profile)]
        
        grad = (q2 - q1)/(r2 - r1)
        
        return grad

#Drawing contours over resolved plots

def get_contour_levels(hist, drawAt = [15, 35, 55, 75, 95]):
    
    highest = np.max(hist)
    levels = []
    for i in range(len(drawAt)):
        levels.append(highest * (drawAt[i]/100))
    
    return levels

    
    
    
