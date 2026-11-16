# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 10:42:28 2026

@author: Riley
"""

import sys
#import matplotlib.pyplot as plt  
import numpy as np               
from astropy.io import fits  
import glob
from photutils.aperture import CircularAperture, CircularAnnulus, aperture_photometry, ApertureStats    
import csv
import os

def norm_flat(flat_file, bias_file):
    raw_flat = fits.getdata(flat_file)
    bias = fits.getdata(bias_file)
    normed_flats = (raw_flat - bias)/np.median(raw_flat - bias) 
    return normed_flats
            
def correct_data(data_file, normed_flat):
    raw_image = fits.getdata(data_file)
    corrected_image = raw_image/normed_flat
    return corrected_image

def extract_aperture_and_annulus(corrected_image,center_x,center_y):
    aperture_radius = 7
    r_in = 10
    r_out = 13
    aperture = CircularAperture((center_x, center_y), r=aperture_radius) 
    annulus = CircularAnnulus((center_x, center_y), r_in, r_out)
    
    phot_table = aperture_photometry(corrected_image, aperture)
    flux = phot_table['aperture_sum']
    
    bkg_stats = ApertureStats(corrected_image, annulus)
    bkg_sigma = bkg_stats.std
    flux_err = np.sqrt(aperture.area)*bkg_sigma
    
    snr = flux/flux_err
    
    return [flux,flux_err,snr]

def write_to_flux_file(output_file,outputs):
    headers = ["Flux","Flux Error","SNR"]
    with open(output_file,"a") as file:
        writer = csv.DictWriter(file, headers)
        if os.path.getsize(output_file) == 0:
            writer.writeheader()
        writer.writerow(outputs)
        
    




def main():
    if len(sys.argv) < 7:
        print("Usage: 'python data_reduction.py"
              " <raw data file> <raw flats file> <bias file>'"
              " <x pixel of star center> <y pixel of star center>"
              " <output file")
    else:
        raw_data_filename = sys.argv[1]
        raw_flat_filename = sys.argv[2]
        bias_filename = sys.argv[3]
        #returns a numpy array
        normed_flat = norm_flat(raw_flat_filename,bias_filename)
        corrected_data = correct_data(raw_data_filename,normed_flat)
        center_x_pixel = float(sys.argv[4])
        center_y_pixel = float(sys.argv[5])
        outputs = extract_aperture_and_annulus(corrected_data, 
                                                         center_x_pixel, 
                                                        center_y_pixel)
        output_file = sys.argv[6]
        write_to_flux_file(output_file, outputs)
        