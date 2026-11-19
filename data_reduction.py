# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 10:42:28 2026

@author: Riley
"""
"""
Run this code by typing:
python data_reduction.py <raw-science-folder> <master-flat> <master-bias> <center-x-on-first-science-image> 
<center-y-on-first-science-image> <name-of-output-file>

For the last three inputs, use these values:
On wasp2.0001.fits: center x = 955, center y = 1277
name of output file: reduced_flux_file

"""

import sys
#import matplotlib.pyplot as plt  
import numpy as np               
from astropy.io import fits  
import glob
from photutils.aperture import CircularAperture, CircularAnnulus, aperture_photometry, ApertureStats    
from photutils.centroids import centroid_sources, centroid_quadratic
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

def find_centroid(corrected_image, x_guess, y_guess, box_size=21):
    xs, ys = centroid_sources(corrected_image, x_guess, y_guess,
                              box_size=box_size, centroid_func=centroid_quadratic)
    x, y = float(xs[0]), float(ys[0])
    half = box_size / 2
    # Accept the fit only if it is finite and inside the search box
    if (np.isfinite(x) and np.isfinite(y)
            and abs(x - x_guess) <= half and abs(y - y_guess) <= half):
        return x, y
    return x_guess, y_guess

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
    
    return [float(flux[0]), float(flux_err), float(snr[0])]

def write_to_flux_file(output_file, outputs):
    headers = ["Flux", "Flux Error", "SNR"]
    with open(output_file, "a", newline="") as file:
        writer = csv.writer(file)
        if os.path.getsize(output_file) == 0:
            writer.writerow(headers)
        writer.writerow(outputs)
    




def main():
    if len(sys.argv) < 7:
        print("Usage: 'python data_reduction.py"
              " <raw data file> <raw flats file> <bias file>'"
              " <x pixel of star center> <y pixel of star center>"
              " <output file")
    else:
        raw_data_folder = sys.argv[1]
        raw_flat_filename = sys.argv[2]
        bias_filename = sys.argv[3]
        #returns a numpy array
        normed_flat = norm_flat(raw_flat_filename,bias_filename)
        center_x_pixel = float(sys.argv[4])
        center_y_pixel = float(sys.argv[5])
        output_file = sys.argv[6]
        for raw_data_filename in sorted(glob.glob(os.path.join(raw_data_folder, "*.fits"))):
            corrected_data = correct_data(raw_data_filename,normed_flat)
            center_x_pixel, center_y_pixel = find_centroid(corrected_data,
                                                           center_x_pixel,
                                                           center_y_pixel)
            #print(os.path.basename(raw_data_filename), center_x_pixel, center_y_pixel)
            outputs = extract_aperture_and_annulus(corrected_data, 
                                                            center_x_pixel, 
                                                            center_y_pixel)
            write_to_flux_file(output_file, outputs)


main()
