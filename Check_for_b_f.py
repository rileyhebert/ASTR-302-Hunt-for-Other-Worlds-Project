import os
import numpy as np
from astropy.io import fits

input_path = "/home/<netID>/hunt_other_worlds/rflat2"
output_path = "/home/<netID>/hunt_other_worlds"

def types(i, f):
    return [os.path.join(input_path, f"r.{n:04d}.fits") for n in range(i, f+1)]

bias = [a for a in types(1, 99)] #bias files in rflat2
flat = [a for a in types(100, 699)] #flat files in R band, in rflat2
print(f"Found {len(bias)} bias and {len(flat)} flat")

#Median of the biases
median_bias = np.median([fits.getdata(a) for a 
                         in bias], axis=0)
fits.writeto(os.path.join(output_path, "median_bias.fits"),
             median_bias, overwrite=True)

#median of the raw flats
median_flat = np.median([fits.getdata(a) for a 
                         in flat], axis=0)
hdr = fits.Header()
hdr["IMAGETYP"] = "MEDIAN FLAT"
hdr["Filter"] = "R"
hdr["NCOMBINE"] = len(flat)
fits.writeto(os.path.join(output_path, "median_flat.fits"),
             median_flat, hdr, overwrite=True)
print(f"Wrote median bias and flat files to {output_path}")
