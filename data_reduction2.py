#Cell 1
from pathlib import Path
import csv

import numpy as np
from astropy.io import fits

from photutils.aperture import (
    CircularAperture,
    CircularAnnulus,
    aperture_photometry,
    ApertureStats,
)
from photutils.centroids import centroid_sources, centroid_quadratic

#Cell 2
def norm_flat(flat_file, bias_file):
    """
    Bias-subtract and normalize a flat-field image.

    Parameters
    ----------
    flat_file : str or pathlib.Path
        Path to the raw flat-field FITS file.
    bias_file : str or pathlib.Path
        Path to the bias FITS file.

    Returns
    -------
    numpy.ndarray
        Normalized flat field.
    """
    raw_flat = fits.getdata(flat_file)
    bias = fits.getdata(bias_file)

    flat_minus_bias = raw_flat - bias
    normalized_flat = flat_minus_bias / np.median(flat_minus_bias)

    return normalized_flat


def correct_data(data_file, flat_value=1.0):
    """
    Load a science image and divide it by a single scalar flat value.

    Parameters
    ----------
    data_file : str or pathlib.Path
        Science FITS file.
    flat_value : float, optional
        Scalar correction factor applied to every pixel.
        Use 1.0 to apply no correction.
    """
    raw_image = fits.getdata(data_file).astype(float)

    if raw_image.ndim != 2:
        raise ValueError(
            f"Expected a 2D image in {data_file}, "
            f"but found shape {raw_image.shape}"
        )

    if not np.isfinite(flat_value) or flat_value == 0:
        raise ValueError(
            f"flat_value must be finite and nonzero; got {flat_value}"
        )

    return raw_image / float(flat_value)


def find_centroid(corrected_image, x_guess, y_guess, box_size=21):
    """
    Refine the centroid near an initial pixel-coordinate guess.
    """
    xs, ys = centroid_sources(
        corrected_image,
        x_guess,
        y_guess,
        box_size=box_size,
        centroid_func=centroid_quadratic,
    )

    x, y = float(xs[0]), float(ys[0])
    half_box = box_size / 2

    # Accept the centroid only if it is finite and remains near the guess.
    if (
        np.isfinite(x)
        and np.isfinite(y)
        and abs(x - x_guess) <= half_box
        and abs(y - y_guess) <= half_box
    ):
        return x, y

    return x_guess, y_guess


def extract_aperture_and_annulus(corrected_image, center_x, center_y):
    """
    Measure aperture flux, background noise, and signal-to-noise ratio.
    """
    aperture_radius = 5.2
    annulus_inner_radius = 10
    annulus_outer_radius = 13

    aperture = CircularAperture(
        (center_x, center_y),
        r=aperture_radius,
    )

    annulus = CircularAnnulus(
        (center_x, center_y),
        r_in=annulus_inner_radius,
        r_out=annulus_outer_radius,
    )

    phot_table = aperture_photometry(corrected_image, aperture)
    flux = float(phot_table["aperture_sum"][0])

    background_stats = ApertureStats(corrected_image, annulus)
    background_sigma = float(background_stats.std)

    flux_error = np.sqrt(aperture.area) * background_sigma
    snr = flux / flux_error if flux_error != 0 else np.nan

    return flux, float(flux_error), float(snr)


def run_reduction(
    raw_data_folder,
    flat_value=1.0,
    initial_x=0,
    initial_y=0,
    output_file=None,
    box_size=21,
):
    """
    Process all FITS science images using a scalar flat value.

    No bias image or 2D flat-field image is used.
    """
    raw_data_folder = Path(raw_data_folder)

    center_x = float(initial_x)
    center_y = float(initial_y)
    results = []

    science_files = sorted(raw_data_folder.glob("*.fits"))

    if not science_files:
        raise FileNotFoundError(
            f"No FITS files found in {raw_data_folder.resolve()}"
        )

    for science_file in science_files:
        corrected_data = correct_data(
            science_file,
            flat_value=flat_value,
        )

        # Update the centroid for this image and use it as the next guess.
        center_x, center_y = find_centroid(
            corrected_data,
            center_x,
            center_y,
            box_size=box_size,
        )

        flux, flux_error, snr = extract_aperture_and_annulus(
            corrected_data,
            center_x,
            center_y,
        )

        result = {
            "Flux": flux,
            "Flux Error": flux_error,
            "SNR": snr,
        }

        results.append(result)

        print(
            f"{science_file.name}: "
            f"x={center_x:.2f}, y={center_y:.2f}, "
            f"flux={flux:.3f}, SNR={snr:.2f}"
        )

    if output_file is not None:
        output_file = Path(output_file)

        fieldnames = [
            "Flux",
            "Flux Error",
            "SNR",
        ]

        with output_file.open("w", newline="") as file:
            writer = csv.DictWriter(file, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(results)

        print(f"\nResults written to: {output_file.resolve()}")

    return results
