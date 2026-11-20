#Cell 1

import csv
import matplotlib.pyplot as plt
import numpy as np


def get_data(filename):
    with open(filename, "r", newline="") as file:
        return list(csv.DictReader(file))


def plot_flux_and_snr(data):
    if not data:
        raise ValueError("The CSV file contains no data.")

    required_columns = {"Flux", "Flux Error", "SNR"}
    actual_columns = set(data[0].keys())
    missing_columns = required_columns - actual_columns

    if missing_columns:
        raise ValueError(
            f"Missing columns: {sorted(missing_columns)}\n"
            f"Found columns: {sorted(actual_columns)}"
        )

    image_number = np.arange(len(data))

    flux = np.array(
        [float(row["Flux"]) for row in data],
        dtype=float,
    )

    flux_error = np.array(
        [float(row["Flux Error"]) for row in data],
        dtype=float,
    )

    snr = np.array(
        [float(row["SNR"]) for row in data],
        dtype=float,
    )

    # Create two vertically stacked plots.
    fig, axes = plt.subplots(
        2,
        1,
        figsize=(10, 8),
        sharex=True,
    )

    # Flux plot with error bars.
    axes[0].errorbar(
        image_number,
        flux,
        yerr=flux_error,
        fmt="o",
        markersize=4,
        capsize=3,
    )

    axes[0].set_ylabel("Flux (counts)")
    axes[0].set_title("Flux of the Star Over Time")
    axes[0].grid(True, alpha=0.3)

    # SNR plot.
    axes[1].scatter(
        image_number,
        snr,
        s=20,
    )

    axes[1].set_xlabel("Image Number")
    axes[1].set_ylabel("SNR")
    axes[1].set_title("Signal-to-Noise Ratio")
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()


def main():
    filename = input("Name of output file:\n").strip()

    data = get_data(filename)
    plot_flux_and_snr(data)


if __name__ == "__main__":
    main()

#Cell 2
def plot_flux_and_snr(data, minimum_flux=5e5):
    image_number = []
    flux = []
    flux_error = []
    snr = []

    for i, row in enumerate(data):
        current_flux = float(row["Flux"])
        current_error = float(row["Flux Error"])
        current_snr = float(row["SNR"])

        if (
            np.isfinite(current_flux)
            and np.isfinite(current_error)
            and np.isfinite(current_snr)
            and current_flux >= minimum_flux
            and current_error >= 0
        ):
            image_number.append(i)
            flux.append(current_flux)
            flux_error.append(current_error)
            snr.append(current_snr)

    image_number = np.array(image_number)
    flux = np.array(flux)
    flux_error = np.array(flux_error)
    snr = np.array(snr)

    fig, axes = plt.subplots(
        2,
        1,
        figsize=(10, 8),
        sharex=True,
    )

    axes[0].errorbar(
        image_number,
        flux,
        yerr=flux_error,
        fmt="o",
        markersize=4,
        capsize=3,
    )

    axes[0].set_ylabel("Flux (counts)")
    axes[0].set_title("Flux of the Star Over Time")
    axes[0].grid(True, alpha=0.3)

    axes[1].scatter(
        image_number,
        snr,
        s=20,
    )

    axes[1].set_xlabel("Image Number")
    axes[1].set_ylabel("SNR")
    axes[1].set_title("Signal-to-Noise Ratio")
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()

#Cell 3
import csv

filename = "reduced_flux_file.csv"

with open(filename, "r", newline="") as file:
    reader = csv.DictReader(file)
    data = list(reader)

#Cell 4
for row in data[:10]:
    print(
        "Flux:", row["Flux"],
        "Flux Error:", row["Flux Error"],
        "SNR:", row["SNR"]
    )

#Cell 5
plot_flux_and_snr(data, minimum_flux=5e5)

#Cell 6
for i, row in enumerate(data):
    flux = float(row["Flux"])
    flux_error = float(row["Flux Error"])
    snr = float(row["SNR"])

    print(
        f"Image {i:3d}: "
        f"Flux={flux:12.2f}, "
        f"Error={flux_error:10.2f}, "
        f"SNR={snr:8.2f}"
    )

#Cell 7
fluxes = np.array([float(row["Flux"]) for row in data])

normal_flux = np.median(fluxes[fluxes > 1e5])

print(f"Typical flux: {normal_flux:.2f}")

for i, row in enumerate(data):
    flux = float(row["Flux"])

    if flux < 0.5 * normal_flux:
        print(
            f"Possible bad measurement: "
            f"image index={i}, flux={flux:.2f}, "
            f"fraction of typical={flux / normal_flux:.3f}"
        )

#Cell 8
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# File names
target_file = "./reduced_flux_file_T_N.csv"

reference_files = [
    "./reduced_flux_file_R1_N.csv",
    "./reduced_flux_file_R2_N.csv",
    "./reduced_flux_file_R3_N.csv",
]


def read_flux_file(filename):
    """
    Read a CSV file containing Flux and Flux Error columns.
    """
    data = pd.read_csv(filename, skipinitialspace=True)

    required_columns = {"Flux", "Flux Error"}

    missing_columns = required_columns - set(data.columns)

    if missing_columns:
        raise ValueError(
            f"{filename} is missing columns: {missing_columns}\n"
            f"Columns found: {list(data.columns)}"
        )

    data["Flux"] = pd.to_numeric(data["Flux"], errors="coerce")
    data["Flux Error"] = pd.to_numeric(
        data["Flux Error"],
        errors="coerce",
    )

    if data[["Flux", "Flux Error"]].isna().any().any():
        raise ValueError(
            f"{filename} contains non-numeric or missing flux values."
        )

    return data


# Read the target and reference-star data.
target = read_flux_file(target_file)
references = [
    read_flux_file(filename)
    for filename in reference_files
]


# Make sure all files have the same number of images.
lengths = [len(target)] + [len(reference) for reference in references]

if len(set(lengths)) != 1:
    raise ValueError(
        f"The files do not contain the same number of rows: {lengths}"
    )


# Convert reference-star data into arrays.
reference_fluxes = np.array([
    reference["Flux"].to_numpy()
    for reference in references
])

reference_errors = np.array([
    reference["Flux Error"].to_numpy()
    for reference in references
])


# Calculate the arithmetic mean reference flux for each image.
mean_reference_flux = np.mean(
    reference_fluxes,
    axis=0,
)


# Propagate the independent uncertainties in the mean.
mean_reference_error = (
    np.sqrt(np.sum(reference_errors ** 2, axis=0))
    / reference_fluxes.shape[0]
)


# Calculate SNR of the mean reference flux.
mean_reference_snr = np.divide(
    mean_reference_flux,
    mean_reference_error,
    out=np.full_like(mean_reference_flux, np.nan),
    where=mean_reference_error != 0,
)


# Target-star measurements.
target_flux = target["Flux"].to_numpy()
target_error = target["Flux Error"].to_numpy()


# Differential flux:
# target flux divided by the mean reference flux.
relative_flux = np.divide(
    target_flux,
    mean_reference_flux,
    out=np.full_like(target_flux, np.nan),
    where=mean_reference_flux != 0,
)


# Propagate uncertainty in target/reference division.
relative_flux_error = relative_flux * np.sqrt(
    (target_error / target_flux) ** 2
    + (mean_reference_error / mean_reference_flux) ** 2
)


# Create the output table.
output = pd.DataFrame({
    "Image Number": np.arange(len(target)),
    "Target Flux": target_flux,
    "Target Flux Error": target_error,
    "Mean Reference Flux": mean_reference_flux,
    "Mean Reference Flux Error": mean_reference_error,
    "Mean Reference SNR": mean_reference_snr,
    "Relative Flux": relative_flux,
    "Relative Flux Error": relative_flux_error,
})


# Save the combined results.
output.to_csv(
    "combined_flux_results.csv",
    index=False,
)

print(output.head())
print("\nSaved as combined_flux_results.csv")

#Cell 9
plt.figure(figsize=(10, 5))

plt.errorbar(
    output["Image Number"],
    output["Mean Reference Flux"],
    yerr=output["Mean Reference Flux Error"],
    fmt="o",
    markersize=3,
    capsize=3,
)

plt.xlabel("Image Number")
plt.ylabel("Mean Reference Flux (counts)")
plt.title("Mean Flux of Three Reference Stars")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

#Cell 10
import numpy as np
import matplotlib.pyplot as plt

target_flux = target["Flux"].to_numpy(dtype=float)
target_error = target["Flux Error"].to_numpy(dtype=float)

reference_fluxes = np.array([
    reference["Flux"].to_numpy(dtype=float)
    for reference in references
])

reference_errors = np.array([
    reference["Flux Error"].to_numpy(dtype=float)
    for reference in references
])

# Mean reference-star flux for each image.
mean_reference_flux = np.mean(
    reference_fluxes,
    axis=0,
)

# Uncertainty in the arithmetic mean of the three reference stars.
mean_reference_error = (
    np.sqrt(np.sum(reference_errors**2, axis=0))
    / len(references)
)

# Calculate mean reference flux divided by target flux.
reference_to_target_flux = np.divide(
    mean_reference_flux,
    target_flux,
    out=np.full_like(mean_reference_flux, np.nan),
    where=target_flux != 0,
)

# Propagate uncertainty for reference_flux / target_flux.
reference_to_target_error = (
    reference_to_target_flux
    * np.sqrt(
        (mean_reference_error / mean_reference_flux) ** 2
        + (target_error / target_flux) ** 2
    )
)

#Cell 11
image_numbers = np.arange(len(target_flux))

valid = (
    np.isfinite(reference_to_target_flux)
    & np.isfinite(reference_to_target_error)
    & (target_flux > 0)
    & (mean_reference_flux > 0)
    & good
)

baseline = np.median(
    reference_to_target_flux[valid]
)

normalized_reference_to_target = (
    reference_to_target_flux / baseline
)

normalized_reference_to_target_error = (
    reference_to_target_error / baseline
)

#Cell 12
plt.figure(figsize=(12, 6))

plt.errorbar(
    image_numbers[valid],
    normalized_reference_to_target[valid],
    yerr=normalized_reference_to_target_error[valid],
    fmt="o",
    markersize=3,
    capsize=2,
    alpha=0.7,
)

plt.axhline(
    1.0,
    color="black",
    linestyle="--",
    label="Normalized baseline",
)

plt.xlabel("Image Number")
plt.ylabel("Mean Reference Flux / Target Flux")
plt.title("Inverse Differential Light Curve")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()
