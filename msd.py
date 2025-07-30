import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress # Import for linear regression

def calculate_msd_and_diffusion(h5_filename="run.h5"):
    """
    Calculates the MSD and Diffusion Coefficient (D) from a simulation trajectory.

    The function reads particle positions, computes the MSD for various time lags,
    identifies the diffusive regime, performs a linear fit to find D, and
    generates a log-log plot of MSD vs. time lag.

    Args:
        h5_filename (str): The path to the HDF5 file containing the trajectory data.
    """
    print(f"Analyzing file: {h5_filename}")

    try:
        with h5py.File(h5_filename, 'r') as f:
            # --- 1. Load Data ---
            frame_group = f['/frames']
            sorted_frame_names = sorted(frame_group.keys())
            positions = np.array([frame_group[name][:] for name in sorted_frame_names])

            if positions.ndim == 4 and positions.shape[2] == 1:
                positions = positions.squeeze(axis=2)

        num_frames, num_particles, _ = positions.shape
        print(f"Loaded {num_frames} frames for {num_particles} particles.")

        time_per_frame = 1.0

        # --- 2. Calculate MSD ---
        max_lag = num_frames // 2
        lags = np.unique(np.logspace(0, np.log10(max_lag), num=50).astype(int))

        msd_values = []

        print(f"Calculating MSD for {len(lags)} time lags...")
        for lag in lags:
            if lag == 0: continue
            displacements = positions[lag:] - positions[:-lag]
            squared_displacements = np.sum(displacements**2, axis=2)
            msd_values.append(np.mean(squared_displacements))

        if lags[0] == 0: lags = lags[1:]

        dt_values = lags * time_per_frame
        msd_values = np.array(msd_values)

        # --- 3. Calculate Diffusion Coefficient (D) --- 🧪
        print("Finding diffusive regime to calculate D...")

        # Calculate the local slope on the log-log plot
        log_slopes = np.diff(np.log(msd_values)) / np.diff(np.log(dt_values))

        # Find indices where the slope is close to 1 (e.g., between 0.9 and 1.1)
        # We also start looking after the initial "ballistic" regime (e.g., after the first 10% of data)
        start_index = len(log_slopes) // 10
        diffusive_indices = np.where((log_slopes[start_index:] > 0.9) & (log_slopes[start_index:] < 1.1))[0] + start_index

        diffusion_coefficient = None
        if len(diffusive_indices) > 5: # Require at least 5 points for a decent fit
            # Get the data from the identified diffusive regime
            fit_lags = dt_values[diffusive_indices]
            fit_msd = msd_values[diffusive_indices]

            # Perform linear regression: msd = slope * t + intercept
            slope, intercept, r_value, _, _ = linregress(fit_lags, fit_msd)

            # The diffusion coefficient in 2D is D = slope / 4
            diffusion_coefficient = slope / 4.0

            print(f"\n----------------------------------------------------")
            print(f"Diffusion Coefficient (D) = {diffusion_coefficient:.4e}")
            print(f"Fit performed on {len(fit_lags)} data points.")
            print(f"R-squared value of the fit: {r_value**2:.4f}")
            print(f"----------------------------------------------------")
        else:
            print("\nCould not find a clear diffusive regime to calculate D.")


        # --- 4. Plotting ---
        print("Plotting results...")
        plt.style.use('seaborn-v0_8-whitegrid')
        fig, ax = plt.subplots(figsize=(8, 6))

        # Plot the main MSD data
        ax.loglog(dt_values, msd_values, 'o', markersize=5, label='MSD Data')

        # If D was found, plot the fitted line
        if diffusion_coefficient is not None:
            fit_line = slope * fit_lags + intercept
            ax.loglog(fit_lags, fit_line, 'r-', linewidth=2.5, label=f'Linear Fit (D={diffusion_coefficient:.2e})')

        ax.set_xlabel('log(dt) [Simulation Time]')
        ax.set_ylabel('log(MSD) [Distance²]')
        ax.set_title('Mean Squared Displacement vs. Time Lag')
        ax.legend()
        ax.set_aspect('equal', adjustable='box')

        plt.tight_layout()
        plt.show()

    except FileNotFoundError:
        print(f"Error: The file '{h5_filename}' was not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == '__main__':
    calculate_msd_and_diffusion("run.h5")
