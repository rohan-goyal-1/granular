import argparse
import numpy as np
import matplotlib.pyplot as plt
from plot_reader import PlotReader


def draw_tangent(ax, x, y, idx, slope, color, label):
    """Draw a tangent line at x[idx], y[idx] with given slope."""
    x0, y0 = x[idx], y[idx]
    x_vals = np.array([x0 - 5, x0 + 5])
    y_vals = y0 + slope * (x_vals - x0)
    ax.plot(x_vals, y_vals, linestyle='--', color=color, label=label)


def plot_with_tangents(plot):
    for dataset in plot.datasets:
        x, y = dataset.x, dataset.y
        dx = np.gradient(x)
        dy = np.gradient(y)
        slopes = dy / dx

        # Find index where slope is closest to 2 and 1
        idx_2 = np.argmin(np.abs(slopes - 2))
        idx_1 = np.argmin(np.abs(slopes - 1))

        fig, ax = plt.subplots(figsize=(8, 6))
        ax.plot(x, y, label=dataset.label, color='blue')

        draw_tangent(ax, x, y, idx_2, 2, 'red', 'slope ≈ 2')
        draw_tangent(ax, x, y, idx_1, 1, 'green', 'slope ≈ 1')

        ax.set_title(f"{plot.title} — {dataset.label}")
        ax.set_xlabel(plot.xlabel)
        ax.set_ylabel(plot.ylabel)
        ax.legend()
        ax.grid(True)
        plt.tight_layout()
        plt.show()


def main():
    parser = argparse.ArgumentParser(description="Plot datasets and draw tangents at slopes ≈ 2 and 1.")
    parser.add_argument("hdf5_file", help="Path to the HDF5 file.")
    parser.add_argument(
        "-p", "--plots", nargs="+", help="Names of specific plots to render. If omitted, all plots are processed."
    )

    args = parser.parse_args()
    reader = PlotReader(args.hdf5_file)

    if args.plots:
        for name in args.plots:
            plot = reader.get_plot(name)
            if plot is None:
                print(f"Plot '{name}' not found in file.")
            else:
                plot_with_tangents(plot)
    else:
        for plot in reader.get_all():
            plot_with_tangents(plot)


if __name__ == "__main__":
    main()
