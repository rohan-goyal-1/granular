import argparse
import matplotlib.pyplot as plt
from plot_reader import PlotReader

def plot_interactively(plot):
    plt.figure(figsize=(8, 6))
    for dataset in plot.datasets:
        plt.plot(dataset.x, dataset.y, label=dataset.label)

    plt.title(plot.title)
    plt.xlabel(plot.xlabel)
    plt.ylabel(plot.ylabel)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def main():
    parser = argparse.ArgumentParser(description="Plot graphs from an HDF5 file.")
    parser.add_argument("hdf5_file", help="Path to the HDF5 file.")
    parser.add_argument(
        "-p", "--plots", nargs="+", help="Names of specific plots to render. If omitted, all plots are shown."
    )

    args = parser.parse_args()
    reader = PlotReader(args.hdf5_file)

    if args.plots:
        for name in args.plots:
            plot = reader.get_plot(name)
            if plot is None:
                print(f"Plot '{name}' not found in file.")
            else:
                plot_interactively(plot)
    else:
        for plot in reader.get_all():
            plot_interactively(plot)

if __name__ == "__main__":
    main()

