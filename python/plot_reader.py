import h5py

class PlotDataset:
    def __init__(self, name, x, y, label):
        self.name = name
        self.x = x
        self.y = y
        self.label = label

class Plot:
    def __init__(self, name, title, xlabel, ylabel, datasets):
        self.name = name
        self.title = title
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.datasets = datasets  # list of PlotDataset

class PlotReader:
    def __init__(self, filename):
        self.filename = filename
        self._plots = {}
        self._load_file()

    def _load_file(self):
        with h5py.File(self.filename, 'r') as f:
            for plot_name in f:
                plot_group = f[plot_name]
                title = plot_group.attrs.get('title', plot_name)
                xlabel = plot_group.attrs.get('xlabel', '')
                ylabel = plot_group.attrs.get('ylabel', '')
                datasets = []

                for ds_name in plot_group:
                    if isinstance(plot_group[ds_name], h5py.Group):
                        ds_group = plot_group[ds_name]
                        x = ds_group['x'][()]
                        y = ds_group['y'][()]
                        label = ds_group.attrs.get('label', ds_name)
                        datasets.append(PlotDataset(ds_name, x, y, label))

                self._plots[plot_name] = Plot(plot_name, title, xlabel, ylabel, datasets)

    def list_plots(self):
        return list(self._plots.keys())

    def get_plot(self, name):
        return self._plots.get(name)

    def get_all(self):
        return list(self._plots.values())
