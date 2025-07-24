import h5py
import numpy as np

class Frame:
    def __init__(self, phi, sigma, vertices, adj_contacts):
        # self.time = time
        # self.num_collisions = num_collisions
        self.phi = phi
        self.sigma = sigma  # scalar value
        self.vertices = vertices  # shape: (num_p, num_v, ndim)
        self.adj_contacts = adj_contacts  # shape: (num_p, num_p)

class System:
    def __init__(self, system_id, num_p, num_v, ndim):
        self.system_id = system_id
        self.num_p = num_p
        self.num_v = num_v
        self.ndim = ndim
        self.frames = []

class SystemReader:
    def __init__(self, file_path):
        self.file_path = file_path
        self.systems = {}
        self._load()

    def _load(self):
        with h5py.File(self.file_path, 'r') as f:
            for system_key in f.keys():
                if not system_key.startswith("system_"):
                    continue
                sys_group = f[system_key]
                system_id = int(system_key.split("_")[1])

                num_p = sys_group.attrs["num_p"]
                num_v = sys_group.attrs["num_v"]
                ndim = sys_group.attrs["ndim"]

                system = System(system_id, num_p, num_v, ndim)

                for frame_key in sys_group.keys():
                    frame_group = sys_group[frame_key]

                    # time = frame_group.attrs["time"]
                    # num_collisions = frame_group.attrs["num_collisions"]
                    phi = frame_group.attrs["phi"]
                    sigma = frame_group.attrs["sigma"]

                    vertices = frame_group["vertices"][()].reshape((num_p, num_v, ndim))
                    adj_contacts = frame_group["adj_contacts"][()]

                    # frame = Frame(phi, time, num_collisions, sigma, vertices, adj_contacts)
                    frame = Frame(phi, sigma, vertices, adj_contacts)
                    system.frames.append(frame)

                self.systems[system_id] = system
