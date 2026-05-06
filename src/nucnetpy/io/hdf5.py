from __future__ import annotations
from pathlib import Path
import numpy as np
from ..core import Network, Zone
from ..species import Species


def write_hdf5(network: Network, path) -> None:
    import h5py
    with h5py.File(path, 'w') as h:
        h.attrs['format'] = 'nucnetpy'
        names = network.species_names()
        h.create_dataset('species', data=np.array(names, dtype='S'))
        if network.zones:
            _, arr = network.abundance_matrix(names)
            h.create_dataset('abundances', data=arr)
            labels = np.array(['|'.join(z.label) for z in network.zones], dtype='S')
            h.create_dataset('zone_labels', data=labels)


def read_hdf5(path) -> Network:
    import h5py
    net = Network()
    with h5py.File(path, 'r') as h:
        names = [n.decode() if isinstance(n, bytes) else str(n) for n in h['species'][()]]
        for n in names:
            try: net.add_species(Species.parse(n))
            except Exception: pass
        if 'abundances' in h:
            arr = h['abundances'][()]
            labels = [x.decode() if isinstance(x, bytes) else str(x) for x in h.get('zone_labels', [])]
            for i, row in enumerate(arr):
                label = tuple((labels[i].split('|') if i < len(labels) else [str(i),'0','0'])[:3])
                while len(label) < 3: label = label + ('0',)
                z = Zone(label=label)
                for n, y in zip(names, row):
                    if y != 0: z.set_abundance(n, float(y))
                net.add_zone(z)
    return net
