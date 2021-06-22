import numpy as np
import porepy as pp


def import_grid(data, cartesian=False):
    mesh_kwargs = {}
    mesh_kwargs = {
        "mesh_size_frac": data["mesh_size"],
        "mesh_size_min": data["mesh_size"] / 20,
    }

    domain = {"xmin": 0, "xmax": 1, "ymin": 0, "ymax": 1}

    file_name = "case1/network.csv"
    network = pp.fracture_importer.network_2d_from_csv(file_name)
    if cartesian:
        pts = network.pts
        fracs= []
        for i in range(0, pts.shape[1], 2):
            fracs.append(pts[:, i:(i+2)])
        domain = [self.domain['xmax'], self.domain['ymax']]
        nx = [10, 10]
        gb = pp.meshing.cart_grid(fracs, nx, physdims=domain)
    else:
        gb = network.mesh(mesh_kwargs, domain=domain)
    gb.compute_geometry()

    return gb

