import pickle

import numpy as np
import porepy as pp


def create_grid(from_file=True, generate_network=False, tol=1e-4):
    """ Obtain domain and grid bucket. Default is to load a pickled bucket;
    alternatively, a .geo file is available.
    """
    if generate_network:
        file_csv = "case4/fracture_network.csv"
        domain = {
            "xmin": -500,
            "xmax": 350,
            "ymin": 100,
            "ymax": 1500,
            "zmin": -100,
            "zmax": 500,
        }

        network = pp.fracture_importer.network_3d_from_csv(
            file_csv, has_domain=False, tol=tol
        )
        network.impose_external_boundary(domain)
        network.find_intersections()
        network.split_intersections()

        pickle.dump(network, open("case4/network_52_fracs", "wb"))

    network = pickle.load(open("case4/network_52_fracs", "rb"))
    domain = network.domain
    if from_file:
#        gb = pickle.load(open("gridbucket_case4.grid", "rb"))
        gb = pp.fracture_importer.dfm_from_gmsh("case4/gmsh.msh", 3)
    else:
        gb = network.mesh({'mesh_size_frac':20, 'mesh_size_min': 20})#20
        # gb = pp.fracture_importer.dfm_from_gmsh(
        #     "gmsh_frac_file.msh", 3
        # )
        pickle.dump(gb, open("case4/gridbucket_case4.grid", "wb"))

    return gb, domain


