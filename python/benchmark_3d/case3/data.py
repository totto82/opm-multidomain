import pickle
import os

import numpy as np
import porepy as pp

# ------------------------------------------------------------------------------#
def create_grid(fn):

    domain = {"xmin": 0, "xmax": 1, "ymin": 0, "ymax": 2.25, "zmin": 0, "zmax": 1}
    _, ext = os.path.splitext(fn)

    if ext == ".geo":
        gb = pp.fracture_importer.dfm_from_gmsh(fn, 3)
    elif ext == ".grid":
        gb = pickle.load(open(fn, "rb"))
    else:
        raise ValueError("Not supported data format")

    return gb, domain
