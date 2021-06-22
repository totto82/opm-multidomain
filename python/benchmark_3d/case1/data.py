import numpy as np
import porepy as pp

# ------------------------------------------------------------------------------#


def import_grid(file_geo, tol):

    domain = {"xmin": 0, "xmax": 100, "ymin": 0, "ymax": 100, "zmin": 0, "zmax": 100}

    gb = pp.fracture_importer.dfm_from_gmsh(file_geo, 3)
    gb.compute_geometry()

    return gb, domain
