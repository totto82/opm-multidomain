import numpy as np
import porepy as pp

def import_grid(file_geo, tol):
    gb = pp.fracture_importer.dfm_from_gmsh(file_geo, 3)
    gb.compute_geometry()

    return gb


