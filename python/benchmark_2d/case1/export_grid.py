import porepy as pp
import os

import case1.data as problem_data

def export_grid(path, cartesian=False):
    mesh_size = 0.025
    data = {"mesh_size": mesh_size}

    gb = problem_data.import_grid(data, cartesian)
    gb = pp.utils.grid_utils.merge_grids_of_equal_dim(gb)

    pp.io.grid_writer.dump_grid_bucket_to_file(gb, os.path.join(path, "benchmark1.txt"))
