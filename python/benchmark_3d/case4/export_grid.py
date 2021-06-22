import os
import porepy as pp

import case4.data as problem_data

def export_grid(path):
    gb, domain = problem_data.create_grid(False, True)

    gb = pp.utils.grid_utils.merge_grids_of_equal_dim(gb)
    pp.io.grid_writer.dump_grid_bucket_to_file(gb, os.path.join(path, "benchmark3d_4.txt"))
