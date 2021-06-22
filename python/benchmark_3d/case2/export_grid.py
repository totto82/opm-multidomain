import porepy as pp
import os

import case2.data as problem_data


def export_grid(path, mesh_id=1):
    case_folder = "case2/"
    mesh_sizes = ["mesh500.geo", "mesh4k.geo", "mesh32k.geo"]

    file_geo = case_folder  + mesh_sizes[mesh_id]

    gb = problem_data.import_grid(file_geo, 1e-8)
    gb = pp.utils.grid_utils.merge_grids_of_equal_dim(gb)

    pp.io.grid_writer.dump_grid_bucket_to_file(gb, os.path.join(path, "benchmark3d_2.txt"))
