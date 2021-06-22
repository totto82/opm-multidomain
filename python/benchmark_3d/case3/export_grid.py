import porepy as pp
import os

import case3.data as problem_data


def export_grid(path, mesh_id=0):
    case_folder = "case3/"
    mesh_sizes = ["mesh30k.geo", "mesh140k.geo", "mesh323k.geo", "mesh454k.geo"]

    file_geo = case_folder  + mesh_sizes[0]
    gb, domain = problem_data.create_grid(file_geo)
    gb = pp.utils.grid_utils.merge_grids_of_equal_dim(gb)
    pp.io.grid_writer.dump_grid_bucket_to_file(gb, os.path.join(path, "benchmark3d_3.txt"))
