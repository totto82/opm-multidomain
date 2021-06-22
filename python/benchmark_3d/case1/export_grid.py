import porepy as pp
import os

import case1.data as problem_data


def export_grid(path, mesh_id=1):
    case_folder = "case1/"
    mesh_sizes = ["mesh1k.geo", "mesh10k.geo", "mesh100k.geo"]

    file_geo = case_folder  + mesh_sizes[mesh_id]
    gb, domain = problem_data.import_grid(file_geo, 1e-8)

    pp.io.grid_writer.dump_grid_bucket_to_file(gb, os.path.join(path, "benchmark3d_1.txt"))
