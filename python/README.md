# fracture_benchmark_grids
Use PorePy (https://github.com/pmgbergen/porepy) to generate grids for the
fracture benchmark papers:

Flemish et.al., Benchmarks for single-phase flow in fractured porous media, doi: 10.1016/j.advwatres.2017.10.036

Berre et.al., Verification benchmarks for single-phase flow in three-dimensional fractured porous media, doi: 10.1016/j.advwatres.2020.103759


# USE

This folder contains two subfolders: benchmark_2d and benchmark_3d. The first contains the
files necessary to create the 2d benchmark grids, the second the 3d benchmark files.

To create the grids of a benchmark, go to the desired subfolder. Both subfolder contains
a file export_grids.py. Run this with

python export_grids.py

The exported grids will be placed in opm-multidomain/tests/data/


The code has been validated for PorePy commit: e9f443e509000b57be2549a71399e668cf568a80