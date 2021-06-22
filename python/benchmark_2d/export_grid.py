import case1.export_grid as case1

import case3.export_grid as case3
import case4.export_grid as case4

path = "../../tests/data"

print("Creating and exporting grid for case 1")
case1.export_grid(path, cartesian=False)
print("Creating and exporting grid for case 3")
case3.export_grid(path)
print("Creating and exporting grid for case 4")
case4.export_grid(path)
print("Done")
