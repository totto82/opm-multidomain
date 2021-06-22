import case1.export_grid as case1
import case2.export_grid as case2
import case3.export_grid as case3
import case4.export_grid as case4

path = "../../tests/data"

print("Creating and exporting grid for case 1")
case1.export_grid(path, mesh_id=1)
print("Creating and exporting grid for case 2")
case2.export_grid(path, mesh_id=1)
print("Creating and exporting grid for case 3")
case3.export_grid(path, mesh_id=0)
print("Creating and exporting grid for case 4")
case4.export_grid(path)
print("Done")
