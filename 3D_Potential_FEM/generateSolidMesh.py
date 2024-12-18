from paraview.vtk.numpy_interface import dataset_adapter as dsa
from vtk.util import numpy_support as nps

import ReadAbaqusMesh as mesh
# import rectModes as modes

POINT = True

nodes, elms, bndry = mesh.readAbq(r'C:\Users\cfsan\Desktop\Random_Stuff\VIP\CQEM-FEM\3D_Potential_FEM\1x0.5x0.75cm.inp')
with open(r'C:\Users\cfsan\Desktop\Random_Stuff\VIP\CQEM-FEM\3D_Potential_FEM\test.npy', 'rb') as f:
    fields = numpy.load(f)

print(fields)

points = vtk.vtkPoints()

# Create a cell array to store the tetrahedron
cellArray = vtk.vtkCellArray()
cellArray.Allocate(len(elms), 1000)

# Placeholder point to ensure indexing aligns with node id
points.InsertNextPoint(0, 0, 0)

for node in nodes:
    points.InsertNextPoint(node[1], node[2], node[3])

for elm in elms:
    tetra = vtk.vtkTetra()
    for i in range(4):
        tetra.GetPointIds().SetId(i, elm[i])
    cellArray.InsertNextCell(tetra)

unstructuredGrid = vtk.vtkUnstructuredGrid()
unstructuredGrid.SetPoints(points)
unstructuredGrid.SetCells(vtk.vtkTetra().GetCellType(), cellArray)

# assuming fields array is in order, no gaps

#for i in range(unstructuredGrid.GetNumberOfCells()):
    #H = numpy.append(field, modes.computeH(centroid[0], centroid[1], centroid[2]))
#H = numpy.column_stack((Hx, Hy, Hz))
vtk_data_array = nps.numpy_to_vtk(fields[:,1].ravel(), deep=True)
vtk_data_array.SetNumberOfComponents(fields.shape[1] - 1)
vtk_data_array.SetName('fields')
print(vtk_data_array)
celldata = unstructuredGrid.GetCellData()
celldata.AddArray(vtk_data_array)
# Create an unstructured grid to hold the cell array


#unstructuredGrid.InsertNextCell(tetra.GetCellType(), tetra.GetPointIds())

#print(unstructuredGrid.GetBounds())

self.GetOutput().ShallowCopy(unstructuredGrid)
