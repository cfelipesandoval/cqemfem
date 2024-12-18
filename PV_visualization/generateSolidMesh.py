from paraview.vtk.numpy_interface import dataset_adapter as dsa
from vtk.util import numpy_support as nps

# load in mesh and data

nodes, elms, bndry = readAbq("MESH_FILE")
with open('DATA_FILE', 'rb') as f:
    fields = numpy.load(f)

# create points
points = vtk.vtkPoints()

# Create a cell array to store the tets
cellArray = vtk.vtkCellArray()
cellArray.Allocate(len(elms), 1000)

# Placeholder point to ensure indexing aligns with node id
points.InsertNextPoint(0, 0, 0)

# add nodes to vtk mesh
for node in nodes:
    points.InsertNextPoint(node[1], node[2], node[3])

# add tets to vtk mesh
for elm in elms:
    tetra = vtk.vtkTetra()
    for i in range(4):
        tetra.GetPointIds().SetId(i, elm[i])
    cellArray.InsertNextCell(tetra)

# set up mesh
unstructuredGrid = vtk.vtkUnstructuredGrid()
unstructuredGrid.SetPoints(points)
unstructuredGrid.SetCells(vtk.vtkTetra().GetCellType(), cellArray)

# assuming fields array is in order, no gaps
num_components = fields.shape[1] - 1

# convert numpy array and assign the vtk data

vtk_data_array = nps.numpy_to_vtk(fields[:,1:(num_components + 1)].ravel(), deep=True)
vtk_data_array.SetNumberOfComponents(num_components)
vtk_data_array.SetName('fields')
celldata = unstructuredGrid.GetCellData()
celldata.AddArray(vtk_data_array)
# Create an unstructured grid to hold the cell array

# output of the programmable source is the created unstructured grid
self.GetOutput().ShallowCopy(unstructuredGrid)
