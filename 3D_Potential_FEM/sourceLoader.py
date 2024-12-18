from paraview.simple import *

POINT = False
GLYPH = False

# Function to read an external Python script
def read_script(file_path):
    with open(file_path, 'r') as file:
        return file.read()

# Create a programmable source
prog_source = ProgrammableSource()

# Set the output type
prog_source.OutputDataSetType = 'vtkUnstructuredGrid'  # or 'vtkUnstructuredGrid', etc.

# Read the external Python script and assign it to the Script property
prog_source.Script = read_script(r'C:\Users\cfsan\Desktop\Random_Stuff\VIP\CQEM-FEM\3D_Potential_FEM\generateSolidMesh.py')
print(prog_source)

# Update the programmable source to generate the output data
prog_source.UpdatePipeline()

renderView1 = GetActiveViewOrCreate('RenderView')

if not POINT:
    programmableSource1Display = GetDisplayProperties(prog_source, view=renderView1)
    #ColorBy(programmableSource1Display, ('CELLS', 'fields'))
    #ColorBy(programmableSource1Display, ('CELLS', 'fields', 'Magnitude'))
    cellDatatoPointData1 = CellDatatoPointData(registrationName='CellDatatoPointData1', Input=prog_source)

    display = Show(cellDatatoPointData1, renderView1, 'UnstructuredGridRepresentation')

    Hide(prog_source, GetActiveView())
        
    display.Representation = 'Volume'
    display.RescaleTransferFunctionToDataRange(False, True)
    ColorBy(display, ('POINTS', 'fields'))

    bounds = prog_source.GetDataInformation().GetBounds()
    #display.ScalarOpacityUnitDistance = 1
    Hide(cellDatatoPointData1, GetActiveView())

    slice1 = Slice(registrationName='Slice1', Input=cellDatatoPointData1)
    slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')
    slice1Display.Representation = 'Surface'
    slice1Display.SetScalarBarVisibility(renderView1, True)

    slice1.SliceType.Origin = [(bounds[0] + bounds[1]) / 2, (bounds[2] + bounds[3]) / 2, (bounds[4] + bounds[5]) / 2]
    slice1.SliceType.Normal = [1.0, 0.0, 0.0]

    slice2 = Slice(registrationName='Slice2', Input=cellDatatoPointData1)
    slice2Display = Show(slice2, renderView1, 'GeometryRepresentation')
    slice2Display.Representation = 'Surface'
    slice2Display.SetScalarBarVisibility(renderView1, True)

    slice2.SliceType.Origin = [(bounds[0] + bounds[1]) / 2, (bounds[2] + bounds[3]) / 2, (bounds[4] + bounds[5]) / 2]
    slice2.SliceType.Normal = [0.0, 0.0, 1.0]

    GetActiveView().OrientationAxesVisibility = 0
    GetActiveView().Background = [1.0, 0.0, 0.4980392156862745]

    """
    slice3 = Slice(registrationName='Slice3', Input=cellDatatoPointData1)
    slice3Display = Show(slice3, renderView1, 'GeometryRepresentation')
    slice3Display.Representation = 'Surface'
    slice3Display.SetScalarBarVisibility(renderView1, True)

    slice3.SliceType.Origin = [(bounds[0] + bounds[1]) / 2, (bounds[2] + bounds[3]) / 2, (bounds[4] + bounds[5]) / 2]
    slice3.SliceType.Normal = [0.0, 1.0, 0.0]
    """
    #ColorBy(display, ('POINTS', 'H', 'Magnitude'))
else:
    display = GetDisplayProperties(prog_source, view=renderView1)
    display.Representation = 'Volume'
    ColorBy(display, ('POINTS', 'H', 'Magnitude'))

if GLYPH:

    glyph1 = Glyph(registrationName='Glyph1', Input=GetActiveSource(),
        GlyphType='Arrow')
    glyph1.ScaleArray = ['POINTS', 'H']
    glyph1.OrientationArray = ['POINTS', 'H']
    glyph1Display = Show(glyph1, renderView1, 'GeometryRepresentation')
    glyph1Display.Representation = 'Surface'
    glyph1Display.SetScalarBarVisibility(renderView1, True)
    renderView1.Update()
#Show()
#Render()
#RenderAllViews()

Interact()
