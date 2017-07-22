#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'XML Partitioned Unstructured Grid Reader'
solution00000pvtu = XMLPartitionedUnstructuredGridReader(FileName=['/Users/ACG/Documents/PhD/Aspect/Data_Aspect_2012/Paper1/indentor/smooth_indentor/solution/solution-00000.pvtu'])
solution00000pvtu.PointArrayStatus = ['velocity', 'p', 'T', 'strain_rate', 'viscosity', 'viscosity_ratio']

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
solution00000pvtuDisplay = Show(solution00000pvtu, renderView1)
# trace defaults for the display properties.
solution00000pvtuDisplay.AmbientColor = [0.0, 0.0, 1.0]
solution00000pvtuDisplay.ColorArrayName = [None, '']
solution00000pvtuDisplay.EdgeColor = [0.0, 0.0, 1.0]
solution00000pvtuDisplay.GlyphType = 'Arrow'
solution00000pvtuDisplay.ScalarOpacityUnitDistance = 0.022097086912079615
solution00000pvtuDisplay.SetScaleArray = ['POINTS', 'T']
solution00000pvtuDisplay.ScaleTransferFunction = 'PiecewiseFunction'
solution00000pvtuDisplay.OpacityArray = ['POINTS', 'T']
solution00000pvtuDisplay.OpacityTransferFunction = 'PiecewiseFunction'

# reset view to fit data
renderView1.ResetCamera()

#changing interaction mode based on data extents
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.5, 0.5, 10000.0]
renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]

# set scalar coloring
ColorBy(solution00000pvtuDisplay, ('POINTS', 'p'))

# rescale color and/or opacity maps used to include current data range
solution00000pvtuDisplay.RescaleTransferFunctionToDataRange(True)

# show color bar/color legend
solution00000pvtuDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'p'
pLUT = GetColorTransferFunction('p')

# get opacity transfer function/opacity map for 'p'
pPWF = GetOpacityTransferFunction('p')

# rescale color and/or opacity maps used to exactly fit the current data range
solution00000pvtuDisplay.RescaleTransferFunctionToDataRange(False)

# create a new 'Plot Over Line'
plotOverLine1 = PlotOverLine(Input=solution00000pvtu,
    Source='High Resolution Line Source')

# init the 'High Resolution Line Source' selected for 'Source'
plotOverLine1.Source.Point2 = [1.0, 1.0, 0.0]

# Properties modified on plotOverLine1
plotOverLine1.Tolerance = 2.22044604925031e-16

# Properties modified on plotOverLine1.Source
plotOverLine1.Source.Point1 = [0.0, 1.0, 0.0]

# get layout
viewLayout1 = GetLayout()

# Create a new 'Line Chart View'
lineChartView1 = CreateView('XYChartView')
lineChartView1.ViewSize = [316, 638]

# place view in the layout
viewLayout1.AssignView(2, lineChartView1)

# show data in view
plotOverLine1Display = Show(plotOverLine1, lineChartView1)
# trace defaults for the display properties.
plotOverLine1Display.CompositeDataSetIndex = [0]
plotOverLine1Display.UseIndexForXAxis = 0
plotOverLine1Display.XArrayName = 'arc_length'
plotOverLine1Display.SeriesVisibility = ['p', 'strain_rate', 'T', 'velocity_Magnitude', 'viscosity', 'viscosity_ratio']
plotOverLine1Display.SeriesLabel = ['arc_length', 'arc_length', 'p', 'p', 'strain_rate', 'strain_rate', 'T', 'T', 'velocity_X', 'velocity_X', 'velocity_Y', 'velocity_Y', 'velocity_Z', 'velocity_Z', 'velocity_Magnitude', 'velocity_Magnitude', 'viscosity', 'viscosity', 'viscosity_ratio', 'viscosity_ratio', 'vtkValidPointMask', 'vtkValidPointMask', 'Points_X', 'Points_X', 'Points_Y', 'Points_Y', 'Points_Z', 'Points_Z', 'Points_Magnitude', 'Points_Magnitude']
plotOverLine1Display.SeriesColor = ['arc_length', '0', '0', '0', 'p', '0.89', '0.1', '0.11', 'strain_rate', '0.22', '0.49', '0.72', 'T', '0.3', '0.69', '0.29', 'velocity_X', '0.6', '0.31', '0.64', 'velocity_Y', '1', '0.5', '0', 'velocity_Z', '0.65', '0.34', '0.16', 'velocity_Magnitude', '0', '0', '0', 'viscosity', '0.89', '0.1', '0.11', 'viscosity_ratio', '0.22', '0.49', '0.72', 'vtkValidPointMask', '0.3', '0.69', '0.29', 'Points_X', '0.6', '0.31', '0.64', 'Points_Y', '1', '0.5', '0', 'Points_Z', '0.65', '0.34', '0.16', 'Points_Magnitude', '0', '0', '0']
plotOverLine1Display.SeriesPlotCorner = ['arc_length', '0', 'p', '0', 'strain_rate', '0', 'T', '0', 'velocity_X', '0', 'velocity_Y', '0', 'velocity_Z', '0', 'velocity_Magnitude', '0', 'viscosity', '0', 'viscosity_ratio', '0', 'vtkValidPointMask', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'Points_Magnitude', '0']
plotOverLine1Display.SeriesLineStyle = ['arc_length', '1', 'p', '1', 'strain_rate', '1', 'T', '1', 'velocity_X', '1', 'velocity_Y', '1', 'velocity_Z', '1', 'velocity_Magnitude', '1', 'viscosity', '1', 'viscosity_ratio', '1', 'vtkValidPointMask', '1', 'Points_X', '1', 'Points_Y', '1', 'Points_Z', '1', 'Points_Magnitude', '1']
plotOverLine1Display.SeriesLineThickness = ['arc_length', '2', 'p', '2', 'strain_rate', '2', 'T', '2', 'velocity_X', '2', 'velocity_Y', '2', 'velocity_Z', '2', 'velocity_Magnitude', '2', 'viscosity', '2', 'viscosity_ratio', '2', 'vtkValidPointMask', '2', 'Points_X', '2', 'Points_Y', '2', 'Points_Z', '2', 'Points_Magnitude', '2']
plotOverLine1Display.SeriesMarkerStyle = ['arc_length', '0', 'p', '0', 'strain_rate', '0', 'T', '0', 'velocity_X', '0', 'velocity_Y', '0', 'velocity_Z', '0', 'velocity_Magnitude', '0', 'viscosity', '0', 'viscosity_ratio', '0', 'vtkValidPointMask', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'Points_Magnitude', '0']

# Properties modified on plotOverLine1Display
plotOverLine1Display.SeriesVisibility = []
plotOverLine1Display.SeriesColor = ['arc_length', '0', '0', '0', 'p', '0.889998', '0.100008', '0.110002', 'strain_rate', '0.220005', '0.489998', '0.719997', 'T', '0.300008', '0.689998', '0.289998', 'velocity_X', '0.6', '0.310002', '0.639994', 'velocity_Y', '1', '0.500008', '0', 'velocity_Z', '0.650004', '0.340002', '0.160006', 'velocity_Magnitude', '0', '0', '0', 'viscosity', '0.889998', '0.100008', '0.110002', 'viscosity_ratio', '0.220005', '0.489998', '0.719997', 'vtkValidPointMask', '0.300008', '0.689998', '0.289998', 'Points_X', '0.6', '0.310002', '0.639994', 'Points_Y', '1', '0.500008', '0', 'Points_Z', '0.650004', '0.340002', '0.160006', 'Points_Magnitude', '0', '0', '0']
plotOverLine1Display.SeriesPlotCorner = ['Points_Magnitude', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'T', '0', 'arc_length', '0', 'p', '0', 'strain_rate', '0', 'velocity_Magnitude', '0', 'velocity_X', '0', 'velocity_Y', '0', 'velocity_Z', '0', 'viscosity', '0', 'viscosity_ratio', '0', 'vtkValidPointMask', '0']
plotOverLine1Display.SeriesLineStyle = ['Points_Magnitude', '1', 'Points_X', '1', 'Points_Y', '1', 'Points_Z', '1', 'T', '1', 'arc_length', '1', 'p', '1', 'strain_rate', '1', 'velocity_Magnitude', '1', 'velocity_X', '1', 'velocity_Y', '1', 'velocity_Z', '1', 'viscosity', '1', 'viscosity_ratio', '1', 'vtkValidPointMask', '1']
plotOverLine1Display.SeriesLineThickness = ['Points_Magnitude', '2', 'Points_X', '2', 'Points_Y', '2', 'Points_Z', '2', 'T', '2', 'arc_length', '2', 'p', '2', 'strain_rate', '2', 'velocity_Magnitude', '2', 'velocity_X', '2', 'velocity_Y', '2', 'velocity_Z', '2', 'viscosity', '2', 'viscosity_ratio', '2', 'vtkValidPointMask', '2']
plotOverLine1Display.SeriesMarkerStyle = ['Points_Magnitude', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'T', '0', 'arc_length', '0', 'p', '0', 'strain_rate', '0', 'velocity_Magnitude', '0', 'velocity_X', '0', 'velocity_Y', '0', 'velocity_Z', '0', 'viscosity', '0', 'viscosity_ratio', '0', 'vtkValidPointMask', '0']

# Properties modified on plotOverLine1Display
plotOverLine1Display.SeriesVisibility = ['p']

# Properties modified on plotOverLine1.Source
plotOverLine1.Source.Resolution = 513

# save data
SaveData('/Users/ACG/Documents/PhD/Aspect/Data_Aspect_2012/Paper1/indentor/smooth_indentor/surface_p.csv', proxy=plotOverLine1)
