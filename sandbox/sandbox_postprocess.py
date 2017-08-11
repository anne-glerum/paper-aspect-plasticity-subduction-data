#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'PVD Reader'
solutionpvd = PVDReader(FileName='sandbox/restart_sandbox/solution.pvd')

# find source
solutionpvd_1 = FindSource('solution.pvd')

LoadState("sandbox/Correct_color_scales.pvsm")

# find source
annotateTimeFilter1 = FindSource('AnnotateTimeFilter1')

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
renderView1.ViewSize = [1087, 758]

# show data in view
solutionpvdDisplay = Show(solutionpvd, renderView1)
# trace defaults for the display properties.
solutionpvdDisplay.CubeAxesVisibility = 0
solutionpvdDisplay.Representation = 'Surface'
solutionpvdDisplay.AmbientColor = [0.0, 0.0, 1.0]
solutionpvdDisplay.ColorArrayName = [None, '']
solutionpvdDisplay.DiffuseColor = [1.0, 1.0, 1.0]
solutionpvdDisplay.LookupTable = None
solutionpvdDisplay.MapScalars = 1
solutionpvdDisplay.InterpolateScalarsBeforeMapping = 1
solutionpvdDisplay.Opacity = 1.0
solutionpvdDisplay.PointSize = 2.0
solutionpvdDisplay.LineWidth = 1.0
solutionpvdDisplay.Interpolation = 'Gouraud'
solutionpvdDisplay.Specular = 0.0
solutionpvdDisplay.SpecularColor = [1.0, 1.0, 1.0]
solutionpvdDisplay.SpecularPower = 100.0
solutionpvdDisplay.Ambient = 0.0
solutionpvdDisplay.Diffuse = 1.0
solutionpvdDisplay.EdgeColor = [0.0, 0.0, 1.0]
solutionpvdDisplay.BackfaceRepresentation = 'Follow Frontface'
solutionpvdDisplay.BackfaceAmbientColor = [1.0, 1.0, 1.0]
solutionpvdDisplay.BackfaceDiffuseColor = [1.0, 1.0, 1.0]
solutionpvdDisplay.BackfaceOpacity = 1.0
solutionpvdDisplay.Position = [0.0, 0.0, 0.0]
solutionpvdDisplay.Scale = [1.0, 1.0, 1.0]
solutionpvdDisplay.Orientation = [0.0, 0.0, 0.0]
solutionpvdDisplay.Origin = [0.0, 0.0, 0.0]
solutionpvdDisplay.Pickable = 1
solutionpvdDisplay.Texture = None
solutionpvdDisplay.Triangulate = 0
solutionpvdDisplay.NonlinearSubdivisionLevel = 1
solutionpvdDisplay.GlyphType = 'Arrow'
solutionpvdDisplay.CubeAxesColor = [0.0, 0.0, 1.0]
solutionpvdDisplay.CubeAxesCornerOffset = 0.0
solutionpvdDisplay.CubeAxesFlyMode = 'Closest Triad'
solutionpvdDisplay.CubeAxesInertia = 1
solutionpvdDisplay.CubeAxesTickLocation = 'Inside'
solutionpvdDisplay.CubeAxesXAxisMinorTickVisibility = 1
solutionpvdDisplay.CubeAxesXAxisTickVisibility = 1
solutionpvdDisplay.CubeAxesXAxisVisibility = 1
solutionpvdDisplay.CubeAxesXGridLines = 0
solutionpvdDisplay.CubeAxesXTitle = 'X-Axis'
solutionpvdDisplay.CubeAxesUseDefaultXTitle = 1
solutionpvdDisplay.CubeAxesYAxisMinorTickVisibility = 1
solutionpvdDisplay.CubeAxesYAxisTickVisibility = 1
solutionpvdDisplay.CubeAxesYAxisVisibility = 1
solutionpvdDisplay.CubeAxesYGridLines = 0
solutionpvdDisplay.CubeAxesYTitle = 'Y-Axis'
solutionpvdDisplay.CubeAxesUseDefaultYTitle = 1
solutionpvdDisplay.CubeAxesZAxisMinorTickVisibility = 1
solutionpvdDisplay.CubeAxesZAxisTickVisibility = 1
solutionpvdDisplay.CubeAxesZAxisVisibility = 1
solutionpvdDisplay.CubeAxesZGridLines = 0
solutionpvdDisplay.CubeAxesZTitle = 'Z-Axis'
solutionpvdDisplay.CubeAxesUseDefaultZTitle = 1
solutionpvdDisplay.CubeAxesGridLineLocation = 'All Faces'
solutionpvdDisplay.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
solutionpvdDisplay.CustomBoundsActive = [0, 0, 0]
solutionpvdDisplay.OriginalBoundsRangeActive = [0, 0, 0]
solutionpvdDisplay.CustomRange = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
solutionpvdDisplay.CustomRangeActive = [0, 0, 0]
solutionpvdDisplay.UseAxesOrigin = 0
solutionpvdDisplay.AxesOrigin = [0.0, 0.0, 0.0]
solutionpvdDisplay.CubeAxesXLabelFormat = '%-#6.3g'
solutionpvdDisplay.CubeAxesYLabelFormat = '%-#6.3g'
solutionpvdDisplay.CubeAxesZLabelFormat = '%-#6.3g'
solutionpvdDisplay.StickyAxes = 0
solutionpvdDisplay.CenterStickyAxes = 0
solutionpvdDisplay.SelectionCellLabelBold = 0
solutionpvdDisplay.SelectionCellLabelColor = [0.0, 1.0, 0.0]
solutionpvdDisplay.SelectionCellLabelFontFamily = 'Arial'
solutionpvdDisplay.SelectionCellLabelFontSize = 18
solutionpvdDisplay.SelectionCellLabelItalic = 0
solutionpvdDisplay.SelectionCellLabelJustification = 'Left'
solutionpvdDisplay.SelectionCellLabelOpacity = 1.0
solutionpvdDisplay.SelectionCellLabelShadow = 0
solutionpvdDisplay.SelectionPointLabelBold = 0
solutionpvdDisplay.SelectionPointLabelColor = [1.0, 1.0, 0.0]
solutionpvdDisplay.SelectionPointLabelFontFamily = 'Arial'
solutionpvdDisplay.SelectionPointLabelFontSize = 18
solutionpvdDisplay.SelectionPointLabelItalic = 0
solutionpvdDisplay.SelectionPointLabelJustification = 'Left'
solutionpvdDisplay.SelectionPointLabelOpacity = 1.0
solutionpvdDisplay.SelectionPointLabelShadow = 0
solutionpvdDisplay.ScalarOpacityUnitDistance = 0.00817991466739737
solutionpvdDisplay.SelectMapper = 'Projected tetra'
solutionpvdDisplay.GaussianRadius = 0.0
solutionpvdDisplay.ShaderPreset = 'Sphere'
solutionpvdDisplay.Emissive = 0
solutionpvdDisplay.ScaleByArray = 0
solutionpvdDisplay.SetScaleArray = ['POINTS', 'C_1']
solutionpvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
solutionpvdDisplay.OpacityByArray = 0
solutionpvdDisplay.OpacityArray = ['POINTS', 'C_1']
solutionpvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'Arrow' selected for 'GlyphType'
solutionpvdDisplay.GlyphType.TipResolution = 6
solutionpvdDisplay.GlyphType.TipRadius = 0.1
solutionpvdDisplay.GlyphType.TipLength = 0.35
solutionpvdDisplay.GlyphType.ShaftResolution = 6
solutionpvdDisplay.GlyphType.ShaftRadius = 0.03
solutionpvdDisplay.GlyphType.Invert = 0

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
solutionpvdDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
solutionpvdDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# set scalar coloring
ColorBy(solutionpvdDisplay, ('POINTS', 'density'))

# rescale color and/or opacity maps used to include current data range
solutionpvdDisplay.RescaleTransferFunctionToDataRange(True)

# show color bar/color legend
solutionpvdDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'density'
densityLUT = GetColorTransferFunction('density')

# get opacity transfer function/opacity map for 'density'
densityPWF = GetOpacityTransferFunction('density')

# hide data in view
Hide(solutionpvd_1, renderView1)

# hide data in view
Hide(annotateTimeFilter1, renderView1)

# create a new 'Annotate Time Filter'
annotateTimeFilter2 = AnnotateTimeFilter(Input=solutionpvd)
annotateTimeFilter2.Format = 'Time: %f'
annotateTimeFilter2.Shift = 0.0
annotateTimeFilter2.Scale = 1.0

# set active source
SetActiveSource(annotateTimeFilter1)

# set active source
SetActiveSource(annotateTimeFilter2)

# set active source
SetActiveSource(annotateTimeFilter1)

# set active source
SetActiveSource(annotateTimeFilter2)

# Properties modified on annotateTimeFilter2
annotateTimeFilter2.Format = 'Time: %2.2f h'
annotateTimeFilter2.Scale = 0.00027777777

# show data in view
annotateTimeFilter2Display = Show(annotateTimeFilter2, renderView1)
# trace defaults for the display properties.
annotateTimeFilter2Display.Interactivity = 1
annotateTimeFilter2Display.Color = [1.0, 1.0, 1.0]
annotateTimeFilter2Display.Opacity = 1.0
annotateTimeFilter2Display.FontFamily = 'Arial'
annotateTimeFilter2Display.Bold = 0
annotateTimeFilter2Display.Italic = 0
annotateTimeFilter2Display.Shadow = 0
annotateTimeFilter2Display.FontSize = 18
annotateTimeFilter2Display.Justification = 'Left'
annotateTimeFilter2Display.WindowLocation = 'AnyLocation'
annotateTimeFilter2Display.Position = [0.05, 0.05]

# Properties modified on annotateTimeFilter2Display
annotateTimeFilter2Display.Color = [0.0, 0.0, 0.0]

# set active source
SetActiveSource(annotateTimeFilter1)

# set active source
SetActiveSource(annotateTimeFilter2)

# Properties modified on annotateTimeFilter2Display
annotateTimeFilter2Display.FontSize = 17

# Properties modified on annotateTimeFilter2Display
annotateTimeFilter2Display.FontSize = 16

# Properties modified on annotateTimeFilter2Display
annotateTimeFilter2Display.FontSize = 15

# Properties modified on annotateTimeFilter2Display
annotateTimeFilter2Display.FontSize = 14

# Properties modified on annotateTimeFilter2Display
annotateTimeFilter2Display.FontSize = 13

# Properties modified on annotateTimeFilter2Display
annotateTimeFilter2Display.FontSize = 12

# Properties modified on annotateTimeFilter2Display
annotateTimeFilter2Display.FontSize = 11

# Properties modified on annotateTimeFilter2Display
annotateTimeFilter2Display.FontSize = 10

# set active source
SetActiveSource(solutionpvd)

# hide color bar/color legend
solutionpvdDisplay.SetScalarBarVisibility(renderView1, False)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.1, 0.025, 0.39826142083018456]
renderView1.CameraFocalPoint = [0.1, 0.025, 0.0]
renderView1.CameraParallelScale = 0.07040341550470697

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.1, 0.025, 0.39826142083018456]
renderView1.CameraFocalPoint = [0.1, 0.025, 0.0]
renderView1.CameraParallelScale = 0.07040341550470697

# save animation images/movie
WriteAnimation('sandbox/restart_sandbox/Composition.avi', Magnification=1, FrameRate=10.0, Compression=True)

# set scalar coloring
ColorBy(solutionpvdDisplay, ('POINTS', 'strain_rate'))

# rescale color and/or opacity maps used to include current data range
solutionpvdDisplay.RescaleTransferFunctionToDataRange(True)

# show color bar/color legend
solutionpvdDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'strainrate'
strainrateLUT = GetColorTransferFunction('strainrate')

# get opacity transfer function/opacity map for 'strainrate'
strainratePWF = GetOpacityTransferFunction('strainrate')

# Rescale transfer function
strainrateLUT.RescaleTransferFunction(1e-09, 0.01)

# Rescale transfer function
strainratePWF.RescaleTransferFunction(1e-09, 0.01)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.1, 0.025, 0.39826142083018456]
renderView1.CameraFocalPoint = [0.1, 0.025, 0.0]
renderView1.CameraParallelScale = 0.07040341550470697

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.1, 0.025, 0.39826142083018456]
renderView1.CameraFocalPoint = [0.1, 0.025, 0.0]
renderView1.CameraParallelScale = 0.07040341550470697

# save animation images/movie
WriteAnimation('sandbox/restart_sandbox/Strainrate.avi', Magnification=1, FrameRate=10.0, Compression=True)

# hide color bar/color legend
solutionpvdDisplay.SetScalarBarVisibility(renderView1, False)

# hide data in view
Hide(annotateTimeFilter2, renderView1)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.1, 0.025, 0.39826142083018456]
renderView1.CameraFocalPoint = [0.1, 0.025, 0.0]
renderView1.CameraParallelScale = 0.07040341550470697

# save screenshot
SaveScreenshot('sandbox/restart_sandbox/Strainrate_1cm.png', magnification=1, quality=100, view=renderView1)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.1, 0.025, 0.39826142083018456]
renderView1.CameraFocalPoint = [0.1, 0.025, 0.0]
renderView1.CameraParallelScale = 0.07040341550470697

# save screenshot
SaveScreenshot('sandbox/restart_sandbox/Strainrate_2cm.png', magnification=1, quality=100, view=renderView1)

# set scalar coloring
ColorBy(solutionpvdDisplay, ('POINTS', 'viscosity'))

# rescale color and/or opacity maps used to include current data range
solutionpvdDisplay.RescaleTransferFunctionToDataRange(True)

# show color bar/color legend
solutionpvdDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'viscosity'
viscosityLUT = GetColorTransferFunction('viscosity')

# get opacity transfer function/opacity map for 'viscosity'
viscosityPWF = GetOpacityTransferFunction('viscosity')

# Rescale transfer function
viscosityLUT.RescaleTransferFunction(100.0, 10000000.0)

# Rescale transfer function
viscosityPWF.RescaleTransferFunction(100.0, 10000000.0)

# Properties modified on viscosityLUT
viscosityLUT.NumberOfTableValues = 7

# get color legend/bar for viscosityLUT in view renderView1
viscosityLUTColorBar = GetScalarBar(viscosityLUT, renderView1)

# Properties modified on viscosityLUTColorBar
viscosityLUTColorBar.Title = 'Viscosity [Pas]'
viscosityLUTColorBar.TitleColor = [0.0, 0.0, 0.0]
viscosityLUTColorBar.LabelColor = [0.0, 0.0, 0.0]

# Properties modified on viscosityLUTColorBar
viscosityLUTColorBar.AutoOrient = 0
viscosityLUTColorBar.Orientation = 'Horizontal'

# Properties modified on viscosityLUTColorBar
viscosityLUTColorBar.AutomaticLabelFormat = 0
viscosityLUTColorBar.LabelFormat = '%-#3.3e'

# Properties modified on viscosityLUTColorBar
viscosityLUTColorBar.LabelFormat = '%-#3.0e'

# Properties modified on viscosityLUTColorBar
viscosityLUTColorBar.RangeLabelFormat = '%3.0e'

# Properties modified on viscosityLUTColorBar
viscosityLUTColorBar.LabelFormat = '%3.0e'

# Properties modified on viscosityLUTColorBar
viscosityLUTColorBar.LabelFormat = '3.0e'

# Properties modified on viscosityLUTColorBar
viscosityLUTColorBar.LabelFormat = '%3.0e'

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.1, 0.025, 0.39826142083018456]
renderView1.CameraFocalPoint = [0.1, 0.025, 0.0]
renderView1.CameraParallelScale = 0.07040341550470697

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.1, 0.025, 0.39826142083018456]
renderView1.CameraFocalPoint = [0.1, 0.025, 0.0]
renderView1.CameraParallelScale = 0.07040341550470697

# save animation images/movie
WriteAnimation('sandbox/restart_sandbox/Viscosity.avi', Magnification=1, FrameRate=10.0, Compression=True)

# Properties modified on viscosityLUT
viscosityLUT.NumberOfTableValues = 5

# set active source
SetActiveSource(annotateTimeFilter2)

# show data in view
annotateTimeFilter2Display = Show(annotateTimeFilter2, renderView1)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.1, 0.025, 0.39826142083018456]
renderView1.CameraFocalPoint = [0.1, 0.025, 0.0]
renderView1.CameraParallelScale = 0.07040341550470697

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.1, 0.025, 0.39826142083018456]
renderView1.CameraFocalPoint = [0.1, 0.025, 0.0]
renderView1.CameraParallelScale = 0.07040341550470697

# save animation images/movie
WriteAnimation('sandbox/restart_sandbox/Viscosity.avi', Magnification=1, FrameRate=10.0, Compression=True)

# set active source
SetActiveSource(annotateTimeFilter2)

# set active source
SetActiveSource(solutionpvd)

# set scalar coloring
ColorBy(solutionpvdDisplay, ('POINTS', 'p'))

# rescale color and/or opacity maps used to include current data range
solutionpvdDisplay.RescaleTransferFunctionToDataRange(True)

# show color bar/color legend
solutionpvdDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'p'
pLUT = GetColorTransferFunction('p')

# get opacity transfer function/opacity map for 'p'
pPWF = GetOpacityTransferFunction('p')

# Rescale transfer function
pLUT.RescaleTransferFunction(0.0, 350.0)

# Rescale transfer function
pPWF.RescaleTransferFunction(0.0, 350.0)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.1, 0.025, 0.39826142083018456]
renderView1.CameraFocalPoint = [0.1, 0.025, 0.0]
renderView1.CameraParallelScale = 0.07040341550470697

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.1, 0.025, 0.39826142083018456]
renderView1.CameraFocalPoint = [0.1, 0.025, 0.0]
renderView1.CameraParallelScale = 0.07040341550470697

# save animation images/movie
WriteAnimation('sandbox/restart_sandbox/Pressure.avi', Magnification=1, FrameRate=10.0, Compression=True)

# hide color bar/color legend
solutionpvdDisplay.SetScalarBarVisibility(renderView1, False)

# hide data in view
Hide(annotateTimeFilter2, renderView1)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.1, 0.025, 0.39826142083018456]
renderView1.CameraFocalPoint = [0.1, 0.025, 0.0]
renderView1.CameraParallelScale = 0.07040341550470697

# save screenshot
SaveScreenshot('sandbox/restart_sandbox/Pressure_2cm.png', magnification=1, quality=100, view=renderView1)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.1, 0.025, 0.39826142083018456]
renderView1.CameraFocalPoint = [0.1, 0.025, 0.0]
renderView1.CameraParallelScale = 0.07040341550470697

# save screenshot
SaveScreenshot('sandbox/restart_sandbox/Pressure_1cm.png', magnification=1, quality=100, view=renderView1)

# set scalar coloring
ColorBy(solutionpvdDisplay, ('POINTS', 'viscosity'))

# rescale color and/or opacity maps used to include current data range
solutionpvdDisplay.RescaleTransferFunctionToDataRange(True)

# show color bar/color legend
solutionpvdDisplay.SetScalarBarVisibility(renderView1, True)

# Rescale transfer function
viscosityLUT.RescaleTransferFunction(100.0, 10000000.0)

# Rescale transfer function
viscosityPWF.RescaleTransferFunction(100.0, 10000000.0)

# hide color bar/color legend
solutionpvdDisplay.SetScalarBarVisibility(renderView1, False)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.1, 0.025, 0.39826142083018456]
renderView1.CameraFocalPoint = [0.1, 0.025, 0.0]
renderView1.CameraParallelScale = 0.07040341550470697

# save screenshot
SaveScreenshot('sandbox/restart_sandbox/Viscosity_1cm.png', magnification=1, quality=100, view=renderView1)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.1, 0.025, 0.39826142083018456]
renderView1.CameraFocalPoint = [0.1, 0.025, 0.0]
renderView1.CameraParallelScale = 0.07040341550470697

# save screenshot
SaveScreenshot('sandbox/restart_sandbox/Viscosity_2cm.png', magnification=1, quality=100, view=renderView1)

# set scalar coloring
ColorBy(solutionpvdDisplay, ('POINTS', 'density'))

# rescale color and/or opacity maps used to include current data range
solutionpvdDisplay.RescaleTransferFunctionToDataRange(True)

# show color bar/color legend
solutionpvdDisplay.SetScalarBarVisibility(renderView1, True)

# hide color bar/color legend
solutionpvdDisplay.SetScalarBarVisibility(renderView1, False)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.1, 0.025, 0.39826142083018456]
renderView1.CameraFocalPoint = [0.1, 0.025, 0.0]
renderView1.CameraParallelScale = 0.07040341550470697

# save screenshot
SaveScreenshot('sandbox/restart_sandbox/Compositions_2cm.png', magnification=1, quality=100, view=renderView1)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.1, 0.025, 0.39826142083018456]
renderView1.CameraFocalPoint = [0.1, 0.025, 0.0]
renderView1.CameraParallelScale = 0.07040341550470697

# save screenshot
SaveScreenshot('sandbox/restart_sandbox/Compositions_1cm.png', magnification=1, quality=100, view=renderView1)

# Properties modified on renderView1
renderView1.UseLight = 0

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.1, 0.025, 0.39826142083018456]
renderView1.CameraFocalPoint = [0.1, 0.025, 0.0]
renderView1.CameraParallelScale = 0.07040341550470697

# save screenshot
SaveScreenshot('sandbox/restart_sandbox/Compositions_1cm.png', magnification=1, quality=100, view=renderView1)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.1, 0.025, 0.39826142083018456]
renderView1.CameraFocalPoint = [0.1, 0.025, 0.0]
renderView1.CameraParallelScale = 0.07040341550470697

# save screenshot
SaveScreenshot('sandbox/restart_sandbox/Compositions_2cm.png', magnification=1, quality=100, view=renderView1)

# set active source
SetActiveSource(annotateTimeFilter2)

# show data in view
annotateTimeFilter2Display = Show(annotateTimeFilter2, renderView1)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.1, 0.025, 0.39826142083018456]
renderView1.CameraFocalPoint = [0.1, 0.025, 0.0]
renderView1.CameraParallelScale = 0.07040341550470697

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.1, 0.025, 0.39826142083018456]
renderView1.CameraFocalPoint = [0.1, 0.025, 0.0]
renderView1.CameraParallelScale = 0.07040341550470697

# save animation images/movie
WriteAnimation('sandbox/restart_sandbox/Composition.avi', Magnification=1, FrameRate=10.0, Compression=True)

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.1, 0.025, 0.39826142083018456]
renderView1.CameraFocalPoint = [0.1, 0.025, 0.0]
renderView1.CameraParallelScale = 0.07040341550470697

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
