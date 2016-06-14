import os
import unittest
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging
import numpy, math

#
# USGeometry
#

class USGeometry(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "USGeometry"
    self.parent.categories = ["USBoneSegmentation"]
    self.parent.dependencies = []
    self.parent.contributors = ["Matt Lougheed (Queen's University)"]
    self.parent.helpText = """
    This is an example of scripted loadable module bundled in an extension.
    It performs a simple thresholding on the input volume and optionally captures a screenshot.
    """
    self.parent.acknowledgementText = """
"""

#
# USGeometryWidget
#

class USGeometryWidget(ScriptedLoadableModuleWidget):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setup(self):
    ScriptedLoadableModuleWidget.setup(self)

    # Instantiate and connect widgets ...

    #
    # INPUTS Area
    #
    inputsCollapsibleButton = ctk.ctkCollapsibleButton()
    inputsCollapsibleButton.text = "Inputs"
    self.layout.addWidget(inputsCollapsibleButton)
    inputsFormLayout = qt.QFormLayout(inputsCollapsibleButton)

    #
    # Input volume selector
    #
    self.inputSelector = slicer.qMRMLNodeComboBox()
    self.inputSelector.nodeTypes = ["vtkMRMLScalarVolumeNode"]
    self.inputSelector.selectNodeUponCreation = True
    self.inputSelector.addEnabled = False
    self.inputSelector.removeEnabled = False
    self.inputSelector.noneEnabled = False
    self.inputSelector.showHidden = False
    self.inputSelector.showChildNodeTypes = False
    self.inputSelector.setMRMLScene( slicer.mrmlScene )
    self.inputSelector.setToolTip( "Pick the input volume for geometry viewing." )
    inputsFormLayout.addRow("Input volume: ", self.inputSelector)

    #
    # Configuration file selector
    #
    fileLayout = qt.QHBoxLayout()
    self.configFile = qt.QLineEdit()
    self.configFile.setReadOnly(True)
    self.configFileButton = qt.QPushButton()
    self.configFileButton.setText("Select File")
    fileLayout.addWidget(self.configFile)
    fileLayout.addWidget(self.configFileButton)
    inputsFormLayout.addRow("Configuration file: ", fileLayout)

    #
    # Manual segmentations directory
    #
    directoryLayout = qt.QHBoxLayout()
    self.directory = qt.QLineEdit()
    self.directory.setReadOnly(True)
    self.directoryButton = qt.QPushButton()
    self.directoryButton.setText("Select Directory")
    directoryLayout.addWidget(self.directory)
    directoryLayout.addWidget(self.directoryButton)
    inputsFormLayout.addRow("Manual segmentations directory: ", directoryLayout)

    #
    # Input volume selector
    #
    self.algorithmSegmentationSelector = slicer.qMRMLNodeComboBox()
    self.algorithmSegmentationSelector.nodeTypes = ["vtkMRMLLabelMapVolumeNode"]
    self.algorithmSegmentationSelector.selectNodeUponCreation = True
    self.algorithmSegmentationSelector.addEnabled = False
    self.algorithmSegmentationSelector.removeEnabled = False
    self.algorithmSegmentationSelector.noneEnabled = False
    self.algorithmSegmentationSelector.showHidden = False
    self.algorithmSegmentationSelector.showChildNodeTypes = False
    self.algorithmSegmentationSelector.setMRMLScene( slicer.mrmlScene )
    self.algorithmSegmentationSelector.setToolTip( "Pick the algorithm segmentation to evaluate." )
    inputsFormLayout.addRow("Algorithm segmentation: ", self.algorithmSegmentationSelector)

    #
    # OUTPUT VOLUMES Area
    #
    outputsCollapsibleButton = ctk.ctkCollapsibleButton()
    outputsCollapsibleButton.text = "Output Volumes"
    self.layout.addWidget(outputsCollapsibleButton)
    outputsFormLayout = qt.QFormLayout(outputsCollapsibleButton)

    #
    # Merged manual segmentations
    #
    self.mergedManualSegmentations = slicer.qMRMLNodeComboBox()
    self.mergedManualSegmentations.nodeTypes = ["vtkMRMLLabelMapVolumeNode"]
    self.mergedManualSegmentations.selectNodeUponCreation = True
    self.mergedManualSegmentations.addEnabled = True
    self.mergedManualSegmentations.renameEnabled = True
    self.mergedManualSegmentations.removeEnabled = True
    self.mergedManualSegmentations.noneEnabled = False
    self.mergedManualSegmentations.showHidden = False
    self.mergedManualSegmentations.showChildNodeTypes = False
    self.mergedManualSegmentations.setMRMLScene( slicer.mrmlScene )
    self.mergedManualSegmentations.setToolTip( "Pick the label map for the merged manual segmentations." )
    outputsFormLayout.addRow("Merged manual segmentations: ", self.mergedManualSegmentations)

    #
    # Scanline label map
    #
    self.scanlines = slicer.qMRMLNodeComboBox()
    self.scanlines.nodeTypes = ["vtkMRMLLabelMapVolumeNode"]
    self.scanlines.selectNodeUponCreation = True
    self.scanlines.addEnabled = True
    self.scanlines.renameEnabled = True
    self.scanlines.removeEnabled = True
    self.scanlines.noneEnabled = False
    self.scanlines.showHidden = False
    self.scanlines.showChildNodeTypes = False
    self.scanlines.setMRMLScene( slicer.mrmlScene )
    self.scanlines.setToolTip( "Pick the label map for displaying the input volumes scanlines." )
    outputsFormLayout.addRow("Scanlines: ", self.scanlines)

    #
    # Scanline label map
    #
    self.outputSegmentation = slicer.qMRMLNodeComboBox()
    self.outputSegmentation.nodeTypes = ["vtkMRMLLabelMapVolumeNode"]
    self.outputSegmentation.selectNodeUponCreation = True
    self.outputSegmentation.addEnabled = True
    self.outputSegmentation.renameEnabled = True
    self.outputSegmentation.removeEnabled = True
    self.outputSegmentation.noneEnabled = False
    self.outputSegmentation.showHidden = False
    self.outputSegmentation.showChildNodeTypes = False
    self.outputSegmentation.setMRMLScene( slicer.mrmlScene )
    self.outputSegmentation.setToolTip( "Pick the label map for the output segmentation." )
    outputsFormLayout.addRow("Output segmentation: ", self.outputSegmentation)

    #
    # OUTPUT METRIC Area
    #
    outputMetricsCollapsibleButton = ctk.ctkCollapsibleButton()
    outputMetricsCollapsibleButton.text = "Output Metrics"
    self.layout.addWidget(outputMetricsCollapsibleButton)
    outputMetricsFormLayout = qt.QFormLayout(outputMetricsCollapsibleButton)

    #
    # False negative metric
    #
    self.falseNegativeMetric = qt.QLineEdit()
    self.falseNegativeMetric.setReadOnly(True)
    outputMetricsFormLayout.addRow("False negative percentage: ", self.falseNegativeMetric)

    #
    # False negative metric
    #
    self.truePositiveMetric = qt.QLineEdit()
    self.truePositiveMetric.setReadOnly(True)
    outputMetricsFormLayout.addRow("True positive percentage: ", self.truePositiveMetric)

    #
    # Apply Button
    #
    self.applyButton = qt.QPushButton("Apply")
    self.applyButton.toolTip = "Run the algorithm."
    self.applyButton.enabled = False
    outputMetricsFormLayout.addRow(self.applyButton)

    # Connections
    self.applyButton.connect('clicked(bool)', self.onApplyButton)
    self.inputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
    self.configFileButton.connect('clicked(bool)', self.selectFile)
    self.directoryButton.connect('clicked(bool)', self.selectDirectory)

    # Add vertical spacer
    self.layout.addStretch(1)

    # Refresh Apply button state
    self.onSelect()

  def cleanup(self):
    pass

  def onSelect(self):
    fileName = self.configFile.text
    self.applyButton.enabled = self.inputSelector.currentNode() and os.path.isfile(fileName)

  def selectFile(self):
    fileName = qt.QFileDialog().getOpenFileName()
    self.configFile.setText(fileName)
    self.applyButton.enabled = self.inputSelector.currentNode() and os.path.isfile(fileName)

  def selectDirectory(self):
    directoryName = qt.QFileDialog().getExistingDirectory()
    self.directory.setText(directoryName)

  def onApplyButton(self):
    logic = USGeometryLogic()
    logic.run(self.inputSelector.currentNode(), self.configFile.text, self.directory.text, self.mergedManualSegmentations.currentNode(), self.scanlines.currentNode(), self.outputSegmentation.currentNode())

#
# USGeometryLogic
#

class USGeometryLogic(ScriptedLoadableModuleLogic):
  """This class should implement all the actual
  computation done by your module.  The interface
  should be such that other python code can import
  this class and make use of the functionality without
  requiring an instance of the Widget.
  Uses ScriptedLoadableModuleLogic base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def run(self, inputVolume, configFileName, manualSegmentationsDirectory, mergedVolume, scanlineVolume, outputSegmentationVolume):
    """
    Run the actual algorithm
    """

    # Get the manual segmentations and create a single summed image
    import glob
    manualSegmentationFilenames = glob.glob(manualSegmentationsDirectory+"/*.mha")

    # Get the first image which each successive image will be added to
    reader = vtk.vtkMetaImageReader()
    reader.SetFileName(manualSegmentationFilenames[0])
    reader.Update()
    summedImage = vtk.vtkImageData()
    summedImage.SetExtent(reader.GetOutput().GetExtent())
    summedImage.AllocateScalars(vtk.VTK_UNSIGNED_CHAR,1)
    summedImage.ShallowCopy(reader.GetOutput())

    # Image data for output segmentation
    outputSegmentation = vtk.vtkImageData()
    outputSegmentation.SetExtent(reader.GetOutput().GetExtent())
    outputSegmentation.AllocateScalars(vtk.VTK_UNSIGNED_CHAR,1)


    # Initialize filter to add images together
    mathFilter = vtk.vtkImageMathematics()

    # Iterate list and add each new image
    for currentFile in manualSegmentationFilenames[1:]:
      # Get new image
      reader.SetFileName(currentFile)
      reader.Update()

      # Add it to existing summation
      mathFilter.SetInput1Data(summedImage)
      mathFilter.SetInput2Data(reader.GetOutput())
      mathFilter.Update()

      # Get new summation
      summedImage.ShallowCopy(mathFilter.GetOutput())

    # Add summed image to slicer scene
    rasToIjk = vtk.vtkMatrix4x4()
    ijkToRas = vtk.vtkMatrix4x4()
    inputVolume.GetRASToIJKMatrix(rasToIjk)
    inputVolume.GetIJKToRASMatrix(ijkToRas)
    mergedVolume.SetRASToIJKMatrix(rasToIjk)
    mergedVolume.SetIJKToRASMatrix(ijkToRas)
    mergedVolume.SetAndObserveImageData(summedImage)

    # Input ultrasound volume
    inputImageData = inputVolume.GetImageData()
    inputDimensions = inputImageData.GetDimensions()
    inputExtent = inputImageData.GetExtent()

    # Parse input config file
    US_Geometry = UltrasoundTransducerGeometry(configFileName)

    # Setup the canvas filter to draw scanlines
    imgDim = inputVolume.GetImageData().GetDimensions()
    drawFilter = vtk.vtkImageCanvasSource2D()
    drawFilter.SetExtent(0,imgDim[0]-1,0,imgDim[1]-1,0,0)
    drawFilter.SetDrawColor(1)

    # Points for the output segmentation
    outputSegmentationPoints = vtk.vtkPoints()

    # Iterate over scanlines and draw line
    for i in range(US_Geometry.numberOfScanlines):
      # Compute scanline endpoints
      [startScanline, endScanline] = US_Geometry.scanlineEndPoints(i)
      print("Scanline start: [{} {}".format(startScanline[0],startScanline[1]))
      print("Scanline end: [{} {}".format(endScanline[0],endScanline[1]))

      # Draw the scanline
      drawFilter.FillTube(int(startScanline[0]),int(startScanline[1]),int(endScanline[0]),int(endScanline[1]),1)

      # Compute the manual segmentation points on each scanline
      for z in range(imgDim[2]):
        print("Creating line for slice {}".format(z))
        currentLine = vtk.vtkLineSource()
        currentStartPoint = [startScanline[0], startScanline[1], z]
        currentEndPoint = [endScanline[0], endScanline[1], z]
        currentLine.SetPoint1(currentStartPoint)
        currentLine.SetPoint2(currentEndPoint)
        currentLine.SetResolution(US_Geometry.numberOfSamplesPerScanline)
        currentLine.Update()

        # Get the points on the line
        denom = math.sqrt((currentEndPoint[0] - currentStartPoint[0])**2 + (currentEndPoint[1] - currentStartPoint[1])**2)
        unitVector = [(currentEndPoint[0] - currentStartPoint[0])/denom, (currentEndPoint[1] - currentStartPoint[1])/denom]
        linePoints = currentLine.GetOutput()
        xVals = []
        yVals = []
        for j in range(linePoints.GetNumberOfPoints()):
          currentPoint = linePoints.GetPoint(j)
          summedSegPoint = summedImage.GetScalarComponentAsDouble(int(currentPoint[0]),int(currentPoint[1]),int(currentPoint[2]),0)
          if (summedSegPoint > 0):
            xVals.extend([int(currentPoint[0]) for _ in range(int(summedSegPoint))]) # Add points current X value X number of times where X is overlap count
            yVals.extend([int(currentPoint[1]) for _ in range(int(summedSegPoint))])
            print("The point [{}, {}, {}] has value {}!".format(int(currentPoint[0]), int(currentPoint[1]), int(currentPoint[2]), summedSegPoint))

        # If there are points on the line, compute mean / std
        if (len(xVals) > 0):
          # Compute mean point
          xMean = sum(xVals)/len(xVals)
          yMean = sum(yVals)/len(yVals)
          outputSegmentation.SetScalarComponentFromDouble(xMean,yMean,z,0,255)

          # Compute region edges from standard deviation
          xStd = numpy.std(xVals)
          yStd = numpy.std(yVals)
          stdSum = xStd + yStd
          toleranceFactor = 5 # What should be used as value here?
          acceptableRegionPoint1 = [int(xMean + toleranceFactor*stdSum * unitVector[0]), int(yMean + toleranceFactor*stdSum * unitVector[1])]
          acceptableRegionPoint2 = [int(xMean + -toleranceFactor*stdSum * unitVector[0]), int(yMean + -toleranceFactor*stdSum * unitVector[1])]          
          outputSegmentation.SetScalarComponentFromDouble(acceptableRegionPoint1[0],acceptableRegionPoint1[1],z,0,100)
          outputSegmentation.SetScalarComponentFromDouble(acceptableRegionPoint2[0],acceptableRegionPoint2[1],z,0,100)

          outputSegmentationPoints.InsertNextPoint(xMean,yMean,z)
      print("Finished scanline #{} processing".format(i))

    # Draw scanlines
    drawFilter.Update()

    # Copy scanline slice to match the Z-dimension of input US volume
    numOfSlices = imgDim[2]
    imageAppendFilter = vtk.vtkImageAppend()
    imageAppendFilter.SetAppendAxis(2)
    for _ in range(numOfSlices):
      imageAppendFilter.AddInputData(drawFilter.GetOutput())
    imageAppendFilter.Update()

    # Set scanline imagedata
    scanlineVolume.SetIJKToRASMatrix(ijkToRas)
    scanlineVolume.SetRASToIJKMatrix(rasToIjk)
    scanlineVolume.SetAndObserveImageData(imageAppendFilter.GetOutput())

    # Set output segmentation imagedata
    outputSegmentationVolume.SetIJKToRASMatrix(ijkToRas)
    outputSegmentationVolume.SetRASToIJKMatrix(rasToIjk)
    outputSegmentationVolume.SetAndObserveImageData(outputSegmentation)

    return True

class UltrasoundTransducerGeometry:
  def __init__(self, configFile):
    from xml.dom import minidom
    parser = minidom.parse(configFile)
    scanConversionElement = parser.getElementsByTagName("ScanConversion")
    if (len(scanConversionElement) < 1):
      slicer.util.errorDisplay("Could not find ScanConversion element in configuration file!")
      return False
    elif (len(scanConversionElement) > 1):
      slicer.util.errorDisplay("Found multiple ScanConversion elements in configuration file!")
      return False

    scanConversionElement = scanConversionElement[0]

    # Values common to both linear and curvilinear
    self.transducerGeometry = scanConversionElement.attributes['TransducerGeometry'].value
    self.transducerCenterPixel = scanConversionElement.attributes['TransducerCenterPixel'].value
    self.transducerCenterPixel = map(int, self.transducerCenterPixel.split(" "))
    self.numberOfScanlines = int(scanConversionElement.attributes['NumberOfScanLines'].value)
    self.outputImageSpacing = scanConversionElement.attributes['OutputImageSpacingMmPerPixel'].value
    self.outputImageSpacing = map(float, self.outputImageSpacing.split(" "))
    self.numberOfSamplesPerScanline = int(scanConversionElement.attributes['NumberOfSamplesPerScanLine'].value)

    # Values just for curvilinear
    if (self.transducerGeometry == "CURVILINEAR"):
      self.thetaStartDeg = float(scanConversionElement.attributes['ThetaStartDeg'].value)
      self.thetaStopDeg = float(scanConversionElement.attributes['ThetaStopDeg'].value)
      self.radiusStartMm = float(scanConversionElement.attributes['RadiusStartMm'].value)
      self.radiusStopMm = float(scanConversionElement.attributes['RadiusStopMm'].value)
      self.totalDeg = abs(self.thetaStopDeg - self.thetaStartDeg)
      self.degreesPerScanline = self.totalDeg / self.numberOfScanlines
      self.circleCenter = [self.transducerCenterPixel[0], self.transducerCenterPixel[1] - self.radiusStartMm/self.outputImageSpacing[1]]
    # Values just for linear
    elif (self.transducerGeometry == "LINEAR"):
      self.transducerWidthMm = float(scanConversionElement.attributes['TransducerWidthMm'].value)
      self.imagingDepthMm = float(scanConversionElement.attributes['ImagingDepthMm'].value)
      self.imageWidthPixel = self.transducerWidthMm / float(self.outputImageSpacing[0])
      print("Image width pixel: {}".format(self.imageWidthPixel))
      self.topLeftPixel = [self.transducerCenterPixel[0] - 0.5 * float(self.imageWidthPixel), self.transducerCenterPixel[1]]
      print("Top left pixel: {}".format(self.topLeftPixel))
      self.scanlineSpacingPixels = self.imageWidthPixel / float(self.numberOfScanlines)
      self.scanlineLengthPixels = self.imagingDepthMm / self.outputImageSpacing[1]
      print("Scanline length pixels: {}".format(self.scanlineLengthPixels))

  def scanlineEndPoints(self, scanline):
    import math
    # Compute curvilinear xy values
    if (self.transducerGeometry == "CURVILINEAR"):

      # Compute angle for desired scanline
      angle = self.thetaStartDeg + scanline * self.degreesPerScanline
      angleRadians = math.pi * angle / 180

      # Compute the starting point
      startScanlineX = self.circleCenter[0] + math.sin(angleRadians) * self.radiusStartMm / self.outputImageSpacing[0]
      startScanlineY = self.circleCenter[1] + math.cos(angleRadians) * self.radiusStartMm / self.outputImageSpacing[1]

      # Compute the ending point
      endScanlineX = self.circleCenter[0] + math.sin(angleRadians) * self.radiusStopMm / self.outputImageSpacing[0]
      endScanlineY = self.circleCenter[1] + math.cos(angleRadians) * self.radiusStopMm / self.outputImageSpacing[1]
    # Compute linear xy values
    elif (self.transducerGeometry == "LINEAR"):
      # Compute the starting point
      startScanlineX = self.topLeftPixel[0] + scanline * self.scanlineSpacingPixels
      startScanlineY = self.topLeftPixel[1]

      # Compute the end point
      endScanlineX = startScanlineX # Vertical line so same 'x' value
      endScanlineY = self.topLeftPixel[1] + self.scanlineLengthPixels
    else:
      # ERROR: should not reach here
      print("Error in scanlineEndPoints")

    # Combine XY values and return
    startScanline = [startScanlineX, startScanlineY]
    endScanline = [endScanlineX, endScanlineY]

    return [startScanline, endScanline]

class USGeometryTest(ScriptedLoadableModuleTest):
  """
  This is the test case for your scripted module.
  Uses ScriptedLoadableModuleTest base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setUp(self):
    """ Do whatever is needed to reset the state - typically a scene clear will be enough.
    """
    slicer.mrmlScene.Clear(0)

  def runTest(self):
    """Run as few or as many tests as needed here.
    """
    self.setUp()
    self.test_USGeometry1()

  def test_USGeometry1(self):
    """ Ideally you should have several levels of tests.  At the lowest level
    tests should exercise the functionality of the logic with different inputs
    (both valid and invalid).  At higher levels your tests should emulate the
    way the user would interact with your code and confirm that it still works
    the way you intended.
    One of the most important features of the tests is that it should alert other
    developers when their changes will have an impact on the behavior of your
    module.  For example, if a developer removes a feature that you depend on,
    your test should break so they know that the feature is needed.
    """

    self.delayDisplay("Starting the test")
    #
    # first, get some data
    #
    import urllib
    downloads = (
        ('http://slicer.kitware.com/midas3/download?items=5767', 'FA.nrrd', slicer.util.loadVolume),
        )

    for url,name,loader in downloads:
      filePath = slicer.app.temporaryPath + '/' + name
      if not os.path.exists(filePath) or os.stat(filePath).st_size == 0:
        logging.info('Requesting download %s from %s...\n' % (name, url))
        urllib.urlretrieve(url, filePath)
      if loader:
        logging.info('Loading %s...' % (name,))
        loader(filePath)
    self.delayDisplay('Finished with download and loading')

    volumeNode = slicer.util.getNode(pattern="FA")
    logic = USGeometryLogic()
    self.assertIsNotNone( logic.hasImageData(volumeNode) )
    self.delayDisplay('Test passed!')
