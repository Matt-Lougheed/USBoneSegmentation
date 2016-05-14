import os
import unittest
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging

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
    inputsFormLayout.addRow("Input Volume: ", self.inputSelector)

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
    # OUTPUTS Area
    #
    outputsCollapsibleButton = ctk.ctkCollapsibleButton()
    outputsCollapsibleButton.text = "Outputs"
    self.layout.addWidget(outputsCollapsibleButton)
    outputsFormLayout = qt.QFormLayout(outputsCollapsibleButton)

    #
    # Merged manual segmentations
    #
    self.mergedManualSegmentations = slicer.qMRMLNodeComboBox()
    self.mergedManualSegmentations.nodeTypes = ["vtkMRMLLabelMapVolumeNode"]
    self.mergedManualSegmentations.selectNodeUponCreation = True
    self.mergedManualSegmentations.addEnabled = True
    self.mergedManualSegmentations.removeEnabled = False
    self.mergedManualSegmentations.noneEnabled = False
    self.mergedManualSegmentations.showHidden = False
    self.mergedManualSegmentations.showChildNodeTypes = False
    self.mergedManualSegmentations.setMRMLScene( slicer.mrmlScene )
    self.mergedManualSegmentations.setToolTip( "Pick the label map for the merged manual segmentations." )
    outputsFormLayout.addRow("Merged Manual Segmentations: ", self.mergedManualSegmentations)

    #
    # Scanline label map
    #
    self.scanlines = slicer.qMRMLNodeComboBox()
    self.scanlines.nodeTypes = ["vtkMRMLLabelMapVolumeNode"]
    self.scanlines.selectNodeUponCreation = True
    self.scanlines.addEnabled = True
    self.scanlines.removeEnabled = False
    self.scanlines.noneEnabled = False
    self.scanlines.showHidden = False
    self.scanlines.showChildNodeTypes = False
    self.scanlines.setMRMLScene( slicer.mrmlScene )
    self.scanlines.setToolTip( "Pick the label map for displaying the input volumes scanlines." )
    outputsFormLayout.addRow("Scanlines: ", self.scanlines)

    #
    # Apply Button
    #
    self.applyButton = qt.QPushButton("Apply")
    self.applyButton.toolTip = "Run the algorithm."
    self.applyButton.enabled = False
    outputsFormLayout.addRow(self.applyButton)

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
    logic.run(self.inputSelector.currentNode(), self.configFile.text, self.directory.text, self.mergedManualSegmentations.currentNode(), self.scanlines.currentNode())

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


  def scanlineEndPoints(self, scanlineIndex):
    import math

    # Compute the angle of the scanline
    angle = thetaStartDeg + scanlineIndex*degreesPerScanline
    angleRadians = math.pi * angle / 180

    # Compute the starting point
    startScanlineX = circleCenter[0] + math.sin(angleRadians) * radiusStartMm / outputImageSpacing[0]
    startScanlineY = circleCenter[1] + math.cos(angleRadians) * radiusStartMm / outputImageSpacing[1]
    startScanline = [startScanlineX, startScanlineY]

    # Compute the ending point
    endScanLineX = circleCenter[0] + math.sin(angleRadians) * radiusStopMm / outputImageSpacing[0]
    endScanLineY = circleCenter[1] + math.cos(angleRadians) * radiusStopMm / outputImageSpacing[1]
    endScanline = [endScanLineX, endScanLineY]

    return [startScanline, endScanline]

  def run(self, inputVolume, configFileName, manualSegmentationsDirectory, mergedVolume, scanlineVolume):
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
    from xml.dom import minidom
    parser = minidom.parse(configFileName)
    scanConversionElement = parser.getElementsByTagName("ScanConversion")
    if (len(scanConversionElement) < 1):
      slicer.util.errorDisplay("Could not find ScanConversion element in configuration file!")
      return False
    elif (len(scanConversionElement) > 1):
      slicer.util.errorDisplay("Found multiple ScanConversion elements in configuration file!")
      return False

    scanConversionElement = scanConversionElement[0]

    # Get scan conversion attributes
    global thetaStartDeg, thetaStopDeg, transducerCenterPixel, radiusStartMm, radiusStopMm, outputImageSpacing, circleCenter, degreesPerScanline
    thetaStartDeg = float(scanConversionElement.attributes['ThetaStartDeg'].value)
    thetaStopDeg = float(scanConversionElement.attributes['ThetaStopDeg'].value)
    transducerCenterPixel = scanConversionElement.attributes['TransducerCenterPixel'].value
    transducerCenterPixel = map(int,transducerCenterPixel.split(" "))
    radiusStartMm = float(scanConversionElement.attributes['RadiusStartMm'].value)
    radiusStopMm = float(scanConversionElement.attributes['RadiusStopMm'].value)
    numberOfScanlines = int(scanConversionElement.attributes['NumberOfScanLines'].value)
    outputImageSpacing = scanConversionElement.attributes['OutputImageSpacingMmPerPixel'].value
    outputImageSpacing = map(float,outputImageSpacing.split(" "))

    # Find angle between scanlines
    totalDeg = abs(thetaStopDeg - thetaStartDeg)
    degreesPerScanline = totalDeg / numberOfScanlines

    # Circle center
    circleCenter = [transducerCenterPixel[0], transducerCenterPixel[1] - radiusStartMm/outputImageSpacing[1]]

    # Setup the canvas filter to draw scanlines
    imgDim = inputVolume.GetImageData().GetDimensions()
    drawFilter = vtk.vtkImageCanvasSource2D()
    drawFilter.SetExtent(0,imgDim[0]-1,0,imgDim[1]-1,0,0)
    drawFilter.SetDrawColor(1)

    # Iterate over scanlines and draw line
    for i in range(numberOfScanlines):
      # Compute scanline endpoints
      [startScanline, endScanline] = self.scanlineEndPoints(i)
      print("Scanline start: [{} {}".format(startScanline[0],startScanline[1]))
      print("Scanline end: [{} {}".format(endScanline[0],endScanline[1]))
      drawFilter.FillTube(int(startScanline[0]),int(startScanline[1]),int(endScanline[0]),int(endScanline[1]),1)
    drawFilter.Update()

    numOfSlices = imgDim[2]
    imageAppendFilter = vtk.vtkImageAppend()
    imageAppendFilter.SetAppendAxis(2)
    for _ in range(numOfSlices):
      imageAppendFilter.AddInputData(drawFilter.GetOutput())
    imageAppendFilter.Update()

    # Add scanline labelmap to slicer scene
    scanlineVolume.SetIJKToRASMatrix(ijkToRas)
    scanlineVolume.SetRASToIJKMatrix(rasToIjk)
    scanlineVolume.SetAndObserveImageData(imageAppendFilter.GetOutput())

    return True

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
