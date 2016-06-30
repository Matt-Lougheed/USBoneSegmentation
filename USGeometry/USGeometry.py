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
    # Algorithm segmentation volume selector
    #
    self.algorithmSegmentation = slicer.qMRMLNodeComboBox()
    self.algorithmSegmentation.nodeTypes = ["vtkMRMLLabelMapVolumeNode"]
    self.algorithmSegmentation.selectNodeUponCreation = True
    self.algorithmSegmentation.addEnabled = False
    self.algorithmSegmentation.removeEnabled = False
    self.algorithmSegmentation.noneEnabled = False
    self.algorithmSegmentation.showHidden = False
    self.algorithmSegmentation.showChildNodeTypes = False
    self.algorithmSegmentation.setMRMLScene( slicer.mrmlScene )
    self.algorithmSegmentation.setToolTip( "Pick the algorithm segmentation to evaluate." )
    inputsFormLayout.addRow("Algorithm segmentation: ", self.algorithmSegmentation)

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
    logic = USGeometryLogic(self.configFile.text, self.inputSelector.currentNode())
    logic.sumManualSegmentations(self.directory.text, self.mergedManualSegmentations.currentNode())
    logic.createScanlines(self.scanlines.currentNode())
    logic.computeMergedSegmentationMetrics(self.mergedManualSegmentations.currentNode(), self.outputSegmentation.currentNode(), self.algorithmSegmentation.currentNode(), self.truePositiveMetric)
    '''
    logic.run(self.inputSelector.currentNode(), self.configFile.text, self.directory.text, self.mergedManualSegmentations.currentNode(), self.scanlines.currentNode(), self.outputSegmentation.currentNode(), self.algorithmSegmentation.currentNode(), self.truePositiveMetric)
    '''

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

  def __init__(self, configFile, inputVolume):

    self.inputVolume = inputVolume
    self.rasToIjk = vtk.vtkMatrix4x4()
    self.ijkToRas = vtk.vtkMatrix4x4()
    self.inputVolume.GetRASToIJKMatrix(self.rasToIjk)
    self.inputVolume.GetIJKToRASMatrix(self.ijkToRas)

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
      self.topLeftPixel = [self.transducerCenterPixel[0] - 0.5 * float(self.imageWidthPixel), self.transducerCenterPixel[1]]
      self.scanlineSpacingPixels = self.imageWidthPixel / float(self.numberOfScanlines)
      self.scanlineLengthPixels = self.imagingDepthMm / self.outputImageSpacing[1]

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

  def euclidean_distance(self,point1,point2):
      return math.sqrt((point2[0] - point1[0]) ** 2 + (point2[1] - point1[1]) ** 2 + (point2[2] - point1[2]) ** 2)

  def sumManualSegmentations(self, manualSegmentationsDirectory, mergedVolume):
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
    mergedVolume.SetRASToIJKMatrix(self.rasToIjk)
    mergedVolume.SetIJKToRASMatrix(self.ijkToRas)
    mergedVolume.SetAndObserveImageData(summedImage)

  def createScanlines(self, scanlineVolume):
    imgDim = self.inputVolume.GetImageData().GetDimensions()
    drawFilter = vtk.vtkImageCanvasSource2D()
    drawFilter.SetExtent(0,imgDim[0]-1,0,imgDim[1]-1,0,0)
    drawFilter.SetDrawColor(1)

    for i in range(self.numberOfScanlines):
      [startScanline, endScanline] = self.scanlineEndPoints(i)
      drawFilter.FillTube(int(startScanline[0]),int(startScanline[1]),int(endScanline[0]),int(endScanline[1]),1)

    drawFilter.Update()

    # Copy scanline slice to match the Z-dimension of input US volume
    numOfSlices = imgDim[2]
    imageAppendFilter = vtk.vtkImageAppend()
    imageAppendFilter.SetAppendAxis(2)
    for _ in range(numOfSlices):
      imageAppendFilter.AddInputData(drawFilter.GetOutput())
    imageAppendFilter.Update()

    # Set scanline imagedata
    scanlineVolume.SetIJKToRASMatrix(self.ijkToRas)
    scanlineVolume.SetRASToIJKMatrix(self.rasToIjk)
    scanlineVolume.SetAndObserveImageData(imageAppendFilter.GetOutput())

  def computeMergedSegmentationMetrics(self, summedImage, outputSegmentation, algorithmSegmentation, truePositiveOutput):
    imgDim = self.inputVolume.GetImageData().GetDimensions()
    summedImageData = summedImage.GetImageData()
    algorithmSegmentationImageData = algorithmSegmentation.GetImageData()
    outputSegmentationImageData = vtk.vtkImageData()
    outputSegmentationImageData.SetExtent(summedImageData.GetExtent())
    outputSegmentationImageData.AllocateScalars(vtk.VTK_UNSIGNED_CHAR,1)

    # Values for metric computations
    totalAlgorithmSegmentationPoints = 0
    pointsWithinAcceptableRegion = 0

    # Iterate each scanline
    for i in range(self.numberOfScanlines):
      [startScanline, endScanline] = self.scanlineEndPoints(i)
      for z in range(imgDim[2]):
        # Create the scanline line source
        currentLine = vtk.vtkLineSource()
        currentStartPoint = [startScanline[0], startScanline[1], z]
        currentEndPoint = [endScanline[0], endScanline[1], z]
        currentLine.SetPoint1(currentStartPoint)
        currentLine.SetPoint2(currentEndPoint)
        currentLine.SetResolution(self.numberOfSamplesPerScanline)
        currentLine.Update()

        # Get the points on the line
        denom = math.sqrt((currentEndPoint[0] - currentStartPoint[0])**2 + (currentEndPoint[1] - currentStartPoint[1])**2)
        unitVector = [(currentEndPoint[0] - currentStartPoint[0])/denom, (currentEndPoint[1] - currentStartPoint[1])/denom]
        linePoints = currentLine.GetOutput()
        xVals = []
        yVals = []
        algorithmSegmentationPoints = []

        # Iterate each point on line
        for j in range(linePoints.GetNumberOfPoints()):
          currentPoint = linePoints.GetPoint(j)
          summedSegPoint = summedImageData.GetScalarComponentAsDouble(int(currentPoint[0]),int(currentPoint[1]),int(currentPoint[2]),0) # TODO: summedImage is mrmlNode now so need to get imagedata from it
          if (summedSegPoint > 0):
            # Add points current X/Y value N times where N is overlap count
            xVals.extend([int(currentPoint[0]) for _ in range(int(summedSegPoint))])
            yVals.extend([int(currentPoint[1]) for _ in range(int(summedSegPoint))])

          algorithmSegPoint = algorithmSegmentationImageData.GetScalarComponentAsDouble(int(currentPoint[0]),int(currentPoint[1]),int(currentPoint[2]),0)

          if (algorithmSegPoint > 0):
            algorithmSegmentationPoints.append([int(currentPoint[0]),int(currentPoint[1]),int(currentPoint[2])])
            totalAlgorithmSegmentationPoints += 1

        if (len(xVals) > 0):
          # Compute mean point
          xMean = sum(xVals)/len(xVals)
          yMean = sum(yVals)/len(yVals)
          outputSegmentationImageData.SetScalarComponentFromDouble(xMean,yMean,z,0,255)

          # Compute region edges from standard deviation
          euclideanDistances = []
          for currentX, currentY in zip(xVals, yVals):
            currentDistance = self.euclidean_distance([xMean,yMean,z],[currentX,currentY,z])
            euclideanDistances.append(currentDistance)
          std = numpy.std(euclideanDistances)

          toleranceFactor = 5 # What should be used here?
          acceptableRegionPoint1 = [int(xMean + toleranceFactor*std * unitVector[0]), int(yMean + toleranceFactor*std * unitVector[1])]
          acceptableRegionPoint2 = [int(xMean + -toleranceFactor*std * unitVector[0]), int(yMean + -toleranceFactor*std * unitVector[1])]

          outputSegmentationImageData.SetScalarComponentFromDouble(acceptableRegionPoint1[0],acceptableRegionPoint1[1],z,0,100)
          outputSegmentationImageData.SetScalarComponentFromDouble(acceptableRegionPoint2[0],acceptableRegionPoint2[1],z,0,100)

          # Compute distance to acceptable region edge, can use either point as they are equal just opposite
          acceptableDistance = self.euclidean_distance([xMean,yMean,z],[acceptableRegionPoint1[0],acceptableRegionPoint1[1],z])
          for currentPoint in algorithmSegmentationPoints:
            currentDistance = self.euclidean_distance([xMean,yMean,z], currentPoint)
            if (currentDistance < acceptableDistance):
              pointsWithinAcceptableRegion += 1

    outputSegmentation.SetAndObserveImageData(outputSegmentationImageData)
    outputSegmentation.SetRASToIJKMatrix(self.rasToIjk)
    outputSegmentation.SetIJKToRASMatrix(self.ijkToRas)
    truePositiveValue = float(pointsWithinAcceptableRegion) / float(totalAlgorithmSegmentationPoints) * 100
    truePositiveOutput.setText(truePositiveValue)

  def run(self, inputVolume, configFileName, manualSegmentationsDirectory, mergedVolume, scanlineVolume, outputSegmentationVolume, algorithmSegmentationVolume, truePositiveOutput):
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

    # Algorithm segmentation
    algorithmSegmentationData = algorithmSegmentationVolume.GetImageData()
    totalAlgorithmSegmentationPoints = 0
    pointsWithinAcceptableRegion = 0

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

      # Draw the scanline
      drawFilter.FillTube(int(startScanline[0]),int(startScanline[1]),int(endScanline[0]),int(endScanline[1]),1)

      # Compute the manual segmentation points on each scanline
      for z in range(imgDim[2]):
        # Create line
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
        algorithmSegmentationPoints = []
        for j in range(linePoints.GetNumberOfPoints()):
          # TODO: Get the point from the algorithm segmentation so we can compute if it is within the required region / false negative line
          currentPoint = linePoints.GetPoint(j)
          summedSegPoint = summedImage.GetScalarComponentAsDouble(int(currentPoint[0]),int(currentPoint[1]),int(currentPoint[2]),0)
          if (summedSegPoint > 0):
            xVals.extend([int(currentPoint[0]) for _ in range(int(summedSegPoint))]) # Add points current X value X number of times where X is overlap count
            yVals.extend([int(currentPoint[1]) for _ in range(int(summedSegPoint))])

          algorithmSegPoint = algorithmSegmentationData.GetScalarComponentAsDouble(int(currentPoint[0]),int(currentPoint[1]),int(currentPoint[2]),0)
          if (algorithmSegPoint > 0):
            algorithmSegmentationPoints.append([int(currentPoint[0]),int(currentPoint[1]),int(currentPoint[2])])
            totalAlgorithmSegmentationPoints += 1

        # If there are points on the line, compute mean / std
        if (len(xVals) > 0):
          # Compute mean point
          xMean = sum(xVals)/len(xVals)
          yMean = sum(yVals)/len(yVals)
          outputSegmentation.SetScalarComponentFromDouble(xMean,yMean,z,0,255)

          # Compute region edges from standard deviation
          euclideanDistances = []
          for currentX, currentY in zip(xVals, yVals):
            currentDistance = self.euclidean_distance([xMean,yMean,z],[currentX,currentY,z])
            euclideanDistances.append(currentDistance)
          std = numpy.std(euclideanDistances)

          toleranceFactor = 5 # What should be used as value here?
          acceptableRegionPoint1 = [int(xMean + toleranceFactor*std * unitVector[0]), int(yMean + toleranceFactor*std * unitVector[1])]
          acceptableRegionPoint2 = [int(xMean + -toleranceFactor*std * unitVector[0]), int(yMean + -toleranceFactor*std * unitVector[1])]

          outputSegmentation.SetScalarComponentFromDouble(acceptableRegionPoint1[0],acceptableRegionPoint1[1],z,0,100)
          outputSegmentation.SetScalarComponentFromDouble(acceptableRegionPoint2[0],acceptableRegionPoint2[1],z,0,100)
          outputSegmentationPoints.InsertNextPoint(xMean,yMean,z)

          acceptableDistance = self.euclidean_distance([xMean,yMean,z],[acceptableRegionPoint1[0],acceptableRegionPoint1[1],z])
          for currentPoint in algorithmSegmentationPoints:
            currentDistance = self.euclidean_distance([xMean,yMean,z],currentPoint)
            if (currentDistance < acceptableDistance):
              pointsWithinAcceptableRegion += 1

      print("Finished scanline #{} processing".format(i))

    # Compute metrics and set output text
    truePositive = float(pointsWithinAcceptableRegion) / float(totalAlgorithmSegmentationPoints) * 100
    truePositiveOutput.setText(truePositive)
    
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
  def __init__(self, configFile, inputVolume):

    self.inputVolume = inputVolume

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
      self.topLeftPixel = [self.transducerCenterPixel[0] - 0.5 * float(self.imageWidthPixel), self.transducerCenterPixel[1]]
      self.scanlineSpacingPixels = self.imageWidthPixel / float(self.numberOfScanlines)
      self.scanlineLengthPixels = self.imagingDepthMm / self.outputImageSpacing[1]

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
