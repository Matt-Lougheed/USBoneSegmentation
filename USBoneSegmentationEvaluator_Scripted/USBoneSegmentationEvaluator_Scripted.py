import os
import unittest
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging
import math

#
# USBoneSegmentationEvaluator_Scripted
#

class USBoneSegmentationEvaluator_Scripted(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "USBoneSegmentationEvaluator_Scripted" # TODO make this more human readable by adding spaces
    self.parent.categories = ["USBoneSegmentation"]
    self.parent.dependencies = []
    self.parent.contributors = ["Matt Lougheed (Queen's University)"] # replace with "Firstname Lastname (Organization)"
    self.parent.helpText = """
    This module allows for iterative calls to the USBoneSegmentationEvaluator CLI 
    module in order to compute metrics for multiple dilation values.
    """
    self.parent.acknowledgementText = """
    """

#
# USBoneSegmentationEvaluator_ScriptedWidget
#

class USBoneSegmentationEvaluator_ScriptedWidget(ScriptedLoadableModuleWidget):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setup(self):
    ScriptedLoadableModuleWidget.setup(self)

    # Instantiate and connect widgets ...

    #
    # Parameters Area
    #
    parametersCollapsibleButton = ctk.ctkCollapsibleButton()
    parametersCollapsibleButton.text = "Parameters"
    self.layout.addWidget(parametersCollapsibleButton)

    # Layout within the dummy collapsible button
    parametersFormLayout = qt.QFormLayout(parametersCollapsibleButton)

    #
    # Segmentation volume selector
    #
    self.segmentationSelector = slicer.qMRMLNodeComboBox()
    self.segmentationSelector.nodeTypes = ["vtkMRMLScalarVolumeNode"]
    self.segmentationSelector.selectNodeUponCreation = True
    self.segmentationSelector.addEnabled = False
    self.segmentationSelector.removeEnabled = False
    self.segmentationSelector.noneEnabled = False
    self.segmentationSelector.showHidden = False
    self.segmentationSelector.showChildNodeTypes = False
    self.segmentationSelector.setMRMLScene( slicer.mrmlScene )
    self.segmentationSelector.setToolTip("Choose the segmentation to dilate.")
    parametersFormLayout.addRow("Segmentation: ", self.segmentationSelector)

    #
    # False negative volume selector
    #
    self.falseNegativeSelector = slicer.qMRMLNodeComboBox()
    self.falseNegativeSelector.nodeTypes = ["vtkMRMLScalarVolumeNode"]
    self.falseNegativeSelector.selectNodeUponCreation = True
    self.falseNegativeSelector.addEnabled = False
    self.falseNegativeSelector.removeEnabled = False
    self.falseNegativeSelector.noneEnabled = False
    self.falseNegativeSelector.showHidden = False
    self.falseNegativeSelector.showChildNodeTypes = False
    self.falseNegativeSelector.setMRMLScene( slicer.mrmlScene )
    self.falseNegativeSelector.setToolTip("Choose the false negative test line.")
    parametersFormLayout.addRow("False negative line: ", self.falseNegativeSelector)

    #
    # True positive volume selector
    #
    self.truePositiveSelector = slicer.qMRMLNodeComboBox()
    self.truePositiveSelector.nodeTypes = ["vtkMRMLScalarVolumeNode"]
    self.truePositiveSelector.selectNodeUponCreation = True
    self.truePositiveSelector.addEnabled = False
    self.truePositiveSelector.removeEnabled = False
    self.truePositiveSelector.noneEnabled = False
    self.truePositiveSelector.showHidden = False
    self.truePositiveSelector.showChildNodeTypes = False
    self.truePositiveSelector.setMRMLScene( slicer.mrmlScene )
    self.truePositiveSelector.setToolTip("Choose the true positive test region.")
    parametersFormLayout.addRow("True positive region: ", self.truePositiveSelector)

    #
    # Dilation value
    #
    self.dilationSliderWidget = ctk.ctkSliderWidget()
    self.dilationSliderWidget.singleStep = 1
    self.dilationSliderWidget.minimum = 0
    self.dilationSliderWidget.maximum = 20
    self.dilationSliderWidget.value = 5
    self.dilationSliderWidget.setToolTip("Set threshold value for computing the output image. Voxels that have intensities lower than this value will set to zero.")
    parametersFormLayout.addRow("Dilation factor", self.dilationSliderWidget)

    #
    # Apply Button
    #
    self.applyButton = qt.QPushButton("Apply")
    self.applyButton.toolTip = "Run the algorithm."
    self.applyButton.enabled = False
    parametersFormLayout.addRow(self.applyButton)

    # connections
    self.applyButton.connect('clicked(bool)', self.onApplyButton)
    self.segmentationSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
    self.falseNegativeSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
    self.truePositiveSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)

    # Add vertical spacer
    self.layout.addStretch(1)

    # Refresh Apply button state
    self.onSelect()

  def cleanup(self):
    pass

  def onSelect(self):
    self.applyButton.enabled = self.segmentationSelector.currentNode() and self.falseNegativeSelector.currentNode() and self.truePositiveSelector.currentNode()

  def onApplyButton(self):
    logic = USBoneSegmentationEvaluator_ScriptedLogic()
    dilationFactor = self.dilationSliderWidget.value
    logic.run(self.segmentationSelector.currentNode(),
              self.falseNegativeSelector.currentNode(),
              self.truePositiveSelector.currentNode(),
              dilationFactor)

#
# USBoneSegmentationEvaluator_ScriptedLogic
#

class USBoneSegmentationEvaluator_ScriptedLogic(ScriptedLoadableModuleLogic):
  """This class should implement all the actual
  computation done by your module.  The interface
  should be such that other python code can import
  this class and make use of the functionality without
  requiring an instance of the Widget.
  Uses ScriptedLoadableModuleLogic base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def hasImageData(self,volumeNode):
    """This is an example logic method that
    returns true if the passed in volume
    node has valid image data
    """
    if not volumeNode:
      logging.debug('hasImageData failed: no volume node')
      return False
    if volumeNode.GetImageData() is None:
      logging.debug('hasImageData failed: no image data in volume node')
      return False
    return True

  def isValidInputOutputData(self, inputVolumeNode, outputVolumeNode):
    """Validates if the output is not the same as input
    """
    if not inputVolumeNode:
      logging.debug('isValidInputOutputData failed: no input volume node defined')
      return False
    if not outputVolumeNode:
      logging.debug('isValidInputOutputData failed: no output volume node defined')
      return False
    if inputVolumeNode.GetID()==outputVolumeNode.GetID():
      logging.debug('isValidInputOutputData failed: input and output volume is the same. Create a new volume for output to avoid this error.')
      return False
    return True

  def run(self, segmentationVolume, falseNegativeLine, truePositiveRegion, dilationFactor):
    """
    Run the actual algorithm
    """
    logging.info('Processing started')

    # Iteratively compute the false negative true positive metrics
    falseNegativeResults = slicer.mrmlScene.AddNode(slicer.vtkMRMLDoubleArrayNode())
    truePositiveResults = slicer.mrmlScene.AddNode(slicer.vtkMRMLDoubleArrayNode())
    falseNegativeArray = falseNegativeResults.GetArray()
    truePositiveArray = truePositiveResults.GetArray()
    falseNegativeArray.SetNumberOfTuples(int(dilationFactor))
    truePositiveArray.SetNumberOfTuples(int(dilationFactor))
    for dilateValue in range(int(dilationFactor)):
      cliParams = {'truePositiveGroundTruthVolume': truePositiveRegion.GetID(),
                   'falseNegativeGroundTruthVolume': falseNegativeLine.GetID(),
                   'dilateValue' : dilateValue,
                   'segmentedVolume': segmentationVolume.GetID()}
      cliNode = slicer.cli.run(slicer.modules.usbonesegmentationevaluator, None, cliParams, wait_for_completion=True)
      falseNegativeResult = cliNode.GetParameterAsString("falseNegativePercentage")
      falseNegativeArray.SetComponent(dilateValue, 0, dilateValue)
      falseNegativeArray.SetComponent(dilateValue, 1, float(falseNegativeResult))
      falseNegativeArray.SetComponent(dilateValue, 2, 0)
      truePositiveResult = cliNode.GetParameterAsString("truePositivePercentage")
      truePositiveArray.SetComponent(dilateValue, 0, dilateValue)
      truePositiveArray.SetComponent(dilateValue, 1, float(truePositiveResult))
      truePositiveArray.SetComponent(dilateValue, 2, 0)


    # Change to layout which has a chart view
    layoutNodes = slicer.mrmlScene.GetNodesByClass('vtkMRMLLayoutNode')
    layoutNodes.InitTraversal()
    layoutNode = layoutNodes.GetNextItemAsObject()
    layoutNode.SetViewArrangement(24)

    # Create chart node and chart view node
    chartViewNode = slicer.mrmlScene.AddNode(slicer.vtkMRMLChartViewNode())
    chartNode = slicer.mrmlScene.AddNode(slicer.vtkMRMLChartNode())
    chartNode.AddArray('False Negative Line Coverage', falseNegativeResults.GetID())
    chartNode.AddArray('True Positive Region Coverage', truePositiveResults.GetID())
    chartNode.SetProperty('default', 'xAxisLabel', 'Dilation Value')
    chartNode.SetProperty('default', 'yAxisLabel', 'Percentage')
    chartNode.SetProperty('default', 'type', 'Scatter')
    chartViewNode.SetChartNodeID(chartNode.GetID())

    logging.info('Processing completed')

    return True


class USBoneSegmentationEvaluator_ScriptedTest(ScriptedLoadableModuleTest):
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
    self.test_USBoneSegmentationEvaluator_Scripted1()

  def test_USBoneSegmentationEvaluator_Scripted1(self):
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
    logic = USBoneSegmentationEvaluator_ScriptedLogic()
    self.assertIsNotNone( logic.hasImageData(volumeNode) )
    self.delayDisplay('Test passed!')
