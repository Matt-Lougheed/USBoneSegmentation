#include "USBoneSegmentationEvaluatorCLP.h"
#include "itkPluginUtilities.h"
#include "itkFlatStructuringElement.h"
#include "itkImage.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMaskImageFilter.h"
#include "itkImageToHistogramFilter.h"

#include <time.h>
#include <vector>
#include <cmath>

// Use an anonymous namespace to keep class types and function names
// from colliding when module is used as shared object module.  Every
// thing should be in an anonymous namespace except for the module
// entry point, e.g. main()
//
namespace
{

template <class T>
int DoIt( int argc, char * argv[], T )
{
  PARSE_ARGS;

  typedef itk::Image<T, 3> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typedef itk::ImageFileWriter<ImageType> WriterType;
  typedef itk::ImageRegionIteratorWithIndex<ImageType> IteratorType;

  typename ReaderType::Pointer segmentationReader = ReaderType::New();
  typename ReaderType::Pointer truePositiveReader = ReaderType::New();
  typename ReaderType::Pointer dilatedSegmentationReader = ReaderType::New();
  typename ReaderType::Pointer falseNegativeReader = ReaderType::New();
  typename WriterType::Pointer writer = WriterType::New();

  // Timing
  time_t start,end;
  time(&start);

  // Read in segmented meta image
  segmentationReader->SetFileName(segmentedVolume.c_str());
  try
  {
      segmentationReader->Update();
  }
  catch (itk::ExceptionObject & err)
  {
      std::cerr << "Error reading segmentedVolume!" << std::endl;
      std::cerr << err << std::endl;
  }
  typename ImageType::Pointer segmentedImage = segmentationReader->GetOutput();

  // Read in true positive verification meta image
  truePositiveReader->SetFileName(truePositiveGroundTruthVolume.c_str());
  try
  {
      truePositiveReader->Update();
  }
  catch (itk::ExceptionObject & err)
  {
      std::cerr << "Error reading truePositiveGroundTruthVolume!" << std::endl;
      std::cerr << err << std::endl;
      return EXIT_FAILURE;
  }
  typename ImageType::Pointer truePositiveImage = truePositiveReader->GetOutput();

  // Read in false negative verification meta image
  falseNegativeReader->SetFileName(falseNegativeGroundTruthVolume.c_str());
  try
  {
      falseNegativeReader->Update();
  }
  catch (itk::ExceptionObject & err)
  {
      std::cerr << "Error reading falseNegativeGroundTruthVolume!" << std::endl;
      std::cerr << err << std::endl;
  }
  typename ImageType::Pointer falseNegativeImage = falseNegativeReader->GetOutput();

  // Check that each volume has the same dimensions
  typename ImageType::RegionType segmentedRegion = segmentedImage->GetLargestPossibleRegion();
  typename ImageType::RegionType truePositiveRegion = truePositiveImage->GetLargestPossibleRegion();
  typename ImageType::RegionType falseNegativeRegion = falseNegativeImage->GetLargestPossibleRegion();
  typename ImageType::SizeType segmentedSize = segmentedRegion.GetSize();
  typename ImageType::SizeType truePositiveSize = truePositiveRegion.GetSize();
  typename ImageType::SizeType falseNegativeSize = falseNegativeRegion.GetSize();

  if (!(segmentedSize[0] == truePositiveSize[0] && segmentedSize[0] == falseNegativeSize[0] &&
      segmentedSize[1] == truePositiveSize[1] && segmentedSize[1] == falseNegativeSize[1] &&
      segmentedSize[2] == truePositiveSize[2] && segmentedSize[2] == falseNegativeSize[2]))
  {
    return EXIT_FAILURE;
  }

  // Iterator for segmented meta image dilation
  dilatedSegmentationReader->SetFileName(segmentedVolume.c_str());
  dilatedSegmentationReader->Update();

  // Use cross structure for dilation
  typedef itk::FlatStructuringElement<3> StructElementType;

  // Set cross size to dilation value specified
  StructElementType::SizeType crossSize;
  crossSize[0] = dilateValue;
  crossSize[1] = dilateValue;
  crossSize[2] = 0;
  StructElementType structuringElement = StructElementType::Cross(crossSize);

  // Dilate the algorithms segmentation
  typedef itk::BinaryDilateImageFilter<ImageType, ImageType, StructElementType>
    BinaryDilateImageFilterType;
  typename BinaryDilateImageFilterType::Pointer dilateFilter =
    BinaryDilateImageFilterType::New();
  dilateFilter->SetInput(dilatedSegmentationReader->GetOutput());
  dilateFilter->SetKernel(structuringElement);
  dilateFilter->SetDilateValue(1);
  dilateFilter->Update();
  typename ImageType::Pointer dilatedImage = ImageType::New();
  dilatedImage = dilateFilter->GetOutput();

  // Write the dilated segmentation to output file
  writer->SetInput(dilatedImage);
  writer->SetFileName(outputDilatedSegmentation.c_str());
  try
  {
  writer->Update();
  }
  catch (itk::ExceptionObject & e)
  {
      std::cerr << "Error while writing dilated segmentation" << std::endl;
      std::cerr << e << std::endl;
      return EXIT_FAILURE;
  }

  // Mask false negative line with dilated segmentation line
  typedef itk::MaskImageFilter<ImageType, ImageType> MaskFilterType;
  typename MaskFilterType::Pointer falseNegativeMaskFilter = MaskFilterType::New();
  falseNegativeMaskFilter->SetInput(falseNegativeReader->GetOutput());
  falseNegativeMaskFilter->SetMaskImage(dilateFilter->GetOutput());
  falseNegativeMaskFilter->Update();
  typename ImageType::Pointer falseNegativeMaskImage = falseNegativeMaskFilter->GetOutput();
  falseNegativeMaskImage->DisconnectPipeline();

  // Mask segmentation with true positive region
  typename MaskFilterType::Pointer truePositiveMaskFilter = MaskFilterType::New();
  truePositiveMaskFilter->SetInput(segmentationReader->GetOutput());
  truePositiveMaskFilter->SetMaskImage(truePositiveReader->GetOutput());
  truePositiveMaskFilter->Update();
  typename ImageType::Pointer truePositiveMaskImage = truePositiveMaskFilter->GetOutput();

  // Image list to iterate for computing histogram / statistics
  std::list<typename ImageType::Pointer> imageList;
  imageList.push_back(segmentedImage);
  imageList.push_back(falseNegativeImage);
  imageList.push_back(dilatedImage);
  imageList.push_back(falseNegativeMaskImage);
  imageList.push_back(truePositiveMaskImage);

  // List to hold pixel counts
  std::vector<int> pixelCounts;

  // Setup histogram filter to calculate pixel intensity counts
  typedef itk::Statistics::ImageToHistogramFilter<ImageType>
      ImageToHistogramFilterType;
  const unsigned int MeasurementVectorSize = 1;
  const unsigned int numberOfBins = 256;
  typename ImageToHistogramFilterType::HistogramType::MeasurementVectorType
      lowerBound(numberOfBins);
  lowerBound.Fill(0);
  typename ImageToHistogramFilterType::HistogramType::MeasurementVectorType
      upperBound(numberOfBins);
  upperBound.Fill(255);
  typename ImageToHistogramFilterType::HistogramType::SizeType
      size(MeasurementVectorSize);
  size.Fill(numberOfBins);
  typename ImageToHistogramFilterType::Pointer imageToHistogramFilter =
      ImageToHistogramFilterType::New();
  imageToHistogramFilter->SetHistogramBinMinimum(lowerBound);
  imageToHistogramFilter->SetHistogramBinMaximum(upperBound);
  imageToHistogramFilter->SetHistogramSize(size);

  // Iterate over each image and compute pixel count
  typename ImageToHistogramFilterType::HistogramType* histogram;
  typename ImageType::Pointer currentImage;
  for (typename std::list<typename ImageType::Pointer>::iterator it = imageList.begin();
       it != imageList.end();
       it++)
  {
      currentImage = *it;
      imageToHistogramFilter->SetInput(currentImage);
      imageToHistogramFilter->Update();
      histogram = imageToHistogramFilter->GetOutput();
      pixelCounts.push_back(histogram->GetFrequency(255));
  }

  // Note: these values are dependent on the order of images in previous for loop
  int totalSegmentationVoxelCount = pixelCounts[0];
  int totalFalseNegativeLineCount = pixelCounts[1];
  int totalDilatedVoxelCount = pixelCounts[2];
  int falseNegativeOverlapCount = pixelCounts[3];
  int truePositiveOverlapCount = pixelCounts[4];

    // Calculate metrics for total volume
  falseNegativePercentage =
      (1 - (float)falseNegativeOverlapCount/(float)totalFalseNegativeLineCount) * 100;
  truePositivePercentage =
      truePositiveOverlapCount/(float)totalSegmentationVoxelCount * 100;

  // Write module output values
  std::ofstream rts;
  rts.open(returnParameterFile.c_str());
  rts << "falseNegativePercentage = " << falseNegativePercentage << std::endl;
  rts << "truePositivePercentage = " << truePositivePercentage << std::endl;
  rts.close();

  time(&end);
  double dif = difftime(end,start);
  std::cout << "Time to run: " << dif << std::endl;

  return EXIT_SUCCESS;
}

} // end of anonymous namespace

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  itk::ImageIOBase::IOPixelType     pixelType;
  itk::ImageIOBase::IOComponentType componentType;

  try
    {
    itk::GetImageType(segmentedVolume, pixelType, componentType);

    // This filter handles all types on input, but only produces
    // signed types
    switch( componentType )
      {
      case itk::ImageIOBase::UCHAR:
        return DoIt( argc, argv, static_cast<unsigned char>(0) );
        break;
      case itk::ImageIOBase::CHAR:
        return DoIt( argc, argv, static_cast<char>(0) );
        break;
      case itk::ImageIOBase::USHORT:
        return DoIt( argc, argv, static_cast<unsigned short>(0) );
        break;
      case itk::ImageIOBase::SHORT:
        return DoIt( argc, argv, static_cast<short>(0) );
        break;
      case itk::ImageIOBase::UINT:
        return DoIt( argc, argv, static_cast<unsigned int>(0) );
        break;
      case itk::ImageIOBase::INT:
        return DoIt( argc, argv, static_cast<int>(0) );
        break;
      case itk::ImageIOBase::ULONG:
        return DoIt( argc, argv, static_cast<unsigned long>(0) );
        break;
      case itk::ImageIOBase::LONG:
        return DoIt( argc, argv, static_cast<long>(0) );
        break;
      case itk::ImageIOBase::FLOAT:
        return DoIt( argc, argv, static_cast<float>(0) );
        break;
      case itk::ImageIOBase::DOUBLE:
        return DoIt( argc, argv, static_cast<double>(0) );
        break;
      case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
      default:
        std::cout << "unknown component type" << std::endl;
        break;
      }
    }

  catch( itk::ExceptionObject & excep )
    {
    std::cerr << argv[0] << ": exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
