/*======================================================================
  Attenuation Correction for the ClearPEM scanner
  Suitable both for STIR and List Mode reconstruction algorithms.
  
  Inputs: Initial reconstructed image (can be reconstructed with a big voxel and filter)
          .lm file of the same data

  Outputs: LMF where the weight field has an attenuation correction factor

  Author: C. S. Ferreira

  Modification History: Dec 2011 - Jan 2012
                        May 2012 - Background segmentation and inversion
			Sep 2012 - Background segmentation, fill holes, inversion, fill holes.

  ======================================================================*/

#include "itkImage.h"
#include <iostream>
#include <string>
#include "RawImport.h"
#include <stdlib.h> 

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkExtractImageFilter.h"
#include "itkShiftScaleImageFilter.h"
#include "itkGDCMImageIO.h"
#include "itkMetaDataObject.h"
#include "itkNumericTraits.h"

#include "itkRescaleIntensityImageFilter.h"
#include "itkCastImageFilter.h"

#include <string>
#include <itksys/ios/sstream>
#include <itksys/Base64.h>

#include <vector>
#include <itksys/SystemTools.hxx>
//#include "DicomCLP.h"
#include <math.h>
#include <sys/stat.h>

#include <time.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <assert.h>

#include <itkSize.h>

//Smooth filter
//#include "itkVectorGradientAnisotropicDiffusionImageFilter.h"
#include "itkGradientAnisotropicDiffusionImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"

//Equalization
#include "itkAdaptiveHistogramEqualizationImageFilter.h"

//Flip Image
#include "itkFlipImageFilter.h"

//segmentation
#include "itkConnectedThresholdImageFilter.h"
//#include "itkImage.h"
//#include "itkCastImageFilter.h"
#include <itkMinimumMaximumImageCalculator.h>
#include "itkOtsuThresholdImageFilter.h"
#include "itkOtsuMultipleThresholdsCalculator.h"
#include "itkScalarImageToHistogramGenerator.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkNumericTraits.h"


//Erode
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
// #include "QuickView.h"


//Fill Holes
#include "itkBinaryMorphologicalClosingImageFilter.h"
//#include "itkBinaryBallStructuringElement.h"

#include "itkVotingBinaryIterativeHoleFillingImageFilter.h"

//Invert Intensity Image
#include <itkInvertIntensityImageFilter.h>


//Smooth isosurface before extraction
#include <itkAntiAliasBinaryImageFilter.h>

//Padding
#include "itkConstantPadImageFilter.h"

//MeshExtract
#include "itkBinaryMask3DMeshSource.h"
#include "itkImage.h"
#include "itkMesh.h"

//ITK->VTK
#include "itkMeshTovtkPolyData.h"
#include "itkVTKPolyDataWriter.h"


//Mesh Fill Holes
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkFillHolesFilter.h>

//Mesh Information
#include <vtkMassProperties.h>

//Decimation
#include <vtkDecimatePro.h>
#include <vtkQuadricClustering.h>

//Smooth mesh
#include <vtkWindowedSincPolyDataFilter.h>

//Intersection line->mesh
#include <vtkOBBTree.h>
#include <vtkModifiedBSPTree.h>

//Visualization
// #include "itkImageToVTKImageFilter.h"
#include "vtkImageViewer.h"
#include "vtkRenderWindowInteractor.h"
#include <vtkPolyDataReader.h>


//VTK
#include <vtkSmartPointer.h>
#include <vtkFeatureEdges.h>
#include <vtkPolyDataWriter.h>

//CloseMeshHoles
#include <vtkCleanPolyData.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkAppendPolyData.h>
#include <vtkStripper.h>
#include <vtkTriangleFilter.h>


#include <vtkDataSetSurfaceFilter.h>

#include <vtkClipPolyData.h>
#include <vtkPlane.h>
#include <vtkActor.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
 


#include <time.h>

void ReadInterfile(std::istream& FileHeader, std::string OutVol, std::string path);
//void SmoothGradientAnisotropicDiffusionImageFilter(std::string fileIn, std::string smoothedVol, unsigned int numberOfIterations, double timeStep, double conductance);
//void MeshExtraction(std::string fileIn,std::string MeshVol,int objectValue);
//vtkPolyData * CloseMeshHoles(vtkSmartPointer<vtkPolyData> meshin);

int main(int argc, char *argv[])
{
  
  if( argc < 2 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << "  inputMaskImageFile.hv  auto|low-auto|thresholdValue" << std::endl;
    return EXIT_FAILURE;
    }

   //TIME STAMP
  //clock_t start, end; //, startDecimation,endDecimation, startSmooth,endSmooth, startEvIt, endEvIt; ///clock_t gives wrong results.
  time_t time1,time2;
 
  std::string HeaderFileName = argv[1];
  std::string tempVol = std::string(tempnam(NULL, "itodc")) + ".nii";
  //std::string tempVol = "teste.nii";

  std::ifstream FileHeaderPEM(HeaderFileName.c_str());

  int found;      
  std::string s1,path;

  found=HeaderFileName.find_last_of("/"); //looks for the last / in the full path of header file
  path.assign(HeaderFileName.c_str(), found+1); // assign to path the full path of header file up to the last /
  ReadInterfile(FileHeaderPEM,tempVol,path); //execute the file reading for PEM file, code in Loader.h
  FileHeaderPEM.close();

  /* //read .nrrd file
  typedef itk::Vector<signed short,3>       VectorType;
  typedef itk::Image<VectorType,3>	    DiffusionImageType;
  typedef DiffusionImageType::Pointer	    DiffusionImagePointer;
   
   
  typedef itk::ImageFileReader<DiffusionImageType> FileReaderType;
  typedef itk::ImageToVTKImageFilter<DiffusionImageType>FilterType;
  

  FileReaderType::Pointer reader = FileReaderType::New();
  reader->SetFileName(tempVol);
  //reader->Update();
  */


  //Smooth
  //SmoothGradientAnisotropicDiffusionImageFilter(tempVol,smoothVol,5,0.05,1);  //0.125,1);

  const     unsigned int    Dimension = 3;
  typedef    float    InputPixelType;
  typedef    float    OutputPixelType;
  
  typedef itk::Image< OutputPixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType >  ReaderType;
  
  typedef   float           InternalPixelType;
  typedef itk::Image< InternalPixelType, Dimension >  InternalImageType;
  
  typedef itk::FlipImageFilter< ImageType >  FlipFilterType;
  FlipFilterType::Pointer flipFilter = FlipFilterType::New();


  typedef itk::GradientAnisotropicDiffusionImageFilter< ImageType, ImageType > FilterType;
  FilterType::Pointer smooth = FilterType::New();

  typedef itk::ConnectedThresholdImageFilter< ImageType, ImageType > ConnectedFilterType;
  ConnectedFilterType::Pointer connectedThreshold = ConnectedFilterType::New();

  typedef itk::MinimumMaximumImageCalculator<ImageType> ImageCalculatorFilterType;
  ImageCalculatorFilterType::Pointer maxMinFilter = ImageCalculatorFilterType::New();

  typedef itk::OtsuThresholdImageFilter<ImageType,ImageType >  OtsuFilterType;
  OtsuFilterType::Pointer segment = OtsuFilterType::New();

  
  typedef itk::Statistics::ScalarImageToHistogramGenerator<ImageType > ScalarImageToHistogramGeneratorType;
  typedef ScalarImageToHistogramGeneratorType::HistogramType    HistogramType;
  ScalarImageToHistogramGeneratorType::Pointer scalarImageToHistogramGenerator = ScalarImageToHistogramGeneratorType::New();
  typedef itk::OtsuMultipleThresholdsCalculator< HistogramType >   CalculatorType;
  CalculatorType::Pointer calculator = CalculatorType::New();
  typedef itk::BinaryThresholdImageFilter<ImageType, ImageType >  BinaryFilterType;
  BinaryFilterType::Pointer binaryFilter = BinaryFilterType::New();
  
  //Erode
  //typedef itk::BinaryBallStructuringElement<ImageType::PixelType, 3> StructuringElementType;
  //StructuringElementType structuringElement;
  //typedef itk::BinaryErodeImageFilter <ImageType, ImageType, StructuringElementType> BinaryErodeImageFilterType;
  //BinaryErodeImageFilterType::Pointer erodeFilter = BinaryErodeImageFilterType::New();


  //holes filling
  /*typedef itk::BinaryBallStructuringElement< InternalPixelType, Dimension > StructuringElementType;
  StructuringElementType structuringElement;
  typedef itk::BinaryMorphologicalClosingImageFilter <ImageType, ImageType, StructuringElementType>
          BinaryMorphologicalClosingImageFilterType;
  BinaryMorphologicalClosingImageFilterType::Pointer dilateFilter
          = BinaryMorphologicalClosingImageFilterType::New();
  */
  typedef itk::VotingBinaryIterativeHoleFillingImageFilter<ImageType >  FillHolesFilterType;
  FillHolesFilterType::Pointer dilateFilter = FillHolesFilterType::New();
 
  //Invert Intensity Image
  typedef itk::InvertIntensityImageFilter <ImageType>    InvertIntensityImageFilterType; 
  InvertIntensityImageFilterType::Pointer invertIntensityFilter = InvertIntensityImageFilterType::New();

  //Smooth isosurface before extraction (verify if it makes a good difference)
  typedef itk::AntiAliasBinaryImageFilter <ImageType, ImageType> AntiAliasFilterType;
  AntiAliasFilterType::Pointer antiAliasFilter = AntiAliasFilterType::New ();

  //Padding
  typedef itk::ConstantPadImageFilter <ImageType, ImageType> ConstantPadImageFilterType;
  ConstantPadImageFilterType::Pointer padFilter = ConstantPadImageFilterType::New();

  //Mesh
  //typedef itk::DefaultDynamicMeshTraits<double, 3, 3,double,double> TriangleMeshTraits;
  //typedef itk::Mesh<double,3, TriangleMeshTraits> MeshType;
  typedef itkMeshTovtkPolyData::TriangleMeshType MeshType;
  //typedef itk::Mesh<double>          MeshType;
  typedef itk::BinaryMask3DMeshSource< ImageType, MeshType >   MeshSourceType;
  MeshSourceType::Pointer meshSource = MeshSourceType::New();

  //Writer
  typedef float WritePixelType;
  typedef itk::Image< WritePixelType, Dimension > WriteImageType;
  typedef itk::ImageFileWriter< WriteImageType >  WriterType;


  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( tempVol );
  reader->Update();
  // flipFilter->SetInput(reader->GetOutput() );
  smooth->SetInput( reader->GetOutput() );
  //connectedThreshold->SetInput( smooth->GetOutput() );
  //dilateFilter->SetInput(connectedThreshold->GetOutput());
  //antiAliasFilter->SetInput(dilateFilter->GetOutput());
  //meshSource->SetInput( antiAliasFilter->GetOutput() );
  //padFilter->SetInput(dilateFilter->GetOutput());
  //meshSource->SetInput( imgVol );
  //connector ->SetInput( meshSource -> GetOutput());


  //--------GET IMAGE SIZE---------------//

  ImageType::RegionType inputSize = reader->GetOutput()->GetLargestPossibleRegion();
  ImageType::SizeType size = inputSize.GetSize();
  // get size of the whole 3D image
  float imageSizeX = size[0];
  float imageSizeY = size[1];
  float imageSizeZ = size[2];
  std::cerr << "Input size (voxels): " << size << std::endl;
  //std::cerr<< "sizes:  " <<imageSizeX<<"  "<<imageSizeY<<"  "<<imageSizeZ<<std::endl;
  std::cerr << "Number of Voxels in image" << imageSizeX*imageSizeY*imageSizeZ << std::endl;

  ImageType::SpacingType imageSpacing = reader->GetOutput()->GetSpacing();
  std::cerr << "Voxel Size: " << imageSpacing << std::endl;
  //std::cerr<< "spacing:  " <<imageSpacing[0]<<"  "<<imageSpacing[1]<<"  "<<imageSpacing[2]<<std::endl;


  float axialDim = imageSizeZ*imageSpacing[2];  //mm ---->For check purposes
  float xDim = imageSizeX*imageSpacing[0];  //mm  ----> For check purposes
  float yDim = imageSizeY*imageSpacing[1];  //mm  ----> For check purposes
  std::cerr<< "Dimensions (mm):  " <<xDim<<"  "<<yDim<<"  "<<axialDim<<std::endl;


  //-------------------FLIP AXES----------------//
  /* typedef FlipFilterType::FlipAxesArrayType     FlipAxesArrayType;
  FlipAxesArrayType flipArray;

  flipArray[0] = 1;
  flipArray[1] = 1;
  flipArray[2] = 1;

  flipFilter->SetFlipAxes( flipArray );
  flipFilter->Update();
  */
  
  //------------------- SMOOTH -------------------//
  //const unsigned int numberOfIterations = atoi( 5 );
  //const double       timeStep = atof( 0.05 );
  //const double       conductance = atof( 1 );

  int numberOfIterations = 10;
  double timeStep = 0.05;
  double conductance = 1;

  smooth->SetNumberOfIterations( numberOfIterations );
  smooth->SetTimeStep( timeStep );
  smooth->SetConductanceParameter( conductance );
  smooth->Update();



  //-------------NEW!!! EQUALIZATION-------------//
  typedef  itk::AdaptiveHistogramEqualizationImageFilter< ImageType > AdaptiveHistogramEqualizationImageFilterType;
  AdaptiveHistogramEqualizationImageFilterType::Pointer adaptiveHistogramEqualizationImageFilter = AdaptiveHistogramEqualizationImageFilterType::New();
  // ImageType::SizeType eqRadius; eqRadius[0]=5;eqRadius[1]=5;eqRadius[2]=5;
  //RadiusValueType eqRadius1;

  /*adaptiveHistogramEqualizationImageFilter->SetInput(smooth->GetOutput());
  adaptiveHistogramEqualizationImageFilter->SetBeta(1);
  adaptiveHistogramEqualizationImageFilter->SetAlpha(0.9);
  AdaptiveHistogramEqualizationImageFilterType::ImageSizeType radius;
  radius.Fill(1);
  adaptiveHistogramEqualizationImageFilter->SetRadius(radius);
  adaptiveHistogramEqualizationImageFilter->Update();
  */

  //----------------ZEROING BORDERS-------------------//
  ImageType::Pointer imgVol = ImageType::New();
  imgVol = smooth->GetOutput();
  //imgVol = adaptiveHistogramEqualizationImageFilter->GetOutput();
  InternalImageType::IndexType  voxelIndex;
  
  for (int lineVoxel = 0; lineVoxel < imageSizeX; lineVoxel++)
    {
      for (int columnVoxel = 0; columnVoxel < imageSizeY; columnVoxel++)
	{
	  voxelIndex[0]=lineVoxel;
	  voxelIndex[1]=columnVoxel;
	  voxelIndex[2]=imageSizeZ-1;
	  imgVol->SetPixel(voxelIndex,0);
	  voxelIndex[2]=0;
	  imgVol->SetPixel(voxelIndex,0);	  
	}
      for (int axialVoxel = 0; axialVoxel < imageSizeZ; axialVoxel++)
	{
	  voxelIndex[0]=lineVoxel;
	  voxelIndex[1]=0;
	  voxelIndex[2]=axialVoxel;
	  imgVol->SetPixel(voxelIndex,0);
	  voxelIndex[1]=imageSizeY-1;
	  imgVol->SetPixel(voxelIndex,0);
	}
    }
  
  for (int columnVoxel = 0; columnVoxel < imageSizeY; columnVoxel++)
    {
      for (int axialVoxel = 0; axialVoxel < imageSizeZ; axialVoxel++)
	{
	  voxelIndex[0]=imageSizeX-1;
	  voxelIndex[1]=columnVoxel;
	  voxelIndex[2]=axialVoxel;
	  imgVol->SetPixel(voxelIndex,0);
	  voxelIndex[0]=0;
	  imgVol->SetPixel(voxelIndex,0);
	}
    }
  
  WriterType::Pointer zeroingBordersWriter = WriterType::New();
  zeroingBordersWriter->SetFileName("teste_zeroingBordersWriter.nii" );
  zeroingBordersWriter->SetInput(imgVol); 
  zeroingBordersWriter->Update(); 
  
  /*
  //--------------------PADDING----------------// 
  ImageType::SizeType lowerExtendRegion;
  lowerExtendRegion[0] = 1;
  lowerExtendRegion[1] = 1;
  lowerExtendRegion[2] = 1;
 
  ImageType::SizeType upperExtendRegion;
  upperExtendRegion[0] = 1;
  upperExtendRegion[1] = 1;
  upperExtendRegion[2] = 1;
 
  ImageType::PixelType constantPixel = 0;
 
  //padFilter->SetInput(image);
  //padFilter->SetPadBound(outputRegion); // Calls SetPadLowerBound(region) and SetPadUpperBound(region)
  padFilter->SetInput(smooth->GetOutput());
  padFilter->SetPadLowerBound(lowerExtendRegion);
  padFilter->SetPadUpperBound(upperExtendRegion);
  padFilter->SetConstant(constantPixel);
  std::cerr << "antes pad update" << std::endl;
  padFilter->Update();
  std::cerr << "depois pad update" << std::endl;

  WriterType::Pointer paddingWriter = WriterType::New();
  paddingWriter->SetFileName("teste_paddingWriter.nii" );
  paddingWriter->SetInput(padFilter->GetOutput()); 
  paddingWriter->Update(); 
  */

  //---------------SEGMENTATION------------------//
  //Min and max
  maxMinFilter->SetImage(smooth->GetOutput());
  maxMinFilter->Compute();

  float maximumImage = maxMinFilter->GetMaximum();
  float minimumImage = maxMinFilter->GetMinimum();

  std::cerr<<"Image Maximum counts : "<<maximumImage<<std::endl;
  std::cerr<<"Image Minimum counts : "<<minimumImage<<std::endl;
  
  ImageType::IndexType maximumIndex = maxMinFilter->GetIndexOfMaximum();
  ImageType::IndexType minimumIndex = maxMinFilter->GetIndexOfMinimum();

  std::cerr<<"Image Maximum position : "<<maximumIndex<<std::endl;
  std::cerr<<"Image Minimum position : "<<minimumIndex<<std::endl;
  
  
  //Otsu Segmentation
  segment->SetInput(imgVol);
  //const OutputPixelType otsuOutsideValue = 1.0;
  //const OutputPixelType otsuInsideValue  = 0.0;
  //segment->SetOutsideValue( otsuOutsideValue );
  //segment->SetInsideValue(  otsuInsideValue  );
  //segment->SetNumberOfHistogramBins(128);
  segment->Update();
  int otsuThreshold = segment->GetThreshold();
  std::cout << "Threshold = " << otsuThreshold << std::endl;
  
  //dilateFilter->SetInput(segment->GetOutput());
  


  //Otsu Multiple threshold
  float otsuMultiThreshold = 0;
  if ((strcmp(argv[2], "low-auto") == 0)) //&& (atoi(argv[2]) == 0))
    {  
      std::cerr<<"Multiple Otsu"<<std::endl;
      scalarImageToHistogramGenerator->SetNumberOfBins(180 );
      calculator->SetNumberOfThresholds( 4 );
      
      // const OutputPixelType outsideValue = 0;
      //const OutputPixelType insideValue = 1;
      
      //binaryFilter->SetOutsideValue( outsideValue );
      //binaryFilter->SetInsideValue(  insideValue  );
      
      scalarImageToHistogramGenerator->SetInput( imgVol );
      calculator->SetInputHistogram( scalarImageToHistogramGenerator->GetOutput() );
      binaryFilter->SetInput( imgVol );
      
      scalarImageToHistogramGenerator->Compute();
      
      try
	{calculator->Update();}
      catch( itk::ExceptionObject & excp )
	{std::cerr << "Exception thrown calculator" << excp << std::endl;}
            
      const CalculatorType::OutputType &thresholdVector = calculator->GetOutput();
      CalculatorType::OutputType::const_iterator itNum = thresholdVector.begin();
      
      int num=0;
      //float otsuLowerThreshold;
      for(; itNum < thresholdVector.end(); itNum++)
	{
	  num++;
	  std::cout << "OtsuThreshold["
		    << (int)(itNum - thresholdVector.begin())
		    << "] = " <<
	    static_cast<itk::NumericTraits<CalculatorType::MeasurementType>::PrintType>
	    (*itNum) << std::endl;
	  if (num == 1)
	    {otsuMultiThreshold = static_cast<itk::NumericTraits<CalculatorType::MeasurementType>::PrintType>(*itNum);} 
	}
      
      //uncoment next lines in case you want to use this segmentation process alone
      /*
      binaryFilter->SetLowerThreshold( otsuLowerThreshold );
      std::cerr<< "otsuLowerThreshold   " <<otsuLowerThreshold<< std::endl;
      binaryFilter->Update();
      */
      //dilateFilter->SetInput(binaryFilter->GetOutput());
      //erodeFilter->SetInput(binaryFilter->GetOutput());
    }
  
   //else
   //{

  //-Connected Threshold
  float thresholdValue;
  if ((strcmp(argv[2], "auto")) == 0)
    {thresholdValue = otsuThreshold;}
  else if ((strcmp(argv[2], "low-auto") == 0)) //&& (atoi(argv[2]) == 0)) 
    {thresholdValue = otsuMultiThreshold;}
  else //(atoi(argv[2]) != 0){ 
    {thresholdValue = atoi( argv[2] );
      if ((atoi( argv[2])) == 0) {std::cerr << "Error on threshold definition!! Exiting!!\n"; exit(0);}
      std::cerr<<"Manual "<<std::endl;}
 
  std::cerr<<"Connected"<<std::endl;
  float lowerThreshold = minimumImage;   //atoi( argv[5] );     //(cilindrov2x2) 2589; //(clinico1)900;  //100.0;  //0.0033*maximumImage
  std::cerr << "lowerThreshold: " << lowerThreshold << std::endl;
  float upperThreshold = thresholdValue;   //30000.0;
  std::cerr << "upperThreshold: " << upperThreshold << std::endl;
  
  connectedThreshold->SetLower(  lowerThreshold  );
  connectedThreshold->SetUpper(  upperThreshold  );
  connectedThreshold->SetReplaceValue( 1 );
  connectedThreshold->SetSeed(minimumIndex); //maximumIndex);  //index );
  connectedThreshold->SetInput(imgVol);
  //connectedThreshold->SetInput(padFilter->GetOutput());
  
  dilateFilter->SetInput(connectedThreshold->GetOutput());
  //}
  
  /*//----------ERODE-----------------//
    unsigned int radius = 5;
    
    //typedef itk::Image<unsigned char, 3>    ImageType;
    //typedef itk::ImageFileReader<ImageType> ReaderType;
    //ReaderType::Pointer reader = ReaderType::New();
    //reader->SetFileName(argv[1]);
    
    structuringElement.SetRadius(radius);
    structuringElement.CreateStructuringElement();
    
    //erodeFilter->SetInput(binaryFilter->GetOutput());
    erodeFilter->SetKernel(structuringElement);
    erodeFilter->Update();
  */
  /*QuickView viewer;
    viewer.AddImage(smooth->GetOutput());
    viewer.AddImage(erodeFilter->GetOutput());
    viewer.Visualize();
  */
  // dilateFilter->SetInput(erodeFilter->GetOutput());
  
  
  
  //------------------FILL HOLES DILATE FILTER------------------//
  //Fill Holes 1
  //itkBinaryMorphologicalClosingImageFilter option
  /* unsigned int radius = 1;
     structuringElement.SetRadius(radius);
     structuringElement.CreateStructuringElement();
     dilateFilter->SetKernel(structuringElement);
     dilateFilter->SetForegroundValue(1);
     dilateFilter->Update();
     //std::cerr<< "dilateFilter foreground value: "<<(dilateFilter->GetForegroundValue())<<std::endl;
     */
  
  //Fill Holes2
  //itkVotingBinaryIterativeHoleFillingImageFilter
  //dilateFilter->SetInput(segment->GetOutput());
  const unsigned int radiusX = 1;
  const unsigned int radiusY = 1;
  const unsigned int radiusZ = 1;
  
  ImageType::SizeType indexRadius;
  indexRadius[0] = radiusX; // radius along x
  indexRadius[1] = radiusY; // radius along y
  indexRadius[2] = radiusZ; // radius along z
  dilateFilter->SetRadius( indexRadius );
  
  dilateFilter->SetBackgroundValue( 0 );
  dilateFilter->SetForegroundValue( 1 ); //Set new foreground voxels to this value. Foreground here is the background of the image.
  
  dilateFilter->SetMajorityThreshold( 1 );
  
  //const unsigned int numberOfIterations = atoi( argv[5] );
  //DilateFilter->SetMaximumNumberOfIterations( numberOfIterations );
  
  try
    {
      dilateFilter->Update();
      //const unsigned int iterationsUsed = dilateFilter->GetCurrentNumberOfIterations();
      //std::cerr << "The filter used " << iterationsUsed << " iterations " << std::endl;
      //const unsigned int numberOfPixelsChanged = dilateFilter->GetNumberOfPixelsChanged();
      //std::cerr << "and changed a total of " << numberOfPixelsChanged << " pixels" << std::endl;
    }
  catch(itk::ExceptionObject & excp )
    {
      std::cerr << "Problem updating dilate filter" << std::endl;
      std::cerr << excp << std::endl;
      return EXIT_FAILURE;
    }
  
  /*QuickView viewer;
    viewer.AddImage(smooth->GetOutput());
    viewer.AddImage(dilateFilter->GetOutput());
    viewer.Visualize();
  */
  
  //Invert Intensity Image
  invertIntensityFilter->SetInput(dilateFilter->GetOutput());
  invertIntensityFilter->SetMaximum(1);
  invertIntensityFilter->Update();

  //FILL HOLES DILATE FILTER2 
  dilateFilter->SetInput(invertIntensityFilter->GetOutput());
  dilateFilter->Update();


  
  //-------------STATISTICS----------------//
  int nObjectVoxels = 0;
  ImageType::PixelType pixValue;
  ImageType::IndexType pixIndex; 
  for (int xdir = 0; xdir < imageSizeX ; xdir++)
    {
      for (int ydir = 0; ydir < imageSizeY ; ydir++)
	{
	  for (int zdir = 0; zdir < imageSizeZ; zdir++)
	    {
	      pixIndex[0] = xdir;
	      pixIndex[1] = ydir;
	      pixIndex[2] = zdir;
	      pixValue =  dilateFilter->GetOutput()->GetPixel(pixIndex);
	      if (pixValue > 0){nObjectVoxels++;}
	    }
	}
    }
  std::cerr << "Number of voxels belonging to the object: " <<  nObjectVoxels << std::endl;
  std::cerr << "Percentage of voxels belonging to the object: " <<  (nObjectVoxels/(imageSizeX*imageSizeY*imageSizeZ))*100 << std::endl;
  
  //----------------OUTPUT---------------------//  
  WriterType::Pointer smoothWriter = WriterType::New();
  smoothWriter->SetFileName( "teste_smooth.nii" );
  smoothWriter->SetInput(smooth->GetOutput());
  smoothWriter->Update();
  
  /*WriterType::Pointer segmentWriter = WriterType::New();
  segmentWriter->SetFileName("teste_segmented.nii" );
  //segmentWriter->SetInput(binaryFilter->GetOutput() );
  segmentWriter->SetInput(connectedThreshold->GetOutput() );
  segmentWriter->Update();
  */ 

  WriterType::Pointer fillHolesWriter = WriterType::New();
  fillHolesWriter->SetFileName("teste_segmented_holesFilled.nii" );
  // fillHolesWriter->SetInput(invertIntensityFilter->GetOutput() );
  fillHolesWriter->SetInput(dilateFilter->GetOutput() );
  fillHolesWriter->Update(); 

  
  //----------------VTK viewer-----------------//
  /* typedef itk::ImageToVTKImageFilter<ImageType>ConnectFilterType;
  
  vtkImageViewer * viewer = vtkImageViewer::New();
  vtkRenderWindowInteractor * renderWindowInteractor = vtkRenderWindowInteractor::New();

  ConnectFilterType::Pointer connectorItkVtk = ConnectFilterType::New();
  connectorItkVtk->SetInput(dilateFilter->GetOutput());

  viewer->SetupInteractor(renderWindowInteractor);
  viewer->SetInput( connectorItkVtk->GetOutput());
  viewer->Render();
  viewer->SetColorWindow(225);
  viewer->SetColorLevel(128);
  renderWindowInteractor->Start();
  */
  




  /*
  //Smooth isosurface before extraction (verify if it makes a good difference)
  double maximumRMSError = 0.01;     //0.01;
  unsigned int numberOfIterationsSmoothIsosurface = 10; //50;

  antiAliasFilter->SetInput( dilateFilter->GetOutput() );
  antiAliasFilter->SetMaximumRMSError( maximumRMSError );
  antiAliasFilter->SetNumberOfIterations( numberOfIterationsSmoothIsosurface );
  antiAliasFilter->SetNumberOfLayers( 2 );
  antiAliasFilter->Update();
  */


  //-----------MESH EXTRACTION-------------------//
  //MeshExtraction(smoothVol,MeshVol,1);
  
  struct stat st; // Unused variable, just to call stat
  
  if (!stat("betta_override.nii", &st)) {
    std::cout << ">>>>>>>>> Anche MeshCreator usa betta_override.nii !!! " << std::endl;
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( "betta_override.nii" );
    reader->Update();
    meshSource->SetInput( reader->GetOutput() );
  } else {
    meshSource->SetInput( dilateFilter->GetOutput() );  //dilateFilter->GetOutput() );
  }
  
  float objectValue = 1;
  meshSource->SetObjectValue( objectValue );
  try
    {
      meshSource->Update();
      std::cerr << "Nodes = " << meshSource->GetNumberOfNodes() << std::endl;
      std::cerr << "Cells = " << meshSource->GetNumberOfCells() << std::endl;
    }
  catch(itk::ExceptionObject & excp )
    {
      std::cerr << "Problem updating mesh" << std::endl;
      std::cerr << excp << std::endl;
      return EXIT_FAILURE;
    }

  
  //------------ITK->VTK : write and read file option---------------//
  typedef itk::VTKPolyDataWriter<MeshType>VtkWriterType;
  VtkWriterType::Pointer meshWriter = VtkWriterType::New();
  meshWriter->SetFileName("tmesh.vtk");
  meshWriter->SetInput( meshSource->GetOutput() );
  meshWriter->Update();


   vtkSmartPointer<vtkPolyDataReader> p_reader = vtkSmartPointer<vtkPolyDataReader>::New();
  p_reader->SetFileName("tmesh.vtk");
  p_reader->Update();
  

  //------------------FILL MESH HOLES-----------------//
  //vtkSmartPointer<vtkPolyData> meshInput = vtkSmartPointer<vtkPolyData>::New();
  //meshInput->ShallowCopy(p_reader->GetOutput());
  
  vtkSmartPointer<vtkFillHolesFilter> fillMeshHolesFilter = vtkSmartPointer<vtkFillHolesFilter>::New();
  //fillMeshHolesFilter->SetInputConnection(meshInput->GetProducerPort());
  fillMeshHolesFilter->SetInput(p_reader->GetOutput());
  fillMeshHolesFilter->SetHoleSize(80.0);
  fillMeshHolesFilter->Update();

  vtkSmartPointer<vtkCleanPolyData> cleanMesh =  vtkSmartPointer<vtkCleanPolyData>::New();
  //cleanMesh->SetInput(fillMeshHolesFilter->GetOutput());
  cleanMesh->SetInput(p_reader->GetOutput());
  cleanMesh->Update();


  /*
  vtkSmartPointer<vtkFeatureEdges> featureEdges = vtkSmartPointer<vtkFeatureEdges>::New();
  featureEdges->SetInput(p_reader->GetOutput());
  
  vtkSmartPointer<vtkStripper> stripper = vtkSmartPointer<vtkStripper>::New();
  stripper->SetInput(featureEdges->GetOutput());
  
  vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
  polyData->SetPoints(stripper->GetOutput()->GetPoints());
  polyData->SetPolys(stripper->GetOutput()->GetLines());
  
  vtkSmartPointer<vtkTriangleFilter> triangleFilter = vtkSmartPointer<vtkTriangleFilter>::New();
  triangleFilter->SetInput(polyData);
    
  vtkSmartPointer<vtkAppendPolyData> append = vtkSmartPointer<vtkAppendPolyData>::New();
  append->AddInput(p_reader->GetOutput());
  append->AddInput(triangleFilter->GetOutput());
  append->Update();
  */
  
  
  //Check output mesh
  /*  vtkSmartPointer<vtkPolyDataWriter> p_writer = vtkSmartPointer<vtkPolyDataWriter>::New();
  p_writer->SetFileName("teste_mesh_filled.vtk");
  p_writer->SetInput(fillMeshHolesFilter->GetOutput());
  p_writer->Update();
  */ 

 //-------------CHECK IF SURFACE IS CLOSED----------------//
  vtkSmartPointer<vtkFeatureEdges> featureEdges2 = vtkSmartPointer<vtkFeatureEdges>::New();
  featureEdges2->FeatureEdgesOff();
  featureEdges2->BoundaryEdgesOn();
  featureEdges2->NonManifoldEdgesOn();
  //featureEdges->SetInputConnection(ClosedMeshHoles->GetOutputPort());
  //featureEdges->SetInput(ClosedMeshHoles);
  featureEdges2->SetInputConnection(cleanMesh ->GetOutputPort());
  featureEdges2->Update();
 
  int numberOfOpenEdges = featureEdges2->GetOutput()->GetNumberOfCells();
 
  if(numberOfOpenEdges > 0)
    {
    std::cerr << "Surface is not closed" << std::endl;
    }
  else
    {
    std::cerr << "Surface is closed" << std::endl;
    }
  
   
  //------------MESH INFORMATION---------------------//
  vtkSmartPointer<vtkMassProperties> meshInfo = vtkSmartPointer<vtkMassProperties>::New();
  meshInfo->SetInput(cleanMesh->GetOutput());
  meshInfo->Update();
  std::cerr << "Mesh Volume: " << meshInfo->GetVolume() << " mm3" << std::endl;
  std::cerr << "Mesh Volume Projected :" << meshInfo->GetVolumeProjected() << " mm3" <<std::endl;
  std::cerr << "check (greater than Mesh Volume? this should identify a problem):" << ((meshInfo->GetVolume()-meshInfo->GetVolumeProjected())*10000) <<std::endl;
  std::cerr << "Mesh Surface Area: " << meshInfo->GetSurfaceArea() << " mm3" << std::endl;


 //---------------------SMOOTH------------------------------//
  //startSmooth = clock();
  vtkSmartPointer<vtkWindowedSincPolyDataFilter> meshSmoother = vtkSmartPointer<vtkWindowedSincPolyDataFilter>::New();
  // meshSmoother->SetInputConnection(decimate->GetOutputPort());
  meshSmoother->SetInput(cleanMesh->GetOutput());
  meshSmoother->SetNumberOfIterations(20);
  // meshSmoother->BoundarySmoothingOff(); //nao se ve diferença entre o on e off
  meshSmoother->FeatureEdgeSmoothingOff(); //O Off retira picos (se usado depois do decimate é melhor o on)
  // meshSmoother->SetEdgeAngle(90); //Nao vi diferenca
  meshSmoother->SetFeatureAngle(120.0);   //quanto maior, menor o efeito lego
  meshSmoother->SetPassBand(0.01);    //1 e 0.1 da efeito lego. 0.01 works. 0.001 da para 20it e da lego para menos. 
  meshSmoother->NonManifoldSmoothingOn();
  meshSmoother->NormalizeCoordinatesOn();
  meshSmoother->Update();

  //Notas: O smooth aumenta com: +it, < passband, > FeatureAngle.
  
  /*  //TIME STAMP--------//
  endSmooth = clock();
  cout << "Time required for execution Smooth: "
       << (double)(endSmooth-startSmooth)/CLOCKS_PER_SEC
       << " seconds." << "\n\n";
  //------------------//
  
  //Check output smoothed mesh
  vtkSmartPointer<vtkPolyDataWriter> p_writerSmoothedMesh = vtkSmartPointer<vtkPolyDataWriter>::New();
  p_writerSmoothedMesh->SetFileName("teste_mesh_smoothed.vtk");
  p_writerSmoothedMesh->SetInput(meshSmoother->GetOutput());
  p_writerSmoothedMesh->Update();
  */


  //---------------DECIMATION-----------------//
  vtkSmartPointer<vtkPolyData> beforeDecimation = vtkSmartPointer<vtkPolyData>::New();
  beforeDecimation->ShallowCopy(meshSmoother->GetOutput());
 
  std::cout << "Before decimation" << std::endl << "------------" << std::endl;
  std::cout << "There are " << beforeDecimation->GetNumberOfPoints() << " points." << std::endl;
  std::cout << "There are " << beforeDecimation->GetNumberOfPolys() << " polygons." << std::endl;
  
  //startDecimation = clock();
  vtkSmartPointer<vtkDecimatePro> decimate = vtkSmartPointer<vtkDecimatePro>::New();
  //#if VTK_MAJOR_VERSION <= 5
  //  decimate->SetInputConnection(input->GetProducerPort());
  //#else
  //decimate->SetInputData();
  decimate->SetInput(meshSmoother->GetOutput());
  //#endif
  decimate->SetTargetReduction(.80); //80% reduction (if there was 100 triangles, now there will be 90)
  //decimate->PreserveTopologyOff();
  //decimate->SplitEdgesOff();
  decimate->Update();
 

  /* //Quadric Clustering
  vtkQuadricClustering * decimate = vtkQuadricClustering::New();
  decimate->SetNumberOfXDivisions(32);
  decimate->SetNumberOfYDivisions(32);
  decimate->SetNumberOfZDivisions(32);
  decimate->SetInput(fillMeshHolesFilter->GetOutput());
  decimate->Update();
  */
  
  /* //TIME STAMP--------//
  endDecimation = clock();
  cout << "Time required for execution Decimation: "
       << (double)(endDecimation-startDecimation)/CLOCKS_PER_SEC
       << " seconds." << "\n\n";
  //------------------//
  */

  vtkSmartPointer<vtkPolyData> decimated = vtkSmartPointer<vtkPolyData>::New();
  decimated->ShallowCopy(decimate->GetOutput());
 
  std::cout << "After decimation" << std::endl << "------------" << std::endl;
 
  std::cout << "There are " << decimated->GetNumberOfPoints() << " points." << std::endl;
  std::cout << "There are " << decimated->GetNumberOfPolys() << " polygons." << std::endl;

  //Check output decimated mesh
  vtkSmartPointer<vtkPolyDataWriter> p_writerDecimated = vtkSmartPointer<vtkPolyDataWriter>::New();
  // p_writerDecimated->SetFileName("teste_mesh_decimated.vtk");
  p_writerDecimated->SetFileName("tfinal_mesh.vtk");
  p_writerDecimated->SetInput(decimate->GetOutput());
  p_writerDecimated->Update();
  


 //TIME STAMP--------//
  time(&time2);
  double diff_sec = difftime (time2,time1);
  std::cout << "Time required for execution: " << diff_sec << " seconds."<<std::endl;
  //------------------//

}



void ReadInterfile(std::istream& FileHeader, std::string OutVol, std::string path)
{
  //By Marco Pizzichemi
  
      int key = 1; //key variable to distinguish pem (1), echo (2) and elasto (3)
      std::cout << "-------- PEM --------" << std::endl;
      
      unsigned int i,h;      
      std::string s1,s2;
      
      // defines the strings to look for
      std::string FileName = "name of data file := ";       
      std::string NumberFormat = "!number format := ";             
      std::string SizeX = "!matrix size [1] := ";         
      std::string SizeY = "!matrix size [2] := ";
      std::string SizeZ = "!matrix size [3] := "; 
      std::string SpacingX = "scaling factor (mm/pixel) [1] := ";	 
      std::string SpacingY = "scaling factor (mm/pixel) [2] := ";	 
      std::string SpacingZ = "scaling factor (mm/pixel) [3] := ";

      std::string RawFile,Format, Tempor;            
      float spacing[3],size[3];

      //looks for the strings defined above, assing the rest of the line to the correct variable      
      while(getline(FileHeader,s1))	
      {
	i = s1.find(FileName); 	 
	if (i==0)  
	{
	  RawFile.assign(s1, FileName.size() , s1.size());
          //int found;
          //found = s1.find_last_of("/");
	  //RawFile.assign(s1.substr(found+1));   
		if (RawFile[0] != '/') RawFile = path + RawFile;
	}	
	i = s1.find(NumberFormat);        
	if (i==0)  
	{
	  Format.assign(s1, NumberFormat.size() , s1.size());	 
	}  	
	i = s1.find(SizeX); 
	if (i==0)  
	{ 
	  s2.assign(s1, SizeX.size() , s1.size());	 
	  size[0] = atof(s2.c_str());
	}	
	i = s1.find(SizeY); 
	if (i==0)  
	{	 
	  s2.assign(s1, SizeY.size() , s1.size());	 
	  size[1] = atof(s2.c_str());
	}
	i = s1.find(SizeZ); 
	if (i==0)  
	{	 
	  s2.assign(s1, SizeZ.size() , s1.size());	 
	  size[2] = atof(s2.c_str());  
	} 
	i = s1.find(SpacingX); 
	if (i==0)         
	{
	  s2.assign(s1, SpacingX.size() , s1.size()); 
	  spacing[0] = atof(s2.c_str());
	}    
	i = s1.find(SpacingY); 
	if (i==0)          
	{ 
	  s2.assign(s1, SpacingY.size() , s1.size());	 
	  spacing[1] = atof(s2.c_str()); 
	} 
	i = s1.find(SpacingZ); 
	if (i==0) 
	{
	  s2.assign(s1, SpacingZ.size() , s1.size());
	  spacing[2] = atof(s2.c_str());
	}	
      }
      
      //creating and filling SpatialInfo
      std::vector<float> SpatialInfo;      
      for(h=0;h<9;h++) 
      {	 
        //for PEM files the origin is moved in order to center the PEM image in (0,0,0)
	if(h<3)	
	{
	  SpatialInfo.push_back(-(size[h]*spacing[h])/2.0);  
	}
	if((h>2) && (h<6)) 
	{	
	  SpatialInfo.push_back(spacing[h-3]); 
	}
	if((h>5) && (h<9)) 
	{	 
	  SpatialInfo.push_back(size[h-6]);	 
	} 
      }      
      
      // output the SpatialInfo vector, for control and debugging
      /*  for(h=0;h<9;h++)
	{
	std::cout << "Spatial Info vector for PEM image = " << SpatialInfo[h] << std::endl;
	}
      */
      // output the path, for control and debugging
      // std::cout << "Absolute path of PEM file = " << path.c_str() << std::endl;
  
      //executing RawImport for float data. PEM files should always be float, so RawImport<float> is called         
      RawImport<float>(RawFile,SpatialInfo,OutVol,key);
     
}
