#include "itkRawImageIO.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
//#include "itkOrientedImage.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkSquareImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"


template <class T> void RawImport(std::string rawFile, std::vector<float> SpatialInfo, std::string outVolume,int key)
{
// define reader
typedef itk::Image<T,3> ImageType;
typedef itk::ImageFileReader<ImageType> ReaderType;

typename ReaderType::Pointer reader = ReaderType::New();

typename itk::RawImageIO<T,3>::Pointer io;
io = itk::RawImageIO<T,3>::New();

// define a cast filter. Will be applied to all images, even if unuseful with pem uses
typedef itk::Image<float,3> ImageTypeAfterCast;
typedef itk::CastImageFilter<ImageType,ImageTypeAfterCast> castFilterType;
typename castFilterType::Pointer castFilter = castFilterType::New();

// define a 2 rescale filters
typedef itk::RescaleIntensityImageFilter<ImageTypeAfterCast,ImageTypeAfterCast>  RescaleFilterType;
typename RescaleFilterType::Pointer    rescaleFilter    = RescaleFilterType::New();
typename RescaleFilterType::Pointer    rescaleFilter2    = RescaleFilterType::New();

// define square filter
typedef itk::SquareImageFilter<ImageTypeAfterCast,ImageTypeAfterCast>  squareFilterType;
typename squareFilterType::Pointer squareFilter = squareFilterType::New();

//define a image calculator filter
typedef itk::MinimumMaximumImageCalculator <ImageType> ImageCalculatorFilterType;
typename ImageCalculatorFilterType::Pointer imageCalculatorFilter = ImageCalculatorFilterType::New();

typedef itk::ImageFileWriter<ImageTypeAfterCast> WriterType;
typename WriterType::Pointer writer = WriterType::New();


io->SetFileName(rawFile);
io->SetHeaderSize(0);
io->SetByteOrderToLittleEndian();
io->SetPixelType(itk::ImageIOBase::SCALAR);
io->SetNumberOfComponents(1);

for(unsigned int i=0; i<3; i++){
	io->SetDimensions(i, SpatialInfo[i+6]);
	io->SetSpacing(i, SpatialInfo[i+3]);
	io->SetOrigin(i, SpatialInfo[i]);
}

reader->SetFileName(rawFile);
reader->SetImageIO(io);

reader->Update();

imageCalculatorFilter->SetImage(reader->GetOutput());
imageCalculatorFilter->Compute();


std::cout << "Min pixel position: " << imageCalculatorFilter->GetIndexOfMinimum() << std::endl;
std::cout << "Max pixel position: " << imageCalculatorFilter->GetIndexOfMaximum()  << std::endl;

std::cout << "Minimum pixel value: " << imageCalculatorFilter->GetMinimum() << std::endl;
std::cout << "Maximum pixel value: " << imageCalculatorFilter->GetMaximum()  << std::endl;

T maximum = imageCalculatorFilter->GetMaximum();
T minimum = imageCalculatorFilter->GetMinimum();


//first, cast image to float

castFilter->SetInput(reader->GetOutput());


//TO DO (BUT ONLY FOR ELASTO!) here = setup a filter to pass from m/s (input) to kpascal -> normalize image dividing all input pixels by maximum (which means also convert them to float)  and multiply by 10 (maximum speed always 10 m/s, to investigate more) the final_pixel = 3*(norm_pixel ^2). Furthermore, image will have to be threshold cutting above 100 Kpascal when displaying

rescaleFilter->SetInput(castFilter->GetOutput());
writer->SetFileName(outVolume.c_str());

if(key == 1) //pem image
{ 
  rescaleFilter->SetOutputMinimum(0);
  rescaleFilter->SetOutputMaximum(30000);
  rescaleFilter->Update();
  writer->SetInput(rescaleFilter->GetOutput());//
}
else if(key == 2) // echo image
{
  rescaleFilter->SetOutputMinimum(0);
  rescaleFilter->SetOutputMaximum(30000);
  rescaleFilter->Update();
  writer->SetInput(rescaleFilter->GetOutput());//
}
else if(key == 3) //elasto image
{
  rescaleFilter->SetOutputMinimum(0);
  rescaleFilter->SetOutputMaximum(10); //rescale between 0 and 10 m/s
  rescaleFilter->Update();
  squareFilter->SetInput(rescaleFilter->GetOutput());
  squareFilter->Update();
  rescaleFilter2->SetInput(squareFilter->GetOutput());
  rescaleFilter2->SetOutputMinimum(0);
  rescaleFilter2->SetOutputMaximum(300); //rescale between 0 and 300 KPa
  rescaleFilter2->Update();
  writer->SetInput(rescaleFilter2->GetOutput());
}

writer->Update();

//rescale the image to fit the 1-30000 range

//



// 4 - innesca la chain



}

