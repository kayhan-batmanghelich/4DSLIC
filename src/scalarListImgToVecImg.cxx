
#include "itkImage.h"
#include "itkComposeImageFilter.h"
#include "itkVectorImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
//#include "itkRescaleIntensityImageFilter.h"
#include "itkVectorMagnitudeImageFilter.h"
#include "itkVectorRescaleIntensityImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkMultiplyImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"

 
typedef itk::VectorImage<double, 3>             VectorImageType;
typedef itk::Image<double, 3>                   ScalarImageType;
typedef itk::ImageFileReader<ScalarImageType>   ImageReaderType ;
typedef itk::ImageFileWriter<VectorImageType>   VectorImageWriterType ;
typedef itk::ImageFileWriter<ScalarImageType>   ScalarImageWriterType ;


typedef itk::ComposeImageFilter<ScalarImageType>                                    ImageToVectorImageFilterType;
typedef itk::VectorMagnitudeImageFilter< VectorImageType, ScalarImageType >         VectorMagnitudeFilterType ;
//typedef itk::RescaleIntensityImageFilter< ScalarImageType, ScalarImageType >        rescaleFilterType;
typedef itk::VectorRescaleIntensityImageFilter<VectorImageType, VectorImageType>    VectorRescaleFilterType;
typedef itk::MinimumMaximumImageCalculator<ScalarImageType>                         ImageCalculatorFilterType;
typedef itk::MultiplyImageFilter<ScalarImageType, ScalarImageType, ScalarImageType> MultiplyImageFilterType;
typedef itk::VectorIndexSelectionCastImageFilter<VectorImageType, ScalarImageType>   IndexSelectionType;

 
//static void CreateImage(ScalarImageType::Pointer image);
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

int main(int argc, char *argv[])
{
   
   string   inListFn;
   string   outVecImgFn;
   string   outMagImgFn ;
   float    normMaxVal ;
   if(argc > 4)
   {
      inListFn = string(argv[1]);
      outVecImgFn = string(argv[2]);
      outMagImgFn = string(argv[3]) ;
      normMaxVal = atof(argv[4]) ;
   }
   else
   {
      printf("Usage: %s inputListFilename   outputVolumeImg  outputMagnitudeImg   maxMagValue\n",argv[0]);
      cout << endl ;
      cout << "inputListFilename : a file  containing list of input images (input)"  << endl ;
      cout << "outputVolumeImg : vector image file for output (output)" << endl ;
      cout << "outputMagnitudeImg : magnitude of the vector image (output) " << endl ;
      cout << "maxMagValue : nomalizer value for the scalar image (input) " << endl ;
      cout << endl ;
      return -1;
   }

    ifstream fileList (inListFn.c_str()); 
   
    int i = 0 ;
    ImageToVectorImageFilterType::Pointer imageToVectorImageFilter = ImageToVectorImageFilterType::New();

    // reading in the scale files in a vector image 
    while (! fileList.eof())
    {
        string   fn ;
        getline(fileList,fn) ;
        if (fn != "")
        {
          cout <<  "fn : "  << fn << endl ;
          ImageReaderType::Pointer  reader = ImageReaderType::New() ;
          reader->SetFileName(fn) ;
          reader->Update() ;
          imageToVectorImageFilter->SetInput(i, reader->GetOutput()) ;
          imageToVectorImageFilter->Update() ;
          i++ ;    
        }
    }
    int numImgs = i ;
    VectorImageType::Pointer  vectorImage = imageToVectorImageFilter->GetOutput();


    
    // creating magnitude image
    VectorMagnitudeFilterType::Pointer magnitudeFilter = VectorMagnitudeFilterType::New();
    magnitudeFilter->SetInput( vectorImage );
    magnitudeFilter->Update() ;
 
    // compute the maximume norm of the vector
    ImageCalculatorFilterType::Pointer imageCalculatorFilter   = ImageCalculatorFilterType::New ();
    imageCalculatorFilter->SetImage(magnitudeFilter->GetOutput()) ;
    imageCalculatorFilter->Compute() ;
    double   maxMagValue = imageCalculatorFilter->GetMaximum() ;
    cout << "maximum magnitude image is : " << maxMagValue << endl ;

    // normalizing entries of the vector image
    ImageToVectorImageFilterType::Pointer imageToVectorImageFilter2 = ImageToVectorImageFilterType::New();
    for (i=0;i<numImgs; i++)
    {
           MultiplyImageFilterType::Pointer multiplyImageFilter = MultiplyImageFilterType::New();
           IndexSelectionType::Pointer indexSelectionFilter = IndexSelectionType::New(); 

           indexSelectionFilter->SetIndex(i);
           indexSelectionFilter->SetInput(vectorImage);
           indexSelectionFilter->Update() ;  // get the  i'th entry

           multiplyImageFilter->SetInput( indexSelectionFilter->GetOutput()  ) ;
           multiplyImageFilter->SetConstant(normMaxVal/maxMagValue) ;
           multiplyImageFilter->Update() ; // normalize it

           imageToVectorImageFilter2->SetInput(i, multiplyImageFilter->GetOutput()) ;   // put it back
    }
    VectorImageType::Pointer  normalized_vectorImage = imageToVectorImageFilter2->GetOutput();

    // re-compute the magnitude and re-compute the maximue just to make sure that the rescaling is working
    magnitudeFilter->SetInput( normalized_vectorImage );
    magnitudeFilter->Update() ;
    imageCalculatorFilter->SetImage(magnitudeFilter->GetOutput()) ;
    imageCalculatorFilter->Compute() ;
    cout << "maximum magnitude image after normalization is : " << imageCalculatorFilter->GetMaximum()  << endl ;


    // writing the normalized image to the disk
    cout << "writing the normalized vector image to : " << outVecImgFn << endl ;
    VectorImageWriterType::Pointer  vecImgWriter  = VectorImageWriterType::New() ;
    vecImgWriter->SetInput(normalized_vectorImage) ;
    vecImgWriter->SetFileName(outVecImgFn) ;
    vecImgWriter->Update() ;


    // To write the magnitude image file, we should rescale the gradient values
    // to a reasonable range
    //rescaleFilterType::Pointer rescaler = rescaleFilterType::New();
    //rescaler->SetOutputMinimum(0);
    //rescaler->SetOutputMaximum(normMaxVal);
    //rescaler->SetInput( magnitudeFilter->GetOutput() );
    //rescaler->Update() ;
    //cout << "finish creating magnitude image " << endl ;
    //cout << "writing to " << outMagImgFn << endl ;

    // writing the normalized magnitude image to the dist
    cout << "Writing the normalized magnitude image to the disk :" <<  outMagImgFn << endl ;
    ScalarImageWriterType::Pointer  magImgWriter  = ScalarImageWriterType::New() ;
    magImgWriter->SetInput(magnitudeFilter->GetOutput()) ;
    magImgWriter->SetFileName(outMagImgFn) ;
    magImgWriter->Update() ;

    // re-scaling vector image itself to have the the pre-specified magnitude
    //VectorRescaleFilterType::Pointer rescaleFilter = VectorRescaleFilterType::New();
    //rescaleFilter->SetInput(vectorImage);
    //rescaleFilter->SetOutputMaximumMagnitude(normMaxVal);
    //rescaleFilter->Update();    
    // write the results
    //VectorImageWriterType::Pointer  vecImgWriter2  = VectorImageWriterType::New() ;
    //vecImgWriter->SetInput(rescaleFilter->GetOutput()) ;
    //vecImgWriter->SetFileName("rescaleVecImg.nii,gz") ;
    //vecImgWriter->Update() ;
   
    //return 0 ;
    return EXIT_SUCCESS;
}
 
