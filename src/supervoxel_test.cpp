#include "itkImage.h"
#include "itkVectorImage.h"
#include "itkImageFileReader.h"
 
typedef itk::VectorImage<double, 3>             VectorImageType;
typedef itk::Image<double, 3>                   ScalarImageType;
typedef itk::ImageFileReader<VectorImageType>   VectorImageReaderType ;



 
//static void CreateImage(ScalarImageType::Pointer image);
#include <iostream>
#include <fstream>
#include <string>

#include <LKM.h>
using namespace std;

int main(int argc, char *argv[])
{
   
   string   inVecFn;
   string   outLabFn;
   int      stepSize ;
   int      cubeness ;
   if(argc > 4)
   {
      inVecFn = string(argv[1]);
      outLabFn = string(argv[2]);
      stepSize = atoi(argv[3]) ;
      cubeness = atoi(argv[4]) ;
   }
   else if (argc == 4)
   {
      inVecFn = string(argv[1]);
      outLabFn = string(argv[2]);
      stepSize = atoi(argv[3]) ;
      cubeness = 20 ;
   }
   else
   {
      printf("Usage: %s   inputVectorImg    outputLabelImg  stepSize cubeness \n",argv[0]);
      cout << endl ;
      cout << "inputVoctorImg : input vector image "  << endl ;
      cout << "outputLabelImg : output label image " << endl ;
      cout << "stepSize : Step size, we suggest 10 "  << endl ;
      cout << "cubeness : cubeness, default value is 20 " << endl ;
      cout << " " << endl ;
      return -1;
   }

   VectorImageReaderType::Pointer    vecImg_reader = VectorImageReaderType::New() ;
   vecImg_reader->SetFileName(inVecFn) ;
   vecImg_reader->Update() ;

   // setting up LKM 
   int numlabels ;
   LKM*     lkm = new LKM ;
   sidType** klabels;
   int vstep = 10 ;

   cout << "performing supervoxelization ..." << endl ;
   lkm->DoSupervoxelSegmentationForMultichannelVolume( vecImg_reader->GetOutput(), 
                                                       klabels, 
                                                       numlabels, 
                                                       stepSize,
                                                       cubeness) ;
   cout << "Supervoxelization is Done !" << endl ;
   cout << "Number of Supervoxel is : " << numlabels << endl ;
 
   lkm->WriteLabelVolume(vecImg_reader->GetOutput(),
                         klabels,
                         outLabFn) ;


   // cleaning up
   delete lkm;
   int depth ;
   VectorImageType::SizeType  size ;
   size = vecImg_reader->GetOutput()->GetLargestPossibleRegion().GetSize() ;
   depth = size[2] ;
   for(int d=0;d<depth;d++)
   {
      delete[] klabels[d];
   }
   delete[] klabels;
   
   return EXIT_SUCCESS;
}
 
