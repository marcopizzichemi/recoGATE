/*======================================================================
  Attenuation Correction for the ClearPEM scanner
  Suitable both for STIR and List Mode reconstruction algorithms.
  
  Inputs: Initial reconstructed image (can be reconstructed with a big voxel and filter)
          .lm file of the same data

  Outputs: LMF where the weight field has an attenuation correction factor

  Author: C. S. Ferreira

  Modification History: Dec 2011 - Jan 2012 (final version)

  ======================================================================*/ 

#include <time.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <assert.h>


//Intersection line->mesh
#include <vtkPolyData.h>
#include <vtkOBBTree.h>
#include <vtkModifiedBSPTree.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>


//Decimation
#include <vtkDecimatePro.h>

 
int main(int argc, char *argv[])
{
  
  if( argc < 4 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " inputMeshFile.vtk  inputLmfFile.lm  outputAttenCorrFile.lm  ReconAlgorithm(STIR,LM)" << std::endl;
    return EXIT_FAILURE;
    }


   //TIME STAMP
  time_t time1,time2;
  time(&time1);

  //xDim= atof(argv[1]);
  //yDim = atof(argv[2]);
  //zdim = atof(argv[3]);



//--------------LINE-MESH INTERSECTION-----------------//


  vtkSmartPointer<vtkPolyDataReader> p_reader = vtkSmartPointer<vtkPolyDataReader>::New();
  p_reader->SetFileName(argv[1]);   //"tfinal_mesh.vtk");
  p_reader->Update();



  /* //Option1: vtkOBBTree (TOO SLOW!!!!)
  vtkOBBTree* Tree = vtkOBBTree::New();
  vtkSmartPointer<vtkOBBTree> tree = vtkSmartPointer<vtkOBBTree>::New();
  Tree->SetDataSet(p_reader->GetOutput());
  Tree->BuildLocator();
  */ 
  
  //Option2 (known as faster): ModifiedBSPTree IntersectWithLine
  // Create the tree
  vtkSmartPointer<vtkModifiedBSPTree> bspTree = vtkSmartPointer<vtkModifiedBSPTree>::New();
  bspTree->SetDataSet(p_reader->GetOutput());
  bspTree->BuildLocator();
  double tolerance = 0.001;
  
//......Outputs: 
  //Parametric coordinate of intersection (0 (corresponding to p1) to 1 (corresponding to p2))
  //points // The coordinates of the intersections
  //cellIds
  

  
  ///////////////////////////////////////////////////
  //Read .lm for LOR intersection with mesh  

  struct LmSTIREvent{    //STIRFormat
    float x1; 
    float y1;
    float z1;
    float x2;
    float y2;
    float z2;
    float weight;   //known as increment on lmmerger2
    unsigned char replace; }  __attribute__((__packed__));
  
 
  struct LmLMEvent{   //DKFZFormat
    float x1;
    float y1;
    float z1;
    float x2;
    float y2;
    float z2;
    float angle;
    float weight; } __attribute__((__packed__));
  
  
  struct EventCoords{
    float x1;
    float y1;
    float z1;
    float x2;
    float y2;
    float z2; };
  
  
  //Input
  FILE * infile = NULL;
  if(strcmp(argv[2], "-") == 0) 
    infile = stdin;
  else
    infile = fopen(argv[2], "r");
				
  if (infile == NULL) 
    {
    fprintf(stderr, "File %s does not exist\n", argv[2]);
    return -1;
    }
  
 
  //Output
  //lmEvent outLMF;
  FILE * out = NULL;
  if(strcmp(argv[3], "-") == 0) 
    out = stdout;
  else
    out = fopen(argv[3], "w");
				
  if (out == NULL) 
    {
    fprintf(stderr, "Can't open output file \n");
    return -1;
    }
  
  LmSTIREvent eventSTIR;
  LmLMEvent eventLM;
  EventCoords event, eventCoord;

  //lmEvent EventCoord; 
  int cont=-1;

  vtkPoints * points = vtkPoints::New();
  vtkIdList * cellIds = vtkIdList::New(); 
  float lorLength;
  double intersectionPoint[3];
  double intersectionPoint1[3];
  double intersectionPoint2[3];

  float niu_h2o;
  float ACF;
  //Medium parameters
  niu_h2o = 0.0958;        //(cm-1) niu = densH2O * sigma
  //niu_h2o = 0.0974;        //(cm-1) niu = densGe68cyl * sigma
  niu_h2o = niu_h2o * (1/10.);  //(mm-1)

  int tangentialLORs = 0;
  int noIntersection = 0;
  int doubleIntersection = 0;
  int tripleIntersection = 0;
  int quadIntersection = 0;
  int multipleIntersection = 0;
  int multipleIntersectionOdd = 0;
  int multipleIntersectionEven = 0;
  vtkIdType lastPoint;
  double weightCountInitial = 0;
  double weightCountFinal = 0;

  
  void * buffer;
  int sizeEvent;
  if (strcmp(argv[4], "STIR") == 0) 
    {
      buffer = (void*)&eventSTIR; sizeEvent = sizeof(eventSTIR); 
    } 
  else 
    { 
      buffer = (void*)&eventLM; 
      sizeEvent = sizeof(eventLM); 
    }
  while(fread(buffer, sizeEvent, 1, infile) == 1)
    //while(fread((void*)&Event, sizeof(Event), 1, infile) == 1)
    { 
      
      ACF = 1.0;
      float norm = 0;
      float normTotal=0;
      cont++;
      // if(cont>5){break;}
      if (cont % 100000 == 0) 	{std::cerr << "Event number: " << cont << std::endl;}

      if (strcmp(argv[4], "STIR") == 0)
	{
	  event.x1 = eventSTIR.x1;
	  event.y1 = eventSTIR.y1;
	  event.z1 = eventSTIR.z1;
	  event.x2 = eventSTIR.x2;
	  event.y2 = eventSTIR.y2;
	  event.z2 = eventSTIR.z2;
	  weightCountInitial += eventSTIR.weight;
	}
      
      if (strcmp(argv[4], "LM") == 0)
	{
	  event.x1 = eventLM.x1;
	  event.y1 = eventLM.y1;
	  event.z1 = eventLM.z1;
	  event.x2 = eventLM.x2;
	  event.y2 = eventLM.y2;
	  event.z2 = eventLM.z2;
	  weightCountInitial += eventLM.weight;
	  
	  //Rotation for N angles////////
	  //Rotation Matrix//
	  // 1       0             0
	  // 0   cos(rotation_deg)    sin(rotation_deg)
	  // 0   -sin(rotation_deg)   cos(rotation_deg)
	  //////////////////////////////////
	  float y1 = event.y1 * cosf(eventLM.angle) + event.z1 * sinf(eventLM.angle);
	  float z1 = - event.y1 * sinf(eventLM.angle) + event.z1 * cosf(eventLM.angle);
	  float y2 = event.y2 * cosf(eventLM.angle) + event.z2 * sinf(eventLM.angle);
	  float z2 = - event.y2 * sinf(eventLM.angle ) + event.z2 * cosf(eventLM.angle);

	  //std::cerr<<"angle   "<< eventLM.angle << std::endl;
	  //std::cerr<<"weight   "<< eventLM.weight << std::endl;

	  event.y1 = y1;
	  event.y2 = y2;
	  event.z1 = z1;
	  event.z2 = z2;

	}  


      // ClearPEM -> STIR coordinates
      // vertical -> YYstir = -YYcp
      // horizontal -> XXstir = ZZcp
      // scanner axis -> ZZstir = -XXcp
      eventCoord.x1 = event.z1;
      eventCoord.x2 = event.z2;
      eventCoord.y1 = -event.y1;
      eventCoord.y2 = -event.y2;
      eventCoord.z1 = -event.x1;
      eventCoord.z2 = -event.x2;
      //std::cerr << Event.x1 << ", " << Event.y1 <<  ", " << Event.z1<< "\n" <<std::endl;
      //std::cerr << Event.x2 << ", " << Event.y2 <<  ", " << Event.z2<< "\n" <<std::endl;

      
      double LineP0[3] = {eventCoord.x1,eventCoord.y1,eventCoord.z1};
      double LineP1[3] = {eventCoord.x2,eventCoord.y2,eventCoord.z2};
      //std::cerr << EventCoord.x1 << ", " << EventCoord.y1 <<  ", " << EventCoord.z1<< "\n" <<std::endl;
      //std::cerr << EventCoord.x2 << ", " << EventCoord.y2 <<  ", " << EventCoord.z2<< "\n" <<std::endl;
      
      lorLength = sqrt(pow(eventCoord.x1-eventCoord.x2,2)+(pow(eventCoord.y1-eventCoord.y2,2))+(pow(eventCoord.z1-eventCoord.z2,2)));
      //if (lorLength > (sqrt(pow(axialDim,2)+pow(DH_distance,2)+pow(DH_distance,2)))) {std::cerr << "Error on lor length" << std::endl;}


      //////////////////////////////////////////////////////////////////

      //-----------OBBTree---------------//
      //intersect the locator with the line
      // Tree->IntersectWithLine(LineP0, LineP1, points, cellIds);


      //------------MODIFIED BSPTree--------------------//
      vtkIdType iD = bspTree->IntersectWithLine(LineP0, LineP1, tolerance, points, cellIds); 
      //std::cerr << "iD: " << iD << std::endl;
      //std::cerr << "Number of intersection points : " << points->GetNumberOfPoints() << std::endl;
      //std::cerr << "Number of IDs: " << cellIds->GetNumberOfIds() << std::endl;
      
      if (iD == 0) 
	{
	  ACF = 1.0;
	  noIntersection++;
	}
      else if (points->GetNumberOfPoints() == 1)
	{
	  ACF = 1.0;
	  tangentialLORs++;
	}
      else
	{
	  if (points->GetNumberOfPoints() == 2) {doubleIntersection++;}
	  if (points->GetNumberOfPoints() == 3) {tripleIntersection++;}
	  if (points->GetNumberOfPoints() == 4) {quadIntersection++;}
	  if ((points->GetNumberOfPoints() > 4)) {multipleIntersection++;}
	  if ((points->GetNumberOfPoints() > 4) && ((points->GetNumberOfPoints()) % 2 != 0)) {multipleIntersectionOdd++;}
	  if ((points->GetNumberOfPoints() > 4) && ((points->GetNumberOfPoints()) % 2 == 0)) {multipleIntersectionEven++;}
	  
	  if ((points->GetNumberOfPoints()) % 2 == 0)
	    {
	      for (vtkIdType pt = 0; pt < points->GetNumberOfPoints() ; pt=pt+2)
		{
		  norm = 0;
		  //Way IN
		  points->GetPoint(pt,intersectionPoint);
		  intersectionPoint1[0] = intersectionPoint[0]; 
		  intersectionPoint1[1] = intersectionPoint[1]; 
		  intersectionPoint1[2] = intersectionPoint[2];
		  //std::cerr << "Intersection: " << intersectionPoint1[0] << ", " << intersectionPoint1[1] << ", " << intersectionPoint1[2] << std::endl;
		  
		  //Check coordinates
		  /*	  if ((intersectionPoint[0] > (xDim/2)) || (intersectionPoint[0] < -(xDim/2)) || (intersectionPoint[1] > (yDim/2)) || (intersectionPoint[1] < -(yDim/2)) || (intersectionPoint[2] > (axialDim/2)) || (intersectionPoint[2] < -(axialDim/2))) 
		    {std::cerr << "out of bounds intersection coordinates! " << std::endl;
		      std::cerr << intersectionPoint[0] << " , " << intersectionPoint[1] << " , " << intersectionPoint[2] << std::endl;}
		  */
		  //Way OUT
		  points->GetPoint(pt+1,intersectionPoint);
		  intersectionPoint2[0] = intersectionPoint[0]; 
		  intersectionPoint2[1] = intersectionPoint[1]; 
		  intersectionPoint2[2] = intersectionPoint[2];
		  //std::cerr << "Intersection: " << intersectionPoint2[0] << ", " << intersectionPoint2[1] << ", " << intersectionPoint2[2] << std::endl;
		  
		  //Check coordinates
		  /*	  if ((intersectionPoint[0] > (xDim/2)) || (intersectionPoint[0] < -(xDim/2)) || (intersectionPoint[1] > (yDim/2)) || (intersectionPoint[1] < -(yDim/2)) || (intersectionPoint[2] > (axialDim/2)) || (intersectionPoint[2] < -(axialDim/2)))
		    {std::cerr << "out of bounds intersection coordinates! " << std::endl;
		      std::cerr << intersectionPoint[0] << " , " << intersectionPoint[1] << " , " << intersectionPoint[2] << std::endl;}
		  */
		  
		  //-----------ACF Calculation---------------//	       		  
		  norm = sqrt(pow(intersectionPoint1[0]-intersectionPoint2[0],2)+(pow(intersectionPoint1[1]-intersectionPoint2[1],2))+(pow(intersectionPoint1[2]-intersectionPoint2[2],2)));
		  normTotal = normTotal+norm;
		} 
	    }

	  lastPoint = (points->GetNumberOfPoints())-1;
	  if ((points->GetNumberOfPoints()) % 2 != 0)
	    {
	      //for (vtkIdType pt = 0; pt < points->GetNumberOfPoints() ; pt=pt+2)
	      // {
	      norm = 0;
	      points->GetPoint(0,intersectionPoint);
	      intersectionPoint1[0] = intersectionPoint[0]; 
	      intersectionPoint1[1] = intersectionPoint[1]; 
	      intersectionPoint1[2] = intersectionPoint[2];
	      //Check coordinates
	      /*  if ((intersectionPoint[0] > (xDim/2)) || (intersectionPoint[0] < -(xDim/2)) || (intersectionPoint[1] > (yDim/2)) || (intersectionPoint[1] < -(yDim/2)) || (intersectionPoint[2] > (axialDim/2)) || (intersectionPoint[2] < -(axialDim/2)))
		{std::cerr << "out of bounds intersection coordinates! " << std::endl;
		  std::cerr << intersectionPoint[0] << " , " << intersectionPoint[1] << " , " << intersectionPoint[2] << std::endl;}	      
	      */
	      points->GetPoint(lastPoint,intersectionPoint);
	      intersectionPoint2[0] = intersectionPoint[0]; 
	      intersectionPoint2[1] = intersectionPoint[1]; 
	      intersectionPoint2[2] = intersectionPoint[2];
	      //Check coordinates
	      /*  if ((intersectionPoint[0] > (xDim/2)) || (intersectionPoint[0] < -(xDim/2)) || (intersectionPoint[1] > (yDim/2)) || (intersectionPoint[1] < -(yDim/2)) || (intersectionPoint[2] > (axialDim/2)) || (intersectionPoint[2] < -(axialDim/2))) 
		{std::cerr << "out of bounds intersection coordinates! " << std::endl;
		  std::cerr << intersectionPoint[0] << " , " << intersectionPoint[1] << " , " << intersectionPoint[2] << std::endl;}	
	      */
	      norm = sqrt(pow(intersectionPoint1[0]-intersectionPoint2[0],2)+(pow(intersectionPoint1[1]-intersectionPoint2[1],2))+(pow(intersectionPoint1[2]-intersectionPoint2[2],2)));
	      normTotal = norm;
	      
	    }
	  
	  // std::cerr<<"normTotal : " <<normTotal<<std::endl;
	  ACF = exp(niu_h2o * normTotal);
	  //std::cerr<<"niu_h2o : " <<niu_h2o<<std::endl;
	  points->Reset();
	  cellIds->Reset();
	}
    
      //-------------Write weighed LOR to file----------------//
      if (strcmp(argv[4], "STIR") == 0)
	{
	  LmSTIREvent outLMF;
	  outLMF = eventSTIR;
	  outLMF.weight = eventSTIR.weight * ACF;
	  outLMF.replace = 1 ;   //set this on to do a sinogram with the correction coefficients only
	  fwrite((void*) &outLMF, sizeof(outLMF), 1, out);
	  weightCountFinal += outLMF.weight;
	  // std::cerr << outLMF.x1 << " , "<<outLMF.y1<< " , " <<outLMF.z1 << " , " << outLMF.x2 << " , "<<outLMF.y2<< " , " <<outLMF.z2 << " , "<< outLMF.increment<< " , " << int(outLMF.replace)<< std::endl;
	}

      if (strcmp(argv[4], "LM") == 0)
	{
	  LmLMEvent outLMF;
	  outLMF = eventLM;
	  outLMF.weight = eventLM.weight * ACF;
	  fwrite((void*) &outLMF, sizeof(outLMF), 1, out);
	  weightCountFinal += outLMF.weight;
	  // std::cerr << outLMF.x1 << " , "<<outLMF.y1<< " , " <<outLMF.z1 << " , " << outLMF.x2 << " , "<<outLMF.y2<< " , " <<outLMF.z2 << " , "<< outLMF.increment<< " , " << int(outLMF.replace)<< std::endl;
	}

    }

  
  fclose(infile);
  std::cerr << "Initial events total weight: " << weightCountInitial << std::endl;
  std::cerr << "Final events total weight: " << weightCountFinal << std::endl;
  std::cerr << "No intersection: " << noIntersection << std::endl;
  std::cerr << "Tangential LORs: " << tangentialLORs << std::endl; 
  std::cerr << "double intersection: " << doubleIntersection << std::endl;
  std::cerr << "triple intersection: " << tripleIntersection << std::endl;
  std::cerr << "quad intersection: " << quadIntersection << std::endl;
  std::cerr << "multiple intersection: " << multipleIntersection << std::endl;
  std::cerr << "multiple intersection odd: " << multipleIntersectionOdd << std::endl;
  std::cerr << "multiple intersection even: " << multipleIntersectionEven << std::endl;
   

  //TIME STAMP--------//
  time(&time2);
  double diff_sec = difftime (time2,time1);
  std::cout << "Time required for execution: " << diff_sec << " seconds."<<std::endl;  
  //------------------//
 
  return 0;

}

