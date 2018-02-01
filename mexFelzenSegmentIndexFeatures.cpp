/**
 * mexFelzenSegmentIndexFeatures.cpp
 *
 * A modified version of the modified version of the Felzenszwalb algorithm
 * found in Jasper Uijlings Selective Search code to interface with Matlab.
 *
 * Generalized to handle images where the third dimension is an arbitrary
 * length feature vector for each pixel.
 % Copyright (c) 2017 The Johns Hopkins University/Applie?d Physics Laboratory
% This software was developed at The Johns Hopkins University/Applied Physics 
% Laboratory (“JHU/APL”) that is the author thereof under the “work made for hire” 
% provisions of the copyright law. Permission is hereby granted, free of charge, 
% to any person obtaining a copy of this software and associated documentation 
% (the “Software”), to use the Software without restriction, including without 
% limitation the rights to copy, modify, merge, publish, distribute, sublicense, 
% and/or sell copies of the Software, and to permit others to do so, subject to 
% the following conditions:
% 1. This LICENSE AND DISCLAIMER, including the copyright notice, 
% shall be included in all copies of the Software, 
% including copies of substantial portions of the Software;
% 
% 2. JHU/APL assumes no obligation to provide support of any kind with regard to the Software. 
% This includes no obligation to provide assistance in using the Software nor to 
% provide updated versions of the Software; and
% 
% 3. THE SOFTWARE AND ITS DOCUMENTATION ARE PROVIDED AS IS AND WITHOUT ANY 
% EXPRESS OR IMPLIED WARRANTIES WHATSOEVER. ALL WARRANTIES INCLUDING, BUT NOT 
% LIMITED TO, PERFORMANCE, MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE,
% AND NONINFRINGEMENT ARE HEREBY DISCLAIMED. USERS ASSUME THE ENTIRE RISK AND 
% LIABILITY OF USING THE SOFTWARE. USERS ARE ADVISED TO TEST THE SOFTWARE 
% THOROUGHLY BEFORE RELYING ON IT. IN NO EVENT SHALL THE JOHNS HOPKINS UNIVERSITY 
% BE LIABLE FOR ANY DAMAGES WHATSOEVER, INCLUDING, WITHOUT LIMITATION, ANY 
% LOST PROFITS, LOST SAVINGS OR OTHER INCIDENTAL OR CONSEQUENTIAL DAMAGES, 
% ARISING OUT OF THE USE OR INABILITY TO USE THE SOFTWARE.”
 */
/*This code has been modified to allow more multidimensional features and initial guidances
--Marc Bosch 09/13/2016*/

#include <cmath>
#include "mex.h"
#include "segment-image.h"

#define UInt8 char

#define CENTROID_SIZE_LIMIT 1000


//#define median3(a,b,c) ((a)>(b) ? (((b)>(c)) ? (b) : (a)>(c) ? (c) : (a)) : ((a)) )

inline float med3(float a, float b, float c)
{
    if (a>b)
    {
        if(b>c)
            return b;
        else
        {
            if(a>c)
                return c;
            else
                return a;
        }
    }
    else
    {
        if(a>c)
            return a;
        else
        {
            if(b>c)
                return c;
            else
                b;            
        }
    }
}

/*
 * Find the difference between feature vectors at different x,y coordinates
 * in the image.
 */
float diff_features(image_features *im_feat, int x1, int y1, int x2, int y2) {

  const int width  = im_feat->width;
  const int height = im_feat->height;
  const int size   = width * height;

  float val = 0;

  for (int i = 0; i <  mmin2(3,im_feat->num_features-1); i++) {
    val += square(
      im_feat->features[ x1*height + y1 + i*size ] -
      im_feat->features[ x2*height + y2 + i*size ] );
  }
  /*if ((im_feat->num_features>1) && (im_feat->features[ x1*height + y1 + (im_feat->num_features-1)*size ]!=
          im_feat->features[ x2*height + y2 + (im_feat->num_features-1)*size ]))
      val *= 2;
  if ((im_feat->num_features>1) && (im_feat->features[ x1*height + y1 + (im_feat->num_features-1)*size ]==
          im_feat->features[ x2*height + y2 + (im_feat->num_features-1)*size ]) && 
          (im_feat->features[ x1*height + y1 + (im_feat->num_features-1)*size ] != 0))
      val *= 2;*/
  
  /*if(im_feat->num_features>3)
  {
      //angle between vectors
      float val1=0;
      float val2=0;
      float val3=0;
      for (int i = 0; i < im_feat->num_features-1; i++) {
        val1 += square(im_feat->features[ x1*height + y1 + i*size ]*
                im_feat->features[ x1*height + y1 + i*size ]);
        val2 +=square(im_feat->features[ x2*height + y2 + i*size ]*
                im_feat->features[ x2*height + y2 + i*size ]);
        val3 +=square(im_feat->features[ x2*height + y2 + i*size ]*
                im_feat->features[ x1*height + y1 + i*size ]);

      }
      float den=sqrt(val1)*sqrt(val2);
      den = (den==0)? 1 : den;
      val = val3/den;
   }*/
  return sqrt(val);
}



double *segment_image_index(image_features *im_feat, float c,
        int min_size, int *num_ccs, float threshold, bool check_centroid) {

  int width = im_feat->width, height = im_feat->height;

  // build graph
  edge *edges = new edge[width*height*4];
  int num = 0;
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      if (x < width-1) {
        edges[num].a = y * width + x;
        edges[num].b = y * width + (x+1);
        edges[num].w = 0;        
        edges[num].w += diff_features(im_feat, x, y, x+1, y);
        /*edges[num].w += diff_features(im_feat, x, mmin2(y+1,height-1), x+1, mmin2(y+1,height-1));
        edges[num].w += diff_features(im_feat, x, mMax2(y-1,1), x+1, mMax2(y-1,1));
        edges[num].w /= 3;
        edges[num].w = med3(diff_features(im_feat, x, y, x+1, y),
         diff_features(im_feat, x, mmin2(y+1,height-1), x+1, mmin2(y+1,height-1)),
         diff_features(im_feat, x, mMax2(y-1,1), x+1, mMax2(y-1,1)));*/
        edges[num].w *=im_feat->features[ x*height + y + 3*width*height ];
        num++;
      }

      if (y < height-1) {
        edges[num].a = y * width + x;
        edges[num].b = (y+1) * width + x;        
        edges[num].w = 0;        
        edges[num].w += diff_features(im_feat, x, y, x, y+1);
        /*edges[num].w += diff_features(im_feat, mmin2(x+1,width-1), y, mmin2(x+1,width-1), y+1);
        edges[num].w += diff_features(im_feat, mMax2(x-1,1), y, mMax2(x-1,1), y+1);
        edges[num].w /= 3;
        edges[num].w = med3(diff_features(im_feat, x, y, x, y+1),
                diff_features(im_feat, mmin2(x+1,width-1), y, mmin2(x+1,width-1), y+1),
                diff_features(im_feat, mMax2(x-1,1), y, mMax2(x-1,1), y+1));*/
        edges[num].w *=im_feat->features[ x*height + y + 6*width*height ];
        num++;
      }

      if ((x < width-1) && (y < height-1)) {
        edges[num].a = y * width + x;
        edges[num].b = (y+1) * width + (x+1);
        edges[num].w = 0;        
        edges[num].w += diff_features(im_feat, x, y, x+1, y+1);
       /* edges[num].w += diff_features(im_feat, x, mMax2(y-1,1), x+1, y);
        edges[num].w += diff_features(im_feat, mMax2(x-1,1), y, x, y+1);
        edges[num].w /= 3;
        edges[num].w = med3(diff_features(im_feat, x, y, x+1, y+1),
         diff_features(im_feat, x, mMax2(y-1,1), x+1, y),
         diff_features(im_feat, mMax2(x-1,1), y, x, y+1));*/
        edges[num].w *=im_feat->features[ x*height + y + 4*width*height ];
        num++;
      }

      if ((x < width-1) && (y > 0)) {
        edges[num].a = y * width + x;
        edges[num].b = (y-1) * width + (x+1);
        edges[num].w = 0;        
        edges[num].w += diff_features(im_feat, x, y, x+1, y-1);
       /* edges[num].w += diff_features(im_feat, x, mmin2(y+1,height-1), x+1, y);
        edges[num].w += diff_features(im_feat, mMax2(x-1,1), y, x, y-1);
        edges[num].w /= 3;
        edges[num].w = med3(diff_features(im_feat, x, y, x+1, y-1),
         diff_features(im_feat, x, mmin2(y+1,height-1), x+1, y),
         diff_features(im_feat, mMax2(x-1,1), y, x, y-1));*/
        edges[num].w *=im_feat->features[ x*height + y + 5*width*height ];
        num++;
      }
    }
  }

  // segment
  universe *u = segment_graph(width*height, num, edges, c, im_feat);

  // post process small components
  for (int i = 0; i < num; i++) {
    int a = u->find(edges[i].a);
    int b = u->find(edges[i].b);
    if ((a != b) && ((u->size(a) < min_size) || (u->size(b) < min_size)))
      u->join(a, b);
  }


  // 2nd pass: merge similar and irregular segments
  centroid cent;
  for (int i = 0; i < num; i++) {
    int a = u->find(edges[i].a);
    int b = u->find(edges[i].b);

    if (a != b) {
      // Join segments with similar feature vectors
      bool below_thresh = u->features_diff(a, b) < threshold;

      // Join segments if the centroid is not part of the segment
      bool no_centroid = false;
      if (check_centroid && u->size(a) < CENTROID_SIZE_LIMIT) {
        u->cent_pos(&cent, a);
        int c = u->find( floor(cent.y) * width + floor(cent.x) );
        no_centroid = (a != c);
      }

      if (below_thresh || no_centroid) {
        u->join(a, b);
      }
    }
  }

  // Clean up!
  delete [] edges;
  *num_ccs = u->num_sets();

  // pick random colors for each component
  double *colors = new double[width*height];
  for (int i = 0; i < width*height; i++)
    colors[i] = 0;

  int idx = 1;
  double* indexmap = new double[width * height];
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      int comp = u->find(y * width + x);
      if (!(colors[comp])){
          colors[comp] = idx;
          idx = idx + 1;
      }
      indexmap[x * height + y] = colors[comp];
    }
  }

  delete [] colors;
  delete u;

  return indexmap;
}

void mexFunction(int nlhs, mxArray *out[], int nrhs, const mxArray *input[])
{
    // Checking number of arguments
    if(nlhs > 5){
        mexErrMsgTxt("Function has five return values");
        return;
    }

    if(nrhs != 3 && nrhs != 4 && nrhs != 5){
        mexErrMsgTxt("Usage: mexFelzenSegment(UINT8 im, double c, int minSize, float threshold)");
        return;
    }

    if(!mxIsClass(input[0], "single")){
        mexErrMsgTxt("Only image arrays of the single class are allowed.");
        return;
    }

    // Load in arrays and parameters
    float* matIm = (float*) mxGetPr(input[0]);
    int nrDims = (int) mxGetNumberOfDimensions(input[0]);
    if (nrDims < 3) {
        mexErrMsgTxt("Input image must have at least 3 dims!");
        return;
    }
    int* dims = (int*) mxGetDimensions(input[0]);
    double* c = mxGetPr(input[1]);
    double* minSize = mxGetPr(input[2]);
    int min_size = (int) *minSize;

    float threshold = 0;
    if (nrhs >= 4)
      threshold = (double)*mxGetPr(input[3]);

    bool check_centroid = false;
    if (nrhs >= 5)
      check_centroid = (bool)*mxGetPr(input[4]);

    // Get information about the input features
    int height       = dims[0];
    int width        = dims[1];
    int num_features = dims[2];

    image_features im_feat;
    im_feat.features     = matIm;
    im_feat.width        = width;
    im_feat.height       = height;
    im_feat.num_features = num_features;

    /*if (threshold <= 0)
      mexPrintf("Warning: Not merging similar segments!\n");
    if (!check_centroid)
      mexPrintf("Warning: Not checking centroids!\n");*/

    // KOEN: Disable randomness of the algorithm
    srand(12345);

    // Call Felzenswalb segmentation algorithm
    int num_css;
    double* segIndices = segment_image_index(&im_feat, *c, min_size,
            &num_css, threshold, check_centroid);

    // The segmentation index image
    out[0] = mxCreateDoubleMatrix(dims[0], dims[1], mxREAL);
    double* outSegInd = mxGetPr(out[0]);

    // Keep track of minimum and maximum of each blob
    out[1] = mxCreateDoubleMatrix(num_css, 4, mxREAL);
    double* minmax = mxGetPr(out[1]);
    for (int i=0; i < 4*num_css; i++)
        minmax[i] = 0;//dims[0];
    /*for (int i= num_css; i < 2 * num_css; i++)
        minmax[i] = dims[1];*/

    // Keep track of neighbouring blobs using square matrix
    out[2] = mxCreateDoubleMatrix(dims[0], dims[1], mxREAL);
    double* nn = mxGetPr(out[2]);
    /*for (int i=0; i < height*width; i++)
        nn[i] = 0;*/
    out[3] = mxCreateDoubleMatrix(dims[0], dims[1], mxREAL);
    double* nn1 = mxGetPr(out[3]);
    /*for (int i=0; i < height*width; i++)
        nn1[i] = 0;*/
    out[4] = mxCreateDoubleMatrix(dims[0], dims[1], mxREAL);
    double* nn2 = mxGetPr(out[4]);
    /*for (int i=0; i < height*width; i++)
        nn2[i] = 0;*/
    // Copy the contents of segIndices
    // Keep track of neighbours
    // Get minimum and maximum
    // These actually comprise of the bounding boxes
    double currDouble;
    int mprev, curr, prevHori, mcurr;
    for(int x = 0; x < width; x++){
        mprev = segIndices[x * height]-1;
        for(int y=0; y < height; y++){
            //mexPrintf("x: %d y: %d\n", x, y);
            int idx = x * height + y;
            //mexPrintf("idx: %d\n", idx);
            //currDouble = segIndices[idx]; 
            //mexPrintf("currDouble: %d\n", currDouble);
            curr = segIndices[idx]; 
            //mexPrintf("curr: %d\n", curr);
            outSegInd[idx] = curr; // copy contents
            //mexPrintf("outSegInd: %f\n", outSegInd[idx]);
            mcurr = curr-1;

            // Get neighbours (vertical)
            //mexPrintf("idx: %d", curr * num_css + mprev);
            //mexPrintf(" %d\n", curr + num_css * mprev);
            //mexPrintf("mprev: %d\n", mprev);
           // nn[(mcurr) * num_css + mprev] = 1;
            //nn[(mcurr) + num_css * mprev] = 1;

            // Get horizontal neighbours
            //mexPrintf("Get horizontal neighbours\n");
            //if (x > 0){
//                prevHori = outSegInd[(x-1) * height + y] - 1;
  //              nn[mcurr * num_css + prevHori] = 1;
    //            nn[mcurr + num_css * prevHori] = 1;
      //      }

            // Keep track of min and maximum index of blobs
            //mexPrintf("Keep track of min and maximum index\n");
            //size of superpixels
            minmax[mcurr] += 1;
            //average of superpixel
            minmax[mcurr + num_css] += (double)im_feat.features[idx];
            //std of superpixel
            minmax[mcurr + 2 * num_css] += (double)(im_feat.features[idx]*im_feat.features[idx]);
            /*if (minmax[mcurr] > y)
                minmax[mcurr] = y;
            if (minmax[mcurr + num_css] > x)
                minmax[mcurr + num_css] = x;
            if (minmax[mcurr + 2 * num_css] < y)
                minmax[mcurr + 2 * num_css] = y;
            if (minmax[mcurr + 3 * num_css] < x)
                minmax[mcurr + 3 * num_css] = x;*/

            //mexPrintf("Mprev = mcurr");
            mprev = mcurr;
        }
    }

    
    for (int i=0; i < num_css; i++)        
    {
        double den=(minmax[i]<2) ? 2 : minmax[i];
        double mean =  minmax[i + num_css]/den;
        double variance =  (minmax[i + 2 * num_css]-minmax[i + num_css]*mean)/(den-1);
        minmax[i + num_css] = mean;
        minmax[i + 2 * num_css] = sqrt(variance);
    }
    for(int x = 0; x < width; x++)
    {        
        for(int y=0; y < height; y++)
        {
            int idx = x * height + y;
           // int idxlab = x * height + y;            
            curr = segIndices[idx]; 
            nn[idx] = minmax[curr-1];
            nn1[idx] = minmax[curr-1+num_css];
            nn2[idx] = minmax[curr-1+2*num_css];
        }
    }
    /*for (int i=0; i < 4 * num_css; i++)
        minmax[i] += 1;*/

    delete [] segIndices;

    return;
}


