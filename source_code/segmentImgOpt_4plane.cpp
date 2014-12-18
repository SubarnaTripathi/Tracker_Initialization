/*Including the mexFunction as a matlab interface
usage:
a = segmentImgOpt_4plane( sigm, k, min, img, OltputFolder, PPMoption);
 Input--
 sigm
 k
 min
 img: rgb 3 color channel 2d matrix
 OltputFolder : path of the  output folder
 PPMoption : 1 if storaging PPM file
 
 Output--
 a : unsorted index superpixel 2d matrix*/

#include <stdio.h>
#include <cstdlib>
#include <cstring>
#include "image.h"
#include "misc.h"
#include "pnmfile.h"
#include "filter.h"
#include "segment-graph.h"
//#include "segment-image-opt.h"
#include "mex.h"
//#include <cstdio>


/* */
#define PPM_PATH        prhs[4]
#define PPM_OPTION      prhs[5]



// random color
rgb random_rgb(){ 
  rgb c;
  double r;
  
  c.r = (uchar)rand(); //random();
  c.g = (uchar)rand(); //random();
  c.b = (uchar)rand(); //random();

  return c;
}

// dissimilarity measure between pixels
static inline float diff(image<float> *r, image<float> *g, image<float> *b,
			 int x1, int y1, int x2, int y2) {
  return sqrt(square(imRef(r, x1, y1)-imRef(r, x2, y2)) +
	      square(imRef(g, x1, y1)-imRef(g, x2, y2)) +
	      square(imRef(b, x1, y1)-imRef(b, x2, y2)));
}

//# define POS_WEIGHT 0.4    /* default : 0.6 ( 0 to 1 ) */
//# define MOTION_WEIGHT 0.3  /* default : 0.3 ( 0 to 1 ) */

//4-pane segmentation 
// dissimilarity measure between pixels
static inline float diff_4plane(image<float> *r, image<float> *g, image<float> *b, image<float> *m,
			 int x1, int y1, int x2, int y2) 
{
  float val =  0.4*(sqrt(square(imRef(r, x1, y1)-imRef(r, x2, y2)) +
	      square(imRef(g, x1, y1)-imRef(g, x2, y2)) +
	      square(imRef(b, x1, y1)-imRef(b, x2, y2)))
	      + 0.2*square(imRef(m, x1, y1)-imRef(m, x2, y2))
		  + 0.5*(square(x1-x2) + square(y1-y2))
		  );

  return val;
}


/*
 * Segment an image
 *
 * Returns a color image representing the segmentation.
 *
 * im: image to segment.
 * sigma: to smooth the image.
 * c: constant for treshold function.
 * min_size: minimum component size (enforced by post-processing stage).
 * num_ccs: number of connected components in the segmentation.
 */
image<rgb> * segment_image( unsigned char *Imptr, double SegOut[], int height, int width, int depth, float sigma, float c, int min_size, double *PpmOption )
{
  printf("within function segment_image\n" );

  image<float> *r = new image<float>(width, height);
  image<float> *g = new image<float>(width, height);
  image<float> *b = new image<float>(width, height);
  //4-pane segmentation //subarna
  image<float> *m = new image<float>(width, height);


  // smooth each color channel 
  int Area = width* height;
  for (int x = 0; x < width; x++) 
  {
    int HightTimesX = height*x;
    for (int y = 0; y < height; y++) 
	{
      int temp = HightTimesX+y;
      imRef(r, x, y) = Imptr[temp];

	  temp += Area;
      imRef(g, x, y) = Imptr[temp];

	  temp += Area;
      imRef(b, x, y) = Imptr[temp];
	  
	  /**** add another field for motion mask ****/
	  //4-plane segmentation //subarna
	  if ( 4 == depth)
	  {
		temp += Area;
		imRef(m, x, y) = Imptr[temp];
	  }
    }
  }

  
	image<float> *smooth_r = smooth(r, sigma);
	image<float> *smooth_g = smooth(g, sigma);
	image<float> *smooth_b = smooth(b, sigma);

	//4-pane segmentation //subarna
	image<float> *smooth_m = smooth(m, sigma);
  


  delete r;
  delete g;
  delete b;
 
  //4-pane segmentation //subarna
  delete m;

 
  printf("\nstart building graph\n" );

  // build graph
  edge *edges = new edge[width*height*4];
  int num = 0;
  for (int y = 0; y < height; y++) 
  {
    int YtimesWidth = y * width;
    for (int x = 0; x < width; x++) 
	{
      if (x < width-1) 
	  {
		edges[num].a = YtimesWidth + x;
		edges[num].b = YtimesWidth + (x+1);

		// 4_PLANE
		edges[num].w = diff_4plane(smooth_r, smooth_g, smooth_b, smooth_m, x, y, x+1, y);
		//edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x+1, y);

		num++;
      }

      if (y < height-1) 
	  {
		edges[num].a = YtimesWidth + x;
		edges[num].b = YtimesWidth + width + x;
		
		// 4_PLANE
		edges[num].w = diff_4plane(smooth_r, smooth_g, smooth_b, smooth_m, x, y, x, y+1);		
		//edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x, y+1);
		
		num++;
      }

      if ((x < width-1) && (y < height-1)) 
	  {
		edges[num].a = YtimesWidth + x;
		edges[num].b = YtimesWidth + width + (x+1);

		// 4_PLANE
		edges[num].w = diff_4plane(smooth_r, smooth_g, smooth_b, smooth_m, x, y, x+1, y+1);		
		//edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x+1, y+1);

		num++;
      }

      if ((x < width-1) && (y > 0)) 
	  {
		edges[num].a = YtimesWidth + x;
		edges[num].b = YtimesWidth - width + (x+1);

		// 4_PLANE
		edges[num].w = diff_4plane(smooth_r, smooth_g, smooth_b, smooth_m, x, y, x+1, y-1);	    
		//edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x+1, y-1);

		num++;
      }
    }
  }
  delete smooth_r;
  delete smooth_g;
  delete smooth_b;
  // 4_PLANE
  delete smooth_m;


  printf("start segmenting graph\n" );


  // segment
  universe *u = segment_graph(width*height, num, edges, c);

  printf("post processing of small components started\n" );
  
  // post process small components
  for (int i = 0; i < num; i++) 
  {
    int a = u->find(edges[i].a);
    int b = u->find(edges[i].b);
    if ((a != b) && ((u->size(a) < min_size) || (u->size(b) < min_size)))
      u->join(a, b);
  }
  delete [] edges;

  image<rgb> *output = new image<rgb>(width, height);
  rgb *colors = new rgb[width*height];
/*  image<int> *output = new image<int>(width, height);*/

  // pick random colors for each component
 if (*PpmOption)
 {
  for (int i = 0; i < width*height; i++)
    colors[i] = random_rgb();
 }
 

  for (int y = 0; y < height-1; y++) 
  {
    int YtimesWidth = y * width;
    int temp = y;
    for (int x = 0; x < width-1; x++) 
	{
      int comp = u->find(YtimesWidth + x);
      if (*PpmOption)
	  {
         imRef(output, x, y) = colors[comp];
      }
      SegOut[ temp] = (double) comp;
      temp += height;
    }
  }  

  delete u;

  return output;
/*  return;*/
}




void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
  if(nrhs!=6) 
  {
    mexErrMsgTxt("usage: segmentedImage = SEGMENT(sigma, k, min, inputimage, OltputFolder, PPMoption)");
  } else if(nlhs>1) 
  {
    mexErrMsgTxt("Too many output arguments");
  }
  if (!mxIsChar( PPM_PATH) || mxIsComplex( PPM_PATH) ) 
  {
        mexErrMsgTxt("WRL path is not a string.");
  }

  //const int *dims;  
  int height, width, ndim;
  char *PpmFolder;
  float sigma = (float)*mxGetPr(prhs[0]);
  float k = (float)*mxGetPr(prhs[1]);
  int min_size = (int)*mxGetPr(prhs[2]);
  unsigned char *Imptr;
  double *SegOut;
  PpmFolder = mxArrayToString(PPM_PATH);
  double *PpmOption = mxGetPr( PPM_OPTION);
/*  float *Imptr;*/
  unsigned short int *output;
  
  height=(mxGetDimensions(prhs[3]))[0];
  width=(mxGetDimensions(prhs[3]))[1];
  ndim=(mxGetDimensions(prhs[3]))[2];

  printf("height= %d, width = %d, depth = %d\n", height, width, ndim);

  Imptr = (unsigned char*) mxGetPr(prhs[3]);

  plhs[0] = mxCreateDoubleMatrix ( height, width, mxREAL);
  SegOut =  mxGetPr(plhs[0]);

  printf("starting image_segmentation ......\n");

  image<rgb> *seg = segment_image( Imptr, SegOut, height, width, ndim, sigma, k, min_size, PpmOption);

  printf("image segmentation done......\n");

  if (*PpmOption){
      savePPM( seg, PpmFolder);
  }

  /*delete [] r;
  delete [] g;
  delete [] b;*/
  return;
}


