
#include "bitmap.h"

#include <pthread.h>
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <sys/time.h>

#define NUM_THREADS 2

int iteration_to_color( int i, int max );
int iterations_at_point( double x, double y, int max );
void compute_image( struct bitmap *bm, double xmin, double xmax, double ymin, double ymax, int max, int hstart, int hend );

const char *outfile = "mandel.bmp";
double xcenter = 0;
double ycenter = 0;
double scale = 4;
int    image_width = 500;
int    image_height = 500;
int    max = 1000;

void show_help()
{
	printf("Use: mandel [options]\n");
	printf("Where options are:\n");
	printf("-m <max>    The maximum number of iterations per point. (default=1000)\n");
	printf("-x <coord>  X coordinate of image center point. (default=0)\n");
	printf("-y <coord>  Y coordinate of image center point. (default=0)\n");
	printf("-s <scale>  Scale of the image in Mandlebrot coordinates. (default=4)\n");
	printf("-W <pixels> Width of the image in pixels. (default=500)\n");
	printf("-H <pixels> Height of the image in pixels. (default=500)\n");
	printf("-o <file>   Set output file. (default=mandel.bmp)\n");
	printf("-h          Show this help text.\n");
	printf("\nSome examples are:\n");
	printf("mandel -x -0.5 -y -0.5 -s 0.2\n");
	printf("mandel -x -.38 -y -.665 -s .05 -m 100\n");
	printf("mandel -x 0.286932 -y 0.014287 -s .0005 -m 1000\n\n");
}

typedef struct parameters
{
  int image_height_start;
  int image_height_end;
  int thread_id;
  int xcenter;
  int ycenter;
  int scale;
  int max;
  struct bitmap *bm;
}thrstruct;

void * run_me( void * arg )
{
 
  image_height = * (int*)arg;
 
  thrstruct* paramt = (thrstruct *) arg;
  
  int hstart = paramt->image_height_start;
  int hend = paramt->image_height_end;
  int xcenter = paramt->xcenter;
  int scale = paramt->scale;
  int ycenter = paramt->ycenter;
  int max = paramt->max;
  struct bitmap *bm = paramt->bm;
 
	// Compute the Mandelbrot image
	compute_image(bm,xcenter-scale,xcenter+scale,ycenter-scale,ycenter+scale,max,hstart,hend);

	// Save the image in the stated file.
	if(!bitmap_save(bm,outfile)) {
		fprintf(stderr,"mandel: couldn't write to %s: %s\n",outfile,strerror(errno));
		
	}

  return NULL;
}

int main( int argc, char *argv[] )
{
	char c;

	// These are the default configuration values used
	// if no command line arguments are given.

	
        
    pthread_t tid[NUM_THREADS];
    struct parameters params[NUM_THREADS];
	
    struct timeval begin_time;
           struct timeval end_time;

           gettimeofday( &begin_time, NULL );
	// For each command line argument given,
	// override the appropriate configuration value.

	while((c = getopt(argc,argv,"x:y:s:W:H:m:o:h"))!=-1) {
		switch(c) {
			case 'x':
				xcenter = atof(optarg);
				break;
			case 'y':
				ycenter = atof(optarg);
				break;
			case 's':
				scale = atof(optarg);
				break;
			case 'W':
				image_width = atoi(optarg);
				break;
			case 'H':
				image_height = atoi(optarg);
				break;
			case 'm':
				max = atoi(optarg);
				break;
			case 'o':
				outfile = optarg;
				break;
			case 'h':
				show_help();
				exit(1);
				break;
		}
	}

	// Display the configuration of the image.
	printf("mandel: x=%lf y=%lf scale=%lf max=%d outfile=%s\n",xcenter,ycenter,scale,max,outfile);
        int i;
        
        
        // Create a bitmap of the appropriate size.
	struct bitmap *bm = bitmap_create(image_width,image_height);

        
	// Fill it with a dark blue, for debugging
	bitmap_reset(bm,MAKE_RGBA(0,0,255,0));


        for( i = 0; i < NUM_THREADS; i++ )
          { 
           
            params[i].thread_id = i;
            params[i].xcenter = xcenter;
            params[i].ycenter = ycenter;
            params[i].scale = scale;
            params[i].max = max;
            params[i].bm = bm;
             if (i==0)
             {
              params[i].image_height_start = 0;
              params[i].image_height_end = 99;   
             }
             else
             {
              params[i].image_height_start = params[i-1].image_height_end + 1;
              params[i].image_height_end = params[i-1].image_height_end + 99;   
                 
             }
            
            pthread_create( &tid[i], NULL, run_me, (void *) &params[i] );
            
          }
	
        
        
       
        for( i = 0; i < NUM_THREADS; i++ )
        { 
            pthread_join( tid[i], NULL );

        }  
		
	    gettimeofday( &end_time, NULL );

	           long time_to_execute = ( end_time.tv_sec * 1000000 + end_time.tv_usec ) -
	                                  ( begin_time.tv_sec * 1000000 + begin_time.tv_usec );

	           printf("This code took %ld microseconds to execute\n", time_to_execute);
        
       
	return 0;
}

/*
Compute an entire Mandelbrot image, writing each point to the given bitmap.
Scale the image to the range (xmin-xmax,ymin-ymax), limiting iterations to "max"
*/

void compute_image( struct bitmap *bm, double xmin, double xmax, double ymin, double ymax, int max, int hstart, int hend )
{
	int i,j;

	int width = bitmap_width(bm);
	int height = bitmap_height(bm);

	// For every pixel in the image...

	for(j=hstart;j<hend;j++) {

		for(i=0;i<width;i++) {

			// Determine the point in x,y space for that pixel.
			double x = xmin + i*(xmax-xmin)/width;
			double y = ymin + j*(ymax-ymin)/height;

			// Compute the iterations at that point.
			int iters = iterations_at_point(x,y,max);

			// Set the pixel in the bitmap.
			bitmap_set(bm,i,j,iters);
		}
	}
}

/*
Return the number of iterations at point x, y
in the Mandelbrot space, up to a maximum of max.
*/

int iterations_at_point( double x, double y, int max )
{
	double x0 = x;
	double y0 = y;

	int iter = 0;

	while( (x*x + y*y <= 4) && iter < max ) {

		double xt = x*x - y*y + x0;
		double yt = 2*x*y + y0;

		x = xt;
		y = yt;

		iter++;
	}

	return iteration_to_color(iter,max);
}

/*
Convert a iteration number to an RGBA color.
Here, we just scale to gray with a maximum of imax.
Modify this function to make more interesting colors.
*/

int iteration_to_color( int i, int max )
{
	int gray = 255*i/max;
	return MAKE_RGBA(gray,gray,gray,0);
}