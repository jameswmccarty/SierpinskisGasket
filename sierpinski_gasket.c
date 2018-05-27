#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <pthread.h>
#include <tiffio.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include "sierpinski.h"

int
random_bit (void)
{
  return rand () & 0x01;
}


void
gasket_init (gasket * sier)
{

  /* set default values for structure */
  /* reference sierpinski.h */

  sier->xres = XRES;
  sier->yres = YRES;
  sier->xmax = XMAX;
  sier->ymax = YMAX;
  sier->xmin = XMIN;
  sier->ymin = YMIN;
  sier->xdiv = XDIV;
  sier->ydiv = YDIV;
  sier->n = NUMV;
  sier->seed = SEED;
  sier->samples = SAMPLES;
  sier->gamma = GAMMA;
  sier->symmetry = 1;
  sier->increment = 1;
  sier->invert = 0;
  sier->file = PATH;
  sier->sup = SUPER;
  sier->palfile = NULL;
  sier->iterations = ITT;
  sier->num_threads = THREADGROUPSIZE;
  sier->R = -1;
  sier->G = -1;
  sier->B = -1;
}

/* initialize the color palette with user inputs
 * or a default state if no input is provided. */
void
pal_init (palette * pal, char *infile)
{
  FILE *palette;

  if (infile == NULL)
    {
      /* a nice light blue default */
      pal->a_r = 0.39;
      pal->a_g = 0.55;
      pal->a_b = 0.5;

      pal->b_r = 0.55;
      pal->b_g = 0.26;
      pal->b_b = 0.68;

      pal->c_r = 0.5;
      pal->c_g = 1.5;
      pal->c_b = 0.0;

      pal->d_r = 0.26;
      pal->d_g = 0.11;
      pal->d_b = 0.24;
    }
  else
    {
      if ((palette = fopen (infile, "r")) == NULL)
	{
	  printf ("Error reading input file %s.\n", infile);
	  exit (EXIT_FAILURE);
	}
      printf ("Using provided color palette.\n");
      /* WARNING -- poor checks for malformed input here. */
      assert (fscanf
	      (palette, "%lf %lf %lf\n", &(pal->a_r), &(pal->a_g),
	       &(pal->a_b)) != EOF);
      assert (fscanf
	      (palette, "%lf %lf %lf\n", &(pal->b_r), &(pal->b_g),
	       &(pal->b_b)) != EOF);
      assert (fscanf
	      (palette, "%lf %lf %lf\n", &(pal->c_r), &(pal->c_g),
	       &(pal->c_b)) != EOF);
      assert (fscanf
	      (palette, "%lf %lf %lf\n", &(pal->d_r), &(pal->d_g),
	       &(pal->d_b)) != EOF);
      (void) fclose (palette);
    }				/* end else */
}				/* end pal_init */

/* continuous color range based on normalized value t */
void
color_pxl (double t, palette * pal, double *r_out, double *g_out,
	   double *b_out)
{
  *r_out =
    255. * (pal->a_r +
	    pal->b_r * cos (M_PI * 2. * (pal->c_r * t + pal->d_r)));
  *g_out =
    255. * (pal->a_g +
	    pal->b_g * cos (M_PI * 2. * (pal->c_g * t + pal->d_g)));
  *b_out =
    255. * (pal->a_b +
	    pal->b_b * cos (M_PI * 2. * (pal->c_b * t + pal->d_b)));
}

/* allocate memory for the image buffer */
void
buffer_init (gasket * sier)
{

  int y;
  /* last minute sanity checks */
  sier->xres = sier->xres <= 0 ? XRES : sier->xres;
  sier->yres = sier->yres <= 0 ? YRES : sier->yres;

  if (sier->xmin > sier->xmax)
    {
      double t = sier->xmin;
      sier->xmin = sier->xmax;
      sier->xmax = t;
    }

  if (sier->ymin > sier->ymax)
    {
      double t = sier->ymin;
      sier->ymin = sier->ymax;
      sier->ymax = t;
    }

  sier->ranx = sier->xmax - sier->xmin;
  sier->rany = sier->ymax - sier->ymin;

  sier->xres *= sier->sup;
  sier->yres *= sier->sup;

  /* malloc new memory array for the image */
  /* then memory can be freed after being written to disk */

  printf ("Size of pixel structure is %ld bytes.\n", sizeof (pixel));
  printf ("Attempting to allocate %ld MiB of RAM.\n",
	  ((sier->yres * sier->xres * sizeof (pixel))) / (1024 * 1024));

  sier->pixels = malloc (sier->yres * sizeof (pixel *));
  if (sier->pixels == NULL)
    {
      printf ("malloc() in buffer_init () failed.\n");
      exit (EXIT_FAILURE);
    }
  sier->lock = malloc (sier->yres * sizeof (pthread_mutex_t));
  if (sier->lock == NULL)
    {
      printf ("malloc() in buffer_init () failed.\n");
      exit (EXIT_FAILURE);
    }
  for (y = 0; y < sier->yres; y++)
    {
      sier->pixels[y] = malloc (sier->xres * sizeof (pixel));
      if (sier->pixels[y] == NULL)
	{
	  printf ("malloc() failed\n");
	  exit (EXIT_FAILURE);
	}
      memset (sier->pixels[y], '\0', sier->xres * sizeof (pixel));
      if (0 != pthread_mutex_init (&(sier->lock[y]), NULL))
	{
	  printf ("Error: mutex init failed.\n");
	}
    }

  printf ("Done!\n");

}

void
print_usage ()
{
  /* print program use */

  printf ("Sierpinski's gasket program usage:\n");
  printf ("sierpinski [-options ...]\n\n");
  printf ("options include:\n");

  printf ("\t-h\t\t\tprint this screen\n");
  printf ("\t-R NUM\t\t\tseed randomizer with NUM\n");
  printf ("\t-f NAME [%s]\tfile to write\n", PATH);
  printf ("\t-S NUM>1 [1]\t\tenable rotational symmetry axis\n");
  printf ("\t-I\t\t\tinvert colors in final image\n");
  printf ("\t-x XRES [%d]\t\timage x resolution\n", XRES);
  printf ("\t-y YRES [%d]\t\timage y resolution\n", YRES);
  printf ("\t-xd XDIV [%f]\t\tx division scale\n", (float) XDIV);
  printf ("\t-yd YDIV [%f]\t\ty division scale\n", (float) YDIV);
  printf ("\t-m XMIN [%f]\tgraph x axis minimum\n", (float) XMIN);
  printf ("\t-M XMAX [%f]\tgraph x axis maximum\n", (float) XMAX);
  printf ("\t-l YMIN [%f]\tgraph y axis minimum\n", (float) YMIN);
  printf ("\t-L YMAX [%f]\tgraph y axis maximum\n", (float) YMAX);
  printf ("\t-n NUMV [%d]\t\tnumber of points around circle to start with\n",
	  NUMV);
  printf ("\t-s SAMPLES [%d]\tnumber of image samples\n", SAMPLES);
  printf ("\t-i NUM>20 [%d]\tnumber of iterations to run per sample\n", ITT);
  printf ("\t-c NUM [%d]\t\tamount to increment counter per hit\n", INC);
  printf ("\t-r 0<=NUM<=255\t\tset static RED channel value\n");
  printf ("\t-g 0<=NUM<=255\t\tset static GREEN channel value\n");
  printf ("\t-b 0<=NUM<=255\t\tset static BLUE channel value\n");
  printf ("\t-sup NUM [%d] \t\tsuper sample NUM^2 bit buckets\n", SUPER);
  printf ("\t-G NUM [%f]\tcorrectional gamma factor\n", (float) GAMMA);
  printf ("\t-p FILE\t\t\tcolor palette file\n");
  printf ("\t-T THREADS [%d]\tnumber of threads to run\n", THREADGROUPSIZE);

  fflush (stdout);
}

void
parse_args (int argc, char **argv, gasket * sier)
{
  int i = 1;
  while (i < argc)
    {
      if (!strcmp (argv[i], "-h"))
	{
	  print_usage ();
	  exit (EXIT_SUCCESS);
	}
      else if (!strcmp (argv[i], "-I"))
	{
	  sier->invert = 1;
	  i++;
	}
      else if (!strcmp (argv[i], "-R"))
	{
	  sier->seed = atoi (argv[i + 1]);
	  i += 2;
	}
      else if (!strcmp (argv[i], "-T"))
	{
	  sier->num_threads = atoi (argv[i + 1]);
	  i += 2;
	}
      else if (!strcmp (argv[i], "-f"))
	{
	  sier->file = argv[i + 1];
	  i += 2;
	}
      else if (!strcmp (argv[i], "-p"))
	{
	  sier->palfile = argv[i + 1];
	  i += 2;
	}
      else if (!strcmp (argv[i], "-x"))
	{
	  sier->xres = atoi (argv[i + 1]);
	  i += 2;
	}
      else if (!strcmp (argv[i], "-y"))
	{
	  sier->yres = atoi (argv[i + 1]);
	  i += 2;
	}
      else if (!strcmp (argv[i], "-m"))
	{
	  sier->xmin = (double) atof (argv[i + 1]);
	  i += 2;
	}
      else if (!strcmp (argv[i], "-M"))
	{
	  sier->xmax = (double) atof (argv[i + 1]);
	  i += 2;
	}
      else if (!strcmp (argv[i], "-l"))
	{
	  sier->ymin = (double) atof (argv[i + 1]);
	  i += 2;
	}
      else if (!strcmp (argv[i], "-L"))
	{
	  sier->ymax = (double) atof (argv[i + 1]);
	  i += 2;
	}
      else if (!strcmp (argv[i], "-S"))
	{
	  sier->symmetry = atoi (argv[i + 1]);
	  i += 2;
	  if (sier->symmetry <= 0)
	    sier->symmetry = 1;
	}
      else if (!strcmp (argv[i], "-n"))
	{
	  sier->n = atoi (argv[i + 1]);
	  i += 2;
	}
      else if (!strcmp (argv[i], "-c"))
	{
	  sier->increment = atoi (argv[i + 1]);
	  i += 2;
	}
      else if (!strcmp (argv[i], "-xd"))
	{
	  sier->xdiv = (double) atof (argv[i + 1]);
	  i += 2;
	}
      else if (!strcmp (argv[i], "-yd"))
	{
	  sier->ydiv = (double) atof (argv[i + 1]);
	  i += 2;
	}
      else if (!strcmp (argv[i], "-s"))
	{
	  sier->samples = atoi (argv[i + 1]);
	  i += 2;
	}
      else if (!strcmp (argv[i], "-i"))
	{
	  sier->iterations =
	    atol (argv[i + 1]) < 20 ? 1000 : atol (argv[i + 1]);
	  i += 2;
	}
      else if (!strcmp (argv[i], "-r"))
	{
	  sier->R = atoi (argv[i + 1]) % 256;
	  i += 2;
	}
      else if (!strcmp (argv[i], "-g"))
	{
	  sier->G = atoi (argv[i + 1]) % 256;
	  i += 2;
	}
      else if (!strcmp (argv[i], "-b"))
	{
	  sier->B = atoi (argv[i + 1]) % 256;
	  i += 2;
	}
      else if (!strcmp (argv[i], "-G"))
	{
	  sier->gamma =
	    (((double) atof (argv[i + 1])) ==
	     0.0) ? 1.0 : (double) atof (argv[i + 1]);
	  i += 2;
	}
      else if (!strcmp (argv[i], "-sup"))
	{
	  sier->sup = atoi (argv[i + 1]) == 0 ? 1 : atoi (argv[i + 1]);
	  i += 2;
	}
      else
	{
	  print_usage ();
	  exit (EXIT_FAILURE);
	}
    }
}

void
write_to_tiff (gasket * fractal)
{
  int row, col;
  TIFF *output;
  char *raster;
  int invert = fractal->invert;

  printf ("Opening output image...\n");
  /* Open the output image */
  if ((output = TIFFOpen (fractal->file, "w")) == NULL)
    {
      fprintf (stderr, "Could not open outgoing image.\n");
      exit (EXIT_FAILURE);
    }

  printf ("Success!\n");
  /* malloc space for the image lines */
  raster = malloc (fractal->xres * 3 * sizeof (char));
  if (raster == NULL)
    {
      printf ("malloc() failed in write_to_tiff.\n");
      exit (EXIT_FAILURE);
    }
  printf ("Space allocated.\n");
  /* Write the tiff tags to the file */

  TIFFSetField (output, TIFFTAG_IMAGEWIDTH, (*fractal).xres);
  TIFFSetField (output, TIFFTAG_IMAGELENGTH, (*fractal).yres);
  TIFFSetField (output, TIFFTAG_COMPRESSION, COMPRESSION_DEFLATE);
  TIFFSetField (output, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
  TIFFSetField (output, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
  TIFFSetField (output, TIFFTAG_BITSPERSAMPLE, 8);
  TIFFSetField (output, TIFFTAG_SAMPLESPERPIXEL, 3);

  printf ("Tags written.\n");
  for (row = 0; row < (*fractal).yres; row++)
    {
      for (col = 0; col < (*fractal).xres; col++)
	{
	  raster[(col * 3)] = invert == 1 ? ~((*fractal).pixels[row][col].r) :
	    (*fractal).pixels[row][col].r;
	  raster[(col * 3) + 1] =
	    invert ==
	    1 ? ~((*fractal).pixels[row][col].
		  g) : (*fractal).pixels[row][col].g;
	  raster[(col * 3) + 2] =
	    invert ==
	    1 ? ~((*fractal).pixels[row][col].
		  b) : (*fractal).pixels[row][col].b;
	}
      if (TIFFWriteScanline (output, raster, row, (*fractal).xres * 3) != 1)
	{
	  fprintf (stderr, "Could not write image\n");
	  exit (EXIT_FAILURE);
	}
      free (fractal->pixels[row]);
    }

  free (raster);
  /* close the file */
  TIFFClose (output);
}

point
next_point (point * start, point * target, double xdiv, double ydiv)
{
  point result;
  result.x = start->x + target->x;
  result.y = start->y + target->y;
  result.x /= xdiv;
  result.y /= ydiv;
  return result;
}


void *
render (void *fract)
{
  double newx, newy;
  int num, s;
  long int step;
  point p, k;
  gasket *sier = (gasket *) fract;
  newx = 0.0;
  newy = 0.0;

  for (num = 0; num < sier->samples; num++)
    {
      p = sier->points[rand () % (sier->n)];

      for (step = 0; step < sier->iterations; step++)
	{

	  k =
	    next_point (&p, &(sier->points[rand () % (sier->n)]), sier->xdiv,
			sier->ydiv);
	  newx = k.x;
	  newy = k.y;

	  if (step > 0)
	    {
	      unsigned int x1, y1;
	      pixel *pix;
	      double theta2, x_rot, y_rot;

	      theta2 = 0.0;

	      for (s = 0; s < sier->symmetry; s++)
		{

		  theta2 += ((2 * M_PI) / (sier->symmetry));
		  x_rot = newx * cos (theta2) - newy * sin (theta2);
		  y_rot = newx * sin (theta2) + newy * cos (theta2);

		  if (x_rot >= sier->xmin && x_rot <= sier->xmax
		      && y_rot >= sier->ymin && y_rot <= sier->ymax)
		    {
		      x1 =
			sier->xres -
			(unsigned int) (((sier->xmax - x_rot) / sier->ranx) *
					sier->xres);
		      y1 =
			sier->yres -
			(unsigned int) (((sier->ymax - y_rot) / sier->rany) *
					sier->yres);

		      if (x1 >= 0 && x1 < sier->xres && y1 >= 0
			  && y1 < sier->yres)
			{
			  pthread_mutex_lock (&(sier->lock[y1]));
			  pix = &sier->pixels[y1][x1];
			  pix->value.counter += sier->increment;
			  pthread_mutex_unlock (&(sier->lock[y1]));
			}
		    }
		}
	    }
	  p.x = newx;
	  p.y = newy;
	}
    }
  return NULL;
}

/* apply gamma color correction and log correction */
void
gamma_log (gasket * fractal)
{

  double max;
  int row, col;
  double gamma = fractal->gamma;
  double rp, gp, bp;

  max = 0.0;
  for (row = 0; row < fractal->yres; row++)
    {
      for (col = 0; col < fractal->xres; col++)
	{
	  if (fractal->pixels[row][col].value.counter != 0)
	    {
	      fractal->pixels[row][col].value.normal =
		log ((double) fractal->pixels[row][col].value.counter);
	      max = max > fractal->pixels[row][col].value.normal ? max :
		fractal->pixels[row][col].value.normal;
	    }
	  else
	    {
	      fractal->pixels[row][col].value.normal = 0.0;
	    }
	}
    }

  /* avoid divide by zero */
  if (max == 0.0)
    {
      max = 1.0;
    }

  for (row = 0; row < fractal->yres; row++)
    {
      for (col = 0; col < fractal->xres; col++)
	{
	  fractal->pixels[row][col].value.normal /= max;
	  color_pxl (fractal->pixels[row][col].value.normal,
		     &(fractal->pal), &rp, &gp, &bp);

	  fractal->pixels[row][col].r =
	    (unsigned char) ((float) (rp)) *
	    pow (fractal->pixels[row][col].value.normal, (1.0 / gamma));
	  fractal->pixels[row][col].g =
	    (unsigned char) ((float) (gp)) *
	    pow (fractal->pixels[row][col].value.normal, (1.0 / gamma));
	  fractal->pixels[row][col].b =
	    (unsigned char) ((float) (bp)) *
	    pow (fractal->pixels[row][col].value.normal, (1.0 / gamma));
	}

    }
}


/* take a fractal rendered with larger bit bucket and shrink it */
void
reduce (gasket * fractal)
{
  unsigned int R, G, B, count, y, x;
  int sx, sy;
  int old_yres;
  int sample = fractal->sup;
  pixel **reduction;
  old_yres = (*fractal).yres;
  (*fractal).xres = (*fractal).xres / sample;
  (*fractal).yres = (*fractal).yres / sample;
  reduction = malloc ((*fractal).yres * sizeof (pixel *));
  if (reduction == NULL)
    {
      printf ("malloc() failed\n");
      exit (EXIT_FAILURE);
    }
  for (y = 0; y < (*fractal).yres; y++)
    {
      reduction[y] = malloc ((*fractal).xres * sizeof (pixel));
      if (reduction[y] == NULL)
	{
	  printf ("malloc() failed\n");
	  exit (EXIT_FAILURE);
	}
    }

  /* simple grid algorithm anti-aliasing */
  /* numerical average of colors in a square region */

  for (y = 0; y < (*fractal).yres; y++)
    {
      for (x = 0; x < (*fractal).xres; x++)
	{
	  R = 0;
	  G = 0;
	  B = 0;
	  count = 0;
	  for (sy = 0; sy < sample; sy++)
	    {
	      for (sx = 0; sx < sample; sx++)
		{
		  R += (*fractal).pixels[y * sample + sy][x * sample + sx].r;
		  G += (*fractal).pixels[y * sample + sy][x * sample + sx].g;
		  B += (*fractal).pixels[y * sample + sy][x * sample + sx].b;
		  count +=
		    (*fractal).pixels[y * sample + sy][x * sample +
						       sx].value.counter;
		}
	    }
	  reduction[y][x].r = (unsigned char) (R / (sample * sample));
	  reduction[y][x].g = (unsigned char) (G / (sample * sample));
	  reduction[y][x].b = (unsigned char) (B / (sample * sample));
	  reduction[y][x].value.counter = count;
	}
    }

  /* replace pixel array with new, smaller array */

  for (y = 0; y < old_yres; y++)
    {
      free ((*fractal).pixels[y]);
    }
  free ((*fractal).pixels);
  (*fractal).pixels = reduction;
}


int
main (int argc, char **argv)
{

  gasket sier;
  int i;
  pthread_t *threads;
  double theta = 0.0;

  /* initialize our Flame Fractal */
  printf ("Initializing...\n");
  gasket_init (&sier);
  printf ("Initialized!\n");
  /* parse arguments from the command line */
  printf ("Parsing user arguments...\n");
  parse_args (argc, argv, &sier);
  pal_init (&(sier.pal), sier.palfile);
  printf ("Done!\n");
  /* seed the randomizer */
  srandom (sier.seed);
  srand48 (random ());
  /* allocate our memory buffer */
  printf ("Allocating memory...\n");
  buffer_init (&sier);

  sier.points = malloc (sier.n * sizeof (point));
  if (sier.points == NULL)
    {
      printf
	("Error: malloc() failed in fractal_init.  Tried to malloc %d * %ld bytes.\n",
	 sier.n, sizeof (point *));
      exit (EXIT_FAILURE);
    }


  for (i = 0; i < sier.n; i++)
    {
      sier.points[i].x = cos (theta);
      sier.points[i].y = sin (theta);
      /* printf("%f, %f, %f\n", (float)theta, (float)sier.points[i].x, (float)sier.points[i].y); */
      theta += ((2 * M_PI) / (sier.n));
    }

  /* correct for threads */
  if (sier.num_threads <= 0)
    {
      sier.num_threads = 1;
    }
  threads = (pthread_t *) malloc (sier.num_threads * sizeof (pthread_t));
  if (threads == NULL)
    {
      printf ("Error: malloc() failed in main.\n");
      exit (EXIT_FAILURE);
    }
  sier.samples /= sier.num_threads;
  /* render the image */
  for (i = 0; i < sier.num_threads; i++)
    {
      printf ("Spawning thread %d\n", i);
      if (0 != pthread_create (&threads[i], NULL, render, (void *) &sier))
	exit (EXIT_FAILURE);
    }

  for (i = 0; i < sier.num_threads; i++)
    {
      printf ("Joining thread %d\n", i);
      if (0 != pthread_join (threads[i], NULL))
	exit (EXIT_FAILURE);
    }
  /* gamma and log correct */
  printf ("Finalizing and writing out...\n");

  printf ("Gamma correction...\n");
  gamma_log (&sier);
  if (sier.sup > 1)
    {
      printf ("Super-sampling...\n");
      reduce (&sier);
    }

  write_to_tiff (&sier);
  /* clean up */
  free (threads);
  free (sier.lock);
  free (sier.points);
  printf ("Done!\n");
  return 0;
}
