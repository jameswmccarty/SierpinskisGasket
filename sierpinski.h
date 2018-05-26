/* default values */
#define XRES 1000
#define YRES 1000
#define XMIN -1.5
#define XMAX 1.5
#define YMIN -1.5
#define YMAX 1.5
#define XDIV 2.0
#define YDIV 2.0
#define SAMPLES 20000
#define ITT 1000
#define SUPER 1
#define GAMMA 1.0
#define SEED 1
#define NUMV 3
#define PATH "/tmp/gasket.tif"
#define THREADGROUPSIZE 16
#define INTERVAL 0.5

/* structs */

typedef struct
{
  double x;
  double y;
} point;

typedef union
{
  unsigned int counter;		/* number of times pixel has been incremented */
  float normal;			/* normalized value at pixel */
} hitcounter;

typedef struct
{
  double a_r, a_g, a_b;
  double b_r, b_g, b_b;
  double c_r, c_g, c_b;
  double d_r, d_g, d_b;
} palette;

typedef struct
{
  hitcounter value;
  unsigned char r, g, b;	/* color content of a pixel: RGB channels */
} pixel;

typedef struct
{
  int xres, yres;		/* x and y resolution of image */
  double xmin, ymin, xmax, ymax;	/* axis bounds */
  double ranx, rany;		/* numerical range of x/y axis */
  double xdiv, ydiv;		/* divisor for point mid dist (default of 2.0) */
  int R, G, B;			/* fixed color channels */
  int sup;			/* super sample value  */
  int samples;			/* number of flame samples */
  long int iterations;		/* number of iterations per sample */
  int invert;			/* use inverse colors? 0 false, else true */
  int symmetry;			/* use symmetrical rotation axis? set to greater than 1 */
  int seed;			/* random seed */
  int num_threads;		/* number of threads to use in render */
  double gamma;			/* gamma correction factor */
  char *file;			/* output file path */
  pixel **pixels;		/* image buffer */
  char *palfile;		/* color palette coefficient variables */
  pthread_mutex_t *lock;	/* lock so that only one section of memory buffer being written at a time */
  int n;			/* number of points considered in chaos game */
  palette pal;
  point *points;		/* all points being considered in chaos game */
} gasket;

/* prototypes */
void print_usage ();
void reduce (gasket * fractal);
void gasket_init (gasket * sier);
void gamma_log (gasket * fractal);
void parse_args (int argc, char **argv, gasket * fractal);
void buffer_init (gasket * fractal);
void write_to_tiff (gasket * fractal);
void *render (void *fract);
