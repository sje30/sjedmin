#include <R.h>
#include <S.h>			/* for seed_in, seed_out */

/* Use the definitions by Venables & Ripley. */
#  define RANDIN  seed_in((long *)NULL)
#  define RANDOUT seed_out((long *)NULL)
#  define UNIF unif_rand()

/* Sfloat already defined in R.h:55
 * typedef double Sfloat;
 */


/* TODO items:
 * Maybe pass all parameters in one vector, rather than using
 * many different arguments.  makes the C/R interface shorter.
 *
 * Birth & death routine: return the dmin value used at each iteration?
 * To do this, I need to know in advance how many iterations to run
 * for -- maybe have default mm as 40 * npts in R code?
 */


#define MAXREJECTS 999999
/* This is the maximum number of cells that we can reject before
 * abandoning this simulation.  (This is the total through the simulation,
 * not the total per cell.)
 */

#define MAXDISTANCE 9999999	/* some large number used in nnd_n. */
void h_amac(Sfloat *t, Sfloat *xd1, Sfloat *p, Sfloat *b, Sfloat *h);

void nnd_n(Sfloat *xs, Sfloat *ys, int n, Sfloat a, Sfloat b,
	   int *idx, Sfloat *min);

void nnd(Sfloat *xs, Sfloat *ys, int n, Sfloat a, Sfloat b, int ignore,
	 int *idx, Sfloat *min);

void nnd_n_3d(Sfloat *xs, Sfloat *ys, Sfloat *zs,
	      int n, Sfloat a, Sfloat b, Sfloat c,
	      int *idx, Sfloat *min);

void nnd_3d(Sfloat *xs, Sfloat *ys, Sfloat *zs, int n,
	    Sfloat a, Sfloat b, Sfloat c, int ignore,
	    int *idx, Sfloat *min);

void bdmin_check(Sfloat *xpts, Sfloat *ypts, int n1, int n2,
		 int i, Sfloat d, Sfloat d12, int *okay, int *id);

void hlookup(Sfloat *h, Sfloat *ds, int *n, Sfloat *d, Sfloat *res);
     
void pairwise_amac(Sfloat *wid, Sfloat *ht, int *numcells,
		   Sfloat *xd1, Sfloat *p, Sfloat *b,
		   Sfloat *xpts, Sfloat *ypts)

{
  /* Create NUMCELLS cells distributed in an area of size wid*ht.
   * Minimum distance between cells should be at least dmin.
   * On exit, we return an array of cells stored in (xs, ys).
   * NJRECTS stores the number of rejected cells.
   */

  /* This routine adapted from Brian Ripley's code for generating
   * Strauss process, from his "spatial" library.  This code
   * implements the pairwise amacrine function suggested by Diggle and
   * Gratton(1984) for modelling amacrine distribution.  All the birth
   * and death processes in this file follow the same structure as
   * this.
   */
  
  int i,j, mm, id;
  Sfloat x1, y1, h, dist, u, h1;
  int n = *numcells;

  RANDIN;

  /* Create initial set of points. */
  for (i=0; i<n; i++) {
    xpts[i] = UNIF * (*wid); ypts[i] = UNIF * (*ht);
  }

  mm = 4 * n;			/* Number of repeats */
  mm *= 10;			/* if starting from Poisson distribution. */
  for (i=0; i<= mm; i++) {
    id = (int)(n*UNIF);	/* pick one point at random. */
    /*Rprintf("select new point for unit %d\n", id);*/
    xpts[id] = xpts[0]; ypts[id] = ypts[0];

    do {

      /* generate a new point. */
      xpts[0] = UNIF * (*wid); ypts[0] = UNIF * (*ht);
      u = UNIF;
      h = 1.0;

      for(j=1; j<n; j++) {
	x1 = xpts[j] - xpts[0]; y1 = ypts[j] - ypts[0];
	dist = (sqrt((x1*x1) + (y1*y1)));
	h_amac(&dist, xd1, p, b, &h1);
	/*printf("xxx valu %lf\n", h1);*/
	h *= h1;
	/*Rprintf("t %lf %lf %lf\n", dist, h, h1); */
	/*h = 1;*/
      }
    } while ( h < u);
  }

  RANDOUT;
}

void h_amac(Sfloat *t, Sfloat *xd1, Sfloat *p, Sfloat *b, Sfloat *h)
{
  /* Pairwise interaction function for amacrines. This is a helper
   * function for pairwise_amac, above. */
  if (*t < *xd1)
    *h = 0.0;
  else if (*t <= *p)
    *h = ( pow( ((*t-*xd1)/(*p-*xd1)), *b));
  else
    *h = 1.0;
  /*printf("ret valu %lf\n", *h);*/
}


void dminlulbd(Sfloat *pw, int *numcells,
	       Sfloat *pdmin, Sfloat *psd,
	       Sfloat *plower, Sfloat *pupper,
	       Sfloat *xpts, Sfloat *ypts,
	       Sfloat *dmins, int *nrejects)
{
  /* Birth and death version of dmin rule. */
  /* Create NUMCELLS cells distributed in an area of size wid*ht.
   * Minimum distance between cells should be at least dmin.
   * On exit, we return an array of cells stored in (xs, ys).
   * NJRECTS stores the number of rejected cells.
   */
  
  int i,j, mm, id, okay, generate_r;
  Sfloat x1, y1, dist2, u, r, rr, lower, upper;
  Sfloat xmin, xmax, ymin, ymax, wid, ht;
  int n = *numcells;

  RANDIN;
  xmin = pw[0]; xmax = pw[1]; ymin = pw[2]; ymax = pw[3];
  wid = xmax - xmin; ht = ymax - ymin;
  lower = *plower; upper = *pupper;

  /* Create initial set of points. */
  for (i=0; i<n; i++) {
    xpts[i] = xmin + (UNIF * wid); ypts[i] = ymin + (UNIF * ht);
  }

  mm = 4 * n;			/* Number of repeats */
  mm *= 10;			/* if starting from Poisson distribution. */
  for (i=0; i<= mm; i++) {
    id = (int)(n*UNIF);	/* pick one point at random. */
    if ( 0 && (i%100)==0)
      Rprintf("it %d:select new point for unit %d\n", i,id);
    xpts[id] = xpts[0]; ypts[id] = ypts[0];

    do {
      
      /* generate a new point. */
      xpts[0] = xmin + (UNIF * wid); ypts[0] = ymin + (UNIF * ht);
/*       Rprintf("generate new trial at %f %f\n", xpts[0], ypts[0]); */
      u = UNIF;
      okay = 1;
      
      /* genereate a new dmin value. */
      generate_r = 1;
      while (generate_r) {
	r = (*psd * norm_rand()) + *pdmin;
	if ( (r > lower ) && ((upper <0) || (r < upper)))
	  generate_r = 0; /* okay r value */
      }
      /*Rprintf("choose %.3f\n", r);*/
      rr = r*r;
      for(j=1; j<n; j++) {
	x1 = xpts[j] - xpts[0]; y1 = ypts[j] - ypts[0];
	dist2 = (x1*x1) + (y1*y1);
	if (dist2 < rr) {
	  okay = 0;
	  break; /* quit for loop and try again. */
	}
      }
    } while ( !okay);

    /* dmins and nrejects not returned here. */
    dmins[id] = r;		/* store r value. */
  }

  RANDOUT;
}

void dminsd(Sfloat *pwid, Sfloat *pht, int *pnumcells,
	    Sfloat *pdmin, Sfloat *psd,
	    Sfloat *xpts, Sfloat *ypts,
	    Sfloat *dmins, int *nrejects)

{
  /* Create NUMCELLS cells distributed in an area of size wid*ht.
   * Minimum distance between cells should be at least dmin.
   * On exit, we return an array of cells stored in (xs, ys).
   * NJRECTS stores the number of rejected cells.
   */
  

  int num_rejects = 0, this_cell_rejects = 0;
  int looking;
  int i;
  Sfloat x,y;
  Sfloat min; int idx;
  Sfloat this_dmin;


  RANDIN;
  for (i=0; i<*pnumcells; i++) {

    this_dmin = *pdmin + (*psd * norm_rand());
    /*Rprintf("this val %lf\n", this_dmin); */
    dmins[i] = this_dmin;
    looking = 1; this_cell_rejects = 0;
    while (looking) {
      x = UNIF * (*pwid); y = UNIF * (*pht);
      
      nnd_n( xpts, ypts, i, x, y, &idx, &min);
      /*min = 1000; */
      if ( min < this_dmin) {
	/* reject cell and try another. */
	num_rejects++;
	this_cell_rejects++;
	if (num_rejects > MAXREJECTS) {
	  /* If we cannot fit more cells in, return and indicate error
	   * by setting first x coordinate to negative value.
	   */
	  printf("dminsd error: num_rejects is too high (%s:%d)\n",
		 __FILE__, __LINE__);
	  xpts[0] = -1;
	  RANDOUT;
	  return;
	}
      }
      else {
	looking = 0;
      }
    }
    /* Accept this cell. */
    xpts[i] = x; ypts[i] = y;
    nrejects[i] = this_cell_rejects;
  }

  Rprintf("#rejects %d\tdminsd packing density %.3f\n", num_rejects,
	  ((*pnumcells * PI * *pdmin * *pdmin)/(4*  *pwid * *pht)));
  
  RANDOUT;
}


void dminl(Sfloat *pwid, Sfloat *pht, int *pnumcells,
	   Sfloat *pdmin, Sfloat *psd,
	   int *quiet,
	   Sfloat *xpts, Sfloat *ypts, 
	   Sfloat *dmins, int *nrejects)

{
  /* Create NUMCELLS cells distributed in an area of size wid*ht.
   * Minimum distance between cells should be at least dmin.  On exit,
   * we return an array of cells stored in (xs, ys).  NJRECTS stores
   * the number of rejected cells.  This is Lucia's version, where a
   * new dmin value is created every time, regardless of whether the
   * trial cell was accepted or not.
   */
  

  int num_rejects = 0, this_cell_rejects = 0;
  int looking;
  int i;
  Sfloat x,y;
  Sfloat min; int idx;
  Sfloat this_dmin;

  Sfloat min_dmin = 0;		/* Just constrain dmin to be >0 */

  RANDIN;
  /* Just an esitmate since we don't know the actual dmin values */

  for (i=0; i<*pnumcells; i++) {

    looking = 1; this_cell_rejects = 0;
    while (looking) {

      this_dmin = *pdmin + (*psd * norm_rand());
      if (this_dmin < min_dmin) this_dmin = min_dmin;
      /*Rprintf("this val %lf\n", this_dmin); */

      x = UNIF * (*pwid); y = UNIF * (*pht);
      
      nnd_n( xpts, ypts, i, x, y, &idx, &min);

      if ( min < this_dmin) {
	/* reject cell and try another. */
	num_rejects++;
	this_cell_rejects++;
	if (num_rejects > MAXREJECTS) {
	  /* If we cannot fit more cells in, return and indicate error
	   * by setting first x coordinate to negative value.
	   */
	  printf("dminl error: num_rejects is too high (%s:%d)\n",
		 __FILE__, __LINE__);
	  xpts[0] = -1;
	  RANDOUT;
	  return;
	}
      }
      else {
	looking = 0;
      }
    }
    /* Accept this cell. */
    dmins[i] = this_dmin;
    xpts[i] = x; ypts[i] = y;
    nrejects[i] = this_cell_rejects;
  }

  if(! *quiet) 
    Rprintf("#rejects %d\tdminl packing density %.3f\n", num_rejects,
	    ((*pnumcells * PI * *pdmin * *pdmin)/(4 * *pwid * *pht)));

  RANDOUT;
}


void dminlul(Sfloat *pw, int *pnumcells,
	     Sfloat *pdmin, Sfloat *psd,
	     Sfloat *plower, Sfloat *pupper,
	     int *quiet,
	     Sfloat *xpts, Sfloat *ypts, 
	     Sfloat *dmins, int *nrejects)

{
  /* Create NUMCELLS cells distributed in an area of size wid*ht.
   * Minimum distance between cells should be at least dmin.  On exit,
   * we return an array of cells stored in (xs, ys).  NJRECTS stores
   * the number of rejected cells.  This is Lucia's version, where a
   * new dmin value is created every time, regardless of whether the
   * trial cell was accepted or not.  Also provide arguments LOWER
   * and UPPER: dmin values outside this range are rejected.
   */
  

  int num_rejects = 0, this_cell_rejects = 0;
  int looking, generate_dmin;
  int i;
  Sfloat x,y;
  Sfloat min; int idx;
  Sfloat this_dmin;
  Sfloat lower, upper;
  Sfloat xmin, xmax, ymin, ymax, wid, ht;


  RANDIN;
  xmin = pw[0]; xmax = pw[1]; ymin = pw[2]; ymax = pw[3];
  wid = xmax - xmin; ht = ymax - ymin;
  lower = *plower; upper = *pupper;
  for (i=0; i<*pnumcells; i++) {

    looking = 1; this_cell_rejects = 0;
    while (looking) {

      generate_dmin = 1;
      while (generate_dmin) {
	this_dmin = *pdmin + (*psd * norm_rand());
	if ( (this_dmin > lower) &&
	     ((upper <0) || (this_dmin < upper)))
	  generate_dmin = 0;
	/*else printf("dminlul: dmin %f outside range\n", this_dmin);*/
      }

      /*Rprintf("this val %lf\n", this_dmin); */
      x = xmin + (UNIF * wid); y = ymin + (UNIF * ht);

      
      nnd_n( xpts, ypts, i, x, y, &idx, &min);
      /*min = 1000; */
      if ( (min < this_dmin) ) {
	/* || (min > upper)*/
	/* Cannot test if the min distance is too big... this won't
	 * work for early cells!  i.e. when putting the first cell in,
	 * the upper distance will be the value returned by nnd_n,
	 * ie.s. some huge number. */

	 
	/* reject cell and try another. */
	num_rejects++;
	this_cell_rejects++;
	if (num_rejects > MAXREJECTS) {
	  /* If we cannot fit more cells in, return and indicate error
	   * by setting first x coordinate to negative value.
	   */
	  printf("dminlul error: num_rejects is too high (%s:%d)\n",
		 __FILE__, __LINE__);
	  xpts[0] = -1;
	  RANDOUT;
	  return;
	}
      }
      else {
	looking = 0;
      }
    }
    /* Accept this cell. */
    dmins[i] = this_dmin;
    xpts[i] = x; ypts[i] = y;
    nrejects[i] = this_cell_rejects;
  }

  if (!*quiet)
    Rprintf("#rejects %d\tdminlul packing density %.3f\n", num_rejects,
	    ((*pnumcells * PI * *pdmin * *pdmin)/( 4 * wid * ht)));

  RANDOUT;
}

void dminlulfix2(Sfloat *pw, int *pnumcells,
		 Sfloat *pdmin, Sfloat *psd,
		 Sfloat *plower, Sfloat *pupper,
		 int *quiet,
		 Sfloat *xpts, Sfloat *ypts, 
		 Sfloat *dmins, int *nrejects,
		 Sfloat *x2, Sfloat *y2, int *n2, Sfloat *d12)

{
  /* Create NUMCELLS cells distributed in an area of size wid*ht.
   * Minimum distance between cells should be at least dmin.  On exit,
   * we return an array of cells stored in (xs, ys).  NJRECTS stores
   * the number of rejected cells.  This is Lucia's version, where a
   * new dmin value is created every time, regardless of whether the
   * trial cell was accepted or not.  Also provide arguments LOWER
   * and UPPER: dmin values outside this range are rejected.
   */
  

  int num_rejects = 0, this_cell_rejects = 0;
  int looking, generate_dmin;
  int i;
  Sfloat x,y;
  Sfloat min; int idx;
  Sfloat this_dmin;
  Sfloat lower, upper;

  Sfloat xmin, xmax, ymin, ymax, wid, ht;


  RANDIN;

  xmin = pw[0]; xmax = pw[1]; ymin = pw[2]; ymax = pw[3];
  /*Rprintf("%f %f %f %f\n", xmin, xmax, ymin, ymax);*/
  wid = xmax - xmin; ht = ymax - ymin;
  lower = *plower; upper = *pupper;
  for (i=0; i<*pnumcells; i++) {

    looking = 1; this_cell_rejects = 0;
    while (looking) {

      generate_dmin = 1;
      while (generate_dmin) {
	this_dmin = *pdmin + (*psd * norm_rand());
	if ( (this_dmin > lower) &&
	     ((upper <0) || (this_dmin < upper)))
	  generate_dmin = 0;
	/*else printf("dminlul: dmin %f outside range\n", this_dmin);*/
      }

      /*Rprintf("this val %lf\n", this_dmin); */

      x = xmin + (UNIF * wid); y = ymin + (UNIF * ht);
      
      nnd_n( xpts, ypts, i, x, y, &idx, &min);
      /*min = 1000; */
      if ( (min < this_dmin) ) {
	/* || (min > upper)*/
	/* Cannot test if the min distance is too big... this won't
	 * work for early cells!  i.e. when putting the first cell in,
	 * the upper distance will be the value returned by nnd_n,
	 * ie.s. some huge number. */

	 
	/* reject cell and try another. */
	num_rejects++;
	this_cell_rejects++;
	if (num_rejects > MAXREJECTS) {
	  /* If we cannot fit more cells in, return and indicate error
	   * by setting first x coordinate to negative value.
	   */
	  printf("dminlul error: num_rejects is too high (%s:%d)\n",
		 __FILE__, __LINE__);
	  xpts[0] = -1;
	  RANDOUT;
	  return;
	}
      }
      else {
	/* min dist for this population okay, but what about 2nd population? */
	nnd_n( x2, y2, *n2, x, y, &idx, &min);
	if ( min > *d12)
	  looking = 0;		/* set if we are okay for 2nd popn */
	else {
	  nrejects[*pnumcells]++; /* keep track of rejects */
	  ;
	  /*Rprintf("point %f %f too close to 2nd pop cell %d\n", x, y, idx);*/
	}
      }	

    }
    /* Accept this cell. */
    dmins[i] = this_dmin;
    xpts[i] = x; ypts[i] = y;
    nrejects[i] = this_cell_rejects;
  }

  if (!*quiet)
    Rprintf("#rejects %d\tdminlul packing density %.3f\n", num_rejects,
	    ((*pnumcells * PI * *pdmin * *pdmin)/( 4 * wid * ht)));

  RANDOUT;
}

void bdmin_bd(Sfloat *pw, int *pn1, int *pn2,
	      Sfloat *params,
	      int *pnsweeps, int *pverbose,
	      Sfloat *xpts, Sfloat *ypts, 
	      int *nrejects)
{
  /* Bivariate DMIN model, with birth&death algorithm.
   * W is the window of the field, [xmin xmax ymin ymax]
   * N1 is number of type 1 cells, N2 number of type 2 cells.
   * PARAMS is a vector storing the parameters for the exclusion zones.
   * NSWEEPS is the number of sweeps to perform.
   * *PVERBOSE is 1 => print verbose output.
   */

  int num_rejects = 0, this_cell_rejects;
  int looking, generate_dmin;
  int i, sweep, constraint, id;
  int n, n1, n2;
  Sfloat x,y;
  Sfloat min; int idx;
  Sfloat this_dmin;
  Sfloat lower, upper;

  Sfloat xmin, xmax, ymin, ymax, wid, ht;
  Sfloat d1, sd1, d2, sd2, d12;
  Sfloat dmin_mu, dmin_sd;
  RANDIN;

  /* Extract relevant exclusion zone parameters. */
  d1  = params[0]; sd1 = params[1];
  d2  = params[2]; sd2 = params[3];
  d12 = params[4];
  lower = params[5]; upper = params[6];
  
  sweep = *pnsweeps;
  n1 = *pn1; n2 = *pn2; n = n1 + n2;
  
  xmin = pw[0]; xmax = pw[1]; ymin = pw[2]; ymax = pw[3];
  if (0 && *pverbose) {
    Rprintf("field %f %f %f %f\n", xmin, xmax, ymin, ymax);
    Rprintf("n1 %d n2 %d\n", n1, n2);
    Rprintf("d1 %f +/- %f d2 %f +/- %f d12 %f\n",
	    d1, sd1, d2, sd2, d12);
  }
  wid = xmax - xmin; ht = ymax - ymin;

  while( sweep-- > 0) {
    /*Perform one complete sweep, updating positions of cells.  */
    if (*pverbose) 
      Rprintf("sweep %d\n", sweep);

    for (i=0; i<n; i++) {

      /* move cell i. */
      if (i < n1) {
	dmin_mu = d1; dmin_sd = sd1;
      } else {
	dmin_mu = d2; dmin_sd = sd2;
      }

      looking = 1; this_cell_rejects = 0;
      while (looking) {
      
	generate_dmin = 1;
	while (generate_dmin) {
	  this_dmin = dmin_mu + (dmin_sd * norm_rand());
	  if ( (this_dmin > lower) &&
	       ((upper <0) || (this_dmin < upper)))
	    generate_dmin = 0;
	  /*else printf("dminlul: dmin %f outside range\n", this_dmin);*/
	}

	/* generate a trial position at random. */
	x = xmin + (UNIF * wid); y = ymin + (UNIF * ht);
	xpts[i] = x; ypts[i] = y;

	bdmin_check(xpts, ypts, n1, n2, i,
		    this_dmin, d12, &constraint, &id);

	if (*pverbose) 
	  Rprintf("cell %d pos %f %f constraint %d id %d\n",
		  i, x, y, constraint, id);
	if (constraint == 0) {
	  /* constraints met fine. */
	  looking = 0;
	} else if ( constraint > 0) {
	  /* homotypic constraint broken. */
	  nrejects[0]++;
	} else {
	  /* heterotypic constraint broken. */
	  nrejects[1]++;
	}
	if (looking) {
	  if (this_cell_rejects++ > 99999) {
	    /* Taking far too long to replace one cell, so error. */
	    nrejects[0]= -1;
	    Rprintf("egg %d\n", this_cell_rejects);
	    sweep = 0; looking=0; i=n; /* end all loops */
	  }
	}
      }
      if(0 && *pverbose) {
	Rprintf("cell %d rejects %d\n", i, this_cell_rejects);
      }
    } /* next cell */
  } /* next sweep */

  if (*pverbose) {
    Rprintf("#rejects: %d %d\n", nrejects[0], nrejects[1]);
  }
  RANDOUT;
  
}

void bdmin_check(Sfloat *xpts, Sfloat *ypts, int n1, int n2,
		int i, Sfloat d, Sfloat d12, int *okay, int *id) {
  /* Check whether bivariate constraints are satisfied between cell I
   * and the neighbouring cells. */
  /* Return values:
   * OKAY - Set to:
   *  0 - constraints okay.
   * +1 -  broke the homotypic constraint.
   * -1 -  broke the heterotypic constraint.
   * In the case of a broken constraint, *ID stores the cell id that
   * broke the constraint.
   */

  Sfloat x, y;			/* coords of cell being checked. */
  Sfloat dist2, mind, dx, dy;
  Sfloat min_sofar =-1;		/* debugging */
  int n;			/* total number of cells */
  int j, homotypic;

  /* By default, everything okay */
  *okay = 0; *id = -1;
  
  x = xpts[i]; y = ypts[i];
  n = n1 + n2;

  /* All distance comparisons done with squared distances. */
  d *= d; d12 *= d12;
  for (j=0; j < n; j++) {
    if (i==j)
      continue;
    
    homotypic = ( (i < n1) == (j < n1)); /* 1 iff (i,j) pair is homotypic */
    dx = x-xpts[j]; dy = y-ypts[j];
    dist2 = (dx*dx) + (dy*dy);
    mind = (homotypic?d:d12);
    if ( (min_sofar <0) || (mind < dist2))
      min_sofar = dist2;
    if (dist2 < mind) {
      *okay = homotypic?1:-1;
      *id = j;
      break;
    }
  }

}

void dminlul3d(Sfloat *pwid, Sfloat *pht, Sfloat *pdepth, int *pnumcells,
	       Sfloat *pdmin, Sfloat *psd,
	       Sfloat *plower, Sfloat *pupper,
	       int *quiet,
	       Sfloat *xpts, Sfloat *ypts, Sfloat *zpts,
	       Sfloat *dmins, int *nrejects)

{
  /* 3-d version of dminlul.
   * Create NUMCELLS cells distributed in an volume of size wid*ht*depth.
   * Minimum distance between cells should be at least dmin.  On exit,
   * we return an array of cells stored in (xs, ys, zs).  NJRECTS stores
   * the number of rejected cells.  This is Lucia's version, where a
   * new dmin value is created every time, regardless of whether the
   * trial cell was accepted or not.  Also provide arguments LOWER
   * and UPPER: dmin values outside this range are rejected.
   */
  

  int num_rejects = 0, this_cell_rejects = 0;
  int looking, generate_dmin;
  int i;
  Sfloat x,y,z;
  Sfloat min; int idx;
  Sfloat this_dmin;
  Sfloat lower, upper;



  RANDIN;

  lower = *plower; upper = *pupper;
  for (i=0; i<*pnumcells; i++) {

    looking = 1; this_cell_rejects = 0;
    while (looking) {

      generate_dmin = 1;
      while (generate_dmin) {
	this_dmin = *pdmin + (*psd * norm_rand());
	if ( (this_dmin > lower) &&
	     ((upper <0) || (this_dmin < upper)))
	  generate_dmin = 0;
	/*else printf("dminlul: dmin %f outside range\n", this_dmin);*/
      }

      /*Rprintf("this val %lf\n", this_dmin); */

      x = UNIF * (*pwid); y = UNIF * (*pht); z = UNIF * (*pdepth); 
      
      nnd_n_3d( xpts, ypts, zpts, i, x, y, z, &idx, &min);
      if ( (min < this_dmin) ) {
	/* || (min > upper)*/
	/* Cannot test if the min distance is too big... this won't
	 * work for early cells!  i.e. when putting the first cell in,
	 * the upper distance will be the value returned by nnd_n,
	 * ie.s. some huge number. */

	 
	/* reject cell and try another. */
	num_rejects++;
	this_cell_rejects++;
	if (num_rejects > MAXREJECTS) {
	  /* If we cannot fit more cells in, return and indicate error
	   * by setting first x coordinate to negative value.
	   */
	  printf("dminlul error: num_rejects is too high (%s:%d)\n",
		 __FILE__, __LINE__);
	  xpts[0] = -1;
	  RANDOUT;
	  return;
	}
      }
      else {
	looking = 0;
      }
    }
    /* Accept this cell. */
    dmins[i] = this_dmin;
    xpts[i] = x; ypts[i] = y; zpts[i] = z;
    nrejects[i] = this_cell_rejects;
  }


  if (!*quiet) {
    Rprintf("#todo: 3-d packing density needs recalculating.\n");
    Rprintf("#rejects %d\tdminlul packing density %.3f\n", num_rejects,
	    ((*pnumcells * PI * *pdmin * *pdmin)/( 4 * *pwid * *pht)));
  }


  RANDOUT;
}


void pipp_lookup(Sfloat *pw, int *pn,
		 Sfloat *ph, Sfloat *pd, int *hlen,
		 int *pnsweeps, int *pverbose,
		 Sfloat *xpts, Sfloat *ypts, int *okay)
{

  /* PIPP: Pairwise interaction point processes.  The h() function is
   * converted into a lookup table for speed and generality.
   * On exit, if *OKAY is 1, the simulation was fine.
   */
  
  int num_rejects = 0, this_cell_rejects;
  int looking, nrejects;
  int i, j, sweep, constraint, id;
  int n;
  Sfloat x,y;
  Sfloat min; int idx;
  Sfloat prob, dx, dy, dist, p;
  Sfloat lower, upper;

  Sfloat xmin, xmax, ymin, ymax, wid, ht;
  RANDIN;

  sweep = *pnsweeps;
  n = *pn;
  
  xmin = pw[0]; xmax = pw[1]; ymin = pw[2]; ymax = pw[3];
  if (1 && *pverbose) {
    Rprintf("field %f %f %f %f\n", xmin, xmax, ymin, ymax);
    Rprintf("npts %d\n", n);
    Rprintf("LUT has %d entries in range (%f, %f) to (%f, %f)\n",
	    *hlen, pd[0], ph[0], pd[(*hlen)-1], ph[(*hlen)-1]);
  }
  wid = xmax - xmin; ht = ymax - ymin;

  *okay = 1;			/* assume simulation will be okay. */
  while( sweep-- > 0) {
    /*Perform one complete sweep, updating positions of cells.  */
    if (*pverbose) 
      Rprintf("sweep %d\n", sweep);

    for (i=0; i<n; i++) {

      /* move cell i. */
   
      looking = 1; this_cell_rejects = 0;
      while (looking) {
      
	/* generate a trial position at random. */
	x = xmin + (UNIF * wid); y = ymin + (UNIF * ht);
	xpts[i] = x; ypts[i] = y;

	prob = 1;
	for (j=0; j<n; j++) {
	  if (i!=j) {
	    dx = xpts[j] - x; dy = ypts[j] - y;
	    dist = sqrt( (dx*dx) + (dy*dy) );
	    hlookup(ph, pd, hlen, &dist, &p);
	    prob *= p;
	  }
	}

	if ( UNIF < prob) {
	  /* Accept this cell position. */
	  looking = 0;
	} else {
	  if (this_cell_rejects++ > 99999) {
	    /* Taking far too long to replace one cell, so error. */
	    Rprintf("pipp_d: trouble converging %d\n", this_cell_rejects);
	    *okay = 0;
	    sweep = 0; looking=0; i=n; /* end all loops */
	  }
	}
      }
    } /* next cell */
  } /* next sweep */

  RANDOUT;
    
}


void hlookup(Sfloat *h, Sfloat *ds, int *n, Sfloat *d, Sfloat *res)
  /* Find h(d) from the lookup table.  Lookup table has *N entries in
   * it.  Return result in *RES.  Do not use return value of function,
   * so that we can interface with R code.
   */
{
  int i;
  Sfloat d1, d2, h1, h2;

  /*
  Rprintf("lookup value of %f in following table\n", *d);
  for (i=0; i<*n; i++) {
    Rprintf("%.3f %.3f\n", h[i], ds[i]);
  }
  Rprintf("max entry: %f\n", ds[(*n)-1]);
  */
  
  if (*d < ds[0])
    *res = 0.0;			/* d lower than smallest entry in LUT */
  else {
    if (*d > ds[(*n)-1] ) {
      *res = 1.0;		/* d higher than entries in LUT */
    } else {
      /* lookup value of d. */
      i=0;
      while( ds[i] < *d) {
	i++;
      }
      /* d1 is smaller than d, d2 is larger than d. 
       * therefore *RES lies between h1 and h2. 
       * Use linear interpolation to guess *RES.*/
      d1 = ds[i-1]; d2 = ds[i];
      h1 = h[i-1]; h2 = h[i];

      *res = h1 + ( (h2-h1)/(d2-d1) * (*d-d1));
    }
  }
}
  


#ifdef unused
Sfloat hlookup(Sfloat *ph, Sfloat *ph_inc, Sfloat *ph_n, Sfloat *d) {
  /* Find h(d) from the lookup table.
   * All args given as pointers so that we can call it from R to check. */
  int index;
  index = (int) (*d/ *ph_inc);
  if (index >= *ph_n)
    return 1.0;			/* Assume end value. */
  else {
    return ph[index];		/* Could interpolate? */
  }
}
#endif

void dminacc(Sfloat *pwid, Sfloat *pht, int *pnumcells,
	     Sfloat *acc, Sfloat *dmax, Sfloat *inc,
	     Sfloat *xpts, Sfloat *ypts,
	     int *nrejects)

{
  /* Create NUMCELLS cells distributed in an area of size wid*ht.
   * Minimum distance between cells should be at least dmin.  On exit,
   * we return an array of cells stored in (xs, ys).  NJRECTS stores
   * the number of rejected cells.  This is Lucia's version, where a
   * new dmin value is created every time, regardless of whether the
   * trial cell was accepted or not.
   */
  

  int this_cell_rejects = 0;
  int looking;
  int i, dlow, dhigh;
  Sfloat x,y;
  Sfloat min; int idx;
  Sfloat p, plow, phigh;

  RANDIN;

  for (i=0; i<*pnumcells; i++) {

    looking = 1; this_cell_rejects = 0;
    while (looking) {

      x = UNIF * (*pwid); y = UNIF * (*pht);
      
      nnd_n( xpts, ypts, i, x, y, &idx, &min);

      if ( min >  *dmax) {
	looking = 0;
      } else {
	dlow = (int) min; dhigh = dlow+1;
	plow = acc[dlow];
	phigh = acc[dhigh];
	p = plow + ( ((phigh-plow)/ *inc) * (min - dlow));
	if (p > UNIF)
	  looking = 0;		/* accept cell */
	else 
	  this_cell_rejects++;
      }
    }

    /* Accept this cell. */
    xpts[i] = x; ypts[i] = y;
    nrejects[i] = this_cell_rejects;
  }
  
  RANDOUT;
}

void dminacc_bd(Sfloat *wid, Sfloat *ht, int *pnumcells,
		Sfloat *acc, Sfloat *dmax, Sfloat *inc,
		Sfloat *xpts, Sfloat *ypts,
		int *nrejects)

{
  /* Birth and death version of the dminacc function.
   */

  int this_cell_rejects = 0;
  int looking;
  int i, dlow, dhigh;
  int n, mm, id;
  Sfloat u;
  Sfloat min; int idx;
  Sfloat p, plow, phigh;

  RANDIN;

  n = *pnumcells;

  /* Create initial set of points. */
  for (i=0; i<n; i++) {
    xpts[i] = UNIF * (*wid); ypts[i] = UNIF * (*ht);
  }

  mm = 4 * n;			/* Number of repeats */
  mm *= 10;			/* if starting from Poisson distribution. */
  for (i=0; i<= mm; i++) {
    id = (int)(n*UNIF);	/* pick one point at random. */
    /*if (i%100==0) Rprintf("it %d:select new point for unit %d\n", i,id);*/
    xpts[id] = xpts[0]; ypts[id] = ypts[0];

    looking = 1;
    do {
      
      /* generate a new point. */
      xpts[0] = UNIF * (*wid); ypts[0] = UNIF * (*ht);
      u = UNIF;

      nnd( xpts, ypts, n, xpts[0], ypts[0], 0, &idx, &min);

      if ( min >  *dmax) {
	looking = 0;
      } else {
	dlow = (int) min; dhigh = dlow+1;
	plow = acc[dlow];
	phigh = acc[dhigh];
	p = plow + ( ((phigh-plow)/ *inc) * (min - dlow));
	if (p > UNIF)
	  looking = 0;		/* accept cell */
	else 
	  this_cell_rejects++;
      }
    } while (looking);

    /* Accept this cell. */
    nrejects[id] = this_cell_rejects;
  }
  
  RANDOUT;
}


void nnd_n(Sfloat *xs, Sfloat *ys, int n, Sfloat a, Sfloat b,
	   int *idx, Sfloat *min)
{
  /* Return the index (IDX) and distance (MIN) of cell closest to point (A,B).
   * The N datapoints are stored in XS[] and YS[].
   */
  if (n==0) {
    *idx = -1; *min = MAXDISTANCE;

    /* If n=0, no other cells have yet been positioned, so there is no
     * nearest neighbour. Hence, we just return some large distance to
     * indicate that fact. May be better to check that idx is less than zero,
     * since MAXDISTANCE may sometime be exceeded. */
  } else {
    nnd(xs, ys, n, a, b, -1, idx, min);
  }
}

void sjenndp(Sfloat *xs, Sfloat *ys, int *n,
	     Sfloat *a, Sfloat *b, int *ignore,
	     int *idx, Sfloat *min)
{
  /* Simple wrapper to call nnd from R -- make sure everything is a pointer. */
  /*Rprintf("point %f %f\n", *a, *b);*/
  nnd(xs, ys, *n,
      *a,  *b,
      *ignore,
      idx, min);
  /*Rprintf("closest id %d distance %f\n", *idx, *min);*/

}


void nnd(Sfloat *xs, Sfloat *ys, int n, Sfloat a, Sfloat b, int ignore,
	int *idx, Sfloat *min)
{
  /* Return the index (IDX) and distance (MIN) of cell closest to point (A,B).
   * The N datapoints are stored in XS[] and YS[].
   * In comparison to nnd_n(), cell number IGNORE is not allowed to be
   * the closest cell.
   */
  int i;
  int index = -1; Sfloat min_dist = 0;

  int first = 1;		/* check for 1st distance comparison. */

  Sfloat x, y, dx, dy, dist;
  
  for (i=0; i<n; i++) {
    if (i==ignore) {
      continue;
    }
    
    x = xs[i]; y = ys[i];
    dx = x - a; dy = y - b;
    dist = (dx*dx) + (dy*dy);
    if ( first || (dist < min_dist)) {
      first = 0;
      index = i;
      min_dist = dist;
    }
    
  }

  *min = sqrt(min_dist);
  *idx = index;
}
  

void nnd_n_3d(Sfloat *xs, Sfloat *ys, Sfloat *zs,
	      int n, Sfloat a, Sfloat b, Sfloat c,
	      int *idx, Sfloat *min)
{
  /* Return index (IDX) and distance (MIN) of cell closest to point (A,B,C).
   * The N datapoints are stored in XS[],  YS[], ZS[]/
   */
  if (n==0) {
    *idx = -1; *min = MAXDISTANCE;

    /* If n=0, no other cells have yet been positioned, so there is no
     * nearest neighbour. Hence, we just return some large distance to
     * indicate that fact. May be better to check that idx is less than zero,
     * since MAXDISTANCE may sometime be exceeded. */
  } else {
    nnd_3d(xs, ys, zs, n, a, b, c, -1, idx, min);
  }
}

void nnd_3d(Sfloat *xs, Sfloat *ys, Sfloat *zs, int n,
	    Sfloat a, Sfloat b, Sfloat c, int ignore,
	    int *idx, Sfloat *min)
{
  /* Return index (IDX) and distance (MIN) of cell closest to point (A,B,C).
   * The N datapoints are stored in XS[], YS[], ZS[] .
   * In comparison to nnd_n(), cell number IGNORE is not allowed to be
   * the closest cell.
   */
  int i;
  int index = -1; Sfloat min_dist = 0;

  int first = 1;		/* check for 1st distance comparison. */

  Sfloat x, y, z, dx, dy, dz, dist;
  
  for (i=0; i<n; i++) {
    if (i==ignore) {
      continue;
    }
    
    x = xs[i]; y = ys[i]; z = zs[i];
    dx = x - a; dy = y - b;  dz = z - c;
    dist = (dx*dx) + (dy*dy) + (dz*dz);
    if ( first || (dist < min_dist)) {
      first = 0;
      index = i;
      min_dist = dist;
    }
    
  }

  *min = sqrt(min_dist);
  *idx = index;
}
