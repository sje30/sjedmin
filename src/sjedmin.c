#include <R.h>
#include <S.h>			/* for seed_in, seed_out */

/* Use the definitions by Venables & Ripley. */
#  define RANDIN  seed_in((long *)NULL)
#  define RANDOUT seed_out((long *)NULL)
#  define UNIF unif_rand()

/* TODO items:
 * Maybe pass all parameters in one vector, rather than using
 * many different arguments.  makes the C/R interface shorter.
 *
 * Birth & death routine: return the dmin value used at each iteration?
 * To do this, I need to know in advance how many iterations to run
 * for -- maybe have default mm as 40 * npts in R code?
 * TODO: add code for periodic checks for interrupts; Brian Ripley uses:
 * ... 	    if(attempts % 1000 == 0) R_CheckUserInterrupt();
 * nested within his loops.  makes sense.
 */

#define MIN(X, Y)  ((X) < (Y) ? (X) : (Y))

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

Sfloat dist2d_tor(int tor, Sfloat x0, Sfloat y0,
		  Sfloat x1, Sfloat y1,
		  Sfloat wid, Sfloat ht);

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
	/*Rprintf("xxx valu %lf\n", h1);*/
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
  /*Rprintf("ret valu %lf\n", *h);*/
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
  
  int i,j, mm, id, okay, generate_r, num_rejects;
  Sfloat x1, y1, dist2, r, rr, lower, upper;
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
    num_rejects = 0;		/* reset to zero. */
    do {
      
      /* generate a new point. */
      xpts[0] = xmin + (UNIF * wid); ypts[0] = ymin + (UNIF * ht);
/*       Rprintf("generate new trial at %f %f\n", xpts[0], ypts[0]); */
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

	  num_rejects++;
	  if (num_rejects > MAXREJECTS) {
	    /* If we cannot fit more cells in, return and indicate error
	     * by setting first x coordinate to negative value.
	     */
	    Rprintf("dminlulbd error: num_rejects is too high (%s:%d)\n",
		   __FILE__, __LINE__);
	    xpts[0] = -1;
	    RANDOUT;
	    return;
	  }

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
	  Rprintf("dminsd error: num_rejects is too high (%s:%d)\n",
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
	  Rprintf("dminl error: num_rejects is too high (%s:%d)\n",
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
	/*else Rprintf("dminlul: dmin %f outside range\n", this_dmin);*/
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
	  Rprintf("dminlul error: num_rejects is too high (%s:%d)\n",
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
	/*else Rprintf("dminlul: dmin %f outside range\n", this_dmin);*/
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
	  Rprintf("dminlul error: num_rejects is too high (%s:%d)\n",
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
  int n, n1, n2, regen_d12;
  Sfloat x,y;
  Sfloat min; int idx;
  Sfloat this_dmin, this_d12;
  Sfloat lower, upper;

  Sfloat xmin, xmax, ymin, ymax, wid, ht;
  Sfloat d1, sd1, d2, sd2, d12, d12sd;
  Sfloat dmin_mu, dmin_sd;
  RANDIN;

  /* Extract relevant exclusion zone parameters. */
  d1  = params[0]; sd1 = params[1];
  d2  = params[2]; sd2 = params[3];
  d12 = params[4]; d12sd = params[5];
  lower = params[6]; upper = params[7];
  
  sweep = *pnsweeps;
  n1 = *pn1; n2 = *pn2; n = n1 + n2;
  
  xmin = pw[0]; xmax = pw[1]; ymin = pw[2]; ymax = pw[3];
  regen_d12 = (d12sd > 0);
  if (1 && *pverbose) {
    Rprintf("field %f %f %f %f\n", xmin, xmax, ymin, ymax);
    Rprintf("n1 %d n2 %d\n", n1, n2);
    Rprintf("d1 %f +/- %f d2 %f +/- %f d12 %f +/- %f\n",
	    d1, sd1, d2, sd2, d12, d12sd);
    Rprintf("regen d12: %d\n", regen_d12);
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

	if (regen_d12) {
	  generate_dmin = 1;
	  while (generate_dmin) {
	    this_d12 = d12 + (d12sd * norm_rand());
	    if  (this_d12 > lower)
	      generate_dmin = 0;
	  }
	} else {
	  this_d12 = d12;
	}
	generate_dmin = 1;
	while (generate_dmin) {
	  this_dmin = dmin_mu + (dmin_sd * norm_rand());
	  if ( (this_dmin > lower) &&
	       ((upper <0) || (this_dmin < upper)))
	    generate_dmin = 0;
	  /*else Rprintf("dminlul: dmin %f outside range\n", this_dmin);*/
	}

	/* generate a trial position at random. */
	x = xmin + (UNIF * wid); y = ymin + (UNIF * ht);
	xpts[i] = x; ypts[i] = y;

	bdmin_check(xpts, ypts, n1, n2, i,
		    this_dmin, this_d12, &constraint, &id);

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
  /* 3-d version of Dmin rule.
   *
   * (P prefix on variable names indicates that variable is a pointer
   * to the value, not the value.)
   *
   * Create NUMCELLS cells distributed in an volume of size WID*HT*DEPTH.
   * This assumes one corner of the volume is the origin.
   *
   * The dmin value is chosen from a Normal distribution with mean
   * DMIN and s.d. of SD.  The dmin value is also bounded to lie
   * between LOWER and UPPER.  (LOWER is normally set to the size of a
   * cell body and UPPER can be ignored by setting it to a negative
   * number.)
   *
   * On exit, we return an array of cells stored in (XPTS, YPTS,
   * ZPTS).  If the simulation could not be performed (e.g. if DMIN
   * too high), the first value of XPTS is set to a negative value.
   * NJRECTS is a vector; NREJECTS[i] stores the number of cells that
   * were rejected when trying to position cell i.
   *
   * A new dmin value is created every time, regardless of whether
   * the trial cell was accepted or not.
   *
   * QUIET is a debugging option.
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
	/*else Rprintf("dminlul: dmin %f outside range\n", this_dmin);*/
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
	  Rprintf("dminlul error: num_rejects is too high\n");
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
    Rprintf("#rejects %d\tdminlul3d packing density %.3f\n", num_rejects,
	    ((*pnumcells * PI * 1/8 * *pdmin * *pdmin * *pdmin)/
	     ( (3.0/4.0) * *pdepth * *pwid * *pht)));
  }

  RANDOUT;
}


void pipp_lookup(Sfloat *pw, int *pn,
		 Sfloat *ph, Sfloat *pd, int *hlen,
		 int *pnsweeps, int *pverbose, int *tor,
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
	    dist = dist2d_tor(*tor, x, y, xpts[j], ypts[j], wid, ht);
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

Sfloat dist2d_tor(int tor, Sfloat x0, Sfloat y0,
		  Sfloat x1, Sfloat y1,
		  Sfloat wid, Sfloat ht) {
  /* Return distance between (x0,y0) and (x1, y1).
   * If TOR is non-zero, we assume toroidal wrap-around.
   */

  Sfloat dx, dy, dxt, dyt;
  if (tor) {
    dx = fabs(x0-x1); dy = fabs(y0-y1);
    dxt = MIN(dx, wid-dx);
    dyt = MIN(dy, ht-dy);
  } else {
    dxt = x0-x1; dyt = y0-y1;
  }
  return ( sqrt( (dxt*dxt) + (dyt*dyt) ));
  
}

Sfloat dist2d_tor_R(int *tor, Sfloat *x0, Sfloat *y0,
		    Sfloat *x1, Sfloat *y1,
		    Sfloat *wid, Sfloat *ht, Sfloat *d) {
  /* Wrapper function for R, so we can test easily. */
  
  *d = dist2d_tor(*tor, *x0, *y0, *x1, *y1, *wid, *ht);

}

void pipp2_lookup(Sfloat *pw, int *pn1, int *pn2,
		  Sfloat *ph1, Sfloat *pd1, int *hlen1,
		  Sfloat *ph2, Sfloat *pd2, int *hlen2,
		  Sfloat *ph12, Sfloat *pd12, int *hlen12,
		  int *pnsweeps, int *pverbose, int *pfix,
		  int *tor,
		  Sfloat *xpts, Sfloat *ypts, int *okay)
{
  /* PIPP2: Bivariate Pairwise interaction point processes.
   * Like pipp_lookup, but bivariate, so we have three h functions():
   * h1 - type 1
   * h2 - type 2
   * h12 - h for cross type.
   * The initial points are stored in XPTS and YPTS, with type 1 points
   * given first, followed by type 2 points.
   * On exit, if *OKAY is 1, the simulation was fine.
   */
  
  int num_rejects = 0, this_cell_rejects;
  int looking, nrejects;
  int i, j, sweep, constraint, id;
  int n, n1, n2;
  Sfloat x,y;
  Sfloat min; int idx;
  Sfloat prob, dx, dy, dist, p;
  Sfloat lower, upper;
  Sfloat *h_homo, *d_homo;
  int *len_homo;
  int i_lo, i_hi;
  
  Sfloat xmin, xmax, ymin, ymax, wid, ht;
  RANDIN;

    
  sweep = *pnsweeps;
  n1 = *pn1; n2 = *pn2; n = n1+ n2;
  xmin = pw[0]; xmax = pw[1]; ymin = pw[2]; ymax = pw[3];

  switch(*pfix) {
  case 1: 			/* fix type 1 cells, just move type 2. */
    i_lo=n1; i_hi=n; break;	/* fix type 2 cells, just move type 1.  */
  case 2:
    i_lo=0; i_hi = n1; break;
  default: /* includes zero */
    i_lo=0; i_hi=n; break;
    /*Rprintf("In trouble, should not get here!");*/
  }


  if (1 && *pverbose) {
    Rprintf("field %f %f %f %f\n", xmin, xmax, ymin, ymax);
    Rprintf("npts %d (%d type 1; %d type 2)\n", n, n1, n2);
    Rprintf("fix %d lo %d hi %d\n", *pfix, i_lo, i_hi);
    Rprintf("LUT h1 has %d entries in range (%f, %f) to (%f, %f)\n",
	    *hlen1, pd1[0], ph1[0], pd1[(*hlen1)-1], ph1[(*hlen1)-1]);
    Rprintf("LUT h2 has %d entries in range (%f, %f) to (%f, %f)\n",
	    *hlen2, pd2[0], ph2[0], pd2[(*hlen2)-1], ph2[(*hlen2)-1]);
    Rprintf("LUT h12 has %d entries in range (%f, %f) to (%f, %f)\n",
	    *hlen12, pd12[0], ph12[0], pd12[(*hlen12)-1], ph12[(*hlen12)-1]);
  }
  wid = xmax - xmin; ht = ymax - ymin;

  *okay = 1;			/* assume simulation will be okay. */
  while( sweep-- > 0) {
    /*Perform one complete sweep, updating positions of cells.  */
    if (*pverbose) 
      Rprintf("sweep %d\n", sweep);

    for (i=i_lo; i< i_hi; i++) {

      /* move cell i. */

      if (i<n1) {
	/* Cell I is type 1. */
	h_homo = ph1; d_homo = pd1; len_homo = hlen1;
      } else {
	/* Cell I is type 2. */
	h_homo = ph2; d_homo = pd2; len_homo = hlen2;
      }
      
      looking = 1; this_cell_rejects = 0;
      while (looking) {
      
	/* generate a trial position at random. */
	x = xmin + (UNIF * wid); y = ymin + (UNIF * ht);
	xpts[i] = x; ypts[i] = y;

	prob = 1;
	for (j=0; j<n; j++) {
	  if (i!=j) {
	    /* dx = xpts[j] - x; dy = ypts[j] - y; */
	    dist = dist2d_tor(*tor, x, y, xpts[j], ypts[j], wid, ht);
	    /* Do we have a homotypic or heterotypic interaction between
	     * cell I and cell J ? */
	    if ( ( i < n1 ) == (j < n1)) {
	      /* Homotypic interaction between (I,J) */
	      hlookup(h_homo, d_homo, len_homo, &dist, &p);
	    } else {
	      /* Homotypic interaction between (I,J) */
	      hlookup(ph12, pd12, hlen12, &dist, &p);
	    }
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
	    sweep = 0; looking=0; i=i_hi; /* end all loops */
	  }
	}
      }
    } /* next cell */
  } /* next sweep */

  RANDOUT;
    
}

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


#ifdef unused

/* For reference, here is Brian Ripley's implementation of Strauss Process
 *  taken from pps.c in VR bundle (VR/spatial/src/pps.c)
 */
void
VR_simpat(Sint *npt, Sfloat *x, Sfloat *y, Sfloat *c,
	  Sfloat *r, Sint *init)
{
    int   i, attempts = 0, id, j, mm, n = *npt;
    Sfloat cc, rr, ax, ay, d, x1, y1, u;

    testinit();
    cc = *c;
    if (cc >= 1.0) {
	VR_pdata(npt, x, y);
	return;
    }
    RANDIN;
    ax = xu0 - xl0;
    ay = yu0 - yl0;
    rr = (*r) * (*r);
    mm = 4 * n;
    if (*init > 0) mm = 10 * mm;
    for (i = 1; i <= mm; i++) {
	id = floor(n * UNIF);
	x[id] = x[0];
	y[id] = y[0];
	do {
	    attempts++;
	    x[0] = xl0 + ax * UNIF;
	    y[0] = yl0 + ay * UNIF;
	    u = UNIF;
	    d = 1.0;
	    for (j = 1; j < n; j++) {
		x1 = x[j] - x[0];
		y1 = y[j] - y[0];
		if (x1 * x1 + y1 * y1 < rr) {
		    d *= cc;
		    if (d < u) continue;
		}
	    }
	    if(attempts % 1000 == 0) R_CheckUserInterrupt();
	}
	while (d < u);
    }
    RANDOUT;
}


#endif

