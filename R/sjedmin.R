## Do the dmin routine.
##
## Load with source("/home/stephen/langs/R/c-code/dminsd.r")

dyn.load("/home/stephen/langs/R/c-code/pairwise_amac.so");

dminmaxattempts <- 5                    #number of attempts before giving up.

dminsd <- function(wid = 1000, ht = 1000, npts = 200,
                   dmin = 20, dminsd = 2) {

  attempt <- 1
  okay <- TRUE
  trying <- TRUE
  while (trying && (attempt < dminmaxattempts)) {
    z <- .C("dminsd",
            as.double(wid),
            as.double(ht),
            as.integer(npts),
            as.double(dmin),
            as.double(dminsd),
            ## create memory to store return values.
            x = double(npts),
            y = double(npts),
            dmins = double(npts),
            nrejects = integer(npts),
            )
    if (z$x[1] > 0)
      trying <- FALSE
    else
      attempt <- attempt + 1
  }
  if (attempt >= dminmaxattempts) {
    warning("dminsd: couldnt do the dmin after many attempts")
    ## just make a random distribution instead of a nice mosaic.
    z$x <- (runif(npts) * wid)
    z$y <- (runif(npts) * ht)
    okay <- FALSE
  }

  ## Even if we couldn't make it, just return something?
  ## Make up the return list.
  note <- paste("dminsd", dmin, dminsd)
  list(x = z$x, y = z$y, dmins = z$dmins, nrejects = z$nrejects, okay = okay)
  
}


dminl <- function(wid = 1000, ht = 1000, npts = 200,
                   dmin = 20, dminsd = 2) {
  ## Lucia's version of dmin -- new dmin value generated every step.
  attempt <- 1
  trying <- TRUE;
  while (trying && (attempt < dminmaxattempts)) {
    z <- .C("dminl",
            as.double(wid),
            as.double(ht),
            as.integer(npts),
            as.double(dmin),
            as.double(dminsd),
            ## create memory to store return values.
            x = double(npts),
            y = double(npts),
            dmins = double(npts),
            nrejects = integer(npts),
            )
    if (z$x[1] > 0)
      trying <- FALSE
    else
      attempt <- attempt + 1
  }
  if (attempt >= dminmaxattempts)
    warning("couldnt do the dmin after many attempts")

  ## Make up the return list.
  note <- paste("dminl", dmin, dminsd)
  list(x = z$x, y = z$y, dmins = z$dmins, nrejects = z$nrejects, note=note)
}

## dyn.unload("/home/stephen/langs/R/c-code/pairwise_amac.so");
