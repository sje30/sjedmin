## Do the dmin routine.
##
## Load with source("/home/stephen/langs/R/c-code/dminsd.r")

dyn.load("/home/stephen/langs/R/c-code/pairwise_amac.so");

dminmaxattempts <- 20                    #number of attempts before giving up.

dminsd <- function(wid = 1000, ht = 1000, npts = 200,
                   dmin = 20, dminsd = 2)
{
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
    cat(paste ("dminsd: ", dmin, dminsd, "fail after",
                dminmaxattempts, "tries\n"))
    ## just make a random distribution instead of a nice mosaic.
    z$x <- (runif(npts) * wid)
    z$y <- (runif(npts) * ht)
    okay <- FALSE
  }

  ## Even if we couldn't make it, just return something?
  ## Make up the return list.
  note <- paste("dminsd", dmin, dminsd)
  list(x = z$x, y = z$y, dmins = z$dmins, nrejects = z$nrejects, okay = okay,
       attempts = attempt, note = note)
  
}


dminl <- function(wid = 1000, ht = 1000, npts = 200,
                   dmin = 20, dminsd = 2)
{
  attempt <- 1
  okay <- TRUE
  trying <- TRUE
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
  if (attempt >= dminmaxattempts) {
    cat(paste ("dminsd: ", dmin, dminsd, "fail after",
                dminmaxattempts, "tries\n"))
    ## just make a random distribution instead of a nice mosaic.
    z$x <- (runif(npts) * wid)
    z$y <- (runif(npts) * ht)
    okay <- FALSE
  }

  ## Even if we couldn't make it, just return something?
  ## Make up the return list.
  note <- paste("dminl", dmin, dminsd)
  list(x = z$x, y = z$y, dmins = z$dmins, nrejects = z$nrejects, okay = okay,
       attempts = attempt, note = note)
  
}


damac <- function(wid = 1000, ht = 1000, npts = 200,
                  d1 = 10, upper = 20, beta = 1.5)
{
  attempt <- 1
  okay <- TRUE
  trying <- TRUE
  while (trying && (attempt < dminmaxattempts)) {
    z <- .C("pairwise_amac",
            as.double(wid),
            as.double(ht),
            as.integer(npts),
            as.double(d1),
            as.double(upper),
            as.double(beta),
            ## create memory to store return values.
            x = double(npts),
            y = double(npts),
            )
    trying <- FALSE
  }

  ## Even if we couldn't make it, just return something?
  ## Make up the return list.
  note <- paste("damac", d1, beta, upper)
  list(x = z$x, y = z$y, note = note)
  
}

hamac <- function(t, lower, upper, beta)
{
  ## Test the hamac function
  ## ts <- seq(0, 100, by=1);
  ## hs <- lapply(ts, function (x) hamac(x, 19, 76, 0.3)); plot(ts, hs)

  z <- .C("h_amac",
          as.double(t),
          as.double(lower),
          as.double(upper),
          as.double(beta),
          ## create memory to store return values.
          h = double(1),
          )
  z$h
}




## dyn.unload("/home/stephen/langs/R/c-code/pairwise_amac.so");
