## Do the dmin routine.
##
## Load with source("/home/stephen/langs/R/c-code/dminsd.r")

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

  ## Make up the return list.
  note <- paste("dminsd", dmin, dminsd)
  list(x = z$x, y = z$y, dmins = z$dmins, nrejects = z$nrejects, okay = okay,
       attempts = attempt, note = note)
  
}


dminl <- function(wid = 1000, ht = 1000, npts = 200,
                   dmin = 20, dminsd = 2)
{
  ## Simple version of Lucia's dmin.
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
    cat(paste ("dminl: ", dmin, dminsd, "fail after",
                dminmaxattempts, "tries\n"))
    ## just make a random distribution instead of a nice mosaic.
    z$x <- (runif(npts) * wid)
    z$y <- (runif(npts) * ht)
    okay <- FALSE
  }

  ## Make up the return list.
  note <- paste("dminl", dmin, dminsd)
  list(x = z$x, y = z$y, dmins = z$dmins, nrejects = z$nrejects, okay = okay,
       attempts = attempt, note = note)
  
}


dminlul <- function(wid = 1000, ht = 1000, npts = 200,
                    dmin = 20, dminsd = 2, lower = 0, upper = 100)
{
  ## Lower and upper bound version of Lucia's dmin.
  attempt <- 1
  okay <- TRUE
  trying <- TRUE
  while (trying && (attempt < dminmaxattempts)) {
    z <- .C("dminlul",
            as.double(wid),
            as.double(ht),
            as.integer(npts),
            as.double(dmin),
            as.double(dminsd),
            as.double(lower), as.double(upper),
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
    cat(paste ("dminlul: ", dmin, dminsd, "fail after",
                dminmaxattempts, "tries\n"))
    ## just make a random distribution instead of a nice mosaic.
    z$x <- (runif(npts) * wid)
    z$y <- (runif(npts) * ht)
    okay <- FALSE
  }

  ## Make up the return list.
  note <- paste("dminlul", dmin, dminsd, lower, upper,
                (if (!okay) "NOT OKAY"))
  res <- list(x = z$x, y = z$y, dmins = z$dmins,
              nrejects = z$nrejects, okay = okay,
              attempts = attempt, note = note)
  class(res) <- "sjedmin"
}

plot.sjedmin <- function(x) {
  ## Show the results of a dmin simulation.
  plot(x$x, x$y, asp=1, main=x$note)
}

dminlulbd <- function(wid = 1000, ht = 1000, npts = 200,
                    dmin = 20, dminsd = 2, lower = 0, upper = 100)
{
  ## Birth and Death version of dminlul
  attempt <- 1
  okay <- TRUE
  trying <- TRUE
  while (trying && (attempt < dminmaxattempts)) {
    z <- .C("dminlulbd",
            as.double(wid),
            as.double(ht),
            as.integer(npts),
            as.double(dmin),
            as.double(dminsd),
            as.double(lower), as.double(upper),
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
    cat(paste ("dminlulbd: ", dmin, dminsd, "fail after",
                dminmaxattempts, "tries\n"))
    ## just make a random distribution instead of a nice mosaic.
    z$x <- (runif(npts) * wid)
    z$y <- (runif(npts) * ht)
    okay <- FALSE
  }

  ## Make up the return list.
  note <- paste("dminlulbd", dmin, dminsd, lower, upper)
  list(x = z$x, y = z$y, dmins = z$dmins, nrejects = z$nrejects, okay = okay,
       attempts = attempt, note = note)
  
}


dminacc <- function(wid = 1000, ht = 1000, npts = 200,
                    acc, dmax = 10, inc=1.0)
{
  ## Shapiro et al. method for using acceptance function created directly
  ## from the K function of the exptl mosaic.  `acc' is the acceptance
  ## function created by `make.acc' and dmax is the maximal distance over
  ## which acceptance function is measured.  Cells with  NND greater than dmax
  ## are automatically accepted.
  attempt <- 1
  okay <- TRUE
  trying <- TRUE
  while (trying && (attempt < dminmaxattempts)) {
    z <- .C("dminacc",
            as.double(wid),
            as.double(ht),
            as.integer(npts),
            as.double(acc),
            as.double(dmax),
            as.double(inc),
            ## create memory to store return values.
            x = double(npts),
            y = double(npts),
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

  ## Make up the return list.
  note <- paste("dmin.acc", dmax)
  list(x = z$x, y = z$y, dmins = z$dmins, nrejects = z$nrejects, okay = okay,
       attempts = attempt, note = note)
  
}

make.acc <- function (l, dmax) {
  ## Make an acceptance function from a L function.
  ## This is useful for dmin.acc() and dminacc.bd()
  k <- l^2 * pi; kmax <- k[length(k)] 
  acc <- k / kmax; acc <- pmin(acc, 1.0)
  acc
}


dminacc.bd <- function(wid = 1000, ht = 1000, npts = 200,
                    acc, dmax = 10, inc=1.0)
{
  ## Birth and death version of dminacc().
  
  attempt <- 1
  okay <- TRUE
  trying <- TRUE
  while (trying && (attempt < dminmaxattempts)) {
    z <- .C("dminacc_bd",
            as.double(wid),
            as.double(ht),
            as.integer(npts),
            as.double(acc),
            as.double(dmax),
            as.double(inc),
            ## create memory to store return values.
            x = double(npts),
            y = double(npts),
            nrejects = integer(npts),
            )
    if (z$x[1] > 0)
      trying <- FALSE
    else
      attempt <- attempt + 1
  }
  if (attempt >= dminmaxattempts) {
    cat(paste ("dminacc.bd: ", dmax, "fail after",
                dminmaxattempts, "tries\n"))
    ## just make a random distribution instead of a nice mosaic.
    z$x <- (runif(npts) * wid)
    z$y <- (runif(npts) * ht)
    okay <- FALSE
  }

  ## Make up the return list.
  note <- paste("dminacc.bd", dmax)
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

  ## Make up the return list.
  note <- paste("damac", d1, beta, upper)
  list(x = z$x, y = z$y, note = note)
  
}

hamac <- function(t, lower, upper, beta)
{
  ## Test the hamac function in C.
  ## ts <- seq(0, 100, by=1);
  ## hs <- lapply(ts, function (x) hamac(x, 19, 76, 0.3)); plot(ts, hs)
  ## This was taken from the Diggle & Gratton (1984) paper.
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



