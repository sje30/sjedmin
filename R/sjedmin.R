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
            PACKAGE="sjedmin")
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
  res <- list(x = z$x, y = z$y, dmins = z$dmins,
              nrejects = z$nrejects, okay = okay,
              attempts = attempt, note = note)
  class(res) <- "sjedmin"
  res
}


dminl <- function(wid = 1000, ht = 1000, npts = 200,
                   dmin = 20, dminsd = 2, quiet=FALSE)
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
            as.integer(quiet),
            ## create memory to store return values.
            x = double(npts),
            y = double(npts),
            dmins = double(npts),
            nrejects = integer(npts), PACKAGE="sjedmin"
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
  res <- list(x = z$x, y = z$y, dmins = z$dmins,
              nrejects = z$nrejects, okay = okay,
              attempts = attempt, note = note)
  class(res) <- "sjedmin"
  res
}


dminlul <- function(wid = 1000, ht = 1000, npts = 200,
                    dmin = 20, dminsd = 2, lower = 0, upper = 100,
                    quiet = FALSE)
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
            as.integer(quiet),
            ## create memory to store return values.
            x = double(npts),
            y = double(npts),
            dmins = double(npts),
            nrejects = integer(npts), PACKAGE="sjedmin"
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
  res
}

dminlulfix2 <- function(wid = 1000, ht = 1000, npts = 200,
                        dmin = 20, dminsd = 2, lower = 0, upper = 100,
                        quiet = FALSE,
                        p2=matrix( c(100, 100, 200, 400), nrow=2),
                        d12=10)
{
  ## Lower and upper bound version of Lucia's dmin.
  attempt <- 1
  okay <- TRUE
  trying <- TRUE
  while (trying && (attempt < dminmaxattempts)) {
    z <- .C("dminlulfix2",
            as.double(wid),
            as.double(ht),
            as.integer(npts),
            as.double(dmin),
            as.double(dminsd),
            as.double(lower), as.double(upper),
            as.integer(quiet),
            ## create memory to store return values.
            x = double(npts),
            y = double(npts),
            dmins = double(npts),
            nrejects = integer(npts+1),
            as.double(p2[,1]), as.double(p2[,2]), as.integer(dim(p2)[1]),
            as.double(d12),
            PACKAGE="sjedmin"
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
  note <- paste("dminlulfix2", dmin, dminsd, lower, upper,
                (if (!okay) "NOT OKAY"))
  res <- list(x = z$x, y = z$y, dmins = z$dmins,
              nrejects = z$nrejects, okay = okay,
              ##args=match.call(),
              args=list(wid=wid, ht=ht, p2=p2, d12=d12),
              attempts = attempt, note = note)
  class(res) <- "sjedmin2"
  res
}

dminlul3d <- function(wid = 1000, ht = 1000, dep=1000, npts = 200,
                      dmin = 20, dminsd = 2, lower = 0, upper = 100,
                      quiet = FALSE)
{
  ## Lower and upper bound version of Lucia's dmin.
  attempt <- 1
  okay <- TRUE
  trying <- TRUE
  while (trying && (attempt < dminmaxattempts)) {
    z <- .C("dminlul3d",
            as.double(wid),
            as.double(ht),
            as.double(dep),
            as.integer(npts),
            as.double(dmin),
            as.double(dminsd),
            as.double(lower), as.double(upper),
            as.integer(quiet),
            ## create memory to store return values.
            x = double(npts),
            y = double(npts),
            z = double(npts),
            dmins = double(npts),
            nrejects = integer(npts), PACKAGE="sjedmin"
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
    z$z <- (runif(npts) * dep)
    okay <- FALSE
  }

  ## Make up the return list.
  note <- paste("dminlul", dmin, dminsd, lower, upper,
                (if (!okay) "NOT OKAY"))
  res <- list(x = z$x, y = z$y, z=z$z, dmins = z$dmins,
              nrejects = z$nrejects, okay = okay,
              attempts = attempt, note = note)
  class(res) <- "sjedmin"
  res
}

plot.sjedmin <- function(x) {
  ## Show the results of a dmin simulation.
  plot(x$x, x$y, asp=1, main=x$note)
}

     
plot.sjedmin2 <- function(x, r1=12, r2=r1) {
  ## Show the results of a dmin simulation.
  ## r1 is radius of cell 1; r2 is radius of cell 2.
  plot(NA, xlim=c(0, x$args$wid), ylim=c(0,x$args$ht), asp=1,main=x$note)
  symbols(x$x, x$y, circles=rep(r1, length(x$x)),
          inch=FALSE, add=TRUE)
  symbols(x$args$p2, circles=rep( r2, dim(x$args$p2)[1]),
          inch=FALSE, add=TRUE, bg="black")
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
            nrejects = integer(npts), PACKAGE="sjedmin"
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
  res <- list(x = z$x, y = z$y, dmins = z$dmins,
              nrejects = z$nrejects, okay = okay,
              attempts = attempt, note = note)
  class(res) <- "sjedmin"
  res
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
            nrejects = integer(npts), PACKAGE="sjedmin"
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
  res <- list(x = z$x, y = z$y, dmins = z$dmins,
              nrejects = z$nrejects, okay = okay,
              attempts = attempt, note = note)
  class(res) <- "sjedmin"
  res
}

make.acc <- function (l, dmax) {
  ## Make an acceptance function from a L function.
  ## This comes from the Shapiro et al. (1985) paper.
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
            nrejects = integer(npts), PACKAGE="sjedmin"
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
  res <- list(x = z$x, y = z$y, dmins = z$dmins,
              nrejects = z$nrejects, okay = okay,
              attempts = attempt, note = note)
  class(res) <- "sjedmin"
  res
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
            y = double(npts), PACKAGE="sjedmin"
            )
    trying <- FALSE
  }

  ## Make up the return list.
  note <- paste("damac", d1, beta, upper)
  res <- list(x = z$x, y = z$y, note = note)

  res
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
          h = double(1), PACKAGE="sjedmin"
          )
  z$h
}

sjennd <- function(x, y, a, b) {
  ## Call the Nearest neighbour distance function.
  ## x,y are  vectors of x,y coords of a group of points.  (a,b) is one point
  ## in space.  Return the id (and distance) of the point closest to (a,b).
  z <- .C("sjenndp",
          as.double(x),
          as.double(y),
          as.integer(length(x)),
          as.double(a), as.double(b),
          as.double(-1),
          ## create memory to store return values.
          id = integer(1),
          dist = double(1), PACKAGE="sjedmin"
          )
  zero.one.offset <- 1                  #C arrays are zero-based.
  res <- list(id=z$id + zero.one.offset, dist=z$dist)
  res
}
