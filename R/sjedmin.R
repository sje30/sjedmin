## Dmin routines.
## Stephen Eglen


.dminmaxattempts <- 5                    #number of attempts before giving up.

dminsd <- function(wid = 1000, ht = 1000, npts = 200,
                   dmin = 20, dminsd = 2)
{
  attempt <- 1
  okay <- TRUE
  trying <- TRUE
  while (trying && (attempt < .dminmaxattempts)) {
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
  if (attempt >= .dminmaxattempts) {
    cat(paste ("dminsd: ", dmin, dminsd, "fail after",
                .dminmaxattempts, "tries\n"))
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
  while (trying && (attempt < .dminmaxattempts)) {
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
  if (attempt >= .dminmaxattempts) {
    cat(paste ("dminl: ", dmin, dminsd, "fail after",
                .dminmaxattempts, "tries\n"))
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


dminlul <- function(w, npts = 200,
                    dmin = 20, dminsd = 2, lower = 0, upper = 100,
                    quiet = TRUE)
{
  ## Lower and upper bound version of Lucia's dmin.
  ## Safety check since a lot of code uses dminlul.
  if( length(w) != 4)
    stop(paste("w (", paste(w, collapse=' '), ") should be of length 4."))
  attempt <- 1
  okay <- TRUE
  trying <- TRUE
  while (trying && (attempt < .dminmaxattempts)) {
    z <- .C("dminlul",
            as.double(w),
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
  if (attempt >= .dminmaxattempts) {
    cat(paste ("dminlul: ", dmin, dminsd, "fail after",
                .dminmaxattempts, "tries\n"))
    ## just make a random distribution instead of a nice mosaic.
    z$x <- (runif(npts) * 1)
    z$y <- (runif(npts) * 1)
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

dminlulfix2 <- function(w,
                        npts = 200,
                        dmin = 20, dminsd = 2, lower = 0, upper = 100,
                        quiet = TRUE,
                        p2=matrix( c(100, 100, 200, 400), nrow=2),
                        d12=10)
{
  ## Lower and upper bound version of Lucia's dmin.
  attempt <- 1
  okay <- TRUE
  trying <- TRUE
  while (trying && (attempt < .dminmaxattempts)) {
    z <- .C("dminlulfix2",
            as.double(w),
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
  if (attempt >= .dminmaxattempts) {
    cat(paste ("dminlul: ", dmin, dminsd, "fail after",
                .dminmaxattempts, "tries\n"))
    ## just make a random distribution instead of a nice mosaic.
    z$x <- w[1] +  (runif(npts) * (w[2] - w[1]))
    z$y <- w[3] + (runif(npts) * (w[4] - w[3]))
    okay <- FALSE
  }

  ## Make up the return list.
  note <- paste("dminlulfix2", dmin, dminsd, lower, upper,
                "d12", d12,
                (if (!okay) "NOT OKAY"))
  res <- list(x = z$x, y = z$y, dmins = z$dmins,
              nrejects = z$nrejects, okay = okay,
              ##args=match.call(),
              args=list(w=w, dmin=dmin, dminsd=dminsd, p2=p2, d12=d12),
              attempts = attempt, note = note)
  class(res) <- "sjedmin2"
  res
}

bdmin.bd <- function(w=c(0, 1000, 0, 1000),
                     pts=NULL,
                     n1=100, n2=100,
                     d1=20, d1.sd=2,
                     d2=20, d2.sd=2,
                     d12=12, d12.sd=-1,
                     lower = 0, upper = -1,
                     nsweeps=10,
                     verbose = FALSE)
{
  ## Bivariate dmin simulation, with birth&death algorithm.

  ## If d12.sd is negative, it means we do not have any variance on d12.
  
  npts <- n1+n2

  ## If initial PTS are not provided, generate some at random.
  if (is.null(pts))
    pts <- cbind(w[1] +  (runif(npts) * (w[2] - w[1])),
                 w[3] + (runif(npts) * (w[4] - w[3])))
  
  ## Could include a check that pts are within window.
  params <- c(d1, d1.sd, d2, d2.sd, d12, d12.sd, lower, upper)
  attempt <- 1
  okay <- TRUE
  trying <- TRUE
  while (trying && (attempt < .dminmaxattempts)) {
    z <- .C("bdmin_bd",
            as.double(w),
            as.integer(n1), as.integer(n2),
            as.double(params),
            as.integer(nsweeps),
            as.integer(verbose),
            ## create memory to store return values.
            x = as.double(pts[,1]),
            y = as.double(pts[,2]),
            nrejects = integer(2),
            PACKAGE="sjedmin")
    if (z$nrejects[1] >= 0)
      trying <- FALSE
    else
      attempt <- attempt + 1
  }
  if (attempt >= .dminmaxattempts) {
    cat(paste ("bdmin: ", paste(params, collapse=' '), "fail after",
                .dminmaxattempts, "tries\n"))
    ## just make a random distribution
    z$x <- w[1] +  (runif(npts) * (w[2] - w[1]))
    z$y <- w[3] + (runif(npts) * (w[4] - w[3]))
    okay <- FALSE
  }

  ## Make up the return list.
  note <- paste("bdmin")
  res <- list(x = z$x, y = z$y, dmins = z$dmins,
              n1=n1, n2=n2,
              nrejects = z$nrejects, okay = okay,
              ##args=match.call(),
              args=list(w=w, params=params),
              attempts = attempt, note = note, w=w)
  class(res) <- "sjebdmin"
  res
}

plot.sjebdmin <- function(x, ...) {
  ## Plotting function for output from bdmin.bd.
  plot(x$x, x$y, asp=1, pch=19,
       main=title(paste("bdmin", paste(x$args$params, collapse='  '),
         ifelse(x$okay, "OK", "!OK")         )),
       col= c( rep("green", x$n1), rep("orangered", x$n2)),
       ...)
  rect( x$w[1], x$w[3], x$w[2], x$w[4], lty=2)
}


pipp.lookup <- function(w=c(0, 1000, 0, 1000),
                 pts=NULL,
                 n1=100,
                 h, d, 
                 nsweeps=10,
                 verbose = FALSE, tor=FALSE)
{
  ## pipp: Pairwise Interaction Point Process
  ## Works using a lookup table idea.

  npts <- n1
  ## If initial PTS are not provided, generate some at random.
  if (is.null(pts))
    pts <- cbind(w[1] +  (runif(npts) * (w[2] - w[1])),
                 w[3] + (runif(npts) * (w[4] - w[3])))
  
  ## Could include a check that pts are within window.

  attempt <- 1
  okay <- TRUE
  trying <- TRUE
  while (trying && (attempt < .dminmaxattempts)) {
    z <- .C("pipp_lookup",
            as.double(w),
            as.integer(n1),
            as.double(h),
            as.double(d),
            as.integer(length(d)),
            as.integer(nsweeps),
            as.integer(verbose),
            as.integer(tor),
            x = as.double(pts[,1]),
            y = as.double(pts[,2]),
            okay = integer(1),
            PACKAGE="sjedmin")
    if (z$okay > 0)
      trying <- FALSE
    else
      attempt <- attempt + 1
  }
  if (attempt >= .dminmaxattempts) {
    cat(paste ("pipp.lookup fail after", 
                .dminmaxattempts, "tries\n"))
    ## just make a random distribution
    z$x <- w[1] +  (runif(npts) * (w[2] - w[1]))
    z$y <- w[3] + (runif(npts) * (w[4] - w[3]))
    okay <- FALSE
  }

  ## Make up the return list.
  note <- paste("pipp.lookup")
  res <- list(x = z$x, y = z$y, 
              okay = okay,
              ##args=match.call(),
              args=list(w=w, h=h, d=d),
              attempts = attempt, note = note, w=w)
  class(res) <- "pipp"
  res
}

pipp2.lookup <- function(w=c(0, 1000, 0, 1000),
                         pts1=NULL, pts2=NULL,
                         n1=100, n2 = 100,
                         h1, d1, h2, d2, h12, d12,
                         nsweeps=10,
                         fix=0,
                         verbose = FALSE, tor=FALSE)
{
  ## pipp: Pairwise Interaction Point Process
  ## Works using a lookup table idea.

  ## FIX set to 0 (move all cells), 1 (fix type 1; move type 2)
  ## or 2 (fix type 2, move type 1)
  npts <- n1 + n2
  ## If initial PTS are not provided, generate some at random.
  if (is.null(pts1))
    pts1 <- cbind(w[1] +  (runif(n1) * (w[2] - w[1])),
                 w[3] + (runif(n1) * (w[4] - w[3])))

  if (is.null(pts2))
    pts2 <- cbind(w[1] +  (runif(n2) * (w[2] - w[1])),
                  w[3] + (runif(n2) * (w[4] - w[3])))

  pts <- rbind(pts1, pts2)
  ## Could include a check that pts are within window.

  attempt <- 1
  okay <- TRUE
  trying <- TRUE
  while (trying && (attempt < .dminmaxattempts)) {
    z <- .C("pipp2_lookup",
            as.double(w),
            as.integer(n1), as.integer(n2),
            as.double(h1), as.double(d1), as.integer(length(d1)),
            as.double(h2), as.double(d2), as.integer(length(d2)),
            as.double(h12), as.double(d12), as.integer(length(d12)),
            as.integer(nsweeps),
            as.integer(verbose),
            as.integer(fix),
            as.integer(tor),
            x = as.double(pts[,1]),
            y = as.double(pts[,2]),
            okay = integer(1),
            PACKAGE="sjedmin")
    if (z$okay > 0)
      trying <- FALSE
    else
      attempt <- attempt + 1
  }
  if (attempt >= .dminmaxattempts) {
    cat(paste ("pipp2.lookup fail after", 
                .dminmaxattempts, "tries\n"))
    ## just make a random distribution
    z$x <- w[1] +  (runif(npts) * (w[2] - w[1]))
    z$y <- w[3] + (runif(npts) * (w[4] - w[3]))
    okay <- FALSE
  }

  ## Make up the return list.
  note <- paste("pipp.lookup")
  res <- list(x = z$x, y = z$y, 
              okay = okay,
              n1=n1, n2=n2,
              ##args=match.call(),
              args=list(w=w, h1=h1, d1=d1, h2=h2, d2=d2, h12=h12, d12=d12),
              attempts = attempt, note = note, w=w)
  class(res) <- "pipp2"
  res
}

plot.pipp2 <- function(x, ...) {
  ## Plotting function for output from pipp2.lookup().
  plot(x$x, x$y, asp=1, pch=19,
       main=title(paste("pipp2"), 
         ifelse(x$okay, "OK", "!OK") ),
       col= c( rep("green", x$n1), rep("orangered", x$n2)),
       ...)
  rect( x$w[1], x$w[3], x$w[2], x$w[4], lty=2)
}

plot.pipp <- function(x, ...) {
  ## Plotting function for output from pipp.lookup().
  main <- NULL #TODO -- fix this to insert properly
  ## https://stackoverflow.com/questions/17390236/r-using-ellipsis-argument
  if (is.null(main)) {
    main <- paste("pipp", ifelse(x$okay, "OK", "!OK"))
  }
  plot(x$x, x$y, asp=1, pch=19, main=main)
  rect( x$w[1], x$w[3], x$w[2], x$w[4], lty=2)
}

hlookup <- function(h, d, r) {
  z <- .C("hlookup", as.double(h), as.double(d),
          as.integer(length(h)),
          as.double(r),
          ## return value.
          res=double(1),
          PACKAGE="sjedmin")
  z$res
}
                    
.dminlul3d <- function(wid = 1000, ht = 1000, dep=1000, npts = 200,
                      dmin = 20, dminsd = 2, lower = 0, upper = 100,
                      quiet = FALSE)
{
  ## Lower and upper bound version of Lucia's dmin.
  attempt <- 1
  okay <- TRUE
  trying <- TRUE
  while (trying && (attempt < .dminmaxattempts)) {
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
  if (attempt >= .dminmaxattempts) {
    cat(paste ("dminlul: ", dmin, dminsd, "fail after",
                .dminmaxattempts, "tries\n"))
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
  class(res) <- "sjedmin3d"
  res
}

plot.sjedmin <- function(x, ...) {
  ## Show the results of a dmin simulation.
  plot(x$x, x$y, asp=1, main=x$note, xlab='', ylab='', ...)
}

     
plot.sjedmin2 <- function(x, ...) {
  ## Show the results of a dmin simulation.
  ## r1 is radius of cell 1; r2 is radius of cell 2.
  r1=12; r2=r1 # TODO -- allow defaults to be changed by user in ...
  plot(NA, xlab='', ylab='',
       xlim=c(x$args$w[1], x$args$w[2]),
       ylim=c(x$args$w[3], x$args$w[4]),
       asp=1,main=x$note)
  symbols(x$x, x$y, circles=rep(r1, length(x$x)),
          inch=FALSE, add=TRUE)
  symbols(x$args$p2, circles=rep( r2, dim(x$args$p2)[1]),
          inch=FALSE, add=TRUE, bg="black")
}

.dminlulbd <- function(w, npts = 200,
                    dmin = 20, dminsd = 2, lower = 0, upper = 100)
{
  ## Birth and Death version of dminlul
  attempt <- 1
  okay <- TRUE
  trying <- TRUE
  while (trying && (attempt < .dminmaxattempts)) {
    z <- .C("dminlulbd",
            as.double(w),
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
  if (attempt >= .dminmaxattempts) {
    cat(paste ("dminlulbd: ", dmin, dminsd, "fail after",
                .dminmaxattempts, "tries\n"))
    ## just make a random distribution instead of a nice mosaic.
    z$x <- (runif(npts) * 1)
    z$y <- (runif(npts) * 1)
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


.dminacc <- function(wid = 1000, ht = 1000, npts = 200,
                    acc, dmax = 10, inc=1.0)
{
  ## Shapiro et al. method for using acceptance function created directly
  ## from the K function of the exptl mosaic.  `acc' is the acceptance
  ## function created by `.make.acc' and dmax is the maximal distance over
  ## which acceptance function is measured.  Cells with  NND greater than dmax
  ## are automatically accepted.
  attempt <- 1
  okay <- TRUE
  trying <- TRUE
  while (trying && (attempt < .dminmaxattempts)) {
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
  if (attempt >= .dminmaxattempts) {
    cat(paste ("dminsd: ", dmin, dminsd, "fail after",
                .dminmaxattempts, "tries\n"))
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

.make.acc <- function (l, dmax) {
  ## Make an acceptance function from a L function.
  ## This comes from the Shapiro et al. (1985) paper.
  ## This is useful for dmin.acc() and dminacc.bd()
  k <- l^2 * pi; kmax <- k[length(k)] 
  acc <- k / kmax; acc <- pmin(acc, 1.0)
  acc
}


.dminacc.bd <- function(wid = 1000, ht = 1000, npts = 200,
                    acc, dmax = 10, inc=1.0)
{
  ## Birth and death version of dminacc().
  attempt <- 1
  okay <- TRUE
  trying <- TRUE
  while (trying && (attempt < .dminmaxattempts)) {
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
  if (attempt >= .dminmaxattempts) {
    cat(paste ("dminacc.bd: ", dmax, "fail after",
                .dminmaxattempts, "tries\n"))
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


.damac <- function(wid = 1000, ht = 1000, npts = 200,
                  d1 = 10, upper = 20, beta = 1.5)
{
  attempt <- 1
  okay <- TRUE
  trying <- TRUE
  while (trying && (attempt < .dminmaxattempts)) {
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

.hamac <- function(t, lower, upper, beta)
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

.sjennd <- function(x, y, a, b) {
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
