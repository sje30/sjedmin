\name{sjedmin}
\alias{dminsd}
\alias{dminl}
\alias{dminlul}
\alias{dminlulfix2}
\alias{bdmin.bd}
\alias{hlookup}
\alias{plot.sjedmin}
\alias{plot.sjedmin2}
\alias{plot.pipp}
\alias{plot.pipp2}
\alias{plot.sjebdmin}
\alias{pipp.lookup}
\alias{pipp2.lookup}
\title{Dmin mosaics}
\description{
  Create mosaics that respect the dmin rule.
}
\usage{
dminl(wid, ht, npts, dmin, dminsd, quiet)
dminlul(w, npts, dmin, dminsd, lower, upper, quiet)
dminsd(wid = 1000, ht = 1000, npts = 200,
       dmin = 20, dminsd = 2)

dminlulfix2(w, npts = 200,
  dmin = 20, dminsd = 2, lower = 0, upper = 100,
  quiet = TRUE,
  p2=matrix( c(100, 100, 200, 400), nrow=2),
  d12=10)
bdmin.bd(w = c(0, 1000, 0, 1000), pts = NULL, n1 = 100, n2 =
                   100, d1 = 20, d1.sd = 2, d2 = 20, d2.sd = 2, d12 = 12,
                   d12.sd = -1, lower = 0, upper = -1, nsweeps = 10,
         verbose = FALSE)
pipp.lookup(w=c(0, 1000, 0, 1000), pts=NULL,
                 n1=100, h, d, nsweeps=10,
                 verbose = FALSE, tor=FALSE)
pipp2.lookup(w=c(0, 1000, 0, 1000), pts1=NULL, pts2=NULL,
             n1=100, n2 = 100,
             h1, d1, h2, d2, h12, d12,
             nsweeps=10,
             fix=0,
             verbose = FALSE, tor=FALSE)
hlookup(h, d, r)
\method{plot}{sjedmin}(x, ...)
\method{plot}{sjedmin}(x, ...)
\method{plot}{sjedmin2}(x, ...)
\method{plot}{pipp}(x, ...)
\method{plot}{pipp2}(x, ...)
\method{plot}{sjebdmin}(x, ...)
}
\arguments{
  \item{wid,ht}{Width and height of the rectangular region to simulate.}
  \item{w}{A vector of length 4 giving: (xmin, xmax, ymin, ymax)}
  \item{p2}{Similar to w, but a matrix giving the bounding box.}
  \item{pts,pts1,pts2}{If not NULL, can be a matrix of N x 2 points to initialise
    the birth and death algorithm.}
  \item{npts}{Number of (2-d) points to generate.}
  \item{dmin,dminsd}{Mean and standard deviation of the dmin rule.}
  \item{lower,upper}{Lowest and highest possible values of dmin to
    generate from the distribution N(dmin,dminsd).  If values outside
    the low and high bounds are generated, they are rejected and another
    value created.  The lower bound should be set to at least the soma diameter
    if you want to take soma size into account.  If upper is negative,
    no upper bound on dmin values is imposed.}
  \item{quiet}{Set to true if you want no printout from routine.}
  \item{n1,n2}{Number of points of type 1 and type 2 to simulate.}
  \item{d1,d1.sd,d2,d2.sd}{Mean and s.d. of the dmin values for type 1
    and type 2 points.}
  \item{d12,d12.sd}{Mean and s.d. of minimal distance between type 1 and
    type 2 points (the heterotypic interaction).}
  \item{h,h1,h2,h12,d,}{Vector of h values at the corresponding
    distances d.  h1 is the type 1 interactions, h2 for type 2
    interactions and h12 for heterotypic interactions.  (Note, d1 and d2
    are used for distance vectors also in the PIPP routines, but their
    name classes with use of d1 and d2 for dmin routines.}
  \item{nsweeps}{Number of sweeps of the birth and death procedure.}
  \item{verbose}{Set to true if you want printout from routine.}
  \item{tor}{Set to TRUE for toroidal conditions on boundary.}
  \item{fix}{For the Bivariate PIPP, fix=0 means move all points; fix=1
    means keep type 1 points fixed; fix=2 means keep type 2 points fixed.}
  \item{x}{An object returned from one of the dmin functions, to be
    plotted.}
  \item{r}{TODO, used by hlookup}
  \item{...}{Extra arguments to pass to the plotting functions.}
}
\value{
  Res is a list with components:
  \item{xs,ys}{A vector of length npts.
    Each stores the x and y coordinates of each point}
  \item{dmins}{A vector of length npts. The value of each dmin created}
  \item{nrejects}{A vector of length npts.the number of rejects created
    when trying to position each cell.}
  \item{okay}{Boolean value saying whether a mosaic was successfully
    made or not.}
  \item{attempts}{The number of attempts made to create this dmin
    mosaic}.
  \item{note}{A short string describing what this mosaic is.}
}
\details{
  
  dminl and dminlul are the main functions for generating dmin mosaics.

  hlookup() returns the h function TODO.
}
\references{
  Lucia's papers?
}
\seealso{
   The spatial library from Brian Ripley, where some of the birth and
   death code was taken.
}
\examples{
wid <- 500; ht <- 500; ncells <- 400; dmin.mean <- 20; dmin.sd <- 4
dmin.lower <- 12; dmin.higher <- 28
w <- c( 100, 100+wid, 200, 200+ht)
par(mfcol=c(3,3))
d1 <- dminl(wid, ht, ncells, dmin.mean, dmin.sd)
plot(d1)
plot(d1$nrejects, main="rejects increase as function of cells placed",
xlab="cell number", ylab="# rejects")
hist(d1$dmins)
d2 <- dminlul(w, ncells, dmin.mean, dmin.sd, dmin.lower,dmin.higher)
plot(d2)
plot(d2$nrejects, xlab="cell number", ylab="# rejects")
hist(d2$dmins)
d3 <- dminlul(w+2000, ncells, dmin.mean, dmin.sd, dmin.lower,-1)
plot(d3)
plot(d3$nrejects, xlab="cell number", ylab="# rejects")
hist(d3$dmins)

## Generate dmin taking into consideration a set of cells already in the layer.
par(mfrow=c(1,1))
r1 <- 12; r2 <- 50 ## radii of cells
d <- dminlulfix2(dmin=2*r1, dminsd=0.01, d12=r1+r2, npts=800,
                 w=c(150, 900, 200, 1200),
                 p2=matrix(c (200,500, 800, 750, 700,500, 500, 1000),ncol=2))
plot.sjedmin2(d, r1=r1, r2=r2)

r1 <- 8; r2 <- 50;  d12 <- r1+r2; dmin <- r1*2
p2 <- matrix(c (200, 500, 630, 400,    700, 500, 800, 760),ncol=2)
d <- dminlulfix2(w = c( 100, 800, 400, 900),
                 npts = 700,
                 dmin = dmin, dminsd = 0.0001,
                 lower=12, upper=-1, d12=d12, p2=p2)
plot.sjedmin2(d, r1=r1, r2=r2)

## Bivariate simulation.
plot(bdmin.bd(n1=50, n2=30, d1=50, d2=30, d12=80, nsweeps=20))
}
\author{Stephen Eglen}
\keyword{spatial}
