\name{sjedmin}
\alias{dminl}
\alias{dminlul}
\alias{plot.sjedmin}
\title{Dmin mosaics}
\description{
  Create mosaics that respect the dmin rule.
}
\usage{
dminl(wid, ht, npts, dmin, dminsd, quiet)
dminlul(wid, ht, npts, dmin, dminsd, lower, upper, quiet)
plot.sjedmin(res)
}
\arguments{
  \item{wid,ht}{Width and height of the rectangular region to simulate.}
  \item{npts}{Number of (2-d) points to generate.}
  \item{dmin,dminsd}{Mean and standard deviation of the dmin rule.}
  \item{lower,upper}{Lowest and highest possible values of dmin to
    generate from the distribution N(dmin,dminsd).  If values outside
    the low and high bounds are generated, they are rejected and another
    value created.  The lower bound should be set to at least the soma diameter
    if you want to take soma size into account.  If upper is negative,
    no upper bound on dmin values is imposed.}
  \item{quiet}{Set to true if you want no printout from routine.}
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
par(mfcol=c(3,3))
d1 <- dminl(wid, ht, ncells, dmin.mean, dmin.sd)
plot(d1)
plot(d1$nrejects, main="rejects increase as function of cells placed",
xlab="cell number", ylab="# rejects")
hist(d1$dmins)
d2 <- dminlul(wid, ht, ncells, dmin.mean, dmin.sd, dmin.lower,dmin.higher)
plot(d2)
plot(d2$nrejects, xlab="cell number", ylab="# rejects")
hist(d2$dmins)
d3 <- dminlul(wid, ht, ncells, dmin.mean, dmin.sd, dmin.lower,-1)
plot(d3)
plot(d3$nrejects, xlab="cell number", ylab="# rejects")
hist(d3$dmins)

## Generate dmin taking into consideration a set of cells already in the layer.
par(mfrow=c(1,1))
r1 <- 12; r2 <- 50 ## radii of cells
d <- dminlulfix2(dmin=2*r1, dminsd=0.01, d12=r1+r2, npts=800,
                 w=c(150, 900, 200, 1200),
                 p2=matrix(c (200,500, 800, 750, 700,500, 500, 1000),ncol=2))
plot(d, r1=r1, r2=r2)

r1 <- 8; r2 <- 50;  d12 <- r1+r2; dmin <- r1*2
p2 <- matrix(c (200, 500, 630, 400,    700, 500, 800, 760),ncol=2)
d <- dminlulfix2(w = c( 100, 800, 400, 900),
                 npts = 700,
                 dmin = dmin, dminsd = 0.0001,
                 lower=12, upper=-1, d12=d12, p2=p2)
plot(d, r1=r1, r2=r2)

}
\author{Stephen Eglen}
\keyword{spatial}
