\name{sjedmin}
\alias{dminl}
\alias{dminlul}
\alias{plot.sjedmin}
\title{Dmin mosaics}
\description{
  Create mosaics that respect the dmin rule.
}
\usage{
dminl(wid, ht, npts, dmin, dminsd)
dminlul(wid, ht, npts, dmin, dminsd, lower, upper)
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
}
\author{Stephen Eglen}
\keyword{spatial}
