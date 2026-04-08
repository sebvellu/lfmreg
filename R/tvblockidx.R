tvblockidx <- function(lgth, mval) {
  lidx <- ceiling(1:lgth/lgth * mval)
  ridx <- lidx/mval
  tidx <- floor(ridx * lgth) #pmin(pmax(floor(ridx * lgth), 1), lgth)
  return(list(lidx = lidx, ridx = ridx, tidx = tidx))
}
