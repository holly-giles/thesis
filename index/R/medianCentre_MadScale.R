
#*medianCenter_MadScale function:*  to scale viability values according to MAD and then center values at zero. Maximum/minimum size of scaling factor set with deckel (above). Accepts a vector `x` to scale.


medianCenter_MadScale <- function(x) {
  s=0
  (x - s) / deckel(mad(x, center = s), lower = 0.05, upper = 0.2)
}



