
#*deckel function:*  to  set limits of scaling factor. Accepts a number `x`, and two numeric limits, `lower` and `upper`.

deckel <- function(x, lower = -Inf, upper = +Inf) ifelse(x<lower, lower, ifelse(x>upper, upper, x))


