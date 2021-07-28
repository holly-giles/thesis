########### Scaling and Centering of Viability Matrix ################

# medianCenter_MadScale: 
#function to  scale with MAD, center to merdian
medianCenter_MadScale05 <- function(x) {
  s <- median(x)
  #s=0
  (x - s) / mad(x, center = s)
}

# scaleCytResp function:  to apply medianCenter_MadScale row wise to viability matrix
scaleCytResp  <- function(x) t(apply(x, 1, medianCenter_MadScale)) 
