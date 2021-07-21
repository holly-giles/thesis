
#*scaleCytResp function:*  to apply `medianCenter_MadScale` row wise to viability matrix. Accepts `x`, a matrix of log(viability) values. 

scaleCytResp  <- function(x) t(apply(x, 1, medianCenter_MadScale)) 


