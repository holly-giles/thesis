
#*Hierarchical clustering function:* To provide cluster function for running `ConsensusClusterPlus`. Accepts `this_dist`, a dissimilarity structure as produced by `dist` and `k`, the number of clusters to assign. 


myclusterfunction = function(this_dist, k){
  #run hierarchical cluster analysis on dissimilarity structure this_dist  
  hc = hclust(this_dist)
  #cut cut tree into k groups 
  assignment = cutree(hc, k)
  return(assignment)
}

