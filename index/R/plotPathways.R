
plotPathways <-function(dat) 
{
  #Set pathway, number of drugs and group
  Pathway = NULL
  No = NULL
  Group = NULL
  #set ordering of drug cateogries (FDA Approved and  Clinical development/tool compound)
  ordM = sort(table(dat$group), decreasing = TRUE)
  
  #set ordering of pathways within categories
  ordS = tapply(dat$target_category, dat$group, function(pth) {
    sort(table(pth), decreasing = TRUE)
  })
  
  #set overall ordering
  ocur = ordS[names(ordM)]
  
  #get a dataframe of number of drugs for each pathway, split by category
  tmp = do.call(rbind, lapply(names(ocur), function(pathgroup) {
    data.frame(Group = pathgroup, Pathway = names(ocur[[pathgroup]]), 
               No = as.vector(unname(ocur[[pathgroup]])))
  }))
  
  #relevel categories
  tmp$Group = factor(tmp$Group, levels = rev(names(ocur)))
  
  #set plotting order
  lev = sort(tapply(tmp$No, tmp$Pathway, function(x) sum(x)), 
             decreasing = TRUE)
  
  #relevel pathways
  tmp$Pathway = factor(tmp$Pathway, levels = names(lev))
  
  #set height of y axis scale
  widthmax = max(lev) + 1
  
  #make plot
  g = ggplot(tmp, aes(x = Pathway, y = No, fill = Group)) + 
    geom_bar(width = 0.6, stat = "identity") + theme_bw() + 
    scale_fill_manual(values = c(drugpal[1], drugpal[2]), name = "Drug type") + 
    scale_y_continuous(breaks = seq(0, 20, 2), 
                       expand = c(0,0), limits = c(0, widthmax)) + 
    xlab("") + 
    ylab("No. of drugs") + 
    t1 + 
    theme(legend.position = c(0.7, 0.65))
  
}

