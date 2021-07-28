lassoPlot_bin <- function(lassoOut, xIn, yIn, freqCut = 1, coefCut = 0.01, symbolMap = NULL, modelIndex = NULL, textWidth =0.8, rowHeight = 0.0005) {
  
  #for the barplot on the left of the heatmap
  if (is.null(modelIndex)) {
    #average all models
    barValue <- rowMeans(lassoOut$coefMat)
    freqValue <- rowMeans(lassoOut$coefMat !=0)
    barValue <- barValue[abs(barValue) > coefCut & freqValue >= freqCut] # a certain threshold
    barValue <- barValue[order(barValue)]
    if(length(barValue) == 0) {
      return(NA)
    }
  } else {
    #using a specificed model
    barValue <- lassoOut$coefMat[,modelIndex]
    barValue <- barValue[abs(barValue) > coefCut] # a certain threshold
    barValue <- barValue[order(barValue)]
    if(length(barValue) == 0) {
      return(NA)
    }
  }
  
  #for the heatmap and scatter plot below the heatmap
  mapData <- xIn
  scatterData <- yIn
  
  ord <- order(scatterData)
  scatterData <- scatterData[ord]
  mapData <- mapData[ord, names(barValue) ,drop=FALSE]
  idOrd <- rownames(mapData) #to record the order of the rows
  mapData <- mapData %>% as_tibble() %>% mutate(row = idOrd) %>%
    gather(key = "Var", value = "Value", -row) %>%
    mutate(row = factor(row, levels = idOrd),
           Var = factor(Var, levels = names(barValue)))
  #add tri12 status
  tri12meta <- dplyr::select(patMeta, PatientID, trisomy12) 
  tri12meta$PatientID <- as.factor(tri12meta$PatientID)
  colnames(mapData) <- c("PatientID", "Var", "Value")
  mapData_tri12 <- left_join(mapData, tri12meta, by =  "PatientID")
  
  #update labeling 
  mapData_tri12$Var <- revalue(mapData_tri12$Var, c("Resiquimod" = "Resiquimod", "TGF-b1" = "TGF-\u03B21", "sCD40L+IL-4" = "sCD40L + IL4")) 
  
  #plot the heatmap
  p1 <- ggplot(mapData_tri12, aes(x=PatientID, y=Var)) + geom_tile(aes(fill=Value),color = "gray") + 
    theme_bw()  + 
    theme(axis.text.x=element_blank(), axis.ticks=element_blank(), 
          panel.border=element_rect(colour="gainsboro"),  
          plot.title=element_text(size=18), 
          panel.background=element_blank(), 
          panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank()) + 
    xlab("Samples") + ylab("Log(Viability)") + 
    ggtitle("")+
    scale_fill_gradient2(low = palblues[1],high = palreds[8], mid = "white",guide="none")+ 
    facet_grid(~trisomy12, space = "free", scales = "free") + theme(strip.background =element_rect(fill="white"))
  
  #Plot the bar plot on the left of the heatmap
  if (is.null(symbolMap)) {
    barDF = data.frame(barValue=barValue, 
                       nm=factor(names(barValue),levels=names(barValue)))
  } else {
    #a id to symbol dataframe is provided
    barDF = data.frame(barValue = barValue, 
                       nm = symbolMap[names(barValue),1])
    barDF$nm <- factor(barDF$nm, levels = unique(barDF$nm))
  }
  
  p2 <- ggplot(data=barDF, aes(x=nm, y=barValue)) + 
    geom_bar(stat="identity", fill=darkergrey, colour="black", position = "identity", width=.66, size=0.2) + 
    theme_bw() + geom_hline(yintercept=0, size=0.3) + scale_x_discrete(expand=c(0,0.5)) + 
    scale_y_continuous(expand=c(0,0)) + 
    coord_flip(ylim=c(min(barValue),max(barValue))) + 
    theme(panel.grid.major=element_blank(), 
          panel.background=element_blank(), axis.ticks.y = element_blank(),
          panel.grid.minor = element_blank(), axis.text.y=element_text(size=8, hjust=0),
          axis.text.x = element_text(size=10),
          panel.border=element_blank()) +
    xlab("") + 
    ylab("Size of predictor") + 
    geom_vline(xintercept=c(0.5), color="black", size=0.6)
  
  
  #Scale bar for continuous variable
  
  Vgg = ggplot(mapData, aes(x=PatientID, y=Var, col = Value)) + 
    geom_point() + 
    scale_color_gradient2(low = palblues[1],high = palreds[8], mid = "white", name = "z-score") + 
    theme(legend.title=element_text(size=10), legend.text=element_text(size=10))
  
  #Assemble all the plots together
  
  # construct the gtable
  
  wdths = c(1.5, 0.2, 1.3*ncol(mapData), textWidth,0.8)
  hghts = c(0.2, rowHeight*nrow(mapData), 0.2, 0.5)
  gt = gtable(widths=unit(wdths, "in"), heights=unit(hghts, "in"))
  
  ## make grobs
  gg1 = ggplotGrob(p1)
  gg2 = ggplotGrob(p2)
  gg4 = ggplotGrob(Vgg)
  
  ## fill in the gtable
  gt = gtable_add_grob(gt, gtable_filter(gg2, "panel"), 2, 1) # add barplot
  gt = gtable_add_grob(gt, gtable_filter(gg2, "xlab-b"), 1, 1) # add  y axis label to barplot 
  gt = gtable_add_grob(gt, gtable_filter(gg1, "panel"), 2, 3) # add heatmap
  gt = gtable_add_grob(gt, gtable_filter(gg1, "strip-t"), 1, 3) # add facetstip to plot
  
  gt = gtable_add_grob(gt, gtable_filter(gg2, "axis-b"), 3, 1) # y axis for barplot
  gt = gtable_add_grob(gt, gtable_filter(gg2, "axis-l"), 2, 4) # variable names
  gt = gtable_add_grob(gt, gtable_filter(gg4, "guide-box"), 2, 5) # scale bar for continous variables
  
  
  #plot
  #grid.draw(gt)
  gt
}