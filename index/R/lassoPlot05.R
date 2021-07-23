lassoPlot05 <- function(lassoOut, geneMatrix, viabMatrix, freqCut = 1, coefCut = 0.01) {
  
  plotList <- list()
  
  for (seaName in names(lassoOut)) { 
    ###FOR THE BAR PLOT
    #extract coefMat for each stimulus
    barValue <- rowMeans(lassoOut[[seaName]]$coefMat)
    #extract proportion of repeats for which each coefficient is significant
    freqValue <- rowMeans(abs(sign(lassoOut[[seaName]]$coefMat)))
    #filter out coefficients that don't meet freqCut and coefCut thresholds
    barValue <- barValue[abs(barValue) >= coefCut & freqValue >= freqCut] 
    #arrange the bar values in numerical order
    barValue <- barValue[order(barValue)]
    #if there are no sig coefficients, don't plot
    if(length(barValue) == 0) {
      plotList[[seaName]] <- NA
      next
    }
    ###FOR THE HEATMAP AND SCATTER PLOT
    #get feature matrix and response matrix to plot
    allData <- geneMatrix
    viabValue <- unlist(viabMatrix[seaName,])
    
    #get feature matrix for sig coefficients only
    tabValue <- allData[, names(barValue),drop=FALSE]
    ord <- order(viabValue)
    viabValue <- viabValue[ord]
    tabValue <- tabValue[ord, ,drop=FALSE]
    sampleIDs <- rownames(tabValue)
    tabValue <- as_tibble(tabValue)
    tabValue$Sample <- sampleIDs
    
    
    #annotate mutations by mutation, methylation or IGHV
    matValue <- gather(tabValue, key = "Var",value = "Value", -Sample)
    matValue$Type <- "mut"
    
    #for Methylation Cluster
    matValue$Type[grep("Methylation",matValue$Var)] <- "meth"
    
    #for IGHV status
    matValue$Type[grep("IGHV",matValue$Var)] <- "ighv"
    
    #change the scale of the value so that IGHV, Methylation and Mutation do not overlap
    matValue[matValue$Type == "mut",]$Value = matValue[matValue$Type == "mut",]$Value + 10
    matValue[matValue$Type == "meth",]$Value = matValue[matValue$Type == "meth",]$Value + 20
    matValue[matValue$Type == "ighv",]$Value = matValue[matValue$Type == "ighv",]$Value + 30
    
    #change continuous to categorical
    matValue$Value <- factor(matValue$Value,levels = sort(unique(matValue$Value)))
    
    #arrange order of heatmap
    matValue$Var <- factor(matValue$Var, levels = names(barValue))
    matValue$Sample <- factor(matValue$Sample, levels = names(viabValue))
    
    #change labels if mutation is a Doehner mutation 
    matValue$Var <- revalue(matValue$Var, c("del11q" = "del(11q)", "del13q" = "del(13q)", "del17p" = "del(17p)", "trisomy12" = "trisomy 12")) 
    
    
    #plot the heatmap 
    #title
    if(seaName=="TGF-b1"){
      thetitle <- "TGF-\u03B21"     
    } else {
      if(seaName=='IL-4'){ 
        thetitle <- "IL4"  
      } else {
        thetitle <- seaName
      }
    }
    
    p1 <- ggplot(matValue, aes(x=Sample, y=Var)) + 
      geom_tile(aes(fill=Value), color = "white") + #ghost white
      theme_bw()+
      scale_y_discrete(expand=c(0,0)) + 
      theme(axis.title.y = element_text( size=14),
            axis.title.x = element_text( size=14),
            axis.text.x=element_text(hjust=0, size=11),
            axis.text.y=element_text(hjust=0.1, size=11),
            axis.ticks=element_blank(),
            panel.border=element_rect(colour="gainsboro"),  
            plot.title=element_text(face="bold", size = 18, margin = margin(t = -5, b = 1)), 
            panel.background=element_blank(),
            panel.grid.major=element_blank(), 
            panel.grid.minor=element_blank()) + 
      xlab("Mutation status for each patient") + 
      ylab("") + 
      scale_fill_manual(name="Mutated", 
                        values=c(`10`= offwhite,  #WT
                                 `11`="#373A36", #Mutant
                                 `20`= offwhite, #LP
                                 `20.5`= "#707372", #IP
                                 `21` = "#A8A99E", #HP
                                 `30` = offwhite, #IGHV-U
                                 `31` = "#707372"), #IGHV-M
                        guide="none") + 
      ggtitle(thetitle)
    
    
    #Plot the bar plot on the left of the heatmap 
    barDF = data.frame(barValue, nm=factor(names(barValue),levels=names(barValue)))
    
    p2 <- ggplot(data=barDF, aes(x=nm, y=barValue)) + 
      geom_bar(stat="identity", 
               fill=ifelse(barValue<0,
                           palblues[6],palreds[8]), 
               colour="black", 
               size=0.3) + 
      scale_x_discrete(expand=c(0,0.5)) + 
      scale_y_continuous(expand=c(0,0)) + 
      coord_flip(ylim=c(-0.3,0.35)) + #changed from min(barValue) and max(barValue)
      theme(panel.grid.major=element_blank(), 
            panel.background=element_blank(), 
            axis.ticks.y = element_blank(),
            panel.grid.minor = element_blank(), 
            axis.text=element_text(size=11, angle = 45, hjust = 1, vjust = 1),
            axis.title = element_text(size=14), 
            panel.border=element_blank()) +
      ylab("Size of predictor") + 
      geom_vline(xintercept=c(0.5), 
                 color="black", 
                 size=0.6)
    
    #Plot the scatter plot under the heatmap
    scatterDF = data.frame(X=factor(names(viabValue), 
                                    levels=names(viabValue)), 
                           Y=unlist(viabValue))
    
    p3 <- 
      ggplot(scatterDF, aes(x=X, y=Y)) + 
      geom_point(shape=21, 
                 fill="dimgrey", 
                 colour="#707372", #dark grey
                 size=1.2) + 
      theme_bw() +
      theme(panel.grid.minor=element_blank(), 
            panel.grid.major.x=element_blank(), 
            #axis.title.x=element_blank(), 
            axis.ticks.x=element_blank(), 
            axis.text.y=element_text(size=11), 
            axis.title.x=element_text(size=14), 
            panel.border=element_rect(colour="dimgrey", size=0.1),
            panel.background=element_rect(fill="white")) +
      xlab("Log(Viability) for each sample")
    
    
    #Assemble all the plots togehter
    # construct the gtable
    wdths = c(0.2, 1.5, 0.2, 1.3*ncol(matValue), 1.5, 0.2)
    hghts = c(0.2, 0.2, 0.0020*nrow(matValue), 0.3, 1, 0.2)
    gt = gtable(widths=unit(wdths, "in"), heights=unit(hghts, "in"))
    
    ## make grobs
    gg1 = ggplotGrob(p1)
    gg2 = ggplotGrob(p2)
    gg3 = ggplotGrob(p3)
    ## fill in the gtable
    
    #HEATMAP
    #5:1 = "PREDICTORS"
    gt = gtable_add_grob(gt, gtable_filter(gg1, "panel"), 3, 4) # add heatmap
    gt = gtable_add_grob(gt, gtable_filter(gg1, "panel"), 3, 4) #add legend
    gt = gtable_add_grob(gt, gtable_filter(gg1, "title"), 1, 4) #add title to plot
    gt = gtable_add_grob(gt, gtable_filter(gg1, "axis-l"), 3, 5) # variable names
    gt = gtable_add_grob(gt, gtable_filter(gg1, "xlab-b"), 2, 4) # axis title
    
    #BARPLOT
    gt = gtable_add_grob(gt, gtable_filter(gg2, "panel"), 3, 2) # add barplot
    gt = gtable_add_grob(gt, gtable_filter(gg2, "axis-b"), 4, 2) # y axis for barplot
    gt = gtable_add_grob(gt, gtable_filter(gg2, "xlab-b"), 2, 2) # y lab for barplot
    
    #SCATTER PLOT
    gt = gtable_add_grob(gt, gtable_filter(gg3, "panel"), 5, 4) # add scatterplot
    gt = gtable_add_grob(gt, gtable_filter(gg3, "xlab-b"), 6, 4) # x label for scatter plot
    gt = gtable_add_grob(gt, gtable_filter(gg3, "axis-l"), 5, 3) #  axis for scatter plot
    
    
    
    #plot
    #grid.draw(gt)
    plotList[[seaName]] <- gt
  }
  return(plotList)
}