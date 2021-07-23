makelegends <- function (legendFor, colors) {
  x = NULL
  y = NULL
  colors = colors[names(colors) %in% legendFor]
  nleg = length(colors)
  
  #set up gtable
  wdths = c(0.4,2,2,2,1.5)
  hghts = c(2)
  
  gtl = gtable(widths=unit(wdths, "in"), heights=unit(hghts, "in"))
  n = 2
  
  #legend for Methylation Cluster
  if ("M" %in% names(colors)) {
    Mgg = ggplot(data = data.frame(x = 1, 
                                   y = factor(c("LP", "IP", "HP"), 
                                              levels = c("LP", "IP", "HP"))), 
                 aes(x = x, y = y, fill = y)) + 
      geom_tile() + 
      scale_fill_manual(name = "Methylation cluster", 
                        values = setNames(colors[["M"]], 
                                          nm = c("LP", "IP","HP"))) + 
      theme(legend.title = element_text(size = 12), 
            legend.text = element_text(size = 12))
    
    gtl = gtable_add_grob(gtl, gtable_filter(ggplotGrob(Mgg), "guide-box"), 1, n)
    n = n + 1
  }
  
  #Legend for IGHV status
  if ("I" %in% names(colors)) {
    Igg = ggplot(data = data.frame(x = 1, y = factor(c("Unmutated", 
                                                       "Mutated"), levels = c("Unmutated", "Mutated"))), 
                 aes(x = x, y = y, fill = y)) + geom_tile() + scale_fill_manual(name = "IGHV", 
                                                                                values = setNames(colors[["I"]], nm = c("Unmutated", "Mutated"))) + 
      theme(legend.title = element_text(size = 12), 
            legend.text = element_text(size = 12))
    gtl = gtable_add_grob(gtl, gtable_filter(ggplotGrob(Igg), 
                                             "guide-box"), 1, n)
    n = n + 1
  }
  
  #legend for gene mutations
  if ("G" %in% names(colors)) {
    Ggg = ggplot(data = data.frame(x = 1, y = factor(c("Wild Type", 
                                                       "Mutated"), levels = c("Wild Type", "Mutated"))), 
                 aes(x = x, y = y, fill = y)) + geom_tile() + scale_fill_manual(name = "Gene", 
                                                                                values = setNames(colors[["G"]], nm = c("Wild Type", "Mutated"))) + 
      theme(legend.title = element_text(size = 12), 
            legend.text = element_text(size = 12))
    gtl = gtable_add_grob(gtl, gtable_filter(ggplotGrob(Ggg), 
                                             "guide-box"), 1, n)
    n = n + 1
  }
  
  return(list(plot = gtl, width = sum(wdths), height = sum(hghts)))
}
