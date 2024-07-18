#Create clonal expansion barplot from VDJ dataframe

#clones = number of most expanded clones to plot
#color.by = VDJ column
#group.by = "none" or VDJ column

VDJ_clonal_expansion <- function(VDJ,
                                        clones,
                                        group.by,
                                        color.by,
                                        text.size,
                                        sub.type = F){
  #add 1 for plotting purposes
  VDJ$cell <- 1
  
  #if group.by is none, all cells will be plotted together
  if (group.by == "none"){group.by <- "cell"}
  
  #Seperate potential groups
  output <- list()
  for (i in unique(VDJ[,group.by])){
    group_VDJ <- VDJ[VDJ[,group.by]==i,]
    
    #Take most expanded clones and create a new VDJ matrix for this specific group
    clonotypes <- names(table(group_VDJ$clonotype_id)[order(table(group_VDJ$clonotype_id), decreasing = TRUE)][1:clones])
    group_VDJ <- group_VDJ[group_VDJ$clonotype_id %in% clonotypes,]
    
    #Replace NA values with "unknown"
    group_VDJ[is.na(group_VDJ[,color.by]),color.by] <- "Unknown"
    
    #create plots
    group_VDJ$cl <- group_VDJ[,color.by]
    if(color.by == "VDJ_cgene"){
      #Annotate as isotype
      group_VDJ$cl <- gsub("H", "", group_VDJ$cl)
      group_VDJ <- group_VDJ %>%
        mutate(cl = ifelse(cl != "Unknown",
                           paste0(substr(cl, 1, 1), tolower(substr(cl, 2, 2)), substr(cl, 3, nchar(cl))),cl))

      if(sub.type == T){
        cl_colors <- c("IgA" = "red",
                       "IgA1" = "#fb6a4a",
                       "IgA2" = "#de2d26",
                   "IgE" = "purple",
                   "IgG1" = "#a1d99b",
                   "IgG2" = "#74c476",
                   "IgG2B" = "#41ab5d",
                   "IgG2C" = "#238b45",
                   "IgG3" = "#005a32",
                   "IgM" = "black",
                   "IgD" = "blue",
                   "Unknown" = "grey")
      }
      if(sub.type == F){
        group_VDJ$cl <- gsub("IgA1", "IgA", group_VDJ$cl)
        group_VDJ$cl <- gsub("IgA2", "IgA", group_VDJ$cl)
        group_VDJ$cl <- gsub("IgG1", "IgG", group_VDJ$cl)
        group_VDJ$cl <- gsub("IgG2", "IgG", group_VDJ$cl)
        group_VDJ$cl <- gsub("IgG2B", "IgG", group_VDJ$cl)
        group_VDJ$cl <- gsub("IgGB", "IgG", group_VDJ$cl)
        group_VDJ$cl <- gsub("IgGC", "IgG", group_VDJ$cl)
        group_VDJ$cl <- gsub("IgG2C", "IgG", group_VDJ$cl)
        group_VDJ$cl <- gsub("IgG3", "IgG", group_VDJ$cl)

        cl_colors <- c("IgA" = "#fb6a4a",
                       "IgE" = "#B452CD",
                       "IgG" = "#74c476",
                       "IgM" = "black",
                       "IgD" = "#1874CD",
                       "Unknown" = "grey")
      }
    }else{
      unique_cl_values <- sort(unique(group_VDJ$cl))
      cl_colors <- setNames(Hue_Pal(length(unique_cl_values)), unique_cl_values)
    }
    
    p <- ggplot(group_VDJ, aes(x=factor(clonotype_id, levels = clonotypes), y=cell, fill=cl)) +
      geom_bar(, stat = "identity") +
      xlab("Clonal Rank") +
      ylab("Number of cells") +
      theme_minimal() +
      theme(text = element_text(size = text.size),
            #axis.text.x=element_blank(),
            #axis.ticks.x=element_blank(),
            # Hide panel borders and remove grid lines
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
            ) +
      scale_fill_manual(name = "Isotype", values = cl_colors) +
      #scale_x_discrete(labels=1:clones) +
      scale_x_discrete(labels=c(rep("",9),10,rep("",9),20,rep("",9),30)) +
      ggtitle(i)
    pdf(file = paste0("PLM_Likelihoods/figures/Figure1_Supplementary1/ClonalExpansion_",i,".pdf"),
        width = 8,
        height = 6)
    print(p)
    dev.off()
  }
}


