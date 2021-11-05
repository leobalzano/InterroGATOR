#############################################
#####     InterroGATOR_functions.R     ######
#############################################
# Author: Leandro Balzano-Nogueira
# Diabetes Institute, University of Florida (Gainesville)
# Last update: November/5/2021

# This script is to gather together the functions used for the InterroGATOR shiny Application.
#############################################
# Functions:

S4_to_dataframe <- function(s4obj) {
  nms <- slotNames(s4obj)
  
  lst <- lapply(nms, function(nm) slot(s4obj, nm))
  as.data.frame(setNames(lst, nms))
}

PlotExprCODEXspatial_BB <-function (stvea_object, name, type = "protein", high_color = "red", 
                                    high_color2 = "green", low_color = "white", pt_size = 0.8) {
  
  if (is.null(stvea_object@codex_spatial)) {
    stop("stvea_object does not contain CODEX spatial information")
  }
  if (type == "protein") {
    if (!is.null(stvea_object@codex_clean)) {
      plotting_data <- stvea_object@codex_clean
    }
    else if (!is.null(stvea_object@codex_protein)) {
      plotting_data <- stvea_object@codex_protein
    }
    else {
      stop("stvea_object must contain CODEX protein data with type=\"protein\"")
    }
  }
  else if (type == "RNA") {
    plotting_data <- stvea_object@codex_mRNA
  }
  else {
    stop("type must be either \"RNA\" or \"protein\"", call. = FALSE)
  }
  if (length(name) > 2) {
    stop("name must be at most length 2", call. = FALSE)
  }
  rbPal1 <- colorRampPalette(c(alpha(low_color, 0), alpha(high_color, 
                                                          1)), alpha = TRUE)
  color <- rbPal1(100)[as.numeric(cut(plotting_data[, name[1]], 
                                      breaks = 100))]
  subtitle <- paste("Expression of ", name[1], " (", high_color, 
                    ")", sep = "")
  if (length(name) == 2) {
    rbPal2 <- colorRampPalette(c(alpha(low_color, 0), alpha(high_color2, 
                                                            1)), alpha = TRUE)
    color2 <- rbPal2(100)[as.numeric(cut(plotting_data[, 
                                                       name[2]], breaks = 100))]
    color <- sapply(1:length(color), function(m) colorRampPalette(c(color[m], 
                                                                    color2[m]), alpha = TRUE)(3)[2])
    subtitle <- paste(subtitle, " and ", name[2], " (", high_color2, 
                      ")", sep = "")
  }
  x_tmp <- stvea_object@codex_spatial[, "x"]
  x_tmp <- x_tmp - min(x_tmp)
  y_tmp <- stvea_object@codex_spatial[, "y"]
  y_tmp <- y_tmp - min(y_tmp)
  spatial_tmp <- as.data.frame(cbind(x = x_tmp, y = y_tmp))
  ggplot(spatial_tmp, aes(x = x, y = y, color = factor(1:length(color)))) + 
    geom_point(size = 0.5, alpha = 0.5) + scale_color_manual(values = alpha(color, 1)) + guides(color = FALSE) + ylim(max(y_tmp), 0) + labs(title = paste("Spatial", type, "expression"), subtitle = subtitle, ) + theme_void() + 
    coord_fixed() + theme(plot.background = element_rect(fill = "black"), plot.title= element_text(colour = "white"),plot.subtitle = element_text(colour="white"))
}

minmaxnorm <- function(x, na.rm = TRUE) {
  return((x- min(x)) /(max(x)-min(x)))
}

PlotClusterCITE_Independently<- function(stvea_object, color_by,pt_size, highlight, Color="dodgerblue2", Selection_on_Top="no") {
  if (any(color_by==0)) {
    color_by<-as.factor(ifelse(as.numeric(as.character(color_by))<0,-1, as.numeric(as.character(color_by)) +1) )
    
  }
  
  if (-1 %in% color_by) {
    colors <- c("gray", colorspace::rainbow_hcl(length(unique(color_by)) - 
                                                  1, c = 80))
  } else {
    colors <- colorspace::rainbow_hcl(length(unique(color_by)), 
                                      c = 80)
  }
  
  fill2<- ifelse(color_by %in% highlight, Color ,"gray95")
  colors2<- ifelse(color_by %in% highlight, Color ,"gray75")
  colors3<-unique(colors2)
  
  if (Selection_on_Top=="no") {
    ggplot(stvea_object@cite_emb, aes_string(x = colnames(stvea_object@cite_emb)[1], 
                                             y = colnames(stvea_object@cite_emb)[2], 
                                             color = factor(fill2)
    )) + 
      geom_point(fill=fill2,colour=colors2,size = pt_size,pch=21) + 
      labs(title = paste("Group",highlight)) +
      theme_void()
  } else {
    datita<-cbind(stvea_object@cite_emb,fill2=fill2, colors2=colors2)
    datita<-datita %>% arrange(factor(fill2, levels = c("gray95",Color)) )
    print("head(datita)" )
    print(head(datita) )
    print("tail(datita)" )
    print(tail(datita) )
    ggplot(datita, aes_string(x = colnames(datita)[1], 
                                             y = colnames(datita)[2], 
                                             color = factor(datita[,3])
    )) + 
      geom_point(fill=datita[,3],colour=datita[,4],size = pt_size,pch=21) + 
      labs(title = paste("Group",highlight)) +
      theme_void()
  }
}


PlotClusterCITE_v2<-function (stvea_object,color_by, pt_size = 0.5) {
  if (any(color_by==0)) {
    color_by<-as.factor(ifelse(as.numeric(as.character(color_by))<0,-1, as.numeric(as.character(color_by)) +1) )
    
  }
  
  if (-1 %in% color_by) {
    colors <- c("gray", colorspace::rainbow_hcl(length(unique(color_by)) - 
                                                  1, c = 80))
  }
  else {
    colors <- colorspace::rainbow_hcl(length(unique(color_by)), 
                                      c = 80)
  }
  ggplot(stvea_object@cite_emb, aes_string(x = colnames(stvea_object@cite_emb)[1], 
                                           y = colnames(stvea_object@cite_emb)[2], 
                                           color = factor(color_by))) + 
    #geom_point(size = pt_size) + 
    geom_point(aes(fill=factor(color_by)),colour="gray5",size = pt_size,pch=21) + 
    scale_color_manual(values = colors, name = "cluster") + 
    #guides(colour = guide_legend(override.aes = list(size = 5))) + 
    guides(fill=guide_legend(title="cluster")) +
    guides(colour = guide_legend(override.aes = list(size = 5))) + 
    #geom_text(aes(label=factor(color_by)), size=3) +
    theme_void()
}

PlotExprCODEXspatial_v2<-function (stvea_object, name, type = "protein", high_color = "red", 
                                   high_color2 = "green", low_color = "white", pt_size = 0.8,Subset="yes", SmoothScatterplot=FALSE, mult=1, maxPercentile=0.9) {
  alfa=1
  if (is.null(stvea_object@codex_spatial)) {
    stop("stvea_object does not contain CODEX spatial information")
  }
  
  if (tolower(Subset) == "yes" ) {
    if (type == "protein") {
      if (!is.null(stvea_object@codex_clean)) {
        plotting_data <- stvea_object@codex_clean
      }
      else if (!is.null(stvea_object@codex_protein)) {
        plotting_data <- stvea_object@codex_protein
      }
      else {
        stop("stvea_object must contain CODEX protein data with type=\"protein\"")
      }
    }
    else if (type == "RNA") {
      plotting_data <- stvea_object@codex_mRNA
    }
    else {
      stop("type must be either \"RNA\" or \"protein\"", call. = FALSE)
    }
    if (length(name) > 2) {
      stop("name must be at most length 2", call. = FALSE)
    }
    rbPal1 <- colorRampPalette(c(alpha(low_color, 0), alpha(high_color, 
                                                            1)), alpha = TRUE)
    color <- rbPal1(100)[as.numeric(cut(plotting_data[, name[1]], 
                                        breaks = 100))]
    subtitle <- paste("Expression of ", name[1], " (", high_color, 
                      ")", sep = "")
    if (length(name) == 2) {
      rbPal2 <- colorRampPalette(c(alpha(low_color, 0), alpha(high_color2, 
                                                              1)), alpha = TRUE)
      color2 <- rbPal2(100)[as.numeric(cut(plotting_data[, 
                                                         name[2]], breaks = 100))]
      color <- sapply(1:length(color), function(m) colorRampPalette(c(color[m], 
                                                                      color2[m]), alpha = TRUE)(3)[2])
      subtitle <- paste(subtitle, " and ", name[2], " (", high_color2, 
                        ")", sep = "")
    }
    x_tmp <- stvea_object@codex_spatial[, "x"]
    x_tmp <- x_tmp - min(x_tmp)
    x_tmpMAX<-max(x_tmp)
    y_tmp <- stvea_object@codex_spatial[, "y"]
    y_tmp <- y_tmp - min(y_tmp)
    y_tmpMAX<-max(y_tmp)
    spatial_tmp <- as.data.frame(cbind(x = x_tmp, y = y_tmp,col=color))
    
    
    print(paste0("Maximum value in X axis is ",x_tmpMAX," while in Y axis is ",y_tmpMAX))
    
    XMINsubset <- suppressWarnings({as.numeric(as.character(readline(prompt="Subset X axis from: ") ))})
    
    while (is.na(XMINsubset)){
      print ("Please enter a valid number for the minimum X axis value to plot")
      XMINsubset<-suppressWarnings({as.numeric(as.character(readline(prompt="Subset X axis from: ") ))})
    };print(paste0("X axis will be plotted from pixel number ",XMINsubset))
    
    XMAXsubset <- suppressWarnings({as.numeric(as.character(readline(prompt="Subset X axis until: ") ))})
    
    while (is.na(XMAXsubset)){
      print ("Please enter a valid number for the maximum X axis value to plot")
      XMAXsubset<-suppressWarnings({as.numeric(as.character(readline(prompt="Subset X axis until: ") ))})
    };print(paste0("X axis will be plotted until pixel number ", XMAXsubset))
    
    
    YMINsubset <- suppressWarnings({as.numeric(as.character(readline(prompt="Subset Y axis from: ") ))})
    
    while (is.na(YMINsubset)){
      print ("Please enter a valid number for the minimum Y axis value to plot")
      YMINsubset<-suppressWarnings({as.numeric(as.character(readline(prompt="Subset Y axis from: ") ))})
    };print(paste0("Y axis will be plotted from pixel number ",YMINsubset))
    
    YMAXsubset <- suppressWarnings({as.numeric(as.character(readline(prompt="Subset Y axis until: ") ))})
    
    while (is.na(YMAXsubset)){
      print ("Please enter a valid number for the maximum Y axis value to plot")
      YMAXsubset<-suppressWarnings({as.numeric(as.character(readline(prompt="Subset Y axis until: ") ))})
    };print(paste0("Y axis will be plotted until pixel number ", YMAXsubset))
    
    
    spatial_tmp_subset<- spatial_tmp[as.numeric(as.character(spatial_tmp[,"x"]))>=XMINsubset &
                                       as.numeric(as.character(spatial_tmp[,"x"]))<=XMAXsubset &
                                       as.numeric(as.character(spatial_tmp[,"y"]))>=YMINsubset &
                                       as.numeric(as.character(spatial_tmp[,"y"]))<=YMAXsubset,] # 
    
    if (nrow(spatial_tmp_subset) <1){
      print("The area selected is too small. Please select a larger area to plot")
    } else {
      print(paste0("The area selected to plot goes from (", XMINsubset," , ",YMINsubset, ") to (",XMAXsubset," , ",YMAXsubset,")") )
    }
    
    ###
    if (SmoothScatterplot) {
      
      ggplot(spatial_tmp_subset) + aes(x = as.numeric(as.character(spatial_tmp_subset$x)), y = as.numeric(as.character(spatial_tmp_subset$y)), color = factor(1:length(spatial_tmp_subset$col))) + 
        geom_point(size = pt_size, alpha = 0.5) + scale_color_manual(values = alpha(spatial_tmp_subset$col, 1)) + guides(color = FALSE) + 
        ylim(YMAXsubset, YMINsubset) + xlim(XMINsubset,XMAXsubset) +
        labs(title = paste("Spatial",  type, "expression"), subtitle = subtitle) + theme_void() + 
        coord_fixed()  
    } else{
      
      ggplot(spatial_tmp_subset, aes(x = as.numeric(as.character(spatial_tmp_subset$x)), y = as.numeric(as.character(spatial_tmp_subset$y)), color = factor(1:length(spatial_tmp_subset$col)))) + 
        geom_point(size = pt_size, alpha = 0.5) + scale_color_manual(values = alpha(spatial_tmp_subset$col, 1)) + guides(color = FALSE) + ylim(max(y_tmp), 0) + xlim(0,max(x_tmp)) +
        labs(title = paste("Spatial",  type, "expression"), subtitle = subtitle) + theme_void() + 
        coord_fixed()
    }
    
    
  }  else {
    ###
    if (type == "protein") {
      if (!is.null(stvea_object@codex_clean)) {
        plotting_data <- stvea_object@codex_clean
      }
      else if (!is.null(stvea_object@codex_protein)) {
        plotting_data <- stvea_object@codex_protein
      }
      else {
        stop("stvea_object must contain CODEX protein data with type=\"protein\"")
      }
    }
    else if (type == "RNA") {
      plotting_data <- stvea_object@codex_mRNA
    }
    else {
      stop("type must be either \"RNA\" or \"protein\"", call. = FALSE)
    }
    if (length(name) > 2) {
      stop("name must be at most length 2", call. = FALSE)
    }
    rbPal1 <- colorRampPalette(c(alpha(low_color, 0), alpha(high_color, 
                                                            1)), alpha = TRUE)
    color <- rbPal1(100)[as.numeric(cut(plotting_data[, name[1]], 
                                        breaks = 100))]
    subtitle <- paste("Expression of ", name[1], " (", high_color, 
                      ")", sep = "")
    if (length(name) == 2) {
      rbPal2 <- colorRampPalette(c(alpha(low_color, 0), alpha(high_color2, 
                                                              1)), alpha = TRUE)
      color2 <- rbPal2(100)[as.numeric(cut(plotting_data[, 
                                                         name[2]], breaks = 100))]
      color <- sapply(1:length(color), function(m) colorRampPalette(c(color[m], 
                                                                      color2[m]), alpha = TRUE)(3)[2])
      subtitle <- paste(subtitle, " and ", name[2], " (", high_color2, 
                        ")", sep = "")
    }
    
    x_tmp <- stvea_object@codex_spatial[, "x"]
    x_tmp <- x_tmp - min(x_tmp)
    y_tmp <- stvea_object@codex_spatial[, "y"]
    y_tmp <- y_tmp - min(y_tmp)
    spatial_tmp <- as.data.frame(cbind(x = x_tmp, y = y_tmp))
    dim(spatial_tmp)
    
    if (low_color == "transparent"){
      sort(table(color))
      
      
      losblancos<- c("#FFFFFF00")
      
      
      pd<-as.data.frame(plotting_data[,name,drop=F]) 
      summary(pd)
      pd <- as.data.frame(lapply(pd, minmaxnorm))
      summary(pd)
      pd$discriminant<-ifelse(rowSums(pd)==0,0,1)
      spatial_tmp$discriminant<-pd$discriminant
      pd<-pd[pd$discriminant == 1, ]
      spatial_tmp<-spatial_tmp[spatial_tmp$discriminant == 1, ]
      summary(pd)
      summary(spatial_tmp)
      
      summary(pd[,1])
      themax<-quantile(pd[,1], maxPercentile) 
      pd[,1][pd[,1]>themax] <- themax
      
      rbPal1 <- colorRampPalette(c(alpha(low_color, 0), alpha(high_color, 1)), alpha = TRUE)
      thecut<-as.numeric(cut(pd[, name[1]], breaks = 100))
      table(thecut)
      
      if (mult == 1){
        colfac<-round(mean(thecut),digits = 0)  
      } else {
        colfac<-mult
      }

      prenewcolfac<-thecut+colfac
      newcolfac<-ifelse( prenewcolfac >100,100,prenewcolfac  )
      
      color <- rbPal1(100)[newcolfac]
      
      subtitle <- paste("Expression of ", name[1], " (", high_color,")", sep = "")
      if (length(name) == 2) {
        rbPal2 <- colorRampPalette(c(alpha(low_color, 0), alpha(high_color2, 1)), alpha = TRUE)
        thecut2<-as.numeric(cut(pd[, name[2]], breaks = 100))
        table(thecut2)
        colfac2<-round(mean(thecut2),digits = 0)
        prenewcolfac2<-thecut2+colfac2
        newcolfac2<-ifelse( prenewcolfac2 >100,100,prenewcolfac2  )
        
        color2 <- rbPal2(100)[newcolfac2]
        
        
        color <- sapply(1:length(color), function(m) colorRampPalette(c(color[m], color2[m]), alpha = TRUE)(3)[2])
        subtitle <- paste(subtitle, " and ", name[2], " (", high_color2, ")", sep = "")
      }
      alfa<-ifelse(color %in% losblancos, 0,0.5)
      
      
      
    }
    
    ggplot(spatial_tmp, aes(x = x, y = y, color = factor(1:length(color)))) + 
      geom_point(size = pt_size, alpha = alfa) + scale_color_manual(values = alpha(color, 
                                                                                   1)) + guides(color = FALSE) + ylim(max(y_tmp), 0) + labs(title = paste("Spatial", 
                                                                                                                                                          type, "expression"), subtitle = subtitle) + theme_void() + 
      coord_fixed()
  }
}

PlotExprCODEXemb_v2<- function (stvea_object, name, type = "protein", high_color = "red", 
                                high_color2 = "green", low_color = "light gray", mult=1, 
                                pt_size = 0.8) {
  print_type <- type
  alfa=1
  if (type == "protein") {
    if (!is.null(stvea_object@codex_clean)) {
      plotting_data <- stvea_object@codex_clean
    }
    else if (!is.null(stvea_object@codex_protein)) {
      plotting_data <- stvea_object@codex_protein
    }
    else {
      stop("stvea_object must contain CODEX protein data with type=\"protein\"")
    }
    print_type <- "Protein"
  }
  else if (type == "RNA") {
    plotting_data <- stvea_object@codex_mRNA
  }
  else {
    stop("type must be either \"RNA\" or \"protein\"", call. = FALSE)
  }
  if (length(name) > 2) {
    stop("name must be at most length 2", call. = FALSE)
  }
  rbPal1 <- colorRampPalette(c(alpha(low_color, 0), alpha(high_color, 
                                                          1)), alpha = TRUE)
  color <- rbPal1(100)[as.numeric(cut(plotting_data[, name[1]], 
                                      breaks = 100))]
  subtitle <- paste("Expression of ", name[1], " (", high_color, 
                    ")", sep = "")
  if (length(name) == 2) {
    rbPal2 <- colorRampPalette(c(alpha(low_color, 0), alpha(high_color2, 
                                                            1)), alpha = TRUE)
    color2 <- rbPal2(100)[as.numeric(cut(plotting_data[, 
                                                       name[2]], breaks = 100))]
    color <- sapply(1:length(color), function(m) colorRampPalette(c(color[m], 
                                                                    color2[m]), alpha = TRUE)(3)[2])
    subtitle <- paste(subtitle, " and ", name[2], " (", high_color2, 
                      ")", sep = "")
  }
  
  if (low_color == "light gray"){
    sort(table(color))
    
    losblancos<- c("#D3D0D002")
    
    pd<-as.data.frame(plotting_data[,name,drop=F])
    summary(pd)
    pd <- as.data.frame(lapply(pd, minmaxnorm))
    summary(pd)
    pd$discriminant<-ifelse(rowSums(pd)==0,0,1)
    
    rbPal1 <- colorRampPalette(c(alpha(low_color, 0), alpha(high_color, 1)), alpha = TRUE)
    thecut<-as.numeric(cut(pd[, name[1]], breaks = 100))
    table(thecut)
    
    if (mult == 1){
      colfac<-round(mean(thecut),digits = 0)  
    } else {
      colfac<-mult
    }
    prenewcolfac<-thecut+colfac
    newcolfac<-ifelse( prenewcolfac >100,100,prenewcolfac  )
    
    color <- rbPal1(100)[newcolfac]
    subtitle <- paste("Expression of ", name[1], " (", high_color,")", sep = "")
    if (length(name) == 2) {
      rbPal2 <- colorRampPalette(c(alpha(low_color, 0), alpha(high_color2, 1)), alpha = TRUE)
      thecut2<-as.numeric(cut(pd[, name[2]], breaks = 100))
      table(thecut2)
      colfac2<-round(mean(thecut2),digits = 0)
      prenewcolfac2<-thecut2+colfac2
      newcolfac2<-ifelse( prenewcolfac2 >100,100,prenewcolfac2  )
      
      color2 <- rbPal2(100)[newcolfac2]
      
      
      color <- sapply(1:length(color), function(m) colorRampPalette(c(color[m], color2[m]), alpha = TRUE)(3)[2])
      subtitle <- paste(subtitle, " and ", name[2], " (", high_color2, ")", sep = "")
    }
    alfa<-ifelse(color %in% losblancos, 0.1,1)
  }
  
  ggplot(stvea_object@codex_emb, aes_string(x = colnames(stvea_object@codex_emb)[1], 
                                            y = colnames(stvea_object@codex_emb)[2], color = factor(1:length(color)))) + 
    geom_point(size = pt_size,alpha = alfa) + labs(title = paste(print_type, 
                                                                 "expression"), subtitle = subtitle) +      
    scale_color_manual(values = alpha(color, alfa), guide = FALSE) +
    theme_void()
}

PlotExprCITE_v2<- function (stvea_object, name, type = "RNA", high_color = "red", 
                            high_color2 = "green", low_color = "light gray", mult=1,
                            pt_size = 0.8) {
  print_type <- type
  alfa=1
  if (type == "RNA") {
    plotting_data <- stvea_object@cite_mRNA_norm
  }
  else if (type == "protein") {
    plotting_data <- stvea_object@codex_clean
    print_type <- "Protein"
  }
  else {
    stop("type must be either \"RNA\" or \"protein\"", call. = FALSE)
  }
  if (length(name) > 2) {
    stop("name must be at most length 2", call. = FALSE)
  }
  rbPal1 <- colorRampPalette(c(alpha(low_color, 0), alpha(high_color, 
                                                          1)), alpha = TRUE)
  color <- rbPal1(100)[as.numeric(cut(plotting_data[, name[1]], 
                                      breaks = 100))]
  subtitle <- paste("Expression of ", name[1], " (", high_color, 
                    ")", sep = "")
  if (length(name) == 2) {
    rbPal2 <- colorRampPalette(c(alpha(low_color, 0), alpha(high_color2, 
                                                            1)), alpha = TRUE)
    color2 <- rbPal2(100)[as.numeric(cut(plotting_data[, 
                                                       name[2]], breaks = 100))]
    color <- sapply(1:length(color), function(m) colorRampPalette(c(color[m], 
                                                                    color2[m]), alpha = TRUE)(3)[2])
    subtitle <- paste(subtitle, " and ", name[2], " (", high_color2, 
                      ")", sep = "")
  }
  if (low_color == "light gray"){
    sort(table(color))
    
    losblancos<- c("#D3D0D002")
    pd<-as.data.frame(plotting_data[,name,drop=F])
    summary(pd)
    pd <- as.data.frame(lapply(pd, minmaxnorm))
    summary(pd)
    pd$discriminant<-ifelse(rowSums(pd)==0,0,1)
    
    rbPal1 <- colorRampPalette(c(alpha(low_color, 0), alpha(high_color, 1)), alpha = TRUE)
    thecut<-as.numeric(cut(pd[, name[1]], breaks = 100))
    table(thecut)
    
    if (mult == 1){
      colfac<-round(mean(thecut),digits = 0)  
    } else {
      colfac<-mult
    }
    prenewcolfac<-thecut+colfac
    newcolfac<-ifelse( prenewcolfac >100,100,prenewcolfac  )
    
    color <- rbPal1(100)[newcolfac]
    subtitle <- paste("Expression of ", name[1], " (", high_color,")", sep = "")
    if (length(name) == 2) {
      rbPal2 <- colorRampPalette(c(alpha(low_color, 0), alpha(high_color2, 1)), alpha = TRUE)
      thecut2<-as.numeric(cut(pd[, name[2]], breaks = 100))
      table(thecut2)
      colfac2<-round(mean(thecut2),digits = 0)
      prenewcolfac2<-thecut2+colfac2
      newcolfac2<-ifelse( prenewcolfac2 >100,100,prenewcolfac2  )
      
      color2 <- rbPal2(100)[newcolfac2]
      
      
      color <- sapply(1:length(color), function(m) colorRampPalette(c(color[m], color2[m]), alpha = TRUE)(3)[2])
      subtitle <- paste(subtitle, " and ", name[2], " (", high_color2, ")", sep = "")
    }
    alfa<-ifelse(color %in% losblancos, 0.1,1)
  }
  
  ggplot(stvea_object@cite_emb, aes_string(x = colnames(stvea_object@cite_emb)[1], 
                                           y = colnames(stvea_object@cite_emb)[2], color = factor(1:length(color)))) + 
    geom_point(size = pt_size, alpha=alfa) + labs(title = paste(print_type, 
                                                                "expression"), subtitle = subtitle) + scale_color_manual(values = alpha(color, 
                                                                                                                                        1), guide = FALSE) + theme_void()
}

PlotSpatialColored<-function (stvea_object, pt_size = 0.5, type="protein", background="black") 
{
  if (is.null(stvea_object@codex_spatial)) {
    stop("stvea_object does not contain CODEX spatial information")
  }
  if (type == "protein") {
    if (!is.null(stvea_object@codex_clean)) {
      plotting_data <- stvea_object@codex_clean
    }
    else if (!is.null(stvea_object@codex_protein)) {
      plotting_data <- stvea_object@codex_protein
    }
    else {
      stop("stvea_object must contain CODEX protein data with type=\"protein\"")
    }
    if (-1 %in% stvea_object@codex_clusters) {
      colors <- c("gray", rainbow_hcl(length(unique(stvea_object@codex_clusters)) - 
                                        1, c = 80))
    }
    else {
      colors <- rainbow_hcl(length(unique(stvea_object@codex_clusters)), 
                            c = 80)
    }
    
  }
  x_tmp <- stvea_object@codex_spatial[, "x"]
  x_tmp <- x_tmp - min(x_tmp)
  y_tmp <- stvea_object@codex_spatial[, "y"]
  y_tmp <- y_tmp - min(y_tmp)
  spatial_tmp <- as.data.frame(cbind(x = x_tmp, y = y_tmp))
  dim(spatial_tmp)
  
  if (background != "black"){
    ggplot(spatial_tmp, aes(x = x,y = y, color = factor(stvea_object@codex_clusters))) + 
      geom_point(size = pt_size) + scale_color_manual(values = colors, 
                                                      name = "cluster") + guides(colour = guide_legend(override.aes = list(size = 5))) + 
      theme_void()
  } else {
    ggplot(spatial_tmp, aes(x = x, y = y, color = factor(stvea_object@codex_clusters))) + 
      geom_point(size = pt_size, alpha = 0.5) + scale_color_manual(values = colors,name = "cluster") + 
      guides(colour = guide_legend(override.aes = list(size = 5))) +
      theme_void() + 
      coord_fixed() + theme(plot.background = element_rect(fill = "black"), 
                            legend.title = element_text(colour="white"),
                            legend.text = element_text(colour="white")
      )
  }
}

HeatmapbyGroup<- function (Data,Genes, Groups = "all", Ncells=30, Breaks=c(-3,3),Heatmap_Color="purpleyellow", group_order=NULL) {
  # Heatmap_Color could be "bluered" or "purpleyellow"
  pretab<-Data[,c(1:5)]
  Data<-Data[,colnames(Data) %in% Genes, drop=F]
  if (ncol(Data)==1){
    Data<-cbind(groups=pretab[,1],Data,pretab[,2:ncol(pretab)])
    Data <- Data[, !duplicated(colnames(Data))]
  } else {
    Data<-cbind(groups=pretab[,1],Data)  
  }
  
  if (any(Data[,1]==0)) {
    Data[,1]<-as.factor(ifelse(as.numeric(as.character(Data[,1]))<0,-1, as.numeric(as.character(Data[,1])) +1) )
    
  }
  suppressWarnings({
    if (Groups == "all"){
      groupstoplot = as.numeric(as.character(as.vector(sort(unique(Data[,1])))))
    } else {
      groupstoplot = Groups
    }
  })
  DF<-NULL
  DFlist<-list()
  n=Ncells
  for (gg in 1: length(groupstoplot)) {
    g<-groupstoplot[gg]
    df<-Data[Data[,1] == g,]
    print(paste ("group",g))
    print(dim(df))
    if (nrow(df)<n) {
      n=nrow(df)
    } else {}
    set.seed(123) 
    df2<-df[sample(nrow(df), n), ]
    DF<-rbind(DF,df2)
    DFlist[g]<-list(df2)
  }
  dim(DF)
  head(DF)
  if (is.null(group_order)) {
    elorder<-groupstoplot
    DF<-DF[order(factor(DF[,1], levels=unique(elorder))),]
  } else {
    print(paste(c("The order to plot the groups will be ", group_order), collapse = " " ) )
    print(group_order)
    elorder<-group_order
    DF<-DF[order(factor(DF[,1], levels=unique(elorder))),]
  }
  DF[,1:3]
  library(data.table)
  count.dups <- function(DF){
    
    DT <- data.table(DF)
    DT[,.N, by = names(DT)]
  }
  conteo<-count.dups(DF[,1])
  colores<-colorRampPalette(rainbow(6))(length(unique (DF[,1])))
  ColSideColores=rep(colores, times=conteo$N)
  
  thebreaks<-seq(Breaks[1],Breaks[2],0.1)

  if (Heatmap_Color == "purpleyellow") {
    thecolors<-colorRampPalette(c("#CC3399","#000000","yellow")) (length(thebreaks)-1)
    thesepcolor="white"
    thena.color="black"
  } else {
    thecolors<-colorRampPalette(c("darkblue","white","darkred")) (length(thebreaks)-1)
    thesepcolor="black"
    thena.color="white"
  }
  
  heatmap.2(t(as.matrix(DF[,2:ncol(DF)]) ),
            labCol=FALSE,
            ylab="Genes",
            xlab="Clusters", 
            ColSideColors=ColSideColores,
            trace="none",
            na.color = thena.color,
            col = thecolors,
            breaks=thebreaks,
            scale="row",
            dendrogram="none",
            Colv = conteo$DF, 
            colsep = c(0,cumsum(conteo$N)),
            sepcolor=thesepcolor,
            Rowv = TRUE
  )
  legend(x=0,y=0.8, pch=c(15,15), cex=0.8, bty="n",
         legend=c(conteo$DF),
         title=c("Cell group"), col=unique(ColSideColores))
  
}


PlotExprCITE_v3<- function (stvea_object, name, type = "RNA", high_color = "red", 
                            high_color2 = "darkred", low_color = "light gray", mult=1,
                            pt_size = 0.8) {
  print_type <- type
  alfa=1
  if (type == "RNA") {
    plotting_data <- stvea_object@cite_mRNA_norm
  } else if (type == "protein") {
    plotting_data <- stvea_object@cite_protein
    print_type <- "Protein"
  } else {
    stop("type must be either \"RNA\" or \"protein\"", call. = FALSE)
  }
  if (length(name) > 2) {
    stop("name must be at most length 2", call. = FALSE)
  }
  rbPal1 <- colorRampPalette(c(alpha(low_color, 0), alpha(high_color, 
                                                          1)), alpha = TRUE)
  color <- rbPal1(100)[as.numeric(cut(plotting_data[, name[1]], 
                                      breaks = 100))]
  subtitle <- paste("Expression of ", name[1], " (", high_color, 
                    ")", sep = "")
  if (length(name) == 2 ) {
    rbPal2 <- colorRampPalette(c(alpha(low_color, 0), alpha(high_color2, 
                                                            1)), alpha = TRUE)
    color2 <- rbPal2(100)[as.numeric(cut(plotting_data[, 
                                                       name[2]], breaks = 100))]
    color <- sapply(1:length(color), function(m) colorRampPalette(c(color[m], 
                                                                    color2[m]), alpha = TRUE)(3)[2])
    subtitle <- paste(subtitle, " and ", name[2], " (", high_color2, 
                      ")", sep = "")
  }
  if (low_color == "light gray"){
    sort(table(color))
    
    losblancos<- c("#D3D0D002")
    pd<-as.data.frame(plotting_data[,name,drop=F])
    summary(pd)
    pd <- as.data.frame(lapply(pd, minmaxnorm))
    summary(pd)
    pd$discriminant<-ifelse(rowSums(pd)==0,0,1)
    
    rbPal1 <- colorRampPalette(c(alpha(low_color, 0), alpha(high_color, 1)), alpha = TRUE)
    thecut<-as.numeric(cut(pd[, name[1]], breaks = 100))
    table(thecut)
    
    if (mult == 1){
      colfac<-round(mean(thecut),digits = 0)  
    } else {
      colfac<-mult
    }
    prenewcolfac<-thecut+colfac
    newcolfac<-ifelse( prenewcolfac >100,100,prenewcolfac  )
    
    color <- rbPal1(100)[newcolfac]
    subtitle <- paste("Expression of ", name[1], " (", high_color,")", sep = "")
    if (length(name) == 2) {
      rbPal2 <- colorRampPalette(c(alpha(low_color, 0), alpha(high_color2, 1)), alpha = TRUE)
      thecut2<-as.numeric(cut(pd[, name[2]], breaks = 100))
      table(thecut2)
      colfac2<-round(mean(thecut2),digits = 0)
      prenewcolfac2<-thecut2+colfac2
      newcolfac2<-ifelse( prenewcolfac2 >100,100,prenewcolfac2  )
      
      color2 <- rbPal2(100)[newcolfac2]
      
      color <- sapply(1:length(color), function(m) colorRampPalette(c(color[m], color2[m]), alpha = TRUE)(3)[2])
      subtitle <- paste(subtitle, " and ", name[2], " (", high_color2, ")", sep = "")
    }
    alfa<-ifelse(color %in% losblancos, 0.1,1)
  }
  
  breaks<-seq(min(plotting_data[, name[1]]),max(plotting_data[, name[1]]),0.1)
  thecolors<-colorRampPalette(c(low_color,high_color)) (length(breaks)-1)
  
  if (type=="protein") {
    ggplot(stvea_object@cite_emb, aes_string(x = colnames(stvea_object@cite_emb)[1], 
                                             y = colnames(stvea_object@cite_emb)[2], 
                                             color = plotting_data[, name[1]]
    )) + 
      geom_point(size = pt_size, alpha=alfa) + 
      labs(title = paste(print_type, "expression"), subtitle = subtitle) + 
      scale_colour_gradient2  (low = high_color2,mid = low_color, high = high_color) +
      
      theme_void() 
  } else {
    ggplot(stvea_object@cite_emb, aes_string(x = colnames(stvea_object@cite_emb)[1], 
                                             y = colnames(stvea_object@cite_emb)[2], 
                                             color = plotting_data[, name[1]]
    )) + 
      geom_point(size = pt_size, alpha=alfa) + 
      labs(title = paste(print_type, "expression"), subtitle = subtitle) + 
      scale_colour_continuous(low=low_color,high=high_color) +
      theme_void()
  }
}


UMAPFeatureExpression<- function (stvea_object,features, type="RNA",low_color="gray75",high_color = "darkred", high_color2 = "dodgerblue2", pt_size=0.8) {
  prevec<-c("p1","p2","p3","p4","p5","p6","p7","p8","p9")
  for (i in 1:length(features)){
    vec<-prevec[1:length(features)]
    assign(vec[i], PlotExprCITE_v3 (stvea_object, name=features[i], type = type,low_color=low_color,high_color=high_color, high_color2=high_color2, pt_size=pt_size) )
    print (paste(vec[i], " calculated!"))
    
  }
  
  if (length(vec) ==1) {
    p<-get("p1")
  } else if (length(vec) ==2) {
    p<-grid.arrange(p1,p2)
  } else if (length(vec) ==3) {
    p<-grid.arrange(p1,p2,p3)
  } else if (length(vec) ==4) {
    p<-grid.arrange(p1,p2,p3,p4)
  } else if (length(vec) ==5) {
    p<-grid.arrange(p1,p2,p3,p4,p5)
  } else if (length(vec) ==6) {
    p<-grid.arrange(p1,p2,p3,p4,p5,p6)
  } else if (length(vec) ==7) {
    p<-grid.arrange(p1,p2,p3,p4,p5,p6,p7)
  } else if (length(vec) ==8) {
    p<-grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8)
  } else if (length(vec) ==9) {
    p<-grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9)
  }  else {
    stop("For convenience, You can only plot a maximum of 9 plots at a time")
  }
  return(p)
}

FactortoDummy<-function (y, Name = NULL) {
  if (is.null(Name)) 
    Name = "C"
  ncat = length(levels(y))
  n = length(y)
  Z = matrix(0, n, ncat)
  for (i in 1:n) Z[i, as.numeric(y[i])] = 1
  colnames(Z) <- paste(Name, levels(y), sep = "-")
  return(Z)
}


DEgenesPairwiseComparison <- function(X, groups){
  #table(X[,1])
  if (any(X[,1]==0)) {
    X[,1]<-as.factor(ifelse(as.numeric(as.character(X[,1]))<0,-1, as.numeric(as.character(X[,1])) +1) )
    
  }
  
  if (length(groups)>2) {
    stop("Error!: You only can insert one or two groups at a time") 
  } else if (length(groups)==2) {
    datpiece<-X[X[,1] %in%groups,]
    dim(datpiece)
    table (datpiece[,1])
  } else {
    datpiece <- X
    vectito<-ifelse(datpiece[,1] == groups, datpiece[,1], "ctrol")
    datpiece[,1]<-vectito
  }
  
  if (is.numeric(groups) ) {
    datpiece[,1]<-as.factor(as.numeric(as.character(datpiece[,1])))
  } else {
    datpiece[,1]<-as.factor(datpiece[,1])
  }
  
  table (datpiece[,1])
  predesign<-datpiece[,1]
  datpiece<-datpiece[,-1]
  datpiece[1:5,1:5]
  
  
  Design<-FactortoDummy (predesign)
  colnames(Design)<-c("Case","Control")
  
  ContrastMat<-limma::makeContrasts(Diff=Case-Control,
                                    levels = colnames(Design) )
  ContrastMat
  
  dim(datpiece)
  dim(Design)
  class(datpiece)
  class(Design)
  datpiece<-as.matrix(datpiece)
  anyNA(Design)
  anyNA(datpiece)
  
  fit<-lmFit(t(datpiece),Design)
  fit2<-contrasts.fit(fit, ContrastMat)
  suppressWarnings(fit3<-eBayes(fit2) )
  
  deg<-topTable(fit3, coef = "Diff",
                p.value = 0.05, adjust.method = "fdr",lfc = log2(1.5),
                number = nrow(datpiece))
  
  degALL<-topTable(fit3, coef = "Diff",
                   p.value = 1, adjust.method = "fdr",lfc = log2(0),
                   number = nrow(datpiece))
  dim(deg);dim (degALL)
  summary (deg$logFC)
  deg2<- deg[order(deg$logFC),]
  deg2<-cbind(deg2,DE_notDE=rep("DE",nrow(deg2)) )
  deg2ALL<- degALL[order(degALL$logFC),]
  deg2ALL<-deg2ALL[setdiff(rownames(deg2ALL),rownames(deg2)), ]
  deg2ALL<-cbind(deg2ALL,DE_notDE=rep("non_DE",nrow(deg2ALL)) )
  dim(deg);dim(deg2);dim (degALL);dim (deg2ALL)
  return(list (DE=deg2,nDE=deg2ALL) )
}


RidgeplotbyGroup<- function (Data,Genes, Groups) {
  pretab<-Data[,c(1:5)]
  Data<-Data[,colnames(Data) %in% Genes, drop=F]
  if (ncol(Data)==1){
    Data<-cbind(groups=pretab[,1],Data,pretab[,2:ncol(pretab)])
    Data <- Data[, !duplicated(colnames(Data))]
  } else {
    Data<-cbind(groups=pretab[,1],Data)  
  }
  
  if (any(Data[,1]==0)) {
    Data[,1]<-as.factor(ifelse(as.numeric(as.character(Data[,1]))<0,-1, as.numeric(as.character(Data[,1])) +1) )
    
  }
  suppressWarnings({
    if (Groups == "all"){
      groupstoplot = as.numeric(as.character(as.vector(sort(unique(Data[,1])))))
      
      
    } else {
      groupstoplot = Groups
    }
  })
  #####
  
  cptoplot2<-reshape2::melt(Data)
  head(cptoplot2)
  p1<-ggplot(cptoplot2, aes(x = value, fill = variable)) +
    geom_density_ridges_gradient(aes(y=variable),scale = 3) +
    theme_ridges(font_size = 13, grid = TRUE) + theme(axis.title.y = element_blank()) + ggtitle("Distribution of the entire data"); p1
  
  dat3<- Data[Data[,1] %in% groupstoplot,]
  cptoplot3<-reshape2::melt(dat3[,-1])
  
  if (length(Groups)>9) {
    eltitle<- paste("More than 9 Groups are represented")  
    
    p1<-ggplot(cptoplot2, aes(x = value, fill = variable)) +
      geom_density_ridges_gradient(aes(y=variable),scale = 3) +
      theme_ridges(font_size = 13, grid = TRUE) + theme(axis.title.y = element_blank()) + ggtitle("Distribution of the entire data"); p1
    
    p2<-ggplot(cptoplot3, aes(x = value, fill = variable)) +
      geom_density_ridges_gradient(aes(y=variable),scale = 3) +
      theme_ridges(font_size = 13, grid = TRUE) + theme(axis.title.y = element_blank(), plot.title = element_text(size = 10, face = "bold")) + ggtitle(eltitle); p2
    
  } else if (length(Groups)>2 & length(Groups)<9) {
    eltitle<- paste(c("Groups:",Groups),collapse = " ")  
    p1<-ggplot(cptoplot2, aes(x = value, fill = variable)) +
      geom_density_ridges_gradient(aes(y=variable),scale = 3) +
      theme_ridges(font_size = 13, grid = TRUE) + theme(axis.title.y = element_blank()) + ggtitle("Distribution of the entire data"); p1
    
    p2<-ggplot(cptoplot3, aes(x = value, fill = variable)) +
      geom_density_ridges_gradient(aes(y=variable),scale = 3) +
      theme_ridges(font_size = 13, grid = TRUE) + theme(axis.title.y = element_blank(), plot.title = element_text(size = 6, face = "bold")) + ggtitle(eltitle); p2
    
    cptoplot4<-reshape2::melt(dat3)
    head(cptoplot4)
    var<-paste(cptoplot4[,2],cptoplot4[,1], sep = "_")
    cptoplot5<-cbind(var,cptoplot4);cptoplot5<-cptoplot5[,-c(2,3)]
    head(cptoplot5)
    p3<-ggplot(cptoplot5, aes(x = value, fill = var)) +
      geom_density_ridges_gradient(aes(y=var),scale = 3) +
      theme_ridges(font_size = 8, grid = TRUE) + theme(axis.title.y = element_blank(),legend.position = "none") + ggtitle(eltitle); p3
    
  } else if (length(Groups)==2) {
    eltitle<- paste("Distribution of the studied groups: ", Groups[1], " and ", Groups[2]) 
    p1<-ggplot(cptoplot2, aes(x = value, fill = variable)) +
      geom_density_ridges_gradient(aes(y=variable),scale = 3) +
      theme_ridges(font_size = 13, grid = TRUE) + theme(axis.title.y = element_blank()) + ggtitle("Distribution of the entire data"); p1
    
    p2<-ggplot(cptoplot3, aes(x = value, fill = variable)) +
      geom_density_ridges_gradient(aes(y=variable),scale = 3) +
      theme_ridges(font_size = 13, grid = TRUE) + theme(axis.title.y = element_blank()) + ggtitle(eltitle); p2
    
    cptoplot4<-reshape2::melt(dat3)
    head(cptoplot4)
    var<-paste(cptoplot4[,2],cptoplot4[,1], sep = "_")
    cptoplot5<-cbind(var,cptoplot4);cptoplot5<-cptoplot5[,-c(2,3)]
    head(cptoplot5)
    p3<-ggplot(cptoplot5, aes(x = value, fill = var)) +
      geom_density_ridges_gradient(aes(y=var),scale = 3) +
      theme_ridges(font_size = 8, grid = TRUE) + theme(axis.title.y = element_blank(),legend.position = "none") + ggtitle(eltitle); p3
  } else {
    eltitle<- paste("Distribution of the studied group: ", Groups[1])  
    p1<-ggplot(cptoplot2, aes(x = value, fill = variable)) +
      geom_density_ridges_gradient(aes(y=variable),scale = 3) +
      theme_ridges(font_size = 13, grid = TRUE) + theme(axis.title.y = element_blank()) + ggtitle("Distribution of the entire data"); p1
    
    p2<-ggplot(cptoplot3, aes(x = value, fill = variable)) +
      geom_density_ridges_gradient(aes(y=variable),scale = 3) +
      theme_ridges(font_size = 13, grid = TRUE) + theme(axis.title.y = element_blank()) + ggtitle(eltitle); p2
    
    cptoplot4<-reshape2::melt(dat3)
    head(cptoplot4)
    var<-paste(cptoplot4[,2],cptoplot4[,1], sep = "_")
    cptoplot5<-cbind(var,cptoplot4);cptoplot5<-cptoplot5[,-c(2,3)]
    head(cptoplot5)
    p3<-ggplot(cptoplot5, aes(x = value, fill = var)) +
      geom_density_ridges_gradient(aes(y=var),scale = 3) +
      theme_ridges(font_size = 8, grid = TRUE) + theme(axis.title.y = element_blank(),legend.position = "none") + ggtitle(eltitle); p3
  }
  
  #####
  return(list(p1,p2,p3))
  
}

ViolinplotbyGroup<- function (Data,Genes, Groups) {
  
  pretab<-Data[,c(1:5)]
  Data<-Data[,colnames(Data) %in% Genes, drop=F]
  if (ncol(Data)==1){
    Data<-cbind(groups=pretab[,1],Data,pretab[,2:ncol(pretab)])
    Data <- Data[, !duplicated(colnames(Data))]
  } else {
    Data<-cbind(groups=pretab[,1],Data)  
  }
  
  if (any(Data[,1]==0)) {
    Data[,1]<-as.factor(ifelse(as.numeric(as.character(Data[,1]))<0,-1, as.numeric(as.character(Data[,1])) +1) )
    
  }
  suppressWarnings({
    if (Groups == "all"){
      groupstoplot = as.numeric(as.character(as.vector(sort(unique(Data[,1])))))
      
      
    } else {
      groupstoplot = Groups
    }
  })
  #####
  
  cptoplot2<-reshape2::melt(Data)
  head(cptoplot2)
  
  if (length(Groups)>9) {
    eltitle<- paste("More than 9 Groups are represented")  
    p1 <- ggplot(cptoplot2, aes(factor(variable), value)) +
      geom_violin(aes(fill = variable)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      ggtitle("Distribution of the entire data");p1#+ geom_jitter(height = 0, width = 0.1)
    
    dat3<- Data[Data[,1] %in% groupstoplot,]
    if (ncol(dat3)>2) {
      cptoplot3<-reshape2::melt(dat3[,-1])  
    } else {
      thevar<-rep(Genes,nrow(dat3))
      dat3<- cbind(var=thevar,dat3)
      cptoplot3<-reshape2::melt(dat3[,-1])  
    }
    
    p2 <- ggplot(cptoplot3, aes(factor(variable), value)) +
      geom_violin(aes(fill = variable)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      ggtitle(eltitle);p2#+ geom_jitter(height = 0, width = 0.1)

  } else if (length(Groups)>2 & length(Groups)<9) {
    eltitle<- paste(c("Groups:",Groups),collapse = " ")
    
    dat3<- Data[Data[,1] %in% groupstoplot,]
    if (ncol(dat3)>2) {
      cptoplot3<-reshape2::melt(dat3[,-1])  
    } else {
      thevar<-rep(Genes,nrow(dat3))
      dat3<- cbind(var=thevar,dat3)
      cptoplot3<-reshape2::melt(dat3[,-1])  
    }
    
    p1 <- ggplot(cptoplot2, aes(factor(variable), value)) +
      geom_violin(aes(fill = variable)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      ggtitle("Distribution of the entire data");p1
    
    p2 <- ggplot(cptoplot3, aes(factor(variable), value)) +
      geom_violin(aes(fill = variable)) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),plot.title = element_text(size = 6, face = "bold")) +
      ggtitle(eltitle);p2#+ geom_jitter(height = 0, width = 0.1)
    
    cptoplot4<-reshape2::melt(dat3)
    head(cptoplot4)
    var<-paste(cptoplot4[,2],cptoplot4[,1], sep = "_")
    cptoplot5<-cbind(var,cptoplot4);cptoplot5<-cptoplot5[,-c(2,3)]
    head(cptoplot5)
    p3 <- ggplot(cptoplot5, aes(factor(var), value)) +
      geom_violin(aes(fill = var)) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = "none") +
      ggtitle(eltitle);p3#+ 
    
  } else if (length(Groups)==2) {
    eltitle<- paste("Distribution of the studied groups: ", Groups[1], " and ", Groups[2])  
    p1 <- ggplot(cptoplot2, aes(factor(variable), value)) +
      geom_violin(aes(fill = variable)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      ggtitle("Distribution of the entire data");p1
    
    dat3<- Data[Data[,1] %in% groupstoplot,]
    if (ncol(dat3)>2) {
      cptoplot3<-reshape2::melt(dat3[,-1])  
    } else {
      thevar<-rep(Genes,nrow(dat3))
      dat3<- cbind(var=thevar,dat3)
      cptoplot3<-reshape2::melt(dat3[,-1])  
    }
    
    p2 <- ggplot(cptoplot3, aes(factor(variable), value)) +
      geom_violin(aes(fill = variable)) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),plot.title = element_text(size = 6, face = "bold")) +
      ggtitle(eltitle);p2#+ geom_jitter(height = 0, width = 0.1)
    
    cptoplot4<-reshape2::melt(dat3)
    head(cptoplot4)
    var<-paste(cptoplot4[,2],cptoplot4[,1], sep = "_")
    cptoplot5<-cbind(var,cptoplot4);cptoplot5<-cptoplot5[,-c(2,3)]
    head(cptoplot5)
    p3 <- ggplot(cptoplot5, aes(factor(var), value)) +
      geom_violin(aes(fill = var)) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = "none") +
      ggtitle(eltitle);p3#+ 
    
  } else {
    eltitle<- paste("Distribution of the studied group: ", Groups[1])  
    p1 <- ggplot(cptoplot2, aes(factor(variable), value)) +
      geom_violin(aes(fill = variable)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      ggtitle("Distribution of the entire data");p1
    
    dat3<- Data[Data[,1] %in% groupstoplot,]
    if (ncol(dat3)>2) {
      cptoplot3<-reshape2::melt(dat3[,-1])  
    } else {
      thevar<-rep(Genes,nrow(dat3))
      dat3<- cbind(var=thevar,dat3)
      cptoplot3<-reshape2::melt(dat3[,-1])  
    }
    
    p2<-ggplot(cptoplot3, aes(x = value, fill = variable)) +
      geom_density_ridges_gradient(aes(y=variable),scale = 3) +
      theme_ridges(font_size = 13, grid = TRUE) + theme(axis.title.y = element_blank()) + ggtitle(eltitle); p2
    
    cptoplot4<-reshape2::melt(dat3)
    head(cptoplot4)
    var<-paste(cptoplot4[,2],cptoplot4[,1], sep = "_")
    cptoplot5<-cbind(var,cptoplot4);cptoplot5<-cptoplot5[,-c(2,3)]
    head(cptoplot5)
    p3<-ggplot(cptoplot5, aes(x = value, fill = var)) +
      geom_density_ridges_gradient(aes(y=var),scale = 3) +
      theme_ridges(font_size = 13, grid = TRUE) + theme(axis.title.y = element_blank(),legend.position = "none") + ggtitle(eltitle); p3
  }
  return(list(p1,p2,p3))
  
}
                      

# END

