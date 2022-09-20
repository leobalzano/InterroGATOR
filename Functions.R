# This is the Beta version
#############################################
#####     InterroGATOR_functions.R     ######
#############################################
# Author: Leandro Balzano-Nogueira
# Diabetes Institute, University of Florida (Gainesville)
# Created: November/5/2021
# Last update: September/20/2022

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

outersect <- function(x, y, ...) {
  big.vec <- c(x, y, ...)
  duplicates <- big.vec[duplicated(big.vec)]
  setdiff(big.vec, unique(duplicates))
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
  } else {
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

PlotClusterCITE_v3<-function (object_table,color_by, pt_size = 0.5) {
  if (any(color_by==0)) {
    color_by<-as.factor(ifelse(as.numeric(as.character(color_by))<0,-1, as.numeric(as.character(color_by)) +1) )
    
  }
  
  if (-1 %in% color_by) {
    # colors <- c("gray95", colorspace::rainbow_hcl(length(unique(color_by)) - 1, c = 80))
    colors <- c("gray95", randomcoloR::distinctColorPalette(length(unique(color_by)) - 1))
    
  } else {
    # colors <- colorspace::rainbow_hcl(length(unique(color_by)), c = 80)
    if (length(unique(color_by)) == 1 ) {
      colors <- randomcoloR::distinctColorPalette(length(unique(color_by)))
    } else {
      vecofcols<-c( "#4FAFB7", "#8034E4", "#ADB56F", "#6C84EC", "#E6B93F", "#CBE7EA",
                    "#F19D38", "#65AD8E","#EB5428", "#DAE650","#7E8C8E", "#F5E4B2","#DC4CA7",
                    "#AB7B90","#6AE68B", "#BBA4E6", "#F9CDC1", "#E7E9D4","#7E85D7","#FBE1C4",
                    "#DC9F86","#877780","#75E6E4","#6F9557","#D7A74E", "#D1D0EE",
                    "#63B8DA","#D3A1D5","#C445DF","#852F02"
      )
      
      ccccol<-rep(vecofcols,length(unique(color_by)))
      colors<-ccccol [1:length(unique(color_by))]
      
      #colors <- randomcoloR::distinctColorPalette(length(unique(color_by)))
    }
    
  }
  if ("Unknown" %in% color_by || "Unknown_or_Out_of_Subset" %in% color_by) {
    colors<-c(colors[1:length(colors)-1],"gray95")
    #color_by<-ifelse (color_by == "Unknown", "gray95",color_by)
    
    ggplot(object_table, aes_string(x = colnames(object_table)[1], 
                                    y = colnames(object_table)[2])) +# ,    color = factor(colors))) + 
      #geom_point(size = pt_size) + 
      
      geom_point(aes(colour=factor(color_by)),size = pt_size,pch=16) + 
      scale_color_manual(values = colors, name = "cluster") + 
      guides(colour = guide_legend(override.aes = list(size = 3))) + 
      #geom_text(aes(label=factor(color_by)), size=3) +
      theme_void()
  } else {
    
    ggplot(object_table, aes_string(x = colnames(object_table)[1], 
                                    y = colnames(object_table )[2])) +      
      geom_point(aes(colour=factor(color_by)),size = pt_size,pch=16) + 
      scale_color_manual(values = colors, name = "cluster") + 
      guides(colour = guide_legend(override.aes = list(size = 3))) + 
      theme_void()
    
  }
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

HeatmapbyGroup<- function (Data,Genes, Groups = "all", Ncells=30, Breaks=c(-3,3),Heatmap_Color="purpleyellow", group_order=NULL,
                           legend_size=0.8,feature_size = 0.8
                           ) {
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
    n=Ncells
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
            Rowv = TRUE,
            cexRow=feature_size
  )
  ifelse(length(conteo$DF)>10,legend(
    x=0,y=0.8, pch=c(15,15), cex=legend_size, bty="n",
         legend=c(conteo$DF), ncol = 2,
         title=c("Cell group"), col=unique(ColSideColores)),
    legend(
      x=0,y=0.8, pch=c(15,15), cex=legend_size, bty="n",
      legend=c(conteo$DF),
      title=c("Cell group"), col=unique(ColSideColores))
    )
  
}


PlotExprCITE_v3<- function (stvea_object, name, type = "RNA", high_color = "red", 
                            high_color2 = "darkred", low_color = "light gray", mult=1,
                            pt_size = 0.8,min.cutoff=0) {
  print_type <- type
  alfa=1
  name<-gsub('\\.', '-', name)
  print(name)
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
    
    if (!is.null(min.cutoff) ){
      cutoffValue<-quantile (plotting_data[,name[1]],probs = min.cutoff)
      plotting_data[,name[1]] <- ifelse(plotting_data[,name[1]] > cutoffValue, plotting_data[,name[1]], 0)  
    } else {}
    
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


UMAPFeatureExpression<- function (stvea_object,features, type="RNA",low_color="gray75",high_color = "darkred", high_color2 = "dodgerblue2", pt_size=0.8,min.cutoff=0) {
  prevec<-c("p1","p2","p3","p4","p5","p6","p7","p8","p9")
  for (i in 1:length(features)){
    vec<-prevec[1:length(features)]
    assign(vec[i], PlotExprCITE_v3 (stvea_object, name=features[i], type = type,low_color=low_color,high_color=high_color, high_color2=high_color2, pt_size=pt_size, min.cutoff=min.cutoff) )
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
      ggtitle(eltitle);p2
    
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
      ggtitle(eltitle);p2
    
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


PlotExprCODEXspatial_v3<-function (stvea_object, name, type = "protein", high_color = "red", 
                                   high_color2 = "green", low_color = "white", pt_size = 0.8,Subset="no", SmoothScatterplot=FALSE, mult=1, maxPercentile=0.9, min.cutoff = 0, maximize_differences="no") {
  # This version is to be able to use the transfer matrix as a list format.
  
  alfa=0.5
  color3 = spatial_tmp= NULL
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
      
      ####### # piece to subset the data to be plotted
      
      if (is.list(plotting_data)) {
        plotting_data_2<-NULL
        for (ll in 1: length(plotting_data)) {
          piece<-subset(plotting_data[[ll]], select = name)
          plotting_data_2<-rbind(plotting_data_2,piece)
        }
        plotting_data_2[,1]<-minmaxnorm(plotting_data_2[,1])
        plotting_data_2[,2]<-minmaxnorm(plotting_data_2[,2])
        plotting_data <- plotting_data_2
        summary(plotting_data)
        hist(plotting_data[,1])
        hist(plotting_data[,2])
      } else { }
      
      #######
    }
    else {
      stop("type must be either \"RNA\" or \"protein\"", call. = FALSE)
    }
    if (length(name) > 2) {
      stop("name must be at most length 2", call. = FALSE)
    }
    
    rbPal1 <- colorRampPalette(c(alpha(low_color, 0), alpha(high_color, 
                                                            1)), alpha = TRUE)
    color <- rbPal1(5)[as.numeric(cut(plotting_data[, name[1]], 
                                      breaks = 5))]
    subtitle <- paste("Expression of ", name[1], " (", high_color, 
                      ")", sep = "")
    if (length(name) == 2) {
      rbPal2 <- colorRampPalette(c(alpha(low_color, 0), alpha(high_color2, 
                                                              1)), alpha = TRUE)
      color2 <- rbPal2(5)[as.numeric(cut(plotting_data[, name[2]], breaks = 5))]
      preprecolor3<-lapply(1:length(color), function(m) colorRampPalette(c(color[m], color2[m]), alpha = TRUE)(3)[])
      precolor3 <- data.frame(matrix(unlist(preprecolor3), ncol=3, byrow=TRUE),stringsAsFactors=FALSE); rm(preprecolor3)
      color3 <- ifelse(precolor3[,1] =="#FFFFFF00" & precolor3[,3]=="#FFFFFF00",precolor3[,2],
                       ifelse(precolor3[,1] !="#FFFFFF00" & precolor3[,3]=="#FFFFFF00",precolor3[,1],
                              ifelse(precolor3[,1] =="#FFFFFF00" & precolor3[,3]!="#FFFFFF00",precolor3[,3],precolor3[,2])))
      testeo<-cbind (precolor3,color3)
      
      subtitle <- paste(subtitle, " and ", name[2], " (", high_color2, 
                        ")", sep = "")
    }
    x_tmp <- stvea_object@codex_spatial[, "x"]
    x_tmp <- x_tmp - min(x_tmp)
    x_tmpMAX<-max(x_tmp)
    y_tmp <- stvea_object@codex_spatial[, "y"]
    y_tmp <- y_tmp - min(y_tmp)
    y_tmpMAX<-max(y_tmp)
    spatial_tmp <- as.data.frame(cbind(x = x_tmp, y = y_tmp,col=color3))
    #spatial_tmp_comparison <- as.data.frame(cbind(plotting_data,x = x_tmp, y = y_tmp,col1=color,col2=color2,col3=color3))
    
    
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
    #spatial_tmp_subset<-spatial_tmp_subset[!spatial_tmp_subset$col == "#FFFFFF00", ]
    ###
    if (SmoothScatterplot) {
      
      ggplot(spatial_tmp_subset) + aes(x = as.numeric(as.character(spatial_tmp_subset$x)), y = as.numeric(as.character(spatial_tmp_subset$y)), color = factor(1:length(spatial_tmp_subset$col))) + 
        geom_point(size = pt_size, alpha = 0.5) + 
        scale_color_manual(values = alpha(spatial_tmp_subset$col, 1)) + 
        #scale_color_manual(values = spatial_tmp_subset$col) + 
        guides(color = FALSE) + 
        ylim(YMAXsubset, YMINsubset) + xlim(XMINsubset,XMAXsubset) +
        labs(title = paste("Spatial",  type, "expression"), subtitle = subtitle) + 
        theme_void() + 
        coord_fixed()  
    } else{
      
      ggplot(spatial_tmp_subset, aes(x = as.numeric(as.character(spatial_tmp_subset$x)), y = as.numeric(as.character(spatial_tmp_subset$y)), color = factor(1:length(spatial_tmp_subset$col)))) + 
        geom_point(size = pt_size, alpha = 0.5) + 
        guides(color = FALSE) +
        scale_color_manual(values = spatial_tmp_subset$col) +  
        ylim(max(y_tmp), 0) + xlim(0,max(x_tmp)) +
        labs(title = paste("Spatial",  type, "expression"), subtitle = subtitle) + theme_void() + 
        coord_fixed()
    }
    
    
  }  else {
    ###
    if (type == "protein") {
      if (!is.null(stvea_object@codex_clean)) {
        plotting_data <- stvea_object@codex_clean
        plotting_data<-subset(plotting_data, select = name)
      }
      else if (!is.null(stvea_object@codex_protein)) {
        plotting_data <- stvea_object@codex_protein
        plotting_data<-subset(plotting_data, select = name)
      }
      else {
        stop("stvea_object must contain CODEX protein data with type=\"protein\"")
      }
    }
    else if (type == "RNA") {
      plotting_data <- stvea_object@codex_mRNA
      
      ####### # piece to subset the data to be plotted
      
      if (is.list(plotting_data)) {
        plotting_data_2<-NULL
        for (ll in 1: length(plotting_data)) {
          piece<-subset(plotting_data[[ll]], select = name)
          plotting_data_2<-rbind(plotting_data_2,piece)
        }
        for (fea in 1:ncol(plotting_data_2)) {
          plotting_data_2[,fea]<-minmaxnorm(plotting_data_2[,fea])
        }
        
        plotting_data <- plotting_data_2
        summary(plotting_data)
      } else { 
        plotting_data_2<-NULL
        plotting_data_2<-subset(plotting_data, select = name)
        for (fea in 1:ncol(plotting_data_2) ) {
          plotting_data_2[,fea]<-minmaxnorm(plotting_data_2[,fea])
        }
        
        plotting_data <- plotting_data_2
        summary(plotting_data)
      }
      ####
      if (!is.null(min.cutoff) ){
        cutoffValue<-quantile (plotting_data,probs = min.cutoff)
        plotting_data <- ifelse(plotting_data > cutoffValue, plotting_data, 0)  
      } else {}
      
        
      ####
      #######
    } else {
      stop("type must be either \"RNA\" or \"protein\"", call. = FALSE)
    }
    if (length(name) > 2) {
      stop("name must be at most length 2", call. = FALSE)
    }
    
    rbPal1 <- colorRampPalette(c(alpha(low_color, 0), alpha(high_color, 
                                                            1)), alpha = TRUE)
    color <- rbPal1(1000)[as.numeric(cut(plotting_data[, name[1]], 
                                        breaks = 1000))] # 10 or 1000 test
    minitab_temporal<-cbind(plotting_data[, name[1]],color)
    maxcolor1<-unique(minitab_temporal[minitab_temporal[,1]==max(minitab_temporal[,1]),2])
    
    subtitle <- paste("Expression of ", name[1], " (", high_color, 
                      ")", sep = "")
    if (length(name) == 2) {
      rbPal2 <- colorRampPalette(c(alpha(low_color, 0), alpha(high_color2, 
                                                              1)), alpha = TRUE)
      color2 <- rbPal2(1000)[as.numeric(cut(plotting_data[, 
                                                         name[2]], breaks = 1000))] # 10 or 1000 test
      minitab_temporal_2<-cbind(plotting_data[, name[2]],color2)
      maxcolor2<-unique(minitab_temporal_2[minitab_temporal_2[,1]==max(minitab_temporal_2[,1]),2])
      
      
      
      
      preprecolor3<-lapply(1:length(color), function(m) colorRampPalette(c(color[m], color2[m]), alpha = TRUE)(3)[])
      precolor3 <- data.frame(matrix(unlist(preprecolor3), ncol=3, byrow=TRUE),stringsAsFactors=FALSE); rm(preprecolor3)
      morecommoncolor<-names(sort(table(precolor3$X3), decreasing = T)[1])
      
      precolor3<- data.frame(plotting_data,precolor3)
      rownames(precolor3)<-seq(1:nrow(precolor3))
      
      
      if (tolower(maximize_differences) == "no"){
        color3<-ifelse(precolor3[,1] == 0 & precolor3[,2] == 0 ,low_color,precolor3[,4])
        subtitle <- paste(subtitle, " and ", name[2], " (", high_color2, 
                          ")", ". cut-off= ",(min.cutoff*100)," %", sep = "")
      } else {
        color3 <- ifelse(precolor3[,1] != 0 & precolor3[,2] == 0 ,maxcolor1,
                  ifelse(precolor3[,1] == 0 & precolor3[,2] != 0 ,maxcolor2,
                  ifelse(precolor3[,1] != 0 & precolor3[,2] != 0 ,precolor3[,4],low_color
                  )))
        subtitle <- paste(subtitle, " and ", name[2], " (", high_color2, 
                          ")", ". cut-off= ",(min.cutoff*100)," %. maximizing expression differences", sep = "")
      }
      testeo<-cbind (precolor3,color3)
      
    }
    
    if (is.null(color3)) {
      color3<-color
    }
    if (type == "protein") {
      color3 <- sapply(1:length(color), function(m) colorRampPalette(c(color[m],
                                        color2[m]), alpha = TRUE)(3)[2])
    }
    
    
    x_tmp <- stvea_object@codex_spatial[, "x"]
    x_tmp <- x_tmp - min(x_tmp)
    y_tmp <- stvea_object@codex_spatial[, "y"]
    y_tmp <- y_tmp - min(y_tmp)
    spatial_tmp <- as.data.frame(cbind(x = x_tmp, y = y_tmp,col=color3))
    dim(spatial_tmp) 
    
    if (low_color == "transparent"){
      
      losblancos<- c("#FFFFFF00")
      spatial_tmp<-spatial_tmp[!spatial_tmp[,3]==losblancos,]
      color3<-spatial_tmp[,3]
    }
    
    spatial_tmp_subset<-spatial_tmp
    ggplot(spatial_tmp_subset) + aes(x = as.numeric(as.character(spatial_tmp_subset$x)), y = as.numeric(as.character(spatial_tmp_subset$y)), color = factor(1:length(color3))) + 
      geom_point(size = pt_size, alpha = alfa, shape=16) + 
      scale_color_manual(values = alpha(color3,1)) + guides(color = FALSE) + ylim(max(y_tmp), 0) + labs(title = paste("Spatial", type, "expression"), subtitle = subtitle) + theme_void() + 
      coord_fixed()
  }
}

PlotExprCODEXspatial_v4<-function (stvea_object, name, type = "protein", high_color = "red", 
                                   high_color2 = "green", low_color = "white", pt_size = 0.8,Subset="no", SmoothScatterplot=FALSE, mult=1, maxPercentile=0.9, min.cutoff = 0, maximize_differences="no", CellsaboveThrs= NULL,CellsofInterest=NULL, IdentificationThreshold = "0") {
  # This version is to be able to use the transfer matrix as a list format.
  
  alfa=0.5
  color3 = spatial_tmp= NULL
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
      
      ####### Nuevo
      if (!is.null(CellsaboveThrs)) {
        plotting_data<-plotting_data[rownames(plotting_data) %in% CellsaboveThrs,]
      } else {}
      ####### Fin de Nuevo
    }
    else if (type == "RNA") {
      plotting_data <- stvea_object@codex_mRNA
      ####### # piece to subset the data to be plotted
      
      if (is.list(plotting_data)) {
        plotting_data_2<-NULL
        for (ll in 1: length(plotting_data)) {
          piece<-subset(plotting_data[[ll]], select = name)
          plotting_data_2<-rbind(plotting_data_2,piece)
        }
        plotting_data_2[,1]<-minmaxnorm(plotting_data_2[,1])
        plotting_data_2[,2]<-minmaxnorm(plotting_data_2[,2])
        plotting_data <- plotting_data_2
        summary(plotting_data)
        hist(plotting_data[,1])
        hist(plotting_data[,2])
      } else { }
      
      #######
      ####### Nuevo
      if (!is.null(CellsaboveThrs)) {
        plotting_data<-plotting_data[rownames(plotting_data) %in% CellsaboveThrs,]
      } else {}
      ####### Fin de Nuevo
    }
    else {
      stop("type must be either \"RNA\" or \"protein\"", call. = FALSE)
    }
    if (length(name) > 2) {
      stop("name must be at most length 2", call. = FALSE)
    }
    
    rbPal1 <- colorRampPalette(c(alpha(low_color, 0), alpha(high_color, 
                                                            1)), alpha = TRUE)
    color <- rbPal1(5)[as.numeric(cut(plotting_data[, name[1]], 
                                      breaks = 5))]
    subtitle <- paste("Expression of ", name[1], " (", high_color, 
                      ")", sep = "")
    if (length(name) == 2) {
      rbPal2 <- colorRampPalette(c(alpha(low_color, 0), alpha(high_color2, 
                                                              1)), alpha = TRUE)
      color2 <- rbPal2(5)[as.numeric(cut(plotting_data[, name[2]], breaks = 5))]
      preprecolor3<-lapply(1:length(color), function(m) colorRampPalette(c(color[m], color2[m]), alpha = TRUE)(3)[])
      precolor3 <- data.frame(matrix(unlist(preprecolor3), ncol=3, byrow=TRUE),stringsAsFactors=FALSE); rm(preprecolor3)
      color3 <- ifelse(precolor3[,1] =="#FFFFFF00" & precolor3[,3]=="#FFFFFF00",precolor3[,2],
                       ifelse(precolor3[,1] !="#FFFFFF00" & precolor3[,3]=="#FFFFFF00",precolor3[,1],
                              ifelse(precolor3[,1] =="#FFFFFF00" & precolor3[,3]!="#FFFFFF00",precolor3[,3],precolor3[,2])))
      testeo<-cbind (precolor3,color3)
      
      subtitle <- paste(subtitle, " and ", name[2], " (", high_color2, 
                        ")", sep = "")
    }
    x_tmp <- stvea_object@codex_spatial[, "x"]
    x_tmp <- x_tmp - min(x_tmp)
    x_tmpMAX<-max(x_tmp)
    y_tmp <- stvea_object@codex_spatial[, "y"]
    y_tmp <- y_tmp - min(y_tmp)
    y_tmpMAX<-max(y_tmp)
    spatial_tmp <- as.data.frame(cbind(x = x_tmp, y = y_tmp,col=color3))
     
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

    if (SmoothScatterplot) {
      
      ggplot(spatial_tmp_subset) + aes(x = as.numeric(as.character(spatial_tmp_subset$x)), y = as.numeric(as.character(spatial_tmp_subset$y)), color = factor(1:length(spatial_tmp_subset$col))) + 
        geom_point(size = pt_size, alpha = 0.5) + 
        scale_color_manual(values = alpha(spatial_tmp_subset$col, 1)) + 
        guides(color = FALSE) + 
        ylim(YMAXsubset, YMINsubset) + xlim(XMINsubset,XMAXsubset) +
        labs(title = paste("Spatial",  type, "expression"), subtitle = subtitle) + 
        theme_void() + 
        coord_fixed()  
    } else{
      
      ggplot(spatial_tmp_subset, aes(x = as.numeric(as.character(spatial_tmp_subset$x)), y = as.numeric(as.character(spatial_tmp_subset$y)), color = factor(1:length(spatial_tmp_subset$col)))) + 
        geom_point(size = pt_size, alpha = 0.5) + 
        guides(color = FALSE) +
        scale_color_manual(values = spatial_tmp_subset$col) +  
        ylim(max(y_tmp), 0) + xlim(0,max(x_tmp)) +
        labs(title = paste("Spatial",  type, "expression"), subtitle = subtitle) + theme_void() + 
        coord_fixed()
    }
    
    
  }  else {
    ###
    if (type == "protein") {
      if (!is.null(stvea_object@codex_clean)) {
        plotting_data <- stvea_object@codex_clean
        plotting_data<-subset(plotting_data, select = name)
      }
      else if (!is.null(stvea_object@codex_protein)) {
        plotting_data <- stvea_object@codex_protein
        plotting_data<-subset(plotting_data, select = name)
      }
      else {
        stop("stvea_object must contain CODEX protein data with type=\"protein\"")
      }
      ####### Nuevo
      if (!is.null(CellsaboveThrs)) {
        plotting_data<-plotting_data[rownames(plotting_data) %in% CellsaboveThrs,,drop=FALSE]
      } else {}
      ####### Fin de Nuevo
    }
    else if (type == "RNA") {
      plotting_data <- stvea_object@codex_mRNA
      
      ####### # piece to subset the data to be plotted
      
      if (is.list(plotting_data)) {
        plotting_data_2<-NULL
        for (ll in 1: length(plotting_data)) {
          piece<-subset(plotting_data[[ll]], select = name)
          plotting_data_2<-rbind(plotting_data_2,piece)
        }
        for (fea in 1:ncol(plotting_data_2)) {
          plotting_data_2[,fea]<-minmaxnorm(plotting_data_2[,fea])
        }
        
        plotting_data <- plotting_data_2
        summary(plotting_data)
      } else { 
        plotting_data_2<-NULL
        plotting_data_2<-subset(plotting_data, select = name)
        for (fea in 1:ncol(plotting_data_2) ) {
          plotting_data_2[,fea]<-minmaxnorm(plotting_data_2[,fea])
        }
        
        plotting_data <- plotting_data_2
        summary(plotting_data)
      }
      
      ####### Nuevo
      if (!is.null(CellsaboveThrs)) {
        plotting_data<-plotting_data[rownames(plotting_data) %in% CellsaboveThrs,,drop=FALSE]
      } else {}
      ####### Fin de Nuevo
      
      ####
      if (!is.null(min.cutoff) ){
        cutoffValue<-quantile (plotting_data,probs = min.cutoff)
        plotting_data <- ifelse(plotting_data > cutoffValue, plotting_data, 0)  
      } else {}
      
      
      ####
      #######
    } else {
      stop("type must be either \"RNA\" or \"protein\"", call. = FALSE)
    }
    if (length(name) > 2) {
      stop("name must be at most length 2", call. = FALSE)
    }
    
    rbPal1 <- colorRampPalette(c(alpha(low_color, 0), alpha(high_color, 
                                                            1)), alpha = TRUE)
    color <- rbPal1(1000)[as.numeric(cut(plotting_data[, name[1]], 
                                         breaks = 1000))] # 10 or 1000 test
    minitab_temporal<-cbind(plotting_data[, name[1]],color)
    maxcolor1<-unique(minitab_temporal[minitab_temporal[,1]==max(minitab_temporal[,1]),2])
    
    subtitle <- paste("Expression of ", name[1], " (", high_color, 
                      ")", sep = "")
    if (length(name) == 2) {
      rbPal2 <- colorRampPalette(c(alpha(low_color, 0), alpha(high_color2, 
                                                              1)), alpha = TRUE)
      color2 <- rbPal2(1000)[as.numeric(cut(plotting_data[, 
                                                          name[2]], breaks = 1000))] # 10 or 1000 test
      minitab_temporal_2<-cbind(plotting_data[, name[2]],color2)
      maxcolor2<-unique(minitab_temporal_2[minitab_temporal_2[,1]==max(minitab_temporal_2[,1]),2])
      
      preprecolor3<-lapply(1:length(color), function(m) colorRampPalette(c(color[m], color2[m]), alpha = TRUE)(3)[])
      precolor3 <- data.frame(matrix(unlist(preprecolor3), ncol=3, byrow=TRUE),stringsAsFactors=FALSE); rm(preprecolor3)
      morecommoncolor<-names(sort(table(precolor3$X3), decreasing = T)[1])
      
      precolor3<- data.frame(plotting_data,precolor3)
      subtitle<-paste(subtitle, " and ", name[2], " (", high_color2, 
            ")", sep="")
      if (is.null(CellsaboveThrs) ) {
        precolor3$above_thrs<- rep ("in", nrow(precolor3))
      } else {
        precolor3$above_thrs<- ifelse(rownames(precolor3) %in% CellsaboveThrs,"in","out")  
      }
      if (is.null(maximize_differences)) {
        color3<-ifelse(precolor3[,6] == "out",low_color,
                       ifelse (precolor3[,1] == 0 & precolor3[,2] == 0 ,low_color,precolor3[,4])
        )
      } else if (tolower(maximize_differences) == "no"){
        color3<-ifelse(precolor3[,6] == "out",low_color,
                       ifelse (precolor3[,1] == 0 & precolor3[,2] == 0 ,low_color,precolor3[,4])
                       )
        
        subtitle <- paste(subtitle, " and ", name[2], " (", high_color2, 
                          ")", ". Threshold= ",as.numeric(as.character(IdentificationThreshold)), ". cut-off= ",(min.cutoff*100)," %", sep = "")
      } else {
        
        color3<-ifelse(precolor3[,6] == "out",low_color,
                       ifelse(precolor3[,1] != 0 & precolor3[,2] == 0 ,maxcolor1,
                              ifelse(precolor3[,1] == 0 & precolor3[,2] != 0 ,maxcolor2,
                                     ifelse(precolor3[,1] != 0 & precolor3[,2] != 0 ,precolor3[,4],low_color
                                               ))) )
        
        subtitle <- paste(subtitle, " and ", name[2], " (", high_color2, 
                          ")",". Threshold= ",as.numeric(as.character(IdentificationThreshold)), ". cut-off= ",(min.cutoff*100)," %. maximizing expression differences", sep = "")
      }
      
      testeo<-cbind (precolor3,color3)
      
    }
    
    if (is.null(color3)) {
      color3<-color
    }
    if (type == "protein" & length (name) == 1) {
      color3 <- sapply(1:length(color), function(m) colorRampPalette(c(color[m],
                                                                       "gray95"[m]), alpha = TRUE)(3)[2])
    }
    
    if (!is.null(CellsofInterest)) {
      pd<-data.frame(cbind(plotting_data,color3=color3))
      color3<-ifelse (rownames(pd) %in% CellsofInterest,pd$color3, "gray95")
    }
    
    therownames<-rownames(stvea_object@codex_spatial)
    x_tmp <- stvea_object@codex_spatial[, "x"]
    x_tmp <- x_tmp - min(x_tmp)
    y_tmp <- stvea_object@codex_spatial[, "y"]
    y_tmp <- y_tmp - min(y_tmp)
    prespatial_temp<-as.data.frame(cbind (x = x_tmp, y = y_tmp))
    rownames(prespatial_temp)<-therownames
    prespatial_temp<-prespatial_temp[rownames(prespatial_temp) %in% CellsaboveThrs,,drop=FALSE]
    dim(prespatial_temp)
    prespatial_temp$y<- -1 * prespatial_temp$y
    spatial_tmp <- as.data.frame(cbind(prespatial_temp,col=color3))
    dim(spatial_tmp) 
    
    if (low_color == "transparent"){
      losblancos<- c("#FFFFFF00")
      spatial_tmp<-spatial_tmp[!spatial_tmp[,3]==losblancos,]
      color3<-spatial_tmp[,3]
    }
    
    spatial_tmp_subset<-spatial_tmp
    ggplot(spatial_tmp_subset) + aes(x = as.numeric(as.character(spatial_tmp_subset$x)), y = as.numeric(as.character(spatial_tmp_subset$y)), color = factor(1:length(color3))) + 
      geom_point(size = pt_size, alpha = alfa, shape=16) + 
      scale_color_manual(values = alpha(color3,1)) + guides(color = FALSE) + ylim(0,min(spatial_tmp_subset$y)) + labs(title = paste("Spatial", type, "expression"), subtitle = subtitle) + theme_void() + 
      coord_fixed()
  }
}



AdjScore_calculator<- function (stvea_object,feature_type,Feature_Pairs) {
  if (feature_type=="RNA") {
    print("Calculating adjacency amongst genes")
    AdjScoreGenes(stvea_object, gene_pairs=Feature_Pairs, k=3, num_cores=1) 
  } else {
    print("Calculating adjacency amongst proteins")
    AdjScoreProteins(stvea_object, protein_pairs=Feature_Pairs, k=3, num_cores=1) 
  } 
}


AdjScoreClustersCODEX_vGator <-function (stvea_object,clusters,
                                         infotable,
                                         k=3, num_cores = 1) 
{
  if (is.null(stvea_object@codex_spatial)) {
    stop("stvea_object does not contain CODEX spatial information")
  }
  cluster_ids<-unique(clusters)
  
  if (length(cluster_ids) < 2) {
    stop("Cannot compute adjacency score of fewer than 2 clusters")
  }
  cluster_ids <- cluster_ids[order(cluster_ids)]
  cluster_matrix <- t(sapply(cluster_ids, function(x) (cluster_ids == 
                                                         x) * 1))
  row.names(cluster_matrix) <- cluster_ids
  colnames(cluster_matrix) <- cluster_ids
  
  
  knn_adj <- knn_graph(stvea_object@codex_spatial, k = k)
  AdjScoreClustersCODEX.internal(knn_adj, clusters, num_cores = num_cores)
}

adjacency_score_v2<-function (adj_matrix, f, f_pairs, c, num_perms = 1000, num_cores = 1, 
                              perm_estimate = F, groupings = F, verbose = T) 
{
  ptm <- proc.time()
  if (!is(f, "matrix") && !is(f, "Matrix")) {
    cat("Converting f to matrix\n")
    f <- as.matrix(f)
  }
  if (is(f_pairs, "list")) {
    f_pairs <- matrix(unlist(f_pairs), ncol = 2, byrow = T)
  } else if (is(f_pairs, "numeric") || is(f_pairs, "character")) {
    f_pairs <- matrix(f_pairs, ncol = 2, byrow = T)
  }
  if (c != 0 && groupings) {
    groupings <- FALSE
    cat("Setting groupings to FALSE since c > 0\n")
  }
  if (groupings && num_perms > 0) {
    num_perms <- 0
    cat("Setting num_perms to 0 since using grouping features\n")
  }
  if (verbose) 
    cat("Creating permutation matrices")
  permutations <- NULL
  if (num_perms > 0) {
    permutations <- t(mcmapply(function(x) sample(1:ncol(f)), 
                               1:num_perms, mc.cores = num_cores))
  }
  permutations <- rbind(1:ncol(f), permutations)
  perm_f <- mclapply(1:nrow(f), function(i) keep_sparse(t(sapply(1:nrow(permutations), 
                                                                 function(j) f[i, ][permutations[j, ]])), f), mc.cores = num_cores)
  names(perm_f) <- row.names(f)
  if (!isSymmetric(adj_matrix)) {
    warning("Adjacency matrix is not symmetrical, computing symmetrical matrix")
    adj_sym <- 1 * ((adj_matrix + t(adj_matrix)) > 0)
  }
  else {
    adj_sym <- adj_matrix
  }
  if (c != 0) {
    expm_adj <- expm::expm(c * adj_sym, method = "Higham08")
    expm_adj <- as.matrix(expm_adj)
    expm_adj <- (expm_adj - diag(nrow(adj_sym)))/c
  }
  cornel <- function(fo) {
    f1 <- perm_f[[fo[1]]]
    f2 <- perm_f[[fo[2]]]
    if (c == 0) {
      qt <- rowSums((f1 %*% adj_sym) * f2)
    }
    else {
      qt <- rowSums((f1 %*% expm_adj) * f2)
    }
    ph <- NULL
    ph$score <- qt[1]
    if (groupings) {
      overlap <- sum(f1 * f2)
      if (overlap == 0) {
        wb <- 2 * sum(f1) * sum(f2)
        edges <- sum(adj_sym)/2
        ph$p <- 1 - phyper(qt[1], m = wb, n = length(f1) * 
                             (length(f1) - 1) - wb, k = edges)
      }
      else {
        wb <- overlap * (overlap - 1)
        edges <- sum(adj_sym)/2
        ph$p <- 1 - phyper(floor(qt[1]/2), m = wb, n = length(f1) * 
                             (length(f1) - 1) - wb, k = edges)
      }
    }
    else if (perm_estimate) {
      nfit <- fitdistr(as.numeric(qt), "normal")
      ph$p <- 1 - pnorm(qt[1], mean = nfit$estimate["mean"], 
                        sd = nfit$estimate["sd"])
    }
    else {
      ph$p <- sum(qt > qt[1])/num_perms
    }
    return(ph)
  }
  worker <- function(fu) {
    qh <- NULL
    qh$score <- NULL
    qh$p <- NULL
    for (i in 1:nrow(fu)) {
      d <- cornel(fu[i, ])
      qh$score <- rbind(qh$score, d$score)
      qh$p <- rbind(qh$p, d$p)
    }
    return(data.frame(qh))
  }
  if (num_cores > nrow(f_pairs)) {
    num_cores <- nrow(f_pairs)
  }
  if (verbose) 
    cat(" -", (proc.time() - ptm)[3], "seconds\n")
  ptm <- proc.time()
  if (verbose) 
    cat("Computing adjacency score for each feature pair")
  if (num_cores == 1 || nrow(f_pairs) == 1) {
    qqh <- worker(f_pairs)
  }
  else {
    wv <- floor(nrow(f_pairs)/num_cores)
    wr <- nrow(f_pairs) - wv * num_cores
    work <- list()
    if (wr > 0) {
      for (m in 1:wr) {
        work[[m]] <- (f_pairs[(1 + (m - 1) * (wv + 1)):(m * 
                                                          (wv + 1)), , drop = FALSE])
      }
      for (m in (wr + 1):num_cores) {
        work[[m]] <- (f_pairs[(1 + wr + (m - 1) * wv):(wr + 
                                                         m * wv), , drop = FALSE])
      }
    }
    else {
      for (m in 1:num_cores) {
        work[[m]] <- (f_pairs[(1 + (m - 1) * wv):(m * 
                                                    wv), , drop = FALSE])
      }
    }
    reul <- mclapply(work, worker, mc.cores = num_cores)
    qqh <- reul[[1]]
    for (m in 2:num_cores) {
      qqh <- rbind(qqh, reul[[m]])
    }
  }
  qqh$q <- p.adjust(qqh$p, method = "BH")
  qqh <- data.frame(f = f_pairs[, 1], g = f_pairs[, 2], qqh, 
                    stringsAsFactors = F)
  if (verbose) 
    cat(" -", (proc.time() - ptm)[3], "seconds\n")
  return(qqh)
}



AdjScoreHeatmap_v2<-function (adj_score_output=protein_adj, low_color="blue",high_color="red", title="Adjacency plot", Subset_Matrix=NULL) {
  vecnames<-unique(c(adj_score_output$f,adj_score_output$g))
  heatmap_matrix <- matrix(rep(0, length(vecnames) * 
                                 length(vecnames)), ncol = length(vecnames) )
  row.names(heatmap_matrix) <- vecnames
  colnames(heatmap_matrix) <- vecnames
  for (i in 1:nrow(adj_score_output)) {
    heatmap_matrix[adj_score_output[i, "f"], adj_score_output[i, 
                                                              "g"]] <- log10(adj_score_output[i, "q"] + 1e-15)
    heatmap_matrix[adj_score_output[i, "g"], adj_score_output[i, 
                                                              "f"]] <- log10(adj_score_output[i, "q"] + 1e-15)

  }

  if (!is.null(Subset_Matrix)) {
    heatmap_matrix<-heatmap_matrix[rownames(heatmap_matrix) %in% Subset_Matrix,]
    heatmap_matrix<-heatmap_matrix[,colnames(heatmap_matrix) %in% Subset_Matrix]

  } else {}
  heatmap_matrix_transf<-heatmap_matrix
  heatmap_matrix_transf<-ifelse(heatmap_matrix_transf< (-4),-4,heatmap_matrix_transf)
  thebreaks<-seq(min(heatmap_matrix_transf),max(heatmap_matrix_transf),length.out = 6)
  thecolors<-colorRampPalette(c(low_color,high_color)) (length(thebreaks)-1)
  heatmap(heatmap_matrix, margins = c(10, 10), symm = T, col=thecolors,main= title)
  legend(x="right", legend=c("Low","", "Medium","", "High"),fill=thecolors)
  
}


AdjScoreHeatmap_cute<- function (adj_score_output,title="Adjacency tuned by Pearson's Correlation",Subset_Matrix=NULL) {
  library(RColorBrewer)
  library(pheatmap)
  library(gplots)
  
  
  vecnames<-unique(c(adj_score_output$f,adj_score_output$g))
  heatmap_matrix <- matrix(rep(0, length(vecnames) * 
                                 length(vecnames)), ncol = length(vecnames) )
  row.names(heatmap_matrix) <- vecnames
  colnames(heatmap_matrix) <- vecnames
  for (i in 1:nrow(adj_score_output)) {
    heatmap_matrix[adj_score_output[i, "f"], adj_score_output[i, 
                                                              "g"]] <- log10(adj_score_output[i, "q"] + 1e-15)
    heatmap_matrix[adj_score_output[i, "g"], adj_score_output[i, 
                                                              "f"]] <- log10(adj_score_output[i, "q"] + 1e-15)
  }
  heatmap_matrix2<-cor(heatmap_matrix,method = "pearson")
  heatmap_matrix2[is.na(heatmap_matrix2)]<-1
  
  if (!is.null(Subset_Matrix)) {
    heatmap_matrix2<-heatmap_matrix2[rownames(heatmap_matrix2) %in% Subset_Matrix,]
    heatmap_matrix2<-heatmap_matrix2[,colnames(heatmap_matrix2) %in% Subset_Matrix]
  } else {}
 

  coul<-colorRampPalette(rev(brewer.pal(9,"RdYlBu")), space="Lab") 
  heatmap.2(heatmap_matrix2, margins = c(10, 10), symm = F, col= coul(100),
            breaks = seq(min(heatmap_matrix2),max(heatmap_matrix2),length.out=101),
            trace="none",
            key = T,
            Rowv=T,
            Colv=T,
            sepwidth=c(0.001, 0.001),  # width of the borders
            main=title,
            sepcolor='white',colsep=1:ncol(heatmap_matrix2),rowsep=1:nrow(heatmap_matrix2)        # color of the separation lines
  )
  
}

legend.col <- function(col, lev) {
  
  opar <- par
  
  n <- length(col)
  
  bx <- par("usr")
  
  box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000,
              bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 50)
  box.cy <- c(bx[3], bx[3])
  box.sy <- (bx[4] - bx[3]) / n
  z
  xx <- rep(box.cx, each = 2)
  
  par(xpd = TRUE)
  for(i in 1:n){
    
    yy <- c(box.cy[1] + (box.sy * (i - 1)),
            box.cy[1] + (box.sy * (i)),
            box.cy[1] + (box.sy * (i)),
            box.cy[1] + (box.sy * (i - 1)))
    polygon(xx, yy, col = col[i], border = col[i])
    
  }
  par(new = TRUE)
  plot(0, 0, type = "n",
       ylim = c(min(lev), max(lev)),
       yaxt = "n", ylab = "",
       xaxt = "n", xlab = "",
       frame.plot = FALSE)
  axis(side = 4, las = 2, tick = FALSE, line = .25)
  par <- opar
}


PlotClusterCODEXemb_vGator<- function (stvea_object, cluster_column, infotable=infotable4(), pt_size = 0.5, Selection_on_Top = "no") 
{
  tab<-merge(stvea_object@codex_emb,infotable, by.x=0, by.y = "CODEXname")
  rownames(tab)<-tab[,1];tab<-tab[,-1]
  tab2<-tab[ order(as.numeric(row.names(tab))), ]
  tab2
  clusters<-(tab2[,cluster_column])
  
  if (any(clusters==0)) {
    clusters<-as.factor(ifelse(as.numeric(as.character(clusters))<0,-1, as.numeric(as.character(clusters)) +1) )
    
  }
  
  
  if (-1 %in% clusters) {
    colors <- c("gray", randomcoloR::distinctColorPalette(length(unique(clusters)) - 1, c = 80))
  } else {
    colors <- randomcoloR::distinctColorPalette(length(unique(clusters)))
  }
  
  
  if (any(grep ("Unknown",clusters))) {
    colors<-c(colors[1:length(colors)-1],"gray95")
  } else {}
  
  if (Selection_on_Top=="no") {
    ggplot(stvea_object@codex_emb, aes_string(x = colnames(stvea_object@codex_emb)[1], 
                                              y = colnames(stvea_object@codex_emb)[2], color = factor(clusters))) + 
      geom_point(size = pt_size) + scale_color_manual(values = colors, 
                                                      name = "cluster") + guides(colour = guide_legend(override.aes = list(size = 5))) + 
      theme_void()
  } else {
    tabelita2<-data.frame(cbind(clusters=sort(unique(clusters)),colors=colors))
    tab2
    
    VecColor<-NULL
    for (c in 1:nrow(tab2)) {
      if (tab2$Prob_threshold[c] == "out") {
        vectitocolor<-"gray95"
        VecColor<-c(VecColor,vectitocolor)
      } else{
        vectitocolor<-tabelita2$colors[tab2[c,cluster_column] == tabelita2$clusters]
        VecColor<-c(VecColor,vectitocolor)
      }
    }
    
    datita<-cbind(stvea_object@codex_emb,clusters=clusters, colors2=VecColor)
    
    print("datita")
    print(datita[1:50,])
    
    datita<-datita %>% arrange(factor(colors2, levels = c("gray95",tabelita2$colors[1:length(tabelita2$colors)-1])) )
    
    print("head(datita) sorted" )
    print(head(datita) )
    print("tail(datita) sorted " )
    print(tail(datita) )
    
    ggplot(datita, aes_string(x = colnames(datita)[1], 
                              y = colnames(datita)[2], 
                              color = factor(datita[,3])
    )) + 
      geom_point(size = pt_size) + scale_color_manual(values = colors, 
                                                      name = "cluster") + guides(colour = guide_legend(override.aes = list(size = 5))) + 
      theme_void()
    
      
  }
}


PlotClusterCODEX_Independently<- function(stvea_object, color_by,pt_size, highlight, cluster_column, Color="dodgerblue2", infotable=infotable4(), Selection_on_Top="no") {
  
  tab<-merge(stvea_object@codex_emb,infotable, by.x=0, by.y = "CODEXname")
  rownames(tab)<-tab[,1];tab<-tab[,-1]
  tab2<-tab[ order(as.numeric(row.names(tab))), ]
  tab2
  clusters<-(tab2[,cluster_column])
  
  
  if (any(clusters==0)) {
    clusters<-as.factor(ifelse(as.numeric(as.character(clusters))<0,-1, as.numeric(as.character(clusters)) +1) )
    
  }
  
  if (-1 %in% clusters) {
    colors <- c("gray", colorspace::rainbow_hcl(length(unique(clusters)) - 
                                                  1, c = 80))
  } else {
    colors <- colorspace::rainbow_hcl(length(unique(clusters)), 
                                      c = 80)
  }
  
  fill2<- ifelse(clusters %in% highlight, Color ,"gray95")
  colors2<- ifelse(clusters %in% highlight, Color ,"gray75")
  colors3<-unique(colors2)
  
  if (Selection_on_Top=="no") {
    ggplot(stvea_object@codex_emb, aes_string(x = colnames(stvea_object@codex_emb)[1], 
                                             y = colnames(stvea_object@codex_emb)[2], 
                                             color = factor(fill2)
    )) + 
      geom_point(fill=fill2,colour=colors2,size = pt_size,pch=21) + 
      labs(title = paste("Group",highlight)) +
      theme_void()
  } else {
    datita<-cbind(stvea_object@codex_emb,fill2=fill2, colors2=colors2)
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


# END
