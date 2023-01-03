#############################################
#######     InterroGATOR_server.R     #######
#############################################
# Author: Leandro Balzano-Nogueira
# Diabetes Institute, University of Florida (Gainesville)
# Created: October/19/2021
# Last update: December/22/2022

# Creating and hiding calculations for posterior analyses

# This is a Shiny application for conducting differential expression and Spatial analyses on different cell types in different human tissues. This application allows the user to analyze cell expression profiles from CITEseq and CODEX data previously integrated.


#############################################

# Libraries:
library(shiny)
library(shinydashboard)
library(STvEA)
library (expm)
library (reticulate)
library(ggplot2)
library (Matrix)
library(data.table)
library(gridExtra)
library(patchwork)
library(ggpointdensity)
library(ggridges)
library (reshape2)
library (hexbin)
library(colorspace)
library (multtest)
library (metap)
library (limma)
library (Cairo)
library (gplots)
library (Seurat)
library (dplyr)
library (RColorBrewer)
library (grDevices)
library(randomcoloR)
library(AdjacencyScore)
library(sparseMatrixStats)


#############################################
## server.R ##
source("Functions.R")

options(shiny.maxRequestSize=120000*1024^2)

shinyServer(function(input, output, session) {
    
    # File upload segment
    ######################
    infotable1 <- reactiveVal(NULL);
    infotable2 <- reactiveVal(NULL);
    infotable3 <- reactiveVal(NULL);
    tabproportions<-reactiveVal(NULL);
    prenamefortab4<-reactiveVal(NULL);
    namefortab4<-reactiveVal(NULL);
    thepath<-reactiveVal(NULL);
    FileLocation<-reactiveVal(NULL);
    
    infotable4  <- reactiveVal(NULL);   # Infotable4()
    infotab4_tmp  <- reactiveVal(NULL);   # Infotable4()
    info_spatial_tmp <- reactiveVal(NULL); # Infotable4()
    tabproportionsCODEX<-reactiveVal(NULL); # Infotable4()
    tabproportionsCODEX2<-reactiveVal(NULL); # Infotable4()
    
    infofileAreaSubsetted<-reactiveVal(NULL); #If area subsetted or not
    Tab4piece2<-reactiveVal(NULL); #If area subsetted or not
    
    observeEvent({
      input$file1
    }, {
      if (is.null(input$file1)) {
        return(NULL)
      }
      
      infotable1(readRDS (input$file1$datapath) )
      infotable3(infotable1()@meta.data[,-which(sapply(infotable1()@meta.data, class) == "numeric")])
      
    })
    
    observeEvent({
      input$file2
    }, {
      if (is.null(input$file2)) {
        return(NULL)
      }
        
      
      elnewenv <- new.env()
      elname<- load(input$file2$datapath, envir = elnewenv)
      infotable2(elnewenv[[elname]] )
      
    })
    
    
    observeEvent({
      infotable2()} , {
        req(infotable2())
        req(infotable3())
        
        if(any(infotable2()@cite_clusters == -1) ) {
          ccai<-infotable2()@cite_clusters
          ccai[ccai== -1] = "Unknown"
          ccai <-factor(ccai)
          
          infotable3( cbind(infotable3(),cite_clusters_after_integration=ccai)   )
        } else {
          infotable3( cbind(infotable3(),cite_clusters_after_integration=infotable2()@cite_clusters)  )  
          
        }
        
        
        coloring_choices <- colnames(infotable3() )# It works
        coloring_choices2 <- c(colnames(infotable3() ), "codex_clusters" )# It works
        updateSelectInput(session, inputId = "colby", choices = coloring_choices)
        updateSelectInput(session, inputId = "clusterUMAPIndSamples", choices = coloring_choices )
        updateSelectInput(session, inputId = "colbyDE", choices = coloring_choices2 )
        
        #######   Piece to save table4

        if ( isFALSE(file_test("-f", FileLocation())) ) {
          print ("Calculating infotable4")

          
          v1<-rownames(infotable2()@transfer_matrix)
          divisor=2000
          vector_of_pieces<-unique(c(seq(from =1, to=nrow(infotable2()@transfer_matrix), by=divisor),nrow(infotable2()@transfer_matrix)));vector_of_pieces

          allv2<-allv3<-allv4<-allv5<-tab<-NULL
          count= 0
          for (n in 1:(length(vector_of_pieces)-1)) {
            count<-count+1
            print ("count")
            print (count)
            print ("n")
            print (n)

            print("vector_of_pieces[count]")
            print(vector_of_pieces[count])
            print(vector_of_pieces[count+1]-1)
            if ((n+1) != length(vector_of_pieces)){
              tabelita<-infotable2()@transfer_matrix[vector_of_pieces[count]:(vector_of_pieces[count+1]-1),]
              dim(tabelita)
              vectito <-colnames(tabelita)[apply(tabelita,1,which.max)]
              vecv2<-rowSums(as.matrix(tabelita))
              elmax <-rowMaxs(tabelita)

              allv2<-c(allv2,vecv2)
              allv3<-c(allv3,vectito)
              allv4<-c(allv4,elmax) # maximos

            } else {
              tabelita<-infotable2()@transfer_matrix[vector_of_pieces[count]:(vector_of_pieces[count+1]),]
              vectito <-colnames(tabelita)[apply(tabelita,1,which.max)]
              vecv2<-rowSums(as.matrix(tabelita))
              elmax <-rowMaxs(tabelita)
              allv2<-c(allv2,vecv2)
              allv3<-c(allv3,vectito)
              allv4<-c(allv4,elmax) # maximos
            }
          }
          v2<-allv2
          v3<-allv3
          v4<-allv4
          v5<- (v4 *100)/v2
          v5[is.nan(v5)]<- 0


          # Retain For each row return the column name of the largest value (It takes time!)
          tab<-cbind(CODEXname=v1,rowsum=v2, maxValue = v4, probCITEcell=v5,CITEname=v3)
          print("dim(tab)")
          print(dim(tab))
          
          print (head (tab))
          print (tail (tab))
          tab[order(tab[,"CODEXname"]),]
          rm(v1,v2,v3,v4,v5)
          sort(table(tab[,5]))
          tab2<-tab[tab[,2]!=0,]
          cellstab<-infotable3()
          
          tab3<-merge(tab2,cellstab, by.x=5, by.y=0)
          tab3<-tab3[order(tab3$CODEXname),]
          pretab4<-cbind(infotable2()@codex_spatial,codex_clusters = infotable2()@codex_clusters)
          
          tab4<-merge(tab3,pretab4, by.x = "CODEXname",by.y=0, all.y = TRUE)
          

          if ("predicted.celltype.l2" %in% colnames(tab4)) {
            print("predicted.celltype.l2 is present in the data")
            tab4$predicted.celltype.l2[is.na(tab4$predicted.celltype.l2)] = "Unknown"
            tab4$condition <- with(tab4, ifelse(tab4$predicted.celltype.l2 =="Unknown", 0, 1))
          } else if ("predicted.celltype" %in% colnames(tab4)) {
            print("predicted.celltype is present in the data")
            tab4$predicted.celltype[is.na(tab4$predicted.celltype)] = "Unknown"
            tab4$condition <- with(tab4, ifelse(tab4$predicted.celltype =="Unknown", 0, 1))
          } else if ("predicted.id" %in% colnames(tab4)) {
            print("predicted.id is present in the data and it was used to calculate the condition variable")
            tab4$predicted.id[is.na(tab4$predicted.id)] = "Unknown"
            tab4$condition <- with(tab4, ifelse(tab4$predicted.id =="Unknown", 0, 1))
          }
          if ("Higher_Hierarchy_grouping" %in% colnames(tab4)) {
            print("Higher_Hierarchy_grouping is present in the data")
            tab4$Higher_Hierarchy_grouping[is.na(tab4$Higher_Hierarchy_grouping)] = "Unknown"
          }
          if ("seurat_clusters" %in% colnames(tab4)) {
            print("seurat_clusters is present in the data")
            tab4$seurat_clusters<-as.numeric(as.character(tab4$seurat_clusters))
            tab4$seurat_clusters[is.na(tab4$seurat_clusters)] = "Unknown"
          }
          if ("seurat_clusters_vLP" %in% colnames(tab4)) {
            print("seurat_clusters is present in the data")
            tab4$seurat_clusters_vLP<-as.numeric(as.character(tab4$seurat_clusters_vLP))
            tab4$seurat_clusters_vLP[is.na(tab4$seurat_clusters_vLP)] = "Unknown"
          }
          
          if ("orig.ident" %in% colnames(tab4)) {
            print("orig.ident is present in the data")
            tab4$orig.ident[is.na(tab4$orig.ident)] = "Unknown"
          }
          
          
          if ("assay_name" %in% colnames(tab4)) {
            print("assay_name is present in the data")
            tab4$assay_name[is.na(tab4$assay_name)] = "Unknown"
          }
          if ("integrated_snn_res.0.5" %in% colnames(tab4)) {
            print("integrated_snn_res.0.5 is present in the data")
            tab4$integrated_snn_res.0.5<-as.numeric(as.character(tab4$integrated_snn_res.0.5))
            tab4$integrated_snn_res.0.5[is.na(tab4$integrated_snn_res.0.5)] = "Unknown"
          }
          if ("CITE_snn_res.1" %in% colnames(tab4)) {
            print("CITE_snn_res.1 is present in the data")
            tab4$CITE_snn_res.1<-as.numeric(as.character(tab4$CITE_snn_res.1))
            tab4$CITE_snn_res.1[is.na(tab4$CITE_snn_res.1)] = "Unknown"
          }
          if ("dsb_knn_res.1.5" %in% colnames(tab4)) {
            print("dsb_knn_res.1.5 is present in the data")
            tab4$dsb_knn_res.1.5<-as.numeric(as.character(tab4$dsb_knn_res.1.5))
            tab4$dsb_knn_res.1.5[is.na(tab4$dsb_knn_res.1.5)] = "Unknown"
          }
          if ("wsnn_res.1.5" %in% colnames(tab4)) {
            print("wsnn_res.1.5 is present in the data")
            tab4$wsnn_res.1.5<-as.numeric(as.character(tab4$wsnn_res.1.5))
            tab4$wsnn_res.1.5[is.na(tab4$wsnn_res.1.5)] = "Unknown"
          }
          if ("predicted.id" %in% colnames(tab4)) {
            print("predicted.id is present in the data")
            tab4$predicted.id[is.na(tab4$predicted.id)] = "Unknown"
          }
          if ("cite_clusters_after_integration" %in% colnames(tab4)) {
            print("cite_clusters_after_integration is present in the data")
            tab4$cite_clusters_after_integration<-as.numeric(as.character(tab4$cite_clusters_after_integration))
            tab4$cite_clusters_after_integration[tab4$cite_clusters_after_integration== -1] = "Unknown"
            tab4$cite_clusters_after_integration[is.na(tab4$cite_clusters_after_integration)] = "Unknown"
          }
          
          tab4 <- tab4[order(tab4$condition),]
          infotable4(tab4)
          write.csv(x=tab4,file=FileLocation())
          rm(tab,tab2,tab3,tab4, cellstab)
        } else {
          infotable4(read.csv(FileLocation(), row.names = 1 ) )
          print (tail(infotable4()))
          print ("infotable4 was retreived from previous calculations")
        }

        x_tmp <- infotable4()[, "x"]
        x_tmp <- x_tmp - min(x_tmp)
        y_tmp <- infotable4()[, "y"]
        y_tmp <- y_tmp - min(y_tmp)
        spatial_tmp <- as.data.frame(cbind(x = x_tmp, y = y_tmp))
        info_spatial_tmp(spatial_tmp)
        rm (spatial_tmp)
        coloring_choicestab4 <- colnames(infotable4()[,-which(sapply(infotable4(), class) == "numeric")] ) # 
        toremove<-c("CODEXname","CITEname","v2","X","x","y","z", "nFeature_RNA","bc","droplet_class","nFeature_CITE","nFeature_refAssay","condition")
        idx<- which (coloring_choicestab4 %in% toremove)
        
        coloring_choicestab4<-coloring_choicestab4[-idx]
        
        var_to_highlightAdjacencybyCluster <- coloring_choicestab4
        updateSelectizeInput(session,inputId="selectclusteringtype", choices = var_to_highlightAdjacencybyCluster, options = list(maxItems = 1),selected=var_to_highlightAdjacencybyCluster[1] ) 
        updateSelectizeInput(session,inputId="selectclusteringtypeUMAPclusters", choices = var_to_highlightAdjacencybyCluster, options = list(maxItems = 1),selected=var_to_highlightAdjacencybyCluster[1] ) 

        print("input$selectclusteringtype")
        print(input$selectclusteringtype)
        

        updateSelectInput(session, inputId = "cluster1", choices = coloring_choicestab4 ) # This is brother of this in UI. If this is on, the piece in UI has to be on as well
        updateSelectInput(session, inputId = "cluster2", choices = coloring_choicestab4 ) # This is brother of this in UI. If this is on, the piece in UI has to be on as well
        updateSelectInput(session, inputId = "clusterSpatialExpression", choices = coloring_choicestab4 ) # This is brother of this in UI. If this is on, the piece in UI has to be on as well
        print("input$clusterSpatialExpression")
        print(input$clusterSpatialExpression)

       
      })
    
   
    ######
    # Creating a unique name to save tab4data
    
    observeEvent({
      input$file2
    }, {
      req(input$file1)
      req(input$file2)
      prenamefortab4(paste(substr(input$file1$name,1,nchar(input$file1$name)-4),substr(input$file2$name,1,nchar(input$file2$name)-6), sep= "_") )
      namefortab4 ( paste0(prenamefortab4(),".csv") )
      print("namefortab4")
      print(namefortab4())
      thepath(getwd())
      print("thepath")
      print(thepath())
      ifelse(!dir.exists(paste0(thepath(),"/",prenamefortab4())), dir.create(paste0(thepath(),"/",prenamefortab4()), showWarnings = TRUE, recursive = FALSE, mode = "0777") , FALSE)
      
      FileLocation (paste0(thepath(),"/",prenamefortab4(),"/",namefortab4()) )
      print("FileLocation")
      print(FileLocation())
    })
    
    ######
    
    
    observeEvent({
      input$fileAreaSubsetted
    }, {
      if (is.null(input$fileAreaSubsetted)) {
        return(NULL)
      }
      
      infofileAreaSubsetted (read.csv(input$fileAreaSubsetted$datapath, row.names = 1))
      pretab<-infofileAreaSubsetted()
      pretab$y<-pretab$y * -1
      pretab$x_y<-paste(pretab$x,pretab$y, sep="_")
      infofileAreaSubsetted(pretab)
      
      x_tmp <- infotable2()@codex_spatial[, "x"]
      x_tmp <- x_tmp - min(x_tmp)
      y_tmp <- infotable2()@codex_spatial[, "y"]
      y_tmp <- y_tmp - min(y_tmp)
      spatial_tmp <- as.data.frame(cbind(x = x_tmp, y = y_tmp) )
      spatial_tmp$x_y<-paste(spatial_tmp$x,spatial_tmp$y,sep="_")
      
      print(table(infofileAreaSubsetted()$x_y %in% spatial_tmp$x_y) ) # All have to be present, YES!!!!!!!!!!!!
      
      CodexCells<-spatial_tmp[spatial_tmp$x_y %in% infofileAreaSubsetted()$x_y, ]
      CodexCells
      
      tab4.2<-infotable4()
      x_tmp <- infotable4()[, "x"]
      x_tmp <- x_tmp - min(x_tmp)
      y_tmp <- infotable4()[, "y"]
      y_tmp <- y_tmp - min(y_tmp)
      tab4.2$x_tmp<-x_tmp
      tab4.2$y_tmp<-y_tmp
      tab4.2$x_y<-paste(tab4.2$x_tmp,tab4.2$y_tmp,sep="_")
      table(tab4.2$x_y %in% CodexCells$x_y)
            
      
      Tab4piece<-tab4.2[tab4.2$x_y %in% CodexCells$x_y,]
      dim(Tab4piece)
      
      Tab4piece$CITEname[is.na(Tab4piece$CITEname)]<-"Unknown"
      Tab4piece
      Tab4piece2(Tab4piece[,c("CODEXname","CITEname","x_tmp","y_tmp")])
      dim(Tab4piece2())
      
      
    }) # Closing the observeEvent of fileAreaSubsetted
    
           
    output$p1<- renderPlot({
      req(infotable2())
      req(input$colby)
      req(infotable3())
      req(input$RetainSubset)
      
      print(input$RetainSubset)
      print("input$colby")
      print(input$colby)
      thecolor <- infotable3()[,input$colby];
      
      tab<-melt(table(thecolor))
      colnames(tab)<- c("Category", "Cell_number")
      tabproportions (tab)
      cite_embpiece<-infotable2()@cite_emb[rownames(infotable2()@cite_emb) %in% rownames(infotable3()), ]
      cite_embpiece<-merge(cite_embpiece,infotable3(), by=0);rownames(cite_embpiece)<-cite_embpiece[,1];cite_embpiece<-cite_embpiece[,-1]
      
      head(cite_embpiece)
      dim(cite_embpiece)
      
      
      if (is.null(input$RetainSubset) || input$RetainSubset == "All") {
        cite_embpiece<-infotable2()@cite_emb[rownames(infotable2()@cite_emb) %in% rownames(infotable3()), ]
        cite_embpiece<-merge(cite_embpiece,infotable3(), by=0);rownames(cite_embpiece)<-cite_embpiece[,1];cite_embpiece<-cite_embpiece[,-1]
        thecolor<- cite_embpiece[,input$colby]
      } else if (input$RetainSubset == "Area") {
        losin<-rownames(cite_embpiece[intersect(rownames(cite_embpiece),Tab4piece2()$CITEname),])
        cite_embpiece$losin<-ifelse(rownames(cite_embpiece) %in% losin,"in","out" )
        prethecolor<- infotable3()[intersect(rownames(infotable3()),rownames(cite_embpiece)),input$colby]
        cite_embpiece$prethecolor<-prethecolor
        cite_embpiece$thecolor<- ifelse(cite_embpiece$losin == "out","Unknown_or_Out_of_Subset",cite_embpiece$prethecolor)
        cite_embpiece<-cite_embpiece %>% arrange(desc(cite_embpiece$thecolor) )
        thecolor<-cite_embpiece$thecolor
        
      } else {
        losin<-rownames(cite_embpiece[outersect(rownames(cite_embpiece),Tab4piece2()$CITEname),])
        cite_embpiece$losin<-ifelse(rownames(cite_embpiece) %in% losin,"in","out" )
        prethecolor<- infotable3()[intersect(rownames(infotable3()),rownames(cite_embpiece)),input$colby]
        cite_embpiece$prethecolor<-prethecolor
        cite_embpiece$thecolor<- ifelse(cite_embpiece$losin == "out","Unknown_or_Out_of_Subset",cite_embpiece$prethecolor)
        cite_embpiece<-cite_embpiece %>% arrange(desc(cite_embpiece$thecolor) )
        thecolor<-cite_embpiece$thecolor
        
      }
      
      p=PlotClusterCITE_v3(cite_embpiece ,color_by = thecolor,pt_size =input$size1)
      plot (p)
      
      }) # Closing the plot
    
    
    output$tproportions<- DT::renderDT({
      tabproportions()
    })# Closing the datatable of proportions of classification categories
    
    
    ######################
    
    ## CITE UMAP of samples plotted individually
    #####################
    infotableumapind<- reactiveVal(NULL);
    Subsetinfotableumapind<- reactiveVal(NULL);
    
    
    observeEvent({
      input$clusterUMAPIndSamples
    }, {
      req(infotable1())
      req(infotable2())
      infotableumapind(cbind(assay_name=infotable1()@meta.data$assay_name,clusters=infotable1()@meta.data[,input$clusterUMAPIndSamples],infotable2()@cite_emb))
      groupstoplot <- sort(unique (infotableumapind()$assay_name)) 
      print("groupstoplot")  
      print(groupstoplot)  
      updateSelectInput(session, inputId = "Umapstoplot", choices = groupstoplot) 
      
      if (any(infotableumapind()$clusters==0)) {
        infotableumapind()$clusters<-as.factor(ifelse(as.numeric(as.character(infotableumapind()$clusters))<0,-1, as.numeric(as.character(infotableumapind()$clusters)) +1) )
      }
      
    })
    
    observeEvent({
      input$ploteaumapindsamples}, {
        req( input$Umapstoplot)
        req(input$RetainSubset2)
        
        colors <- randomcoloR::distinctColorPalette(length(unique(infotableumapind()$clusters)))
        
        Subsetinfotableumapind (infotableumapind()[infotableumapind()$assay_name == input$Umapstoplot,] )
        print("dim (Subsetinfotableumapind)")
        print (dim(Subsetinfotableumapind()))
        seurat_clusters<-Subsetinfotableumapind()$clusters
        
        
        output$p5<- renderPlot({
          if (is.null(input$RetainSubset2) || input$RetainSubset2 == "All") {
            Subsetinfotableumapind (infotableumapind()[infotableumapind()$assay_name == input$Umapstoplot,] )
            
          } else if (input$RetainSubset2 == "Area") {
            pretab<-Subsetinfotableumapind()
            losin<-rownames(pretab[intersect(rownames(pretab),Tab4piece2()$CITEname),])
            pretab$losin<-ifelse(rownames(pretab) %in% losin,"in","out" )
            
            pretab$thecolor<- ifelse(pretab$losin == "out","Unknown_or_Out_of_Subset",pretab$clusters)
            vecofcols<-c( "#4FAFB7", "#8034E4", "#ADB56F", "#6C84EC", "#E6B93F", "#CBE7EA",
                          "#F19D38", "#65AD8E","#EB5428", "#DAE650","#7E8C8E", "#F5E4B2","#DC4CA7",
                          "#AB7B90","#6AE68B", "#BBA4E6", "#F9CDC1", "#E7E9D4","#7E85D7","#FBE1C4",
                          "#DC9F86","#877780","#75E6E4","#6F9557","#D7A74E", "#D1D0EE",
                          "#63B8DA","#D3A1D5","#C445DF","#852F02"
            )
            ccccol<-rep(vecofcols,length(unique(pretab$thecolor)))
            colors<-ccccol [1:length(unique(pretab$thecolor))]
            
            if ("Unknown" %in% pretab$thecolor || "Unknown_or_Out_of_Subset" %in% pretab$thecolor) {
              colors<-c(colors[1:length(colors)-1],"gray95") 
            } else {}
            pretab<-pretab %>% arrange(desc(pretab$thecolor) )
            thecolor<-pretab$thecolor
            Subsetinfotableumapind(pretab)
          } else {
            pretab<-Subsetinfotableumapind()
            losin<-rownames(pretab[outersect(rownames(pretab),Tab4piece2()$CITEname),])
            pretab$losin<-ifelse(rownames(pretab) %in% losin,"in","out" )
            pretab$thecolor<- ifelse(pretab$losin == "out","Unknown_or_Out_of_Subset",pretab$clusters)
            vecofcols<-c( "#4FAFB7", "#8034E4", "#ADB56F", "#6C84EC", "#E6B93F", "#CBE7EA",
                          "#F19D38", "#65AD8E","#EB5428", "#DAE650","#7E8C8E", "#F5E4B2","#DC4CA7",
                          "#AB7B90","#6AE68B", "#BBA4E6", "#F9CDC1", "#E7E9D4","#7E85D7","#FBE1C4",
                          "#DC9F86","#877780","#75E6E4","#6F9557","#D7A74E", "#D1D0EE",
                          "#63B8DA","#D3A1D5","#C445DF","#852F02"
            )
            ccccol<-rep(vecofcols,length(unique(pretab$thecolor)))
            colors<-ccccol [1:length(unique(pretab$thecolor))]
            
            if ("Unknown" %in% pretab$thecolor || "Unknown_or_Out_of_Subset" %in% pretab$thecolor) {
              colors<-c(colors[1:length(colors)-1],"gray95") 
            } else {}
            pretab<-pretab %>% arrange(desc(pretab$thecolor) )
            thecolor<-pretab$thecolor
            Subsetinfotableumapind(pretab)
            
          }
          
          if (input$RetainSubset2 == "Area" || input$RetainSubset2 == "Inverse") {
            ggplot(Subsetinfotableumapind(), aes_string(x = colnames(Subsetinfotableumapind() )[3],
                                                        y = colnames(Subsetinfotableumapind() )[4]#,
                                                        #color = factor(seurat_clusters)
            )) +
              
              geom_point(aes(colour=factor(Subsetinfotableumapind()$thecolor)),size = input$sizeumapindsamples,pch=16) + 
              scale_color_manual(values = colors, name = "cluster") + 
              ggtitle(input$Umapstoplot) +
              guides(colour = guide_legend(title="cluster",override.aes = list(size = 3))) + 
              theme_void() +
              theme(plot.title = element_text(hjust = 0.5, face="bold") )
          } else {
          ggplot(Subsetinfotableumapind(), aes_string(x = colnames(Subsetinfotableumapind() )[3],
                                                      y = colnames(Subsetinfotableumapind() )[4]#,
                                                      #color = factor(seurat_clusters)
                                                      )) +
            
            geom_point(aes(colour=factor(seurat_clusters)),size = input$sizeumapindsamples,pch=16) + 
            scale_color_manual(values = colors, name = "cluster") + 
            ggtitle(input$Umapstoplot) +
            guides(colour = guide_legend(title="cluster",override.aes = list(size = 3))) + 
            
            theme_void() +
            theme(plot.title = element_text(hjust = 0.5, face="bold") )
          }
        }) # Closing the plot or plots 
        
        
      })
    
    #####################
    
    ## CITE independent Tab CITEindep
    ######################
    observeEvent({
      input$colby}, {
       req(input$colby)
       req(infotable3())
       var_to_highlight<- sort(unique (infotable3()[,input$colby]))
       updateSelectInput(session, inputId = "highlight", choices = var_to_highlight )
       
     })
    
    observeEvent({ 
      input$ploteaCITEIndep }, {
        req(infotable2())
        req(input$colby)
        req(input$highlight)
        
        output$p2<- renderPlot({
          if (input$RetainSubset3 == "Area" ) { 
            
            pretab<-infotable2()@cite_emb
            clusters<- infotable3()[,input$colby];
            pretab$clusters<-clusters
            losin<-rownames(pretab[intersect(rownames(pretab),Tab4piece2()$CITEname),])
            pretab$losin<-ifelse(rownames(pretab) %in% losin,"in","out" )
            pretab$thecolor<- ifelse(pretab$losin == "out","Unknown_or_Out_of_Subset",pretab$clusters)
            PlotClusterCITE_Independently(infotable2() ,color_by = pretab$thecolor,pt_size =input$size,highlight =input$highlight, Color = input$color, Selection_on_Top=input$SelectiononTop) 
            
            } else if (input$RetainSubset3 == "Inverse") {
              pretab<-infotable2()@cite_emb
              clusters<- infotable3()[,input$colby];
              pretab$clusters<-clusters
              losin<-rownames(pretab[outersect(rownames(pretab),Tab4piece2()$CITEname),])
              pretab$losin<-ifelse(rownames(pretab) %in% losin,"in","out" )
              pretab$thecolor<- ifelse(pretab$losin == "out","Unknown_or_Out_of_Subset",pretab$clusters)
              PlotClusterCITE_Independently(infotable2() ,color_by = pretab$thecolor,pt_size =input$size,highlight =input$highlight, Color = input$color, Selection_on_Top=input$SelectiononTop)
            
            } else {
            thecolor<- infotable3()[,input$colby];
            PlotClusterCITE_Independently(infotable2() ,color_by = thecolor,pt_size =input$size,highlight =input$highlight, Color = input$color, Selection_on_Top=input$SelectiononTop)  
          }
          
        }) # Closing the plot p2
      })
    
    ######################

    ## CITEtoCODEXspatial extrapolation
    #####################
    IdThrs<- reactiveVal(NULL);
    
    observeEvent({
      input$ploteaExtrapolationTab}, {
        
        req(input$IdentificationThreshold)
        IdThrs ( as.numeric(as.character(input$IdentificationThreshold)) )
        req( input$cluster1)

        colnameVect<-c("orig.ident","assay_name","integrated_snn_res.0.5","seurat_clusters",
                       "CITE_snn_res.1","dsb_knn_res.1.5","wsnn_res.1.5",
                       "predicted.id","predicted.celltype","Higher_Hierarchy_grouping",
                       "cite_clusters_after_integration")
        
        if(input$cluster1 %in% colnameVect & IdThrs() != 0 ){
          preinfotab4_tmp<-infotable4()
          preinfotab4_tmp$Prob_threshold<-ifelse(preinfotab4_tmp$probCITEcell<IdThrs() | is.na(preinfotab4_tmp$probCITEcell) ,"out","in") 
          preinfotab4_tmp[,input$cluster1]<-ifelse(preinfotab4_tmp$Prob_threshold == "in", preinfotab4_tmp[,input$cluster1],"Unknown_or_Below_threshold")
          print(table (preinfotab4_tmp[,input$cluster1]))
          
          preinfotab4_tmp[,input$cluster1]<-ifelse(preinfotab4_tmp[,input$cluster1] == "Unknown", "Unknown_or_Below_threshold",preinfotab4_tmp[,input$cluster1])
          
          infotab4_tmp(preinfotab4_tmp)
        } else {
          infotab4_tmp(infotable4())
        }
        
        
        output$p3<- renderPlot({

            colors <- randomcoloR::distinctColorPalette(length(unique(infotab4_tmp()[,input$cluster1]))) #This is new
            vecofcols<-c( "#4FAFB7", "#8034E4", "#ADB56F", "#6C84EC", "#E6B93F", "#ffdc00",
                          "#65AD8E", "#F19D38","#EB5428", "#DAE650","#7E8C8E", "#F5E4B2","#DC4CA7",
                          "#AB7B90","#6AE68B", "#BBA4E6", "#F9CDC1", "#E7E9D4","#7E85D7","#FBE1C4",
                          "#DC9F86","#877780","#75E6E4","#6F9557","#D7A74E", "#D1D0EE",
                          "#63B8DA","#D3A1D5","#C445DF","#852F02"
            )
    
            ccccol<-rep(vecofcols,length(unique(  infotab4_tmp()[,input$cluster1]   )))
            colors<-ccccol [1:length(unique( infotab4_tmp()[,input$cluster1]  ))]

            if (input$cluster1 == "Higher_Hierarchy_grouping" |input$cluster1 == "predicted.celltype"|input$cluster1 == "predicted.id"  & IdThrs() != 0 ){
              temporalsubsetinfotab4<-infotab4_tmp()
              unique(temporalsubsetinfotab4[,input$cluster1] )
              temporalsubsetinfotab4[,input$cluster1]<-ifelse(temporalsubsetinfotab4[,input$cluster1] == "Unknown","Unknown_or_Below_threshold",temporalsubsetinfotab4[,input$cluster1] )
              
              infotab4_tmp(temporalsubsetinfotab4)
              
              ccccol<-rep(vecofcols,length(unique(  infotab4_tmp()[,input$cluster1]   )))
              colors_temp<-ccccol [1:length(unique( infotab4_tmp()[,input$cluster1]  ))]
              
              pretabinfotab4_total<-infotable4()
              unique(pretabinfotab4_total[,input$cluster1] )
              pretabinfotab4_total[,input$cluster1]<-ifelse(pretabinfotab4_total[,input$cluster1] == "Unknown","Unknown_or_Below_threshold",pretabinfotab4_total[,input$cluster1] )
              theclasses<-sort(unique(pretabinfotab4_total[,input$cluster1] ))
              
              ccccol<-rep(vecofcols,length(unique(  pretabinfotab4_total[,input$cluster1]   )))
              colors_total<-ccccol [1:length(unique( pretabinfotab4_total[,input$cluster1]  ))]
              tabclassesandcolors<-data.frame(cbind(class=theclasses,color=colors_total) )
              print(tabclassesandcolors)
              smallertabclassesandcolors<-tabclassesandcolors[tabclassesandcolors$class %in% unique(infotab4_tmp()[,input$cluster1] ),]
              print(smallertabclassesandcolors)
              colors<-smallertabclassesandcolors$color
              
            }
            
            table (infotab4_tmp()[,input$cluster1] )
            
            if ("Unknown" %in% infotab4_tmp()[,input$cluster1] ) {
              colors<-c(colors[1:length(colors)-1],"gray95")
              fills<-c(colors[1:length(colors)-1],"gray95")  
            } else {}
            if ("Unknown_or_Below_threshold" %in% infotab4_tmp()[,input$cluster1] ) {
              colors<-c(colors[1:length(colors)-1],"transparent")
              fills<-c(colors[1:length(colors)-1],"transparent")  
            } else {}
            
            thecolor2<- infotab4_tmp()[,input$cluster1];
            tab<-melt(table(thecolor2))
            colnames(tab)<- c("Category", "Cell_number")
            tabproportionsCODEX (tab)

            
            if (input$RetainSubset4 == "Area" ) { 
              
              pretab<-infotab4_tmp()[,c(1,2)]
              pretab<-cbind(pretab, prethecolor2=thecolor2)
              losin<-intersect(pretab$CODEXname,Tab4piece2()$CODEXname)
              
              pretab$losin<-ifelse(pretab$CODEXname %in% losin,"in","out" )
              pretab$thecolor2<- ifelse(pretab$losin == "out","Unknown_or_Out_of_Subset",pretab$prethecolor2)
              print(table(pretab$thecolor2))
              thecolor2<-pretab$thecolor2
              vecofcols<-c( "#4FAFB7", "#8034E4", "#ADB56F", "#6C84EC", "#E6B93F", "#ffdc00",
                            "#65AD8E", "#F19D38","#EB5428", "#DAE650","#7E8C8E", "#F5E4B2","#DC4CA7",
                            "#AB7B90","#6AE68B", "#BBA4E6", "#F9CDC1", "#E7E9D4","#7E85D7","#FBE1C4",
                            "#DC9F86","#877780","#75E6E4","#6F9557","#D7A74E", "#D1D0EE",
                            "#63B8DA","#D3A1D5","#C445DF","#852F02"
              )
              
              ccccol<-rep(vecofcols,length(unique(  pretab$thecolor2   )))
              colors<-ccccol [1:length(unique( pretab$thecolor2   ))]
              if ("Unknown" %in% pretab$thecolor2 || "Unknown_or_Out_of_Subset" %in% pretab$thecolor2 && "Unknown_or_Below_threshold" %in% pretab$thecolor2) { 
                print(print(table(pretab$thecolor2)))
                colors<-c(colors[1:(length(colors)-2)],"transparent","gray95") 
                print(colors)
              } else if ("Unknown" %in% pretab$thecolor2 || "Unknown_or_Out_of_Subset" %in% pretab$thecolor2 && !"Unknown_or_Below_threshold" %in% pretab$thecolor2) {
                print(print(table(pretab$thecolor2)))
                colors<-c(colors[1:length(colors)-1],"gray95")
                print(colors)
              } else if (!"Unknown" %in% pretab$thecolor2 || !"Unknown_or_Out_of_Subset" %in% pretab$thecolor2 && "Unknown_or_Below_threshold" %in% pretab$thecolor2) {
                print(print(table(pretab$thecolor2)))
                colors<-c(colors[1:length(colors)-1],"transparent")
                print(colors)
              } else {}
              
              tab<-melt(table(thecolor2))
              colnames(tab)<- c("Category", "Cell_number")
              tabproportionsCODEX (tab)
              
            } else if (input$RetainSubset4 == "Inverse") {
              
              pretab<-infotab4_tmp()[,c(1,2)]
              pretab<-cbind(pretab, prethecolor2=thecolor2)
              losin<-outersect(pretab$CODEXname,Tab4piece2()$CODEXname)
              
              pretab$losin<-ifelse(pretab$CODEXname %in% losin,"in","out" )
              pretab$thecolor2<- ifelse(pretab$losin == "out","Unknown_or_Out_of_Subset",pretab$prethecolor2)
              print(table(pretab$thecolor2))
              thecolor2<-pretab$thecolor2
              vecofcols<-c( "#4FAFB7", "#8034E4", "#ADB56F", "#6C84EC", "#E6B93F", "#ffdc00",
                            "#65AD8E", "#F19D38","#EB5428", "#DAE650","#7E8C8E", "#F5E4B2","#DC4CA7",
                            "#AB7B90","#6AE68B", "#BBA4E6", "#F9CDC1", "#E7E9D4","#7E85D7","#FBE1C4",
                            "#DC9F86","#877780","#75E6E4","#6F9557","#D7A74E", "#D1D0EE",
                            "#63B8DA","#D3A1D5","#C445DF","#852F02"
              )
              
              ccccol<-rep(vecofcols,length(unique(  pretab$thecolor2   )))
              colors<-ccccol [1:length(unique( pretab$thecolor2   ))]
              if ("Unknown" %in% pretab$thecolor2 || "Unknown_or_Out_of_Subset" %in% pretab$thecolor2 && "Unknown_or_Below_threshold" %in% pretab$thecolor2) { 
                print(print(table(pretab$thecolor2)))
                colors<-c(colors[1:(length(colors)-2)],"transparent","gray95") 
                print(colors)
              } else if ("Unknown" %in% pretab$thecolor2 || "Unknown_or_Out_of_Subset" %in% pretab$thecolor2 && !"Unknown_or_Below_threshold" %in% pretab$thecolor2) {
                print(print(table(pretab$thecolor2)))
                colors<-c(colors[1:length(colors)-1],"gray95")
                print(colors)
              } else if (!"Unknown" %in% pretab$thecolor2 || !"Unknown_or_Out_of_Subset" %in% pretab$thecolor2 && "Unknown_or_Below_threshold" %in% pretab$thecolor2) {
                print(print(table(pretab$thecolor2)))
                colors<-c(colors[1:length(colors)-1],"transparent")
                print(colors)
              } else {}
            
              tab<-melt(table(thecolor2))
              colnames(tab)<- c("Category", "Cell_number")
              tabproportionsCODEX (tab)
            } # Closing else if (input$RetainSubset4 == "Inverse")
            
            
            ggplot(info_spatial_tmp()) +
              geom_point(aes(x = x,y = y, color = factor(thecolor2)),size = input$size2) +
              scale_color_manual(values = colors, name = "cluster") +
              guides(colour = guide_legend(override.aes = list(size = 5))) +
              theme_void()
        }) # Closing the plot
        
        # This is new
        output$tproportionsCODEX<- DT::renderDT({
          tabproportionsCODEX()
        })# Closing the datatable of proportions of classification categories for CODEX cells
      }) # Closing ObserveEvent ploteaExtrapolationTab
    
    ######################
    
    ## CITEtoCODEXspatial independent
    #####################
    info_spatial_tmpindep <- reactiveVal(NULL);

    observeEvent({
      input$cluster2}, {
        req( input$cluster2)
        req(input$IdentificationThreshold)
        
        var_to_highlight2<- sort(unique (infotable4()[,input$cluster2]))
        updateSelectInput(session, inputId = "clusterincolor", choices = var_to_highlight2 )
        
        
      })

    observeEvent({
      input$plotea}, {
        req( input$clusterincolor)
        req(input$IdentificationThreshold)
        
        x_tmp <- infotable4()[, "x"]
        x_tmp <- x_tmp - min(x_tmp)
        y_tmp <- infotable4()[, "y"]
        y_tmp <- y_tmp - min(y_tmp)
        spatial_tmp <- as.data.frame(cbind(x = x_tmp, y = y_tmp))
        info_spatial_tmpindep(spatial_tmp)
        rm (spatial_tmp)
        print(input$clusterincolor)
       
        IdThrs = as.numeric(as.character(input$IdentificationThreshold))
        colnameVect<-c("orig.ident","assay_name","integrated_snn_res.0.5","seurat_clusters",
                       "CITE_snn_res.1","dsb_knn_res.1.5","wsnn_res.1.5",
                       "predicted.id","predicted.celltype","Higher_Hierarchy_grouping",
                       "cite_clusters_after_integration")
        
        if(input$cluster2 %in% colnameVect & IdThrs != 0 ){
          preinfotab4_tmp<-infotable4()
          preinfotab4_tmp$Prob_threshold<-ifelse(preinfotab4_tmp$probCITEcell<=IdThrs | is.na(preinfotab4_tmp$probCITEcell) ,"out","in") 
          preinfotab4_tmp[,input$cluster2]<-ifelse(preinfotab4_tmp$Prob_threshold == "in", preinfotab4_tmp[,input$cluster2],"Unknown_or_Below_threshold")
          
          preinfotab4_tmp[,input$cluster2]<-ifelse(preinfotab4_tmp[,input$cluster2] == "Unknown", "Unknown_or_Below_threshold",preinfotab4_tmp[,input$cluster2])
          
          infotab4_tmp(preinfotab4_tmp)
        } else {
          infotab4_tmp(infotable4())
        }
        
        #### RetainSubset5
        if (input$RetainSubset5 == "Area") {
          df <- infotab4_tmp()
          losin<-intersect(df$CODEXname,Tab4piece2()$CODEXname)
          df$losin<-ifelse(df$CODEXname %in% losin,"in","out" )
          df$condition2 <- ifelse (infotable4()[,input$cluster2] %in% input$clusterincolor & df$losin == "in", 1,0)
          print(table(df$condition2))
          df<-df[order(df$condition2),]
          infotable4(df)
          rm(df)

        } else if(input$RetainSubset5 == "Inverse") {
          df <- infotab4_tmp()
          losin<-outersect(df$CODEXname,Tab4piece2()$CODEXname)
          df$losin<-ifelse(df$CODEXname %in% losin,"in","out" )
          df$condition2 <- ifelse (infotable4()[,input$cluster2] %in% input$clusterincolor & df$losin == "in", 1,0)
          df<-df[order(df$condition2),]
          infotable4(df)
          rm(df)
          
        } else {
          df <- infotab4_tmp()
          df$condition2 <- ifelse (infotable4()[,input$cluster2] %in% input$clusterincolor, 1,0)
          df<-df[order(df$condition2),]
          infotable4(df)
          rm(df)
        }
        
        df<-info_spatial_tmpindep()
        df$condition2<-infotable4()$condition2
        df<-df[order(df$condition2),]
        info_spatial_tmpindep(df)
        rm(df)
        
        colorsindep <- randomcoloR::distinctColorPalette(length(unique(infotable4()[,input$cluster2])))         
        vecofcols<-c("#75B5D6","#D5E04F",
                     "#E8645A","#CBE6CC",
                     "#DD8DA3","#E5A645",
                     "#D7A6DF","#AB41E2",
                     "#6CE5C3","#899A5A",
                     "#E7DE94","#716E91",
                     "#D4CBDB","#89E951",
                     "#ADE392","#E44DB0",
                     "#7B99E0","#DDB296",
                     "#6CE892","#CD7BD7",
                     "#7C67D8","#86E2E1",
                     "#75E6E4","#6F9557",
                     "#D7A74E","#D1D0EE",
                     "#63B8DA","#D3A1D5",
                     "#C445DF","#852F02"
        )

        ccccol<-rep(vecofcols,length(unique(  infotab4_tmp()[,input$cluster2]   )))
        colorsindep<-ccccol [1:length(unique( infotab4_tmp()[,input$cluster2]  ))]

        groupsindep<- sort(unique (infotable4()[,input$cluster2]))

        colorsindep2<-data.frame(cbind(groupsindep,colorsindep))
        colorsindep2$theones<-ifelse(colorsindep2$groups %in% input$clusterincolor, colorsindep2$colorsindep,"gray95")
        
        if ("Unknown_or_Below_threshold" %in% infotab4_tmp()[,input$cluster2] ) {
          colorsindep2$theones<-c(colorsindep2$theones[1:length(colorsindep2$theones)-1],"transparent")
        } else {}
        
        if (input$RetainSubset5 == "Area" |input$RetainSubset5 == "Inverse") {
          if (IdThrs== 0) {
            thecolorsRetained<-ifelse(info_spatial_tmpindep()$condition2 ==1,infotable4()[,input$cluster2],"Unknown" )  
          } else {
            thecolorsRetained<-ifelse(info_spatial_tmpindep()$condition2 ==1,infotable4()[,input$cluster2],"Unknown_or_Below_threshold" )  
          }
          
          groupsindep<- sort(unique (infotable4()[,input$cluster2]))
          colorsindep2<-data.frame(cbind(groupsindep,colorsindep))
          
          allactors<-c(sort(unique(thecolorsRetained)),setdiff( groupsindep,   sort(unique(thecolorsRetained))))
          colorsindep2<-colorsindep2[match(allactors,colorsindep2$groupsindep),]
          colorsindep2$theones<-ifelse(colorsindep2$groups %in% sort(unique(thecolorsRetained))[1:length(unique(thecolorsRetained))-1], colorsindep2$colorsindep,"gray95")
        } else{
          thecolorsRetained<-infotable4()[,input$cluster2]
        }

        output$p4<- renderPlot({
          
          ppp<-ggplot(info_spatial_tmpindep()) +
            geom_point(aes(x = x,y = y, color = factor(thecolorsRetained)),size = input$sizespatial2) +
            scale_color_manual(values = colorsindep2$theones, name = "cluster") +
            guides(colour = guide_legend(override.aes = list(size = 5))) +
            theme_void()
          plot (ppp)
          plot (ppp)
        }) # Closing the plot
      })


    ######################

    ## DEsTabandUMaps
    #####################
    dataforDE <- reactiveVal(NULL);
    datadeg <- reactiveVal(NULL);
    datadegwithNonDE <- reactiveVal(NULL);
    datatouseinDEpart2 <- reactiveVal(NULL);
    
    observeEvent({
      input$colbyDE}, {
        req(input$colbyDE)
        req(infotable3())
        req(infotable4())

        if (input$colbyDE =="codex_clusters") {
          var_to_highlightDE<- sort(unique (infotable4()[,input$colbyDE]))
        } else {
          var_to_highlightDE<- sort(unique (infotable3()[,input$colbyDE]))    
          }

        updateSelectizeInput(inputId="selectDE", choices = var_to_highlightDE, options = list(maxItems = 2))
        
      })
    
    observeEvent({
      input$calculaDE}, {
        req(infotable1())
        req(infotable2())
        req(infotable4())
        req(input$colbyDE)
        req(input$cluster1)
        req(input$selectDE)
        req(input$datatypeDE)
      
        if (input$colbyDE =="codex_clusters") {
          tab42<-infotable4()[!(is.na(infotable4()$CITEname)),]
          predat1<-tab42[c('CITEname','codex_clusters')]
          
        } else {
          predat1<-infotable3()[,input$colbyDE]
        }
        

        if (input$datatypeDE =="cite_mRNA_norm" ){
          if (is.list(infotable2()@cite_mRNA_norm)) {
            print("Since there are too many cells in the initial data, the DE analysis will be performed in a sample of the data")
            predat2 <- do.call("rbind", infotable2()@cite_mRNA_norm)
            print("dim(predat2)")
            print(dim(predat2))
            set.seed(123)
            my.ids <- sample(rownames(predat2), replace= F, size= round(0.1 * nrow (predat2), digits = 0))
            length(my.ids)
            piece<- predat2[my.ids, ]; dim(piece)
            predat2<-data.frame(piece)
            
            predat1<-infotable3()[,input$colbyDE, drop=F]
            piecepredat1<-predat1[my.ids,];dim(piecepredat1)
            predat1<-piecepredat1
            
            
          } else {
            predat2<-infotable2()@cite_mRNA_norm  
            print("dim(predat2)")
            print(dim(predat2))  
          }
          
        } else {
          predat2<-infotable2()@cite_protein
        }
        if (is.null(dim(predat1))){
          if (class(predat2) == "dgCMatrix"){
            dataforDE(data.frame(cbind(predat1,data.frame(predat2) )) )  
          } else {
            dataforDE(data.frame(cbind(predat1,predat2)) )    
          }  
        } else {
          dataforDE(merge (predat1,predat2,by.x=1, by.y=0) ); dim(dataforDE());dim(predat1);dim(predat2)
          anyNA(dataforDE()$codex_clusters)
          dat4de<-dataforDE()[,-1]
          dataforDE(dat4de)
        }
        
        
        print(input$selectDE)
        
        
        
        output$tDE<- DT::renderDT({
          datadeg (DEgenesPairwiseComparison (X=dataforDE(), groups=c(input$selectDE)))
          print("datadeg()$DE")
          print(datadeg()$DE)
        })# Closing the datatable
        datadegwithNonDE <- reactiveVal(NULL)
        datatouseinDEpart2(datadeg()$DE)
        print("datatouseinDEpart2()")
        print(datatouseinDEpart2())
    })# Closing the observeEvent of the button
    
    
    observeEvent({
      datadeg()}, {
        req(datadeg())

        if (length(rownames(datadeg()$DE)) == 0) {
          stop (paste("There are no significantly different features between the studied groups"))
        } else { }

        nonDEvars<-sort(unique (rownames(datadeg()$nDE )))
        updateSelectizeInput(session,inputId="NonDEfeatures", choices = nonDEvars, server = TRUE)
        
        var_to_highlightUMAPDE<- sort(unique (rownames(datadeg()$DE )))
        print("var_to_highlightUMAPDE")
        print (var_to_highlightUMAPDE)
        if (length(var_to_highlightUMAPDE)<9 ) {
          maxit=length(var_to_highlightUMAPDE)
        } else {
          maxit=9
        }

        if(is.null(input$NonDEfeatures) ) {
          umapchoices <- var_to_highlightUMAPDE
        } else {
          umapchoices <- c(var_to_highlightUMAPDE,input$NonDEfeatures)
        }
        
        updateSelectizeInput(session, inputId="featureUMAPDE", choices = umapchoices, options = list(maxItems = maxit))
        
    })
    
    observeEvent({
      input$NonDEfeatures } , {
        
        if(is.null(input$NonDEfeatures)) {
        } else {
          posibilities<-c(sort(unique (rownames(datadeg()$DE ))),input$NonDEfeatures)
          preselection<-unique(c(input$featureUMAPDE,input$NonDEfeatures))
          updateSelectizeInput(session, inputId="featureUMAPDE", choices =posibilities , options = list(maxItems = 9),selected =preselection )
          
          addition <-datadeg()$nDE[rownames(datadeg()$nDE) %in% input$NonDEfeatures, ,drop=F]
          newtab <- rbind(datadeg()$DE,addition)
          datadegwithNonDE (newtab)
          
          output$tDE<- DT::renderDT({
            print(datadegwithNonDE())
          })# Closing the datatable
          datatouseinDEpart2(datadegwithNonDE())
          print("datatouseinDEpart2() inside the observeEvent of input$NonDEfeatures")
          print(datatouseinDEpart2())
        }
      })

    observeEvent({
      input$ploteaUMAPDE}, {
        req(infotable2())
        req(input$colbyDE)
        req(input$selectDE)
        req(input$datatypeDE)
        req(input$featureUMAPDE)
        req(input$minco)
        
        output$p6 <- renderPlot({
          
          if (input$datatypeDE =="cite_mRNA_norm" ){
            if (is.list (infotable2()@cite_mRNA_norm)) {
              print("cite_mRNA_norm is a list")
              predat2 <- do.call("rbind", infotable2()@cite_mRNA_norm)
              names<-input$featureUMAPDE
              names<-gsub('\\.', '-', names)
              predat2<-predat2[,colnames(predat2) %in% names]
              print("dim(predat2)")
              print(dim(predat2))
              UMAPFeatureExpression(infotable2(), plottingData=predat2,features=input$featureUMAPDE, low_color=input$colorUMAPDE2,type="RNA", high_color = input$colorUMAPDE, pt_size= input$sizeDE, min.cutoff=input$minco)  
              
            } else {
              UMAPFeatureExpression(infotable2(),features=input$featureUMAPDE, low_color=input$colorUMAPDE2,type="RNA", high_color = input$colorUMAPDE, pt_size= input$sizeDE, min.cutoff=input$minco)  
            }
            
          } else {
            print(input$colorUMAPDE)
            print(input$colorUMAPDE2)
            print(input$colorUMAPDE3)
            
            UMAPFeatureExpression(infotable2(),features=input$featureUMAPDE, type = "protein", high_color = input$colorUMAPDE, low_color=input$colorUMAPDE2,high_color2=input$colorUMAPDE3, pt_size= input$sizeDE,min.cutoff=input$minco)
          }

          }) # Closing the plot
      })# Closing the observeEvent of ploteaUMAPDE
        
    ######################
    
    ## DEsHeatmapRidgeandViolin
    #####################
    
    
    output$tDE2<- DT::renderDT({
      datatouseinDEpart2()
    })# Closing the datatable
    
    observeEvent({
      input$allorsomefeatures }, {
        req(datatouseinDEpart2())
        req(input$allorsomefeatures)
        
        if (input$allorsomefeatures == "All") {
          var_to_plotinHeatmap<- sort(unique (rownames(datatouseinDEpart2() )))
          print("var_to_plotinHeatmap")
          print (var_to_plotinHeatmap)
          updateSelectizeInput(session, inputId="selectfeaturesforheatmap", choices = var_to_plotinHeatmap, selected = var_to_plotinHeatmap )
          updateSelectizeInput(session, inputId="selectfeaturesforridgeandviolin", choices = var_to_plotinHeatmap)
        } else {
          var_to_plotinHeatmap<- sort(unique (rownames(datatouseinDEpart2() )))
          print("var_to_plotinHeatmap")
          print (var_to_plotinHeatmap)
          updateSelectizeInput(session, inputId="selectfeaturesforheatmap", choices = var_to_plotinHeatmap)
          updateSelectizeInput(session, inputId="selectfeaturesforridgeandviolin", choices = var_to_plotinHeatmap)
        }
      })
    
    observeEvent({
      input$allorsomeGroups }, {
        req(input$allorsomeGroups)
        req(input$colbyDE)
        
        if (input$allorsomeGroups == "All") {
          if (input$colbyDE =="codex_clusters") {
            var_to_highlightDE<- sort(unique (infotable4()[,input$colbyDE]))  
          } else {
            var_to_highlightDE<- sort(unique (infotable3()[,input$colbyDE]))  
          }
          
          updateSelectizeInput(inputId="selectDE", choices = var_to_highlightDE, options = list(maxItems = 2))
          print("var_to_highlightDE")
          print (var_to_highlightDE)
          updateSelectizeInput(session, inputId="selectgroupsforheatmap", choices = var_to_highlightDE, selected = var_to_highlightDE )
          updateSelectizeInput(session, inputId="selectgroupsforridgeandviolin", choices = var_to_highlightDE )
        } else {
          if (input$colbyDE =="codex_clusters") {
            var_to_highlightDE<- sort(unique (infotable4()[,input$colbyDE]))  
          } else {
            var_to_highlightDE<- sort(unique (infotable3()[,input$colbyDE]))  
          }
          updateSelectizeInput(inputId="selectDE", choices = var_to_highlightDE, options = list(maxItems = 2))
          print("var_to_highlightDE")
          print (var_to_highlightDE)
          updateSelectizeInput(session, inputId="selectgroupsforheatmap", choices = var_to_highlightDE )
          updateSelectizeInput(session, inputId="selectgroupsforridgeandviolin", choices = var_to_highlightDE )
        }
    })

    observeEvent({
      input$selectgroupsforheatmap}, {
        
        req(input$selectgroupsforheatmap)
        var_to_sort<- input$selectgroupsforheatmap
        updateSelectizeInput(inputId="selectgroupsorder", choices = var_to_sort)
      })
    
    observeEvent({
      input$ploteaHeatmap}, {
        req(dataforDE())
        req(input$selectfeaturesforheatmap)
        req(input$selectgroupsforheatmap)
        req(input$heatmapcolor)
        req(input$lowerbreak)
        req(input$upperbreak)
        req(input$heatmapPlotsize)
        req(input$heatmapLegendsize)
        req(input$heatmapFeaturesize)
        

        output$p7<- renderPlot({
          HeatmapbyGroup(Data = dataforDE(),Genes=input$selectfeaturesforheatmap,Groups = input$selectgroupsforheatmap,Heatmap_Color = input$heatmapcolor,Ncells = input$cellnumber,Breaks = c(input$lowerbreak,input$upperbreak),group_order = c(input$selectgroupsorder),
                         legend_size=input$heatmapLegendsize,feature_size = input$heatmapFeaturesize )
        },height =input$heatmapPlotsize)# Closing the p7
      })# Closing the observeEvent of the button ploteaHeatmap
    
    
    observeEvent({
      input$plotearidgeandviolin}, {
        req(dataforDE())
        req(input$selectfeaturesforridgeandviolin)
        req(input$selectgroupsforridgeandviolin)
        p<-RidgeplotbyGroup(Data=dataforDE(),Genes=input$selectfeaturesforridgeandviolin,Groups=input$selectgroupsforridgeandviolin)
        print("length(p)")
        print(length(p))
        
        if (length(p) == 3) {
            output$p8<- renderPlot({
              p[[1]]
            })
            output$p9<- renderPlot({
              p[[2]]
            })
            output$p10<- renderPlot({
              p[[3]]
            })
          } else if (length(p) == 2) {
            output$p8<- renderPlot({
              p[[1]]
            })
            output$p9<- renderPlot({
              p[[2]]
            })
            output$p10<- renderPlot({
              NULL
            })
          }
        
        q<-ViolinplotbyGroup(Data=dataforDE(),Genes=input$selectfeaturesforridgeandviolin,Groups=input$selectgroupsforridgeandviolin)
        
        if (length(p) == 3) {
          output$p11<- renderPlot({
            q[[1]]
          })
          output$p12<- renderPlot({
            q[[2]]
          })
          output$p13<- renderPlot({
            q[[3]]
          })
        } else if (length(p) == 2) {
          output$p11<- renderPlot({
            q[[1]]
          })
          output$p12<- renderPlot({
            q[[2]]
          })
          output$p13<- renderPlot({
            NULL
          })
        }
        
      })# Closing the observeEvent of the button plotearidgeandviolin

    ######################
    
    ## Spatial Expression (CODEX spatial analysis)
    #####################
    observeEvent({
      
      input$clusterSpatialExpression}, {

        req(input$clusterSpatialExpression)
        req(infotable2())
        req(input$IdentificationThreshold2)
        
        var_to_highlightSpatialExpression<- sort(unique (infotable4()[,input$clusterSpatialExpression]))
        updateSelectInput(session, inputId = "clusterincolorSpatialExpression", choices = var_to_highlightSpatialExpression )
        print("input$clusterSpatialExpression")
        print(input$clusterSpatialExpression)
        print("input$clusterincolorSpatialExpression")
        print(input$clusterincolorSpatialExpression)
        
        
      })
    
    observeEvent({
      input$datatypeCODEXspatial }, {
        req(req(infotable2()))
        req(input$datatypeCODEXspatial)
        print("input$datatypeCODEXspatial")
        print(input$datatypeCODEXspatial)

        if (input$datatypeCODEXspatial =="protein" ){
          if (!is.null(infotable2()@codex_clean)) {
            var_to_highlightCODEXspatial <- colnames(infotable2()@codex_clean)
          }
          else {
            var_to_highlightCODEXspatial <- colnames(infotable2()@codex_protein)
          }
          
        } else{ 
          var_to_highlightCODEXspatial <- colnames(infotable2()@codex_mRNA)
          
        }
        
        updateSelectizeInput(inputId="selectfeaturesCODEXspatial", choices = var_to_highlightCODEXspatial, options = list(maxItems = 2))
        
      })
    
    
    observeEvent({
      input$ploteaCODEXspatial}, {
        req(infotable2())
        req(input$clusterSpatialExpression)
        req(input$clusterincolorSpatialExpression)
        req(input$datatypeCODEXspatial)
        req(input$selectfeaturesCODEXspatial)
        req(input$mincoCODEXspatial)
        req(input$IdentificationThreshold2)
        
        IdThrs = as.numeric(as.character(input$IdentificationThreshold2))
        
        if (is.null (input$clusterincolorSpatialExpression)) {
          preinfotab4_tmp<-infotable4()
        } else{
          preinfotab4_tmp<-infotable4()
          preinfotab4_tmp<-preinfotab4_tmp[preinfotab4_tmp[,input$clusterSpatialExpression] %in% input$clusterincolorSpatialExpression,]
          dim(infotable4());dim(preinfotab4_tmp)
        }
        
        preinfotab4_tmp$Prob_threshold<-ifelse(preinfotab4_tmp$probCITEcell<=IdThrs | is.na(preinfotab4_tmp$probCITEcell) ,"out","in")
        preinfotab4_tmp[,input$cluster2]<-ifelse(preinfotab4_tmp$Prob_threshold == "in", preinfotab4_tmp[,input$cluster2],"Unknown_or_Below_threshold")

        preinfotab4_tmp[,input$cluster2]<-ifelse(preinfotab4_tmp[,input$cluster2] == "Unknown", "Unknown_or_Below_threshold",preinfotab4_tmp[,input$cluster2])
        precodexabovethrs<-preinfotab4_tmp[preinfotab4_tmp$Prob_threshold == "in",]
        codexabovethrs<-precodexabovethrs$CODEXname
        length(codexabovethrs)
        

        output$p14<- renderPlot({
          
          var_to_highlightSpatialExpression<- sort(unique (infotable4()[,input$clusterSpatialExpression]))
          if (identical(sort(var_to_highlightSpatialExpression),sort(input$clusterincolorSpatialExpression) ) ) {
            onlyshowclusters = "no"  
          } else {
            onlyshowclusters = "yes"  
          }
          
          
          
          if (input$TGOFanalysis == "Yes") {
           
            if (input$OnlyShowGroups == "Yes") {
              onlyshowclusters="yes"
            } else {
              onlyshowclusters="no"
            }
            
            if (input$RetainSubset6 == "Area" & input$Context =="Isolated"  ) {
              codexabovethrs<-intersect(codexabovethrs,Tab4piece2()$CODEXname)
              length(codexabovethrs)
              vecgroups<-preinfotab4_tmp[,input$clusterSpatialExpression][intersect(preinfotab4_tmp$CODEXname,Tab4piece2()$CODEXname)]
              length(vecgroups)
              
              withProgress(message = 'Plotting... please wait!', value = 0, {
                  PlotExprCODEXspatial_2groups1feature(infotable2(),infotab= infotable4(),classification_category=input$clusterSpatialExpression ,name= input$selectfeaturesCODEXspatial, clusters=input$clusterincolorSpatialExpression, type=input$datatypeCODEXspatial, high_color = input$colorCODEXspatial, low_color=input$colorCODEXspatial2,high_color2=input$colorCODEXspatial3, pt_size= input$sizefeatCODEXspatial, CellsaboveThrs=codexabovethrs, IdentificationThreshold =input$IdentificationThreshold2,only_show_clusters= onlyshowclusters,vectorofGroups=vecgroups,maximize_differences= input$maximize_differences,min.cutoff = input$mincoCODEXspatial)
                })
            } else if (input$RetainSubset6 == "Inverse" & input$Context =="Isolated" ) {
              codexabovethrs<-outersect(codexabovethrs,Tab4piece2()$CODEXname)
              vecgroups<-preinfotab4_tmp[,input$clusterSpatialExpression][outersect(codexabovethrs,Tab4piece2()$CODEXname)]
              withProgress(message = 'Plotting... please wait!', value = 0, {
                PlotExprCODEXspatial_2groups1feature(infotable2(),infotab= infotable4(), classification_category=input$clusterSpatialExpression ,name= input$selectfeaturesCODEXspatial, clusters=input$clusterincolorSpatialExpression, type=input$datatypeCODEXspatial, high_color = input$colorCODEXspatial, low_color=input$colorCODEXspatial2,high_color2=input$colorCODEXspatial3, pt_size= input$sizefeatCODEXspatial, CellsaboveThrs=codexabovethrs, IdentificationThreshold =input$IdentificationThreshold2,only_show_clusters= onlyshowclusters,vectorofGroups=vecgroups,maximize_differences= input$maximize_differences,min.cutoff = input$mincoCODEXspatial)
              }) 
            } else if (input$RetainSubset6 == "Area" & input$Context =="In_Space" ) {
              vecgroups<-preinfotab4_tmp[,input$clusterSpatialExpression][intersect(codexabovethrs,Tab4piece2()$CODEXname)]
              cellsofinterest<-Tab4piece2()$CODEXname
              withProgress(message = paste0 ('Plotting',input$RetainSubset6,input$Context, '... please wait!'), value = 0, {
                PlotExprCODEXspatial_2groups1feature(infotable2(),infotab= infotable4(), classification_category=input$clusterSpatialExpression , name= input$selectfeaturesCODEXspatial, clusters=input$clusterincolorSpatialExpression, type=input$datatypeCODEXspatial, high_color = input$colorCODEXspatial, low_color=input$colorCODEXspatial2,high_color2=input$colorCODEXspatial3, pt_size= input$sizefeatCODEXspatial, CellsaboveThrs=codexabovethrs, IdentificationThreshold =input$IdentificationThreshold2,only_show_clusters= onlyshowclusters,vectorofGroups=vecgroups,CellsofInterest=cellsofinterest,maximize_differences= input$maximize_differences,min.cutoff = input$mincoCODEXspatial)

              })
            } else if(input$RetainSubset6 == "Inverse" & input$Context =="In_Space" ) {
              df<-infotable4()
              vecgroups<-preinfotab4_tmp[,input$clusterSpatialExpression][intersect(codexabovethrs,Tab4piece2()$CODEXname)]
              cellsofinterest<-outersect(df$CODEXname,Tab4piece2()$CODEXname)
              withProgress(message = 'Plotting... please wait!', value = 0, {
                PlotExprCODEXspatial_2groups1feature(infotable2(),infotab= infotable4(),classification_category=input$clusterSpatialExpression , name= input$selectfeaturesCODEXspatial, clusters=input$clusterincolorSpatialExpression, type=input$datatypeCODEXspatial, high_color = input$colorCODEXspatial, low_color=input$colorCODEXspatial2,high_color2=input$colorCODEXspatial3, pt_size= input$sizefeatCODEXspatial, CellsaboveThrs=codexabovethrs, IdentificationThreshold =input$IdentificationThreshold2,only_show_clusters= onlyshowclusters,vectorofGroups=vecgroups,CellsofInterest=cellsofinterest,maximize_differences= input$maximize_differences,min.cutoff = input$mincoCODEXspatial)
                
              })
              
            }
          

          } else {

            if (input$RetainSubset6 == "Area" & input$Context =="Isolated"  ) {
              codexabovethrs<-intersect(codexabovethrs,Tab4piece2()$CODEXname)
              dim(Tab4piece2())
              length(codexabovethrs)
              vecgroups<-preinfotab4_tmp[,input$clusterSpatialExpression][intersect(codexabovethrs,Tab4piece2()$CODEXname)]
              length(vecgroups)
              
             withProgress(message = 'Plotting... please wait!', value = 0, {
                PlotExprCODEXspatial_v5(infotable2(),infotab= infotable4(),classification_category=input$clusterSpatialExpression , name= input$selectfeaturesCODEXspatial, clusters=input$clusterincolorSpatialExpression, type=input$datatypeCODEXspatial, high_color = input$colorCODEXspatial, low_color=input$colorCODEXspatial2,high_color2=input$colorCODEXspatial3, pt_size= input$sizefeatCODEXspatial, CellsaboveThrs=codexabovethrs, IdentificationThreshold =input$IdentificationThreshold2,only_show_clusters= onlyshowclusters,vectorofGroups=vecgroups,CellsofInterest=cellsofinterest,maximize_differences= input$maximize_differences,min.cutoff = input$mincoCODEXspatial)
                
              })
              
              
            } else if (input$RetainSubset6 == "Inverse" & input$Context =="Isolated" ) {
              codexabovethrs<-outersect(codexabovethrs,Tab4piece2()$CODEXname)
              vecgroups<-preinfotab4_tmp[,input$clusterSpatialExpression][intersect(codexabovethrs,Tab4piece2()$CODEXname)]
              
              withProgress(message = 'Plotting... please wait!', value = 0, {
                PlotExprCODEXspatial_v5(infotable2(),infotab= infotable4(),classification_category=input$clusterSpatialExpression , name= input$selectfeaturesCODEXspatial, clusters=input$clusterincolorSpatialExpression, type=input$datatypeCODEXspatial, high_color = input$colorCODEXspatial, low_color=input$colorCODEXspatial2,high_color2=input$colorCODEXspatial3, pt_size= input$sizefeatCODEXspatial, CellsaboveThrs=codexabovethrs, IdentificationThreshold =input$IdentificationThreshold2,only_show_clusters= onlyshowclusters,vectorofGroups=vecgroups,CellsofInterest=cellsofinterest,maximize_differences= input$maximize_differences,min.cutoff = input$mincoCODEXspatial)
              }) 
            } else if (input$RetainSubset6 == "Area" & input$Context =="In_Space" ) {
              vecgroups<-preinfotab4_tmp[,input$clusterSpatialExpression][intersect(codexabovethrs,Tab4piece2()$CODEXname)]
              cellsofinterest<-Tab4piece2()$CODEXname
              withProgress(message = 'Plotting... please wait!', value = 0, {
                PlotExprCODEXspatial_v5(infotable2(),infotab= infotable4(),classification_category=input$clusterSpatialExpression , name= input$selectfeaturesCODEXspatial, clusters=input$clusterincolorSpatialExpression, type=input$datatypeCODEXspatial, high_color = input$colorCODEXspatial, low_color=input$colorCODEXspatial2,high_color2=input$colorCODEXspatial3, pt_size= input$sizefeatCODEXspatial, CellsaboveThrs=codexabovethrs, IdentificationThreshold =input$IdentificationThreshold2,only_show_clusters= onlyshowclusters,vectorofGroups=vecgroups,CellsofInterest=cellsofinterest,maximize_differences= input$maximize_differences,min.cutoff = input$mincoCODEXspatial)
                
              })
            } else if(input$RetainSubset6 == "Inverse" & input$Context =="In_Space" ) {
              df<-infotable4()
              vecgroups<-preinfotab4_tmp[,input$clusterSpatialExpression][intersect(codexabovethrs,Tab4piece2()$CODEXname)]
              cellsofinterest<-outersect(df$CODEXname,Tab4piece2()$CODEXname)
              withProgress(message = 'Plotting... please wait!', value = 0, {
                PlotExprCODEXspatial_v5(infotable2(),infotab= infotable4(),classification_category=input$clusterSpatialExpression , name= input$selectfeaturesCODEXspatial, clusters=input$clusterincolorSpatialExpression, type=input$datatypeCODEXspatial, high_color = input$colorCODEXspatial, low_color=input$colorCODEXspatial2,high_color2=input$colorCODEXspatial3, pt_size= input$sizefeatCODEXspatial, CellsaboveThrs=codexabovethrs, IdentificationThreshold =input$IdentificationThreshold2,only_show_clusters= onlyshowclusters,vectorofGroups=vecgroups,CellsofInterest=cellsofinterest,maximize_differences= input$maximize_differences,min.cutoff = input$mincoCODEXspatial)
              })
              
            } else {
              withProgress(message = 'Plotting... please wait!', value = 0, {
                PlotExprCODEXspatial_v4(infotable2(), name = input$selectfeaturesCODEXspatial,Subset = "no", type = input$datatypeCODEXspatial, high_color = input$colorCODEXspatial, low_color=input$colorCODEXspatial2,high_color2=input$colorCODEXspatial3, pt_size= input$sizefeatCODEXspatial, SmoothScatterplot=FALSE, mult=1, maxPercentile=0.9,min.cutoff = input$mincoCODEXspatial, maximize_differences= input$maximize_differences, CellsaboveThrs=codexabovethrs, IdentificationThreshold =input$IdentificationThreshold2)
              })
            }
          }
        
          
          
      }) # Closing the plot
    })# Closing the observeEvent of ploteaCODEXspatial
    ######################
    
    # Expression CODEX UMAP analysis
    #####################
    
    observeEvent({
      input$ploteaCODEXUMAPclusters}, {
        req(input$ploteaCODEXUMAPclusters)
        req(infotable2())
        req(infotable4())
        req(input$selectclusteringtypeUMAPclusters)
        req(input$sizefeatCODEXUMAPclusters)
        req(input$IdentificationThreshold3)
        req(input$SelectiononTopCODEXUMAPAnalysis)

        IdThrs = as.numeric(as.character(input$IdentificationThreshold3))
        colnameVect<-c("orig.ident","assay_name","integrated_snn_res.0.5","seurat_clusters",
                       "CITE_snn_res.1","dsb_knn_res.1.5","wsnn_res.1.5",
                       "predicted.id","predicted.celltype","Higher_Hierarchy_grouping",
                       "cite_clusters_after_integration")
        
        if(input$selectclusteringtypeUMAPclusters %in% colnameVect & IdThrs != 0 ){
          preinfotab4_tmp<-infotable4()
          preinfotab4_tmp$Prob_threshold<-ifelse(preinfotab4_tmp$probCITEcell<=IdThrs | is.na(preinfotab4_tmp$probCITEcell) ,"out","in") 
          preinfotab4_tmp[,input$selectclusteringtypeUMAPclusters]<-ifelse(preinfotab4_tmp$Prob_threshold == "in", preinfotab4_tmp[,input$selectclusteringtypeUMAPclusters],"Unknown_or_Below_threshold")
          
          preinfotab4_tmp[,input$selectclusteringtypeUMAPclusters]<-ifelse(preinfotab4_tmp[,input$selectclusteringtypeUMAPclusters] == "Unknown", "Unknown_or_Below_threshold",preinfotab4_tmp[,input$selectclusteringtypeUMAPclusters])
          
          infotab4_tmp(preinfotab4_tmp)
        } else {
          preinfotab4_tmp<-infotable4()
          preinfotab4_tmp$Prob_threshold<-ifelse(is.na(preinfotab4_tmp$probCITEcell) ,"out","in") 
          infotab4_tmp(preinfotab4_tmp)
        }
        
        vecofcols<-c("#75B5D6","#D5E04F",
                     "#E7DE94","#CBE6CC",
                     "#DD8DA3","#E5A645",
                     "#D7A6DF","#AB41E2",
                     "#6CE5C3","#899A5A",
                     "#E8645A","#716E91",
                     "#D4CBDB","#89E951",
                     "#ADE392","#E44DB0",
                     "#7B99E0","#DDB296",
                     "#6CE892","#CD7BD7",
                     "#7C67D8","#86E2E1",
                     "#75E6E4","#6F9557",
                     "#D7A74E","#D1D0EE",
                     "#63B8DA","#D3A1D5",
                     "#C445DF","#852F02"
        )
        
        thecolor2<- infotab4_tmp()[,input$selectclusteringtypeUMAPclusters]
        tab<-melt(table(thecolor2))
        colnames(tab)<- c("Category", "Cell_number")
        tabproportionsCODEX2 (tab)
        
        output$p19<- renderPlot({
          withProgress(message = 'Plotting... please wait!', value = 0, {
            
            
            PlotClusterCODEXemb_vGator(infotable2(), infotable=infotab4_tmp(), cluster_column=input$selectclusteringtypeUMAPclusters, pt_size=input$sizefeatCODEXUMAPclusters, Selection_on_Top=input$SelectiononTopCODEXUMAPAnalysis,wished_colors=vecofcols )
          })
          
        }) # Closing the plot
        output$tproportionsCODEX2<- DT::renderDT({
          tabproportionsCODEX2()
        })# Closing the datatable of proportions of classification categories for CODEX cells
        
      })# Closing the observeEvent of ploteaCODEXUMAPclusters
        
 
    ######################
    
    ## CODEX independent Tab CODEXindep
    ######################
    observeEvent({
      input$selectclusteringtypeUMAPclusters}, {
        req(input$selectclusteringtypeUMAPclusters)
        req(infotable2())
        req(infotable4())
        var_to_highlightCODEXindep<- sort(unique (infotable4()[,input$selectclusteringtypeUMAPclusters]))
        updateSelectInput(session, inputId = "highlightCODEXindep", choices = var_to_highlightCODEXindep )
        
      })
    
    observeEvent({ 
      input$ploteaCODEXindep }, {
        req(infotable2())
        req(input$selectclusteringtypeUMAPclusters)
        req(input$highlightCODEXindep)
        req(input$IdentificationThreshold4)
        
        IdThrs = as.numeric(as.character(input$IdentificationThreshold4))

        colnameVect<-c("orig.ident","assay_name","integrated_snn_res.0.5","seurat_clusters",
                       "CITE_snn_res.1","dsb_knn_res.1.5","wsnn_res.1.5",
                       "predicted.id","predicted.celltype","Higher_Hierarchy_grouping",
                       "cite_clusters_after_integration")
        
        if(input$selectclusteringtypeUMAPclusters %in% colnameVect & IdThrs != 0 ){
          preinfotab4_tmp<-infotable4()
          preinfotab4_tmp$Prob_threshold<-ifelse(preinfotab4_tmp$probCITEcell<=IdThrs | is.na(preinfotab4_tmp$probCITEcell) ,"out","in") 
          preinfotab4_tmp[,input$selectclusteringtypeUMAPclusters]<-ifelse(preinfotab4_tmp$Prob_threshold == "in", preinfotab4_tmp[,input$selectclusteringtypeUMAPclusters],"Unknown_or_Below_threshold")
          
          preinfotab4_tmp[,input$selectclusteringtypeUMAPclusters]<-ifelse(preinfotab4_tmp[,input$selectclusteringtypeUMAPclusters] == "Unknown", "Unknown_or_Below_threshold",preinfotab4_tmp[,input$selectclusteringtypeUMAPclusters])
          
          infotab4_tmp(preinfotab4_tmp)
        } else {
          preinfotab4_tmp<-infotable4()
          preinfotab4_tmp$Prob_threshold<-ifelse(is.na(preinfotab4_tmp$probCITEcell) ,"out","in") 
          infotab4_tmp(preinfotab4_tmp)
        }

        output$p20<- renderPlot({
          thecolor<- infotab4_tmp()[,input$selectclusteringtypeUMAPclusters];
          withProgress(message = 'Plotting... please wait!', value = 0, {
          PlotClusterCODEX_Independently(infotable2() ,infotable=infotab4_tmp(), color_by = thecolor,pt_size =input$sizeCODEX,highlight =input$highlightCODEXindep, Color = input$colorCODEXindep, cluster_column=input$selectclusteringtypeUMAPclusters, Selection_on_Top=input$SelectiononTopCODEX)
          }) # Closing the "withProgress"
        }) # Closing the plot p2
      })
    
    ######################
    
 
    ## AdjacencybyFeature
    #####################
    feature_pairs_original <- reactiveVal(NULL);
    feature_pairs <- reactiveVal(NULL);
    feature_pairsAREA <- reactiveVal(NULL);
    feature_pairsANTIAREA <- reactiveVal(NULL);
    feature_adj <- reactiveVal(NULL);
    feature_adjAREA <- reactiveVal(NULL);
    feature_adjANTIAREA <- reactiveVal(NULL);
    
    adjacency_table <- reactiveVal(NULL);
    prenameforadjtab<-reactiveVal(NULL);
    nameadjtab<-reactiveVal(NULL);
    thepathadj<-reactiveVal(NULL);
    FileLocationadj<-reactiveVal(NULL);
    
    prenameAREAforadjtab<-reactiveVal(NULL);
    nameAREAadjtab<-reactiveVal(NULL);
    nameANTIAREAadjtab<-reactiveVal(NULL);
    FileLocationAREAadj<-reactiveVal(NULL);
    FileLocationANTIAREAadj<-reactiveVal(NULL);
    
    
    observeEvent({
      input$datatypeAdjacencybyFeature}, {
        req(input$datatypeAdjacencybyFeature)
        req(infotable2())
        req(infotable4())

        if (input$datatypeAdjacencybyFeature =="protein" ){
          var_to_highlightAdjacencybyFeature <- colnames(infotable2()@codex_protein)
          updateSelectizeInput(session,inputId="selectfeaturesAdjacencybyFeature", choices = var_to_highlightAdjacencybyFeature, options = list(maxItems = 200))
          
        } else  { 
          var_to_highlightAdjacencybyFeature <- colnames(infotable2()@codex_mRNA)
          updateSelectizeInput(session,inputId="selectfeaturesAdjacencybyFeature", choices = var_to_highlightAdjacencybyFeature, options = list(maxItems = 200)) 
          
        }
        
        print("var_to_highlightAdjacencybyFeature")
        print(var_to_highlightAdjacencybyFeature)
        print("input$selectfeaturesAdjacencybyFeature")
        print(input$selectfeaturesAdjacencybyFeature)

      })

    
    
    observeEvent({
      input$ploteaAdjacencybyFeature}, {
        req(infotable2())
        req(infotable4())
        req(input$datatypeAdjacencybyFeature)
        req(input$selectfeaturesAdjacencybyFeature)
        # ##################
        # Here we will insert the new section in which we will calculate and
        # save the adjacency of several features
        
        prenameforadjtab(paste(substr(input$file1$name,1,nchar(input$file1$name)-4),substr(input$file2$name,1,nchar(input$file2$name)-6), sep= "_") )
    
          nameadjtab <- paste0("Adjacency_",input$datatypeAdjacencybyFeature,"_",prenameforadjtab(),".csv")
        
          print("nameadjtab")
          print(nameadjtab)
          
          thepathadj(getwd())
          print("thepathadj")
          print(thepathadj())
          ifelse(!dir.exists(paste0(thepathadj(),"/",prenameforadjtab())), dir.create(paste0(thepathadj(),"/",prenameforadjtab()), showWarnings = TRUE, recursive = FALSE, mode = "0777") , FALSE)

          FileLocationadj (paste0(thepathadj(),"/",prenameforadjtab(),"/",nameadjtab) )
          print("FileLocationadj")
          print(FileLocationadj())
          
          
          feature_pairs_original(t(combn(input$selectfeaturesAdjacencybyFeature,2)) )
          feature_pairs(t(combn(input$selectfeaturesAdjacencybyFeature,2)) )
          feature_adj_old_calculated<-NULL
          
          if (file.exists(FileLocationadj())) {
            feature_adj_old <- read.csv(FileLocationadj(), row.names=1)
            # Here it comes the step to select the pairs to be calculated!
            feature_adj_old$f_g<-paste(feature_adj_old$f,feature_adj_old$g, sep="_")
            feature_pairs( cbind(feature_pairs(),f_g=as.matrix(paste(feature_pairs()[,1],feature_pairs()[,2], sep="_"))) )

            feature_adj_old_calculated<-feature_adj_old[feature_adj_old$f_g %in% feature_pairs()[,3],c(1:5)]
            
            testeo<-try(feature_pairs()[!is.element(feature_pairs()[,3],feature_adj_old$f_g ),c(1,2),drop=FALSE],silent = TRUE)
            class(testeo)
            if (class(testeo) =="try-error") {
              print("All feature pairs will be calculated and added to the saved list" )
              feature_pairs(feature_pairs()[,c(1,2),drop=FALSE])
            } else{
              feature_pairs(feature_pairs()[!is.element(feature_pairs()[,3],feature_adj_old$f_g ),c(1,2),drop=FALSE ])
            }
            
            
            if (dim(feature_pairs())[1] == 0 ){
              print("All feature pairs were previously calculated and retrieved from saved data" )
              feature_pairs (feature_pairs_original() )
              feature_adj (feature_adj_old_calculated)
            } else{
              print(paste(c(nrow(feature_pairs()),"feature(s) pairs to be calculated:",t(feature_pairs())), collapse = " " ))
              withProgress(message = 'Calculating adjacency', value = 0, {
              feature_adj(AdjScore_calculator(infotable2(), feature_type=input$datatypeAdjacencybyFeature,Feature_Pairs=feature_pairs()) )
              })
              invadj<-feature_adj()[,c(2,1,3:5)]
              colnames(invadj)<-colnames(feature_adj())
              
              feature_adj_old<-feature_adj_old[,c(1:5)]
              
              feature_adj_old<-unique(rbind(feature_adj_old,feature_adj(),invadj) )
              write.csv(x=feature_adj_old,file=FileLocationadj())
              
              print ("Some adjacency features were calculated and saved with previous calculations")
              
            }
            
            
          } else {
            print ("Adjacency feature(s) are being calculated and saved for the first time") 
            withProgress(message = 'Calculating adjacency', value = 0, {
            feature_adj (AdjScore_calculator(infotable2(), feature_type=input$datatypeAdjacencybyFeature,Feature_Pairs=feature_pairs()) )
            })
            
            invadj<-feature_adj()[,c(2,1,3:5)]
            colnames(invadj)<-colnames(feature_adj())
            feature_adj_old<-rbind(feature_adj(),invadj)
            
            write.csv(x=feature_adj_old,file=FileLocationadj())
            print ("Some adjacency features were calculated and saved for the first time")
          }
          
        # ##################
        # Here I will insert the new section to calculate adjacency of certain areas and
        # save the adjacency of that area

          if (input$RetainSubset9 == "Area" || input$RetainSubset9 == "Inverse" || input$RetainSubset9 == "Both"){
            prenameAREAforadjtab(paste(substr(input$fileAreaSubsetted$name,1,nchar (input$fileAreaSubsetted$name)-4),
                                       substr(input$file1$name,1,nchar(input$file1$name)-4),substr(input$file2$name,1,nchar(input$file2$name)-6), sep= "__") )
            nameAREAadjtab <- paste0("Adj_Area_",input$datatypeAdjacencybyFeature,"_",prenameAREAforadjtab(),".csv")
            
            nameANTIAREAadjtab <- paste0("Adj_InverseArea_",input$datatypeAdjacencybyFeature,"_",prenameAREAforadjtab(),".csv")
            
            print("nameAREAadjtab")
            print(nameAREAadjtab)
            print("nameINVERSEAREAadjtab")
            print(nameANTIAREAadjtab)
            
            thepathadj(getwd())
            print("thepathadj")
            print(thepathadj())
            
            FileLocationAREAadj (paste0(thepathadj(),"/",prenameforadjtab(),"/",nameAREAadjtab) )
            FileLocationANTIAREAadj (paste0(thepathadj(),"/",prenameforadjtab(),"/",nameANTIAREAadjtab) )
            
            print("FileLocationAREAadj")
            print(FileLocationAREAadj())
            print("FileLocationANTIAREAadj")
            print(FileLocationANTIAREAadj())
            
            
            feature_pairs_original(t(combn(input$selectfeaturesAdjacencybyFeature,2)) )
            feature_pairs(t(combn(input$selectfeaturesAdjacencybyFeature,2)) )
            feature_pairsAREA(feature_pairs())
            feature_pairsANTIAREA(feature_pairs())
            
            feature_adj_old_calculatedAREA<-feature_adj_old_calculatedANTIAREA<-NULL
            
            
            if (file.exists(FileLocationAREAadj()) & file.exists(FileLocationANTIAREAadj()) ) {

              feature_adj_oldAREA <- read.csv(FileLocationAREAadj(), row.names=1)
              feature_adj_oldANTIAREA <- read.csv(FileLocationANTIAREAadj(), row.names=1)

              # Here it comes the step to select the pairs to be calculated!
              feature_adj_oldAREA$f_g<-paste(feature_adj_oldAREA$f,feature_adj_oldAREA$g, sep="_")
              feature_adj_oldANTIAREA$f_g<-paste(feature_adj_oldANTIAREA$f,feature_adj_oldANTIAREA$g, sep="_")
              feature_pairsAREA( cbind(feature_pairsAREA(),f_g=as.matrix(paste(feature_pairsAREA()[,1],feature_pairsAREA()[,2], sep="_"))) )
              feature_pairsANTIAREA( cbind(feature_pairsANTIAREA(),f_g=as.matrix(paste(feature_pairsANTIAREA()[,1],feature_pairsANTIAREA()[,2], sep="_"))) )
              
              feature_adj_old_calculatedAREA<-feature_adj_oldAREA[feature_adj_oldAREA$f_g %in% feature_pairsAREA()[,3],c(1:5)]
              feature_adj_old_calculatedANTIAREA<-feature_adj_oldANTIAREA[feature_adj_oldANTIAREA$f_g %in% feature_pairsANTIAREA()[,3],c(1:5)]
              
              testeo<-try(feature_pairsAREA()[!is.element(feature_pairsAREA()[,3],feature_adj_oldAREA$f_g ),c(1,2),drop=FALSE],silent = TRUE)
              class(testeo)
              if (class(testeo) =="try-error") {
                print("All feature pairs will be calculated and added to the saved list" )
                feature_pairsAREA(feature_pairsAREA()[,c(1,2),drop=FALSE])
              } else{
                feature_pairsAREA(feature_pairsAREA()[!is.element(feature_pairsAREA()[,3],feature_adj_oldAREA$f_g ),c(1,2),drop=FALSE ])
              }
              
              if (dim(feature_pairsAREA())[1] == 0 ){
                print("All feature pairs were previously calculated and retrieved from saved data" )
                feature_pairsAREA (feature_pairs_original() )
                feature_adjAREA (feature_adj_old_calculatedAREA)
              } else{
                print(paste(c(nrow(feature_pairsAREA()),"feature(s) pairs to be calculated:",t(feature_pairsAREA())), collapse = " " ))
                withProgress(message = 'Calculating adjacency', value = 0, {

                  feature_adjAREA(AdjScore_calculator(infotable2(), feature_type=input$datatypeAdjacencybyFeature,Feature_Pairs=feature_pairsAREA()) )
                })
                feature_adj_oldAREA<-feature_adj_oldAREA[,c(1:5)]
                feature_adj_oldAREA<-rbind(feature_adj_oldAREA,feature_adjAREA())
                write.csv(x=feature_adj_oldAREA,file=FileLocationAREAadj())
                
                print ("Some adjacency features of the area were calculated and added to previous calculations")
                
              }
              
              # ANTIAREA
              testeo<-try(feature_pairsANTIAREA()[!is.element(feature_pairsANTIAREA()[,3],feature_adj_oldANTIAREA$f_g ),c(1,2),drop=FALSE],silent = TRUE)
              class(testeo)
              if (class(testeo) =="try-error") {
                print("All feature pairs will be calculated and added to the saved list" )
                feature_pairsANTIAREA(feature_pairsANTIAREA()[,c(1,2),drop=FALSE])
              } else{
                feature_pairsANTIAREA(feature_pairsANTIAREA()[!is.element(feature_pairsANTIAREA()[,3],feature_adj_oldANTIAREA$f_g ),c(1,2),drop=FALSE ])
              }
              
              
              if (dim(feature_pairsANTIAREA())[1] == 0 ){
                print("All feature pairs were previously calculated and retrieved from saved data" )
                feature_pairsANTIAREA (feature_pairs_original() )
                feature_adjANTIAREA (feature_adj_old_calculatedANTIAREA)
              } else{
                print(paste(c(nrow(feature_pairsANTIAREA()),"feature(s) pairs to be calculated:",t(feature_pairsANTIAREA())), collapse = " " ))
                withProgress(message = 'Calculating adjacency', value = 0, {

                  feature_adjANTIAREA(AdjScore_calculator(infotable2(), feature_type=input$datatypeAdjacencybyFeature,Feature_Pairs=feature_pairsANTIAREA()) )
                })
                feature_adj_oldANTIAREA<-feature_adj_oldANTIAREA[,c(1:5)]
                feature_adj_oldANTIAREA<-rbind(feature_adj_oldANTIAREA,feature_adjANTIAREA())
                write.csv(x=feature_adj_oldANTIAREA,file=FileLocationANTIAREAadj())
                
                print ("Some adjacency features of the inverse area were calculated and added to previous calculations")
                
              }
              
            
              
            } else {
              print ("Adjacency feature(s) are being calculated and saved for the first time") 
              withProgress(message = 'Calculating adjacency of the area', value = 0, {

                infotab2_temp<-infotable2()
                
                if (input$datatypeAdjacencybyFeature =="protein" ) {
                  infotab2_temp@codex_protein<- infotab2_temp@codex_protein[rownames(infotab2_temp@codex_protein) %in%Tab4piece2()$CODEXname,]
                  infotab2_temp@codex_clean<- infotab2_temp@codex_clean[rownames(infotab2_temp@codex_clean) %in% Tab4piece2()$CODEXname,]
                  infotab2_temp@codex_spatial<-infotab2_temp@codex_spatial[rownames(infotab2_temp@codex_spatial) %in%Tab4piece2()$CODEXname,]
                  
                } else {
                  infotab2_temp@codex_mRNA<- infotab2_temp@codex_mRNA[rownames(infotab2_temp@codex_mRNA) %in%Tab4piece2()$CODEXname,]
                  infotab2_temp@codex_spatial<-infotab2_temp@codex_spatial[rownames(infotab2_temp@codex_spatial) %in%Tab4piece2()$CODEXname,]
                }

                print("Calculating adjacency scores of the area selected")
                feature_adjAREA (AdjScore_calculator(infotab2_temp, feature_type=input$datatypeAdjacencybyFeature,Feature_Pairs=feature_pairsAREA()) )
              })
              invadj<-feature_adjAREA()[,c(2,1,3:5)]
              colnames(invadj)<-colnames(feature_adjAREA())
              feature_adj_oldAREA<-rbind(feature_adjAREA(),invadj)
              
              
              write.csv(x=feature_adj_oldAREA,file=FileLocationAREAadj())
              print (paste(dim(feature_adjAREA())[1],"adjacency Area scores amongst feature pairs were calculated and saved for the first time"))
              
              withProgress(message = 'Calculating adjacency of the inverse area', value = 0, {
                infotab2_temp<-infotable2()
                
                if (input$datatypeAdjacencybyFeature =="protein" ) {
                  infotab2_temp@codex_protein<- infotab2_temp@codex_protein[!(rownames(infotab2_temp@codex_protein) %in%Tab4piece2()$CODEXname),]
                  infotab2_temp@codex_clean<- infotab2_temp@codex_clean[!(rownames(infotab2_temp@codex_clean) %in% Tab4piece2()$CODEXname),]
                  infotab2_temp@codex_spatial<-infotab2_temp@codex_spatial[!(rownames(infotab2_temp@codex_spatial) %in%Tab4piece2()$CODEXname),]
                  
                } else {
                  infotab2_temp@codex_mRNA<- infotab2_temp@codex_mRNA[!(rownames(infotab2_temp@codex_mRNA) %in%Tab4piece2()$CODEXname),]
                  infotab2_temp@codex_spatial<-infotab2_temp@codex_spatial[!(rownames(infotab2_temp@codex_spatial) %in%Tab4piece2()$CODEXname),]
                }
                
                print("Calculating adjacency scores of the inverse of the area selected")
                feature_adjANTIAREA (AdjScore_calculator(infotab2_temp, feature_type=input$datatypeAdjacencybyFeature,Feature_Pairs=feature_pairsANTIAREA()) )
              })
              
              invadj<-feature_adjANTIAREA()[,c(2,1,3:5)]
              colnames(invadj)<-colnames(feature_adjANTIAREA())
              feature_adj_oldINVAREA<-rbind(feature_adjANTIAREA(),invadj)
              
              write.csv(x=feature_adj_oldINVAREA,file=FileLocationANTIAREAadj())
              print (paste(dim(feature_adj())[1],"adjacency Area scores amongst feature pairs were calculated and saved for the first time"))
              
            }
          } else {}
        
        #### The plots start here
          if (input$RetainSubset9 == "Area") {
            
            output$p15<- renderPlot({
              if (dim(feature_adj_old_calculated)[1]==0 || is.null(feature_adj_old_calculated)) {
                feature_adj_to_plot<-feature_adj()
              } else {
                feature_adj_to_plot<-unique(rbind(feature_adj_old_calculated,feature_adj()) )
              }
              
              
              AdjScoreHeatmap_v2(adj_score_output= feature_adj_to_plot,low_color=input$color1AdjacencybyFeature,high_color =input$color2AdjacencybyFeature, title = "Adjacency plot of all Data")
              
              
            }) # Closing the plot p15 (adjacency)
            
            
            output$p16<- renderPlot({
              if (dim(feature_adj_old_calculated)[1]==0 || is.null(feature_adj_old_calculated)) {
                feature_adj_to_plot<-feature_adj()
              } else {
                feature_adj_to_plot<-unique(rbind(feature_adj_old_calculated,feature_adj()) )
              }
              
              AdjScoreHeatmap_cute(adj_score_output=feature_adj_to_plot, title= "Adjacency of all Data tuned by Pearson's Correlation")
              
            }) # Closing the plot p16 (correlation)
            
            
            output$AdjacencyTable<- DT::renderDT({
              if (dim(feature_adj_old_calculated)[1]==0 || is.null(feature_adj_old_calculated)) {
                feature_adj_to_plot<-feature_adj()
              } else {
                feature_adj_to_plot<-unique(rbind(feature_adj_old_calculated,feature_adj()) )
              }
              
              feature_adj_to_plot
            })# Closing the datatable of proportions of classification categories
            
            output$p17<- renderPlot({
              if (dim(feature_adj_old_calculatedAREA)[1]==0 || is.null(feature_adj_old_calculatedAREA)) {
                feature_adj_to_plot<-feature_adjAREA()
              } else {
                feature_adj_to_plot<-unique(rbind(feature_adj_old_calculatedAREA,feature_adjAREA()) )
              }
              
              AdjScoreHeatmap_v2(adj_score_output= feature_adj_to_plot,low_color=input$color1AdjacencybyFeature,high_color =input$color2AdjacencybyFeature, title = "Adjacency plot of Data inside Area")
              
            }) # Closing the plot p15 (adjacency)
            
            
            output$p18<- renderPlot({

              if (dim(feature_adj_old_calculatedAREA)[1]==0 || is.null(feature_adj_old_calculatedAREA)) {
                feature_adj_to_plot<-feature_adjAREA()
              } else {
                feature_adj_to_plot<-unique(rbind(feature_adj_old_calculatedAREA,feature_adjAREA()) )
              }
              
              AdjScoreHeatmap_cute(adj_score_output=feature_adj_to_plot, title= "Adjacency of Data inside Area tuned by Pearson's Correlation")
              
            }) # Closing the plot p16 (correlation)
            
            
            output$AdjacencyTableArea<- DT::renderDT({
              if (dim(feature_adj_old_calculatedAREA)[1]==0 || is.null(feature_adj_old_calculatedAREA)) {
                feature_adj_to_plot<-feature_adjAREA()
              } else {
                feature_adj_to_plot<-unique(rbind(feature_adj_old_calculatedAREA,feature_adjAREA()) )
              }
              
              colnames(feature_adj_to_plot)[3]<-"score_inside_area"
              feature_adj_to_plot
            })# Closing the datatable of proportions of classification categories
            
            output$p23<- renderPlot({
              NULL          
            }) # Closing the plot p23 
            
            output$p24<- renderPlot({
              NULL
            }) # Closing the plot p24
            
            output$AdjacencyTableANTIArea<- DT::renderDT({
              NULL
            })# Closing the datatable of proportions of classification categories
            
            
            } else if (input$RetainSubset9 == "Inverse"){
              
              output$p15<- renderPlot({

                if (dim(feature_adj_old_calculated)[1]==0 || is.null(feature_adj_old_calculated)) {
                  feature_adj_to_plot<-feature_adj()
                } else {
                  feature_adj_to_plot<-unique(rbind(feature_adj_old_calculated,feature_adj()) )
                }
                
                AdjScoreHeatmap_v2(adj_score_output= feature_adj_to_plot,low_color=input$color1AdjacencybyFeature,high_color =input$color2AdjacencybyFeature, title = "Adjacency plot of all Data")
                
              }) # Closing the plot p15 (adjacency)
              
              
              output$p16<- renderPlot({
                if (dim(feature_adj_old_calculated)[1]==0 || is.null(feature_adj_old_calculated)) {
                  feature_adj_to_plot<-feature_adj()
                } else {
                  feature_adj_to_plot<-unique(rbind(feature_adj_old_calculated,feature_adj()) )
                }
                
                AdjScoreHeatmap_cute(adj_score_output=feature_adj_to_plot, title= "Adjacency of all Data tuned by Pearson's Correlation")
                
              }) # Closing the plot p16 (correlation)
              
              
              output$AdjacencyTable<- DT::renderDT({
                if (dim(feature_adj_old_calculated)[1]==0 || is.null(feature_adj_old_calculated)) {
                  feature_adj_to_plot<-feature_adj()
                } else {
                  feature_adj_to_plot<-unique(rbind(feature_adj_old_calculated,feature_adj()) )
                }
                
                feature_adj_to_plot
              })# Closing the datatable of proportions of classification categories
              
              output$p17<- renderPlot({
                if (dim(feature_adj_old_calculatedANTIAREA)[1]==0 || is.null(feature_adj_old_calculatedANTIAREA)) {
                  feature_adj_to_plot<-feature_adjANTIAREA()
                } else {
                  feature_adj_to_plot<-unique(rbind(feature_adj_old_calculatedANTIAREA,feature_adjANTIAREA()) )
                }
                
                AdjScoreHeatmap_v2(adj_score_output= feature_adj_to_plot,low_color=input$color1AdjacencybyFeature,high_color =input$color2AdjacencybyFeature, title = "Adjacency plot of Data Outside Area")
                
              }) # Closing the plot p15 (adjacency)
              
              
              output$p18<- renderPlot({
                if (dim(feature_adj_old_calculatedANTIAREA)[1]==0 || is.null(feature_adj_old_calculatedANTIAREA)) {
                  feature_adj_to_plot<-feature_adjANTIAREA()
                } else {
                  feature_adj_to_plot<-unique(rbind(feature_adj_old_calculatedANTIAREA,feature_adjANTIAREA()) )
                }
                
                AdjScoreHeatmap_cute(adj_score_output=feature_adj_to_plot, title= "Adjacency of Data ouside Area tuned by Pearson's Correlation")
                
              }) # Closing the plot p16 (correlation)
              
              output$AdjacencyTableArea<- DT::renderDT({
                if (dim(feature_adj_old_calculatedANTIAREA)[1]==0 || is.null(feature_adj_old_calculatedANTIAREA)) {
                  feature_adj_to_plot<-feature_adjANTIAREA()
                } else {
                  feature_adj_to_plot<-unique(rbind(feature_adj_old_calculatedANTIAREA,feature_adjANTIAREA()) )
                }
                
                colnames(feature_adj_to_plot)[3]<-"score_outside_area"
                feature_adj_to_plot
              })# Closing the datatable of proportions of classification categories
              
              output$p23<- renderPlot({
                NULL          
              }) # Closing the plot p23 
              
              output$p24<- renderPlot({
                NULL
              }) # Closing the plot p24
              
              output$AdjacencyTableANTIArea<- DT::renderDT({
                NULL
              })# Closing the datatable of proportions of classification categories
              
              
            } else if (input$RetainSubset9 == "Both"){
        #if (input$RetainSubset9 == "Area" || input$RetainSubset9 == "Inverse"){
        
        output$p15<- renderPlot({
          if (dim(feature_adj_old_calculated)[1]==0 || is.null(feature_adj_old_calculated)) {
            feature_adj_to_plot<-feature_adj()
          } else {
            feature_adj_to_plot<-unique(rbind(feature_adj_old_calculated,feature_adj()) )
          }
          
          AdjScoreHeatmap_v2(adj_score_output= feature_adj_to_plot,low_color=input$color1AdjacencybyFeature,high_color =input$color2AdjacencybyFeature, title = "Adjacency plot of all Data")
          
        }) # Closing the plot p15 (adjacency)
        
        
        output$p16<- renderPlot({
          if (dim(feature_adj_old_calculated)[1]==0 || is.null(feature_adj_old_calculated)) {
            feature_adj_to_plot<-feature_adj()
          } else {
            feature_adj_to_plot<-unique(rbind(feature_adj_old_calculated,feature_adj()) )
          }
          
          AdjScoreHeatmap_cute(adj_score_output=feature_adj_to_plot, title= "Adjacency of all Data tuned by Pearson's Correlation")
          
        }) # Closing the plot p16 (correlation)
        
        
        output$AdjacencyTable<- DT::renderDT({
          if (dim(feature_adj_old_calculated)[1]==0 || is.null(feature_adj_old_calculated)) {
            feature_adj_to_plot<-feature_adj()
          } else {
            feature_adj_to_plot<-unique(rbind(feature_adj_old_calculated,feature_adj()) )
          }
          
          feature_adj_to_plot
        })# Closing the datatable of proportions of classification categories
        
        output$p17<- renderPlot({
          if (dim(feature_adj_old_calculatedAREA)[1]==0 || is.null(feature_adj_old_calculatedAREA)) {
            feature_adj_to_plot<-feature_adjAREA()
          } else {
            feature_adj_to_plot<-unique(rbind(feature_adj_old_calculatedAREA,feature_adjAREA()) )
          }
          
          AdjScoreHeatmap_v2(adj_score_output= feature_adj_to_plot,low_color=input$color1AdjacencybyFeature,high_color =input$color2AdjacencybyFeature, title = "Adjacency plot of Data inside Area")
          
          
        }) # Closing the plot p17 (adjacency)
        
        
        output$p18<- renderPlot({
          if (dim(feature_adj_old_calculatedAREA)[1]==0 || is.null(feature_adj_old_calculatedAREA)) {
            feature_adj_to_plot<-feature_adjAREA()
          } else {
            feature_adj_to_plot<-unique(rbind(feature_adj_old_calculatedAREA,feature_adjAREA()) )
          }
          
          AdjScoreHeatmap_cute(adj_score_output=feature_adj_to_plot, title= "Adjacency of Data inside Area tuned by Pearson's Correlation")
          
        }) # Closing the plot p18 (correlation)
        
        
        output$AdjacencyTableArea<- DT::renderDT({
          if (dim(feature_adj_old_calculatedAREA)[1]==0 || is.null(feature_adj_old_calculatedAREA)) {
            feature_adj_to_plot<-feature_adjAREA()
          } else {
            feature_adj_to_plot<-unique(rbind(feature_adj_old_calculatedAREA,feature_adjAREA()) )
          }
          
          colnames(feature_adj_to_plot)[3]<-"score_inside_area"
          feature_adj_to_plot
        })# Closing the datatable of proportions of classification categories
        
        
        output$p23<- renderPlot({
          if (dim(feature_adj_old_calculatedANTIAREA)[1]==0 || is.null(feature_adj_old_calculatedANTIAREA)) {
            feature_adj_to_plot<-feature_adjANTIAREA()
          } else {
            feature_adj_to_plot<-unique(rbind(feature_adj_old_calculatedANTIAREA,feature_adjANTIAREA()) )
          }
          
          
          AdjScoreHeatmap_v2(adj_score_output= feature_adj_to_plot,low_color=input$color1AdjacencybyFeature,high_color =input$color2AdjacencybyFeature, title = "Adjacency plot of Data Outside Area")
          
        }) # Closing the plot p23 (adjacency)
        
        
        output$p24<- renderPlot({
          if (dim(feature_adj_old_calculatedANTIAREA)[1]==0 || is.null(feature_adj_old_calculatedANTIAREA)) {
            feature_adj_to_plot<-feature_adjANTIAREA()
          } else {
            feature_adj_to_plot<-unique(rbind(feature_adj_old_calculatedANTIAREA,feature_adjANTIAREA()) )
          }
          
          AdjScoreHeatmap_cute(adj_score_output=feature_adj_to_plot, title= "Adjacency of Data ouside Area tuned by Pearson's Correlation")
          
        }) # Closing the plot p16 (correlation)
        
        
        output$AdjacencyTableANTIArea<- DT::renderDT({
          if (dim(feature_adj_old_calculatedANTIAREA)[1]==0 || is.null(feature_adj_old_calculatedANTIAREA)) {
            feature_adj_to_plot<-feature_adjANTIAREA()
          } else {
            feature_adj_to_plot<-unique(rbind(feature_adj_old_calculatedANTIAREA,feature_adjANTIAREA()) )
          }
          
          colnames(feature_adj_to_plot)[3]<-"score_outside_area"
          feature_adj_to_plot
        })# Closing the datatable of proportions of classification categories
      } 
          else {
            
            output$p15<- renderPlot({
              if (dim(feature_adj_old_calculated)[1]==0 || is.null(feature_adj_old_calculated)) {
                feature_adj_to_plot<-feature_adj()
              } else {
                feature_adj_to_plot<-unique(rbind(feature_adj_old_calculated,feature_adj()) )
              }
              
              AdjScoreHeatmap_v2(adj_score_output= feature_adj_to_plot,low_color=input$color1AdjacencybyFeature,high_color =input$color2AdjacencybyFeature,title = "Adjacency plot of all Data")
              
              
            }) # Closing the plot p15 (adjacency)
            
            
            output$p16<- renderPlot({
              if (dim(feature_adj_old_calculated)[1]==0 || is.null(feature_adj_old_calculated)) {
                feature_adj_to_plot<-feature_adj()
              } else {
                feature_adj_to_plot<-unique(rbind(feature_adj_old_calculated,feature_adj()) )
              }
              
              AdjScoreHeatmap_cute(adj_score_output=feature_adj_to_plot, title= "Adjacency of all Data tuned by Pearson's Correlation")
              
            }) # Closing the plot p16 (correlation)
            
            
            output$AdjacencyTable<- DT::renderDT({
              if (dim(feature_adj_old_calculated)[1]==0 || is.null(feature_adj_old_calculated)) {
                feature_adj_to_plot<-feature_adj()
              } else {
                feature_adj_to_plot<-unique(rbind(feature_adj_old_calculated,feature_adj()) )
              }
              
              feature_adj_to_plot
              
            })# Closing the datatable of proportions of classification categories
            
            output$p17<- renderPlot({
              NULL          
            }) # Closing the plot p17 
            
            output$p18<- renderPlot({
              NULL
            }) # Closing the plot p18
            
            output$AdjacencyTableArea<- DT::renderDT({
              NULL
            })# Closing the datatable of proportions of classification categories
            output$p23<- renderPlot({
              NULL          
            }) # Closing the plot p23 
            
            output$p24<- renderPlot({
              NULL
            }) # Closing the plot p24
            
            output$AdjacencyTableANTIArea<- DT::renderDT({
              NULL
            })# Closing the datatable of proportions of classification categories
            
            
          }
        
        
      })# Closing the observeEvent of ploteaAdjacencybyFeature
    
    ######################
    
    ## AdjacencybyCluster
    #####################
    feature_pairs_original <- reactiveVal(NULL);
    feature_pairs <- reactiveVal(NULL);
    feature_adj <- reactiveVal(NULL);
    feature_adjbyClusterAREA <- reactiveVal(NULL);
    feature_adjbyClusterANTIAREA <- reactiveVal(NULL);
    
    adjacency_table <- reactiveVal(NULL);
    prenameforadjtab<-reactiveVal(NULL);
    nameadjtab<-reactiveVal(NULL);
    thepathadj<-reactiveVal(NULL);
    FileLocationadj<-reactiveVal(NULL);
    var_features<-reactiveVal(NULL);
    
    observeEvent({
      input$selectclusteringtype}, {
        req(input$selectclusteringtype)
        req(infotable4())
        
        var_features(infotable4()[,input$selectclusteringtype])
      })
    
    observeEvent({
      input$selectclusteringtype}, {
        req(input$selectclusteringtype)
        req(infotable2())
        req(infotable4())
        var_to_highlightadjacencybycluster<- sort(unique (infotable4()[,input$selectclusteringtype]))
        updateSelectInput(session, inputId = "clusterincolorbyCluster", choices = var_to_highlightadjacencybycluster )
        
      })
    
    observeEvent({
      input$ploteaAdjacencybyCluster}, {
        req(infotable2())
        req(infotable4())
        req(var_features())
        req(input$selectclusteringtype)
        req(input$clusterincolorbyCluster)

        # ##################
        # Here we will insert the new section in which we will calculate and
        # save the adjacency of several clusters
        
        prenameforadjtab(paste(substr(input$file1$name,1,nchar(input$file1$name)-4),substr(input$file2$name,1,nchar(input$file2$name)-6), sep= "_") )
        
        
        nameadjtab <- paste0("Adjacency_by_Cluster_type","_",input$selectclusteringtype,"_",prenameforadjtab(),".csv") 
        
        
        print("nameadjtab")
        print(nameadjtab)
        
        thepathadj(getwd())
        print("thepathadj")
        print(thepathadj())
        ifelse(!dir.exists(paste0(thepathadj(),"/",prenameforadjtab())), dir.create(paste0(thepathadj(),"/",prenameforadjtab()), showWarnings = TRUE, recursive = FALSE, mode = "0777") , FALSE)
        
        FileLocationadj (paste0(thepathadj(),"/",prenameforadjtab(),"/",nameadjtab) )
        print("FileLocationadj")
        print(FileLocationadj())
        
        if (file.exists(FileLocationadj())) {
          print(paste ("Adjacency of", input$selectclusteringtype, "was previously calculated"))
          
          feature_adj ( read.csv(FileLocationadj(), row.names=1))
          
          } else {
            print(paste ("Calculating Adjacency of", input$selectclusteringtype))
            
            withProgress(message = 'Calculating adjacency', value = 0, {
              feature_adj(AdjScoreClustersCODEX_vGator(infotable2(), infotable= infotable4(),clusters=var_features(), k=3, num_cores=1)  )
            })
            
            write.csv(x=feature_adj(),file=FileLocationadj())
          }
            
         if (input$RetainSubset10 == "Area" || input$RetainSubset10 == "Inverse" || input$RetainSubset10 == "Both") {
              prenameAREAforadjtab(paste(substr(input$fileAreaSubsetted$name,1,nchar (input$fileAreaSubsetted$name)-4),
                                         substr(input$file1$name,1,nchar(input$file1$name)-4),substr(input$file2$name,1,nchar(input$file2$name)-6), sep= "__") )
              
              sctyp<- gsub("\\.", "_", input$selectclusteringtype)
              nameAREAadjtab <- paste0("Adj_by_Cluster_Area_",sctyp,"_",prenameAREAforadjtab(),".csv")
              
              nameANTIAREAadjtab <- paste0("Adj_by_Cluster_InverseArea_",sctyp,"_",prenameAREAforadjtab(),".csv")
              
              print("nameAREAadjtab")
              print(nameAREAadjtab)
              print("nameINVERSEAREAadjtab")
              print(nameANTIAREAadjtab)
              
              thepathadj(getwd())
              print("thepathadj")
              print(thepathadj())
              
              FileLocationAREAadj (paste0(thepathadj(),"/",prenameforadjtab(),"/",nameAREAadjtab) )
              FileLocationANTIAREAadj (paste0(thepathadj(),"/",prenameforadjtab(),"/",nameANTIAREAadjtab) )
              
              print("FileLocationAREAadj")
              print(FileLocationAREAadj())
              print("FileLocationANTIAREAadj")
              print(FileLocationANTIAREAadj())
              
              if (file.exists(FileLocationAREAadj()) & file.exists(FileLocationANTIAREAadj()) ) {
                print(paste ("Adjacency by cluster of", input$selectclusteringtype, "was previously calculated"))
                
                feature_adjbyClusterAREA ( read.csv(FileLocationAREAadj(), row.names=1))
                feature_adjbyClusterANTIAREA ( read.csv(FileLocationANTIAREAadj(), row.names=1))
              } else {
                print(paste ("Calculating Adjacency by cluster of area and inverse area of ", input$selectclusteringtype))
                
                withProgress(message = 'Calculating Area adjacency', value = 0, {
                  
                  infotab2_temp<-infotable2()
                  infotab2_temp@codex_spatial<-infotab2_temp@codex_spatial[rownames(infotab2_temp@codex_spatial) %in%Tab4piece2()$CODEXname,]
                  infotab4_temp<-infotable4()[infotable4()$CODEXname %in%Tab4piece2()$CODEXname, ]
                  thevar_features<-infotab4_temp[,input$selectclusteringtype]
                  
                  feature_adjbyClusterAREA(AdjScoreClustersCODEX_vGator(infotab2_temp, infotable= infotab4_temp,clusters=thevar_features, k=3, num_cores=1))
                  
                })
                
                write.csv(x=feature_adjbyClusterAREA(),file=FileLocationAREAadj())
                
                withProgress(message = 'Calculating Inverse Area adjacency', value = 0, {
                  
                  infotab2_temp<-infotable2()
                  infotab2_temp@codex_spatial<-infotab2_temp@codex_spatial[!(rownames(infotab2_temp@codex_spatial) %in%Tab4piece2()$CODEXname),]
                  infotab4_temp<-infotable4()[!(infotable4()$CODEXname %in%Tab4piece2()$CODEXname), ]
                  thevar_features<-infotab4_temp[,input$selectclusteringtype]
                  
                  feature_adjbyClusterANTIAREA(AdjScoreClustersCODEX_vGator(infotab2_temp, infotable= infotab4_temp,clusters=thevar_features, k=3, num_cores=1))
                  
                })
                
                write.csv(x=feature_adjbyClusterANTIAREA(),file=FileLocationANTIAREAadj())
              } 
            } else {}
          
        #### The plots start here
        if (input$RetainSubset10 == "Area") {
        output$p21<- renderPlot({
          AdjScoreHeatmap_v2(adj_score_output= feature_adj(),low_color=input$color1AdjacencybyCluster,high_color =input$color2AdjacencybyCluster,title = "Adjacency plot of all Data",Subset_Matrix = input$clusterincolorbyCluster)
        }) # Closing the plot p21 (adjacency)
        
        output$p22<- renderPlot({       
          
          AdjScoreHeatmap_cute(adj_score_output=feature_adj(),title= "Adjacency of all Data tuned by Pearson's Correlation",Subset_Matrix = input$clusterincolorbyCluster )
          
        }) # Closing the plot p22 (correlation)
        
        output$AdjacencyTableCluster<- DT::renderDT({
          feature_adj()
  
        })# Closing the datatable of proportions of classification categories
        
        output$p25<- renderPlot({
          AdjScoreHeatmap_v2(adj_score_output= feature_adjbyClusterAREA(),low_color=input$color1AdjacencybyCluster,high_color =input$color2AdjacencybyCluster,title = "Adjacency plot of Data inside Area",Subset_Matrix = input$clusterincolorbyCluster)
        }) # Closing the plot p25 (adjacency)
        
        
        output$p26<- renderPlot({       
          
          AdjScoreHeatmap_cute(adj_score_output=feature_adjbyClusterAREA(),title= "Adjacency of Data inside Area tuned by Pearson's Correlation",Subset_Matrix = input$clusterincolorbyCluster )
          
        }) # Closing the plot p26 (correlation)
        
        output$AdjacencyTableClusterArea<- DT::renderDT({
          tab<-feature_adjbyClusterAREA()
          colnames(tab)[3]<-"score_inside_area"
          tab
          
        })# Closing the datatable of proportions of classification categories
        
        output$p27<- renderPlot({
          NULL          
        }) # Closing the plot p27 
        
        output$p28<- renderPlot({
          NULL
        }) # Closing the plot p28
        
        output$AdjacencyTableClusterANTIArea<- DT::renderDT({
          NULL
        })# Closing the datatable of proportions of classification categories
        
        
        
        } else if (input$RetainSubset10 == "Inverse") { # END of the if
          output$p21<- renderPlot({
            AdjScoreHeatmap_v2(adj_score_output= feature_adj(),low_color=input$color1AdjacencybyCluster,high_color =input$color2AdjacencybyCluster,title = "Adjacency plot of all Data",Subset_Matrix = input$clusterincolorbyCluster)
          }) # Closing the plot p21 (adjacency)
          
          output$p22<- renderPlot({       
            
            AdjScoreHeatmap_cute(adj_score_output=feature_adj(),title= "Adjacency of all Data tuned by Pearson's Correlation",Subset_Matrix = input$clusterincolorbyCluster )
            
          }) # Closing the plot p22 (correlation)
          
          output$AdjacencyTableCluster<- DT::renderDT({
            feature_adj()
            
          })# Closing the datatable of proportions of classification categories
          
          output$p25<- renderPlot({
            AdjScoreHeatmap_v2(adj_score_output= feature_adjbyClusterANTIAREA(),low_color=input$color1AdjacencybyCluster,high_color =input$color2AdjacencybyCluster,title = "Adjacency plot of Data outside Area",Subset_Matrix = input$clusterincolorbyCluster)
          }) # Closing the plot p25 (adjacency)
          
          
          output$p26<- renderPlot({       
            
            AdjScoreHeatmap_cute(adj_score_output=feature_adjbyClusterANTIAREA(),title= "Adjacency of Data outside Area tuned by Pearson's Correlation",Subset_Matrix = input$clusterincolorbyCluster )
            
          }) # Closing the plot p26 (correlation)
          
          output$AdjacencyTableClusterArea<- DT::renderDT({
            tab<-feature_adjbyClusterANTIAREA()
            colnames(tab)[3]<-"score_outside_area"
            tab
            
          })# Closing the datatable of proportions of classification categories
          
          output$p27<- renderPlot({
            NULL          
          }) # Closing the plot p27 
          
          output$p28<- renderPlot({
            NULL
          }) # Closing the plot p28
          
          output$AdjacencyTableClusterANTIArea<- DT::renderDT({
            NULL
          })# Closing the datatable of proportions of classification categories 
        
        } else if (input$RetainSubset10 == "Both") {
          output$p21<- renderPlot({
            AdjScoreHeatmap_v2(adj_score_output= feature_adj(),low_color=input$color1AdjacencybyCluster,high_color =input$color2AdjacencybyCluster,title = "Adjacency plot of all Data",Subset_Matrix = input$clusterincolorbyCluster)
          }) # Closing the plot p21 (adjacency)
          
          output$p22<- renderPlot({       
            
            AdjScoreHeatmap_cute(adj_score_output=feature_adj(),title= "Adjacency of all Data tuned by Pearson's Correlation",Subset_Matrix = input$clusterincolorbyCluster )
            
          }) # Closing the plot p22 (correlation)
          
          output$AdjacencyTableCluster<- DT::renderDT({
            feature_adj()
            
          })# Closing the datatable of proportions of classification categories
          
          output$p25<- renderPlot({
            AdjScoreHeatmap_v2(adj_score_output= feature_adjbyClusterAREA(),low_color=input$color1AdjacencybyCluster,high_color =input$color2AdjacencybyCluster,title = "Adjacency plot of Data inside Area",Subset_Matrix = input$clusterincolorbyCluster)
          }) # Closing the plot p25 (adjacency)
          
          
          output$p26<- renderPlot({       
            
            AdjScoreHeatmap_cute(adj_score_output=feature_adjbyClusterAREA(),title= "Adjacency of Data inside Area tuned by Pearson's Correlation",Subset_Matrix = input$clusterincolorbyCluster )
            
          }) # Closing the plot p26 (correlation)
          
          output$AdjacencyTableClusterArea<- DT::renderDT({
            tab<-feature_adjbyClusterAREA()
            colnames(tab)[3]<-"score_inside_area"
            tab
            
          })# Closing the datatable of proportions of classification categories
          
          output$p27<- renderPlot({
            AdjScoreHeatmap_v2(adj_score_output= feature_adjbyClusterANTIAREA(),low_color=input$color1AdjacencybyCluster,high_color =input$color2AdjacencybyCluster,title = "Adjacency plot of Data outside Area",Subset_Matrix = input$clusterincolorbyCluster)
          }) # Closing the plot p27 (adjacency)
          
          
          output$p28<- renderPlot({       
            
            AdjScoreHeatmap_cute(adj_score_output=feature_adjbyClusterANTIAREA(),title= "Adjacency of Data outside Area tuned by Pearson's Correlation",Subset_Matrix = input$clusterincolorbyCluster )
            
          }) # Closing the plot p28 (correlation)
          
          output$AdjacencyTableClusterANTIArea<- DT::renderDT({
            tab<-feature_adjbyClusterANTIAREA()
            colnames(tab)[3]<-"score_outside_area"
            tab
            
          })# Closing the datatable of proportions of classification categories
          
        } else {
          output$p21<- renderPlot({
            AdjScoreHeatmap_v2(adj_score_output= feature_adj(),low_color=input$color1AdjacencybyCluster,high_color =input$color2AdjacencybyCluster,title = "Adjacency plot of all Data",Subset_Matrix = input$clusterincolorbyCluster)
          }) # Closing the plot p21 (adjacency)
          
          output$p22<- renderPlot({       
            
            AdjScoreHeatmap_cute(adj_score_output=feature_adj(),title= "Adjacency of all Data tuned by Pearson's Correlation",Subset_Matrix = input$clusterincolorbyCluster )
            
          }) # Closing the plot p22 (correlation)
          
          output$AdjacencyTableCluster<- DT::renderDT({
            feature_adj()
            
          })# Closing the datatable of proportions of classification categories
          
          output$p25<- renderPlot({
            NULL
          }) # Closing the plot p25 (adjacency)
          
          
          output$p26<- renderPlot({       
            NULL
          }) # Closing the plot p26 (correlation)
          
          output$AdjacencyTableClusterArea<- DT::renderDT({
            NULL
          })# Closing the datatable of proportions of classification categories
          
          output$p27<- renderPlot({
            NULL
          }) # Closing the plot p27 (adjacency)
          
          
          output$p28<- renderPlot({       
            NULL
          }) # Closing the plot p28 (correlation)
          
          output$AdjacencyTableClusterANTIArea<- DT::renderDT({
            NULL
          })# Closing the datatable of proportions of classification categories
        }
        
        
        
        
      })# Closing the observeEvent of ploteaAdjacencybyFeature
    
    ######################
    
    
    
  } # Closing the Function
)  # Closing shinyServer

# END
