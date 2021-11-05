#############################################
#######     InterroGATOR_server.R     #######
#############################################
# Author: Leandro Balzano-Nogueira
# Diabetes Institute, University of Florida (Gainesville)
# Last update: November/5/2021

# Creating and hiding calculations for posterior analyses

# This is a Shiny application for conducting differential expression and Spatial analyses on different cell types in different human tissues. This application allows the user to analyze cell expression profiles from CITEseq and CODEX data previously integrated.


#############################################

# Libraries:
library(shiny)
library(shinydashboard)
library(STvEA)
library (reticulate)
library(ggplot2)
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
    
    observeEvent({
      input$file1
    }, {
      if (is.null(input$file1)) {
        return(NULL)
      }
      
      infotable1(readRDS (input$file1$datapath) )
      #infotable1(readRDS ("/Users/leobalzano/Dropbox (UFL)/Diabetes2021/STvEA/Results/CITE_CODEX_Integration/Seurat47_54_73_79Merged.rds") )
      infotable3(infotable1()@meta.data[,-which(sapply(infotable1()@meta.data, class) == "numeric")])
      coloring_choices <- colnames(infotable3() )
      print(coloring_choices)
      updateSelectInput(session, inputId = "colby", choices = coloring_choices)
      updateSelectInput(session, inputId = "cluster1", choices = coloring_choices ) # This is brother of this in UI. If this is on, the piece in UI has to be on as well
      updateSelectInput(session, inputId = "cluster2", choices = coloring_choices ) # This is brother of this in UI. If this is on, the piece in UI has to be on as well
      updateSelectInput(session, inputId = "colbyDE", choices = coloring_choices )
    })
    
    observeEvent({
      input$file2
    }, {
      if (is.null(input$file2)) {
        return(NULL)
      }
        
      
      elnewenv <- new.env()
      elname<- load(input$file2$datapath, envir = elnewenv)
      #elname<- load("/Users/leobalzano/Dropbox (UFL)/Diabetes2021/STvEA/Results/CITE_CODEX_Integration/HDL_54_SCALEDSeurat47547379MergedStrategy2Sep212021_cp6.RData", envir = elnewenv)
      #elname<- load ("/Users/leobalzano/Dropbox (UFL)/Diabetes2021/STvEA/Results/CITE_CODEX_Integration/steveamock.RData", envir = elnewenv)
      infotable2(elnewenv[[elname]] )
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
    
    output$p1<- renderPlot({
      req(infotable2())
      req(input$colby)
      print("input$colby")
      print(input$colby)
      thecolor<- infotable3()[,input$colby];
      tab<-melt(table(thecolor))
      colnames(tab)<- c("Category", "Cell_number")
      tabproportions (tab)
      PlotClusterCITE_v2(infotable2() ,color_by = thecolor,pt_size =input$size1)
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
      infotableumapind(cbind(assay_name=infotable1()@meta.data$assay_name,clusters=infotable1()@meta.data[,input$clusterUMAPIndSamples],infotable2()@cite_emb))
      
      if (any(infotableumapind()$clusters==0)) {
        infotableumapind()$clusters<-as.factor(ifelse(as.numeric(as.character(infotableumapind()$clusters))<0,-1, as.numeric(as.character(infotableumapind()$clusters)) +1) )
      }
      groupstoplot <- sort(unique (infotableumapind()$assay_name))
      print("groupstoplot")
      print(groupstoplot)
      updateSelectInput(session, inputId = "Umapstoplot", choices = groupstoplot)
      
    })
    
    observeEvent({
      input$ploteaumapindsamples}, {
        req( input$Umapstoplot)
        
        colors <- colorspace::rainbow_hcl(length(unique(infotableumapind()$clusters)),c = 80)
        Subsetinfotableumapind (infotableumapind()[infotableumapind()$assay_name == input$Umapstoplot,] )
        print("dim (Subsetinfotableumapind)")
        print (dim(Subsetinfotableumapind()))
        seurat_clusters<-Subsetinfotableumapind()$clusters
        
        output$p5<- renderPlot({
          ggplot(Subsetinfotableumapind(), aes_string(x = colnames(Subsetinfotableumapind() )[3],
                                                      y = colnames(Subsetinfotableumapind() )[4],
                                                      color = factor(seurat_clusters))) +
            
            geom_point(aes(fill=factor(seurat_clusters)),colour="gray5",size = input$sizeumapindsamples,pch=21) +
            scale_color_manual(values = colors, name = "cluster") +
            guides(fill=guide_legend(title="cluster")) +
            guides(colour = guide_legend(override.aes = list(size = 5))) +
            theme_void()
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
          thecolor<- infotable3()[,input$colby];
          PlotClusterCITE_Independently(infotable2() ,color_by = thecolor,pt_size =input$size,highlight =input$highlight, Color = input$color, Selection_on_Top=input$SelectiononTop)
        }) # Closing the plot p2
      })
    
    ######################

    ## CITEtoCODEXspatial extrapolation
    #####################
    infotable4  <- reactiveVal(NULL);
    info_spatial_tmp <- reactiveVal(NULL);
    
    observeEvent({
      infotable2()} , {
        req(infotable2())
        req(infotable3())
        #######   Piece to save table4
        if ( isFALSE(file_test("-f", FileLocation())) ) {
          print ("Calculating infotable4")
          v1<-rownames(infotable2()@transfer_matrix)
          v2<-rowSums(infotable2()@transfer_matrix)
          v3<-colnames(infotable2()@transfer_matrix)[apply(infotable2()@transfer_matrix,1,which.max)] # Retain For each row return the column name of the largest value (It takes time!)
          tab<-cbind(CODEXname=v1,v2=v2,CITEname=v3)
          print (head (tab))
          rm(v1,v2,v3)
          sort(table(tab[,3]))
          tab2<-tab[tab[,2]!=0,]
          cellstab<-infotable3()
          tab3<-merge(tab2,cellstab, by.x=3, by.y=0)
          tab4<-merge(tab3,infotable2()@codex_spatial, by.x = 2,by.y=0, all.y = TRUE)
          tab4$predicted.celltype.l2[is.na(tab4$predicted.celltype.l2)] = "Unknown"
          tab4$Higher_Hierarchy_grouping[is.na(tab4$Higher_Hierarchy_grouping)] = "Unknown"
          
          tab4$condition <- with(tab4, ifelse(tab4$predicted.celltype.l2 =="Unknown", 0, 1))
          tab4 <- tab4[order(tab4$condition),]
          
          infotable4(tab4)
          write.csv(x=tab4,file=FileLocation())
          rm(tab,tab2,tab3,tab4, cellstab)
        } else {
          infotable4(read.csv(FileLocation() ) )
          print (infotable4())
          print ("infotable4 was retreived from previous calculations")
        }
        
        #######
        x_tmp <- infotable4()[, "x"]
        x_tmp <- x_tmp - min(x_tmp)
        y_tmp <- infotable4()[, "y"]
        y_tmp <- y_tmp - min(y_tmp)
        spatial_tmp <- as.data.frame(cbind(x = x_tmp, y = y_tmp))
        info_spatial_tmp(spatial_tmp)
        rm (spatial_tmp)
      })
    
    
    observeEvent({
      input$cluster1}, {
        req( input$cluster1)
        
      })  
    
    output$p3<- renderPlot({
      colors <- rainbow_hcl(length(unique(infotable4()[,input$cluster1])),
                            c = 80)
      
      colors<-c(colors[1:length(colors)-1],"gray85")
      fills<-c(colors[1:length(colors)-1],"gray95")
      thecolor2<- infotable4()[,input$cluster1];
      print("thecolor2")
      print(thecolor2)
      
      ggplot(info_spatial_tmp()) +
        geom_point(aes(x = x,y = y, color = factor(thecolor2)),size = input$size2) +
        scale_color_manual(values = colors, name = "cluster") +
        guides(colour = guide_legend(override.aes = list(size = 5))) +
        theme_void()
      
    }) # Closing the plot
    ######################
    
    ## CITEtoCODEXspatial independent
    #####################
    info_spatial_tmpindep <- reactiveVal(NULL);

    observeEvent({
      input$cluster2}, {
        req( input$cluster2)
        var_to_highlight2<- sort(unique (infotable4()[,input$cluster2]))
        updateSelectInput(session, inputId = "clusterincolor", choices = var_to_highlight2 )
      })

    observeEvent({
      input$plotea}, {
        req( input$clusterincolor)
        x_tmp <- infotable4()[, "x"]
        x_tmp <- x_tmp - min(x_tmp)
        y_tmp <- infotable4()[, "y"]
        y_tmp <- y_tmp - min(y_tmp)
        spatial_tmp <- as.data.frame(cbind(x = x_tmp, y = y_tmp))
        info_spatial_tmpindep(spatial_tmp)
        rm (spatial_tmp)
        print(input$clusterincolor)
        df <- infotable4()
        df$condition2 <- ifelse (infotable4()[,input$cluster2] %in% input$clusterincolor, 1,0)
        df<-df[order(df$condition2),]
        infotable4(df)
        rm(df)

        df<-info_spatial_tmpindep()
        df$condition2<-infotable4()$condition2
        df<-df[order(df$condition2),]
        info_spatial_tmpindep(df)
        rm(df)
        
        colorsindep <- rainbow_hcl(length(unique(infotable4()[,input$cluster2])),
                                   c = 80)
        colorsindep<-c(colorsindep[1:length(colorsindep)-1],"gray85")
        fillsindep<-c(colorsindep[1:length(colorsindep)-1],"gray95")

        groupsindep<- sort(unique (infotable4()[,input$cluster2]))

        colorsindep2<-data.frame(cbind(groupsindep,colorsindep))
        colorsindep2$theones<-ifelse(colorsindep2$groups %in% input$clusterincolor, colorsindep2$colorsindep,"gray85")

        output$p4<- renderPlot({
          
          ggplot(info_spatial_tmpindep()) +
            geom_point(aes(x = x,y = y, color = factor(infotable4()[,input$cluster2])),size = input$sizespatial2) +
            scale_color_manual(values = colorsindep2$theones, name = "cluster") +
            guides(colour = guide_legend(override.aes = list(size = 5))) +
            theme_void()
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
        var_to_highlightDE<- sort(unique (infotable3()[,input$colbyDE]))
        updateSelectizeInput(inputId="selectDE", choices = var_to_highlightDE, options = list(maxItems = 2))
        
      })
    
    observeEvent({
      input$calculaDE}, {
        req(infotable1())
        req(infotable2())
        req(input$colbyDE)
        req(input$selectDE)
        req(input$datatypeDE)
      
        predat1<-infotable3()[,input$colbyDE]
        if (input$datatypeDE =="cite_mRNA_norm" ){
          predat2<-infotable2()@cite_mRNA_norm  
          print("dim(predat2)")
          print(dim(predat2))
        } else {
          predat2<-infotable2()@cite_protein
        }
        
        dataforDE(data.frame(cbind(predat1,predat2)) )
        print(input$selectDE)
        
        
        
        output$tDE<- DT::renderDT({
          datadeg (DEgenesPairwiseComparison (X=dataforDE(), groups=c(input$selectDE)))
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
        
        output$p6 <- renderPlot({
          
          if (input$datatypeDE =="cite_mRNA_norm" ){
            UMAPFeatureExpression(infotable2(),features=input$featureUMAPDE, type = "RNA", low_color=input$colorUMAPDE2, high_color = input$colorUMAPDE, pt_size= input$sizeDE)
          } else {
            print(input$colorUMAPDE)
            print(input$colorUMAPDE2)
            print(input$colorUMAPDE3)
            UMAPFeatureExpression(infotable2(),features=input$featureUMAPDE, type = "protein", high_color = input$colorUMAPDE, low_color=input$colorUMAPDE2,high_color2=input$colorUMAPDE3, pt_size= input$sizeDE)
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
          var_to_highlightDE<- sort(unique (infotable3()[,input$colbyDE]))
          updateSelectizeInput(inputId="selectDE", choices = var_to_highlightDE, options = list(maxItems = 2))
          print("var_to_highlightDE")
          print (var_to_highlightDE)
          updateSelectizeInput(session, inputId="selectgroupsforheatmap", choices = var_to_highlightDE, selected = var_to_highlightDE )
          updateSelectizeInput(session, inputId="selectgroupsforridgeandviolin", choices = var_to_highlightDE )
        } else {
          var_to_highlightDE<- sort(unique (infotable3()[,input$colbyDE]))
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
        req(input$selectgroupsorder)
        req(input$heatmapcolor)
        req(input$cellnumber)
        req(input$lowerbreak)
        req(input$upperbreak)
        
        output$p7<- renderPlot({
          
          HeatmapbyGroup(Data = dataforDE(),Genes=input$selectfeaturesforheatmap,Groups = input$selectgroupsforheatmap,Heatmap_Color = input$heatmapcolor,Ncells = input$cellnumber,Breaks = c(input$lowerbreak,input$upperbreak),group_order = c(input$selectgroupsorder) )
        })# Closing the p7
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
            
    #browser()
  } # Closing the Function
)  # Closing shinyServer

# END

