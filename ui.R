#############################################
########     InterroGATOR_ui.R     ##########
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
library(DT)
library (RColorBrewer)
library (grDevices)
library(randomcoloR)
library(AdjacencyScore)
library(sparseMatrixStats)


#############################################
## ui.R ##
options(shiny.maxRequestSize=120000*1024^2)

upload_data <- tags$a(span("InterroGATOR", style = "color: white; font-size: 30px;font-weight: bold"),tags$img(src="Gator.jpg", height="40", width= "50"))

shinyUI(
  dashboardPage(skin = "blue", #skin color
   
    ### HEADER
    dashboardHeader(title=upload_data), ## Header

    
    ### Side Bar
    dashboardSidebar(
      sidebarMenu(id = "sidebarmenu", ### create a menu in the sidebar
          # CITE UMAP Analysis
          menuItem("File Upload", tabName = "fileupload", icon = icon("code-branch")) ,
          conditionalPanel("input.sidebarmenu === 'fileupload'",
          fluidPage(
              # File Input 1

              fileInput("file1", label = h4("Insert Seurat rds file")),
              
                                   
              fileInput("file2", label = h4("Insert RData object")),
              
              selectInput(inputId = "colby", label = "Classification category",choices = list() ), 
              sliderInput(inputId = "size1",label = "Size",min = 0,max = 5,value = 2.4,step = 0.1),
              
              radioButtons(
                inputId="AreaAnalysis",
                label="Do you have an area of interest?",
                choices=list(
                  "Yes",
                  "No"
                ),
                selected="No"),
              
              conditionalPanel(
                condition = "input.AreaAnalysis != 'No'",
                  fileInput("fileAreaSubsetted", label = h4("Insert cell coordinates CSV file")),
                  radioButtons(inputId="RetainSubset", label="Subset?", choices= list("Area","Inverse", "All"), selected = "All")
              ),
              
              ) # Closing fluidPage of input  
            ), # Closing conditionalPanel of Input
          
          # CITE UMAP of samples plotted individually
          menuItem("CITEseq UMAP individual samples", tabName = "UMAPIndSamples", icon = icon("code-branch")),
          conditionalPanel("input.sidebarmenu === 'UMAPIndSamples'",
                           fluidPage(
                             # File Input 1
                             selectInput(inputId = "clusterUMAPIndSamples", label = "Classification category",choices = list()),
                             selectInput(inputId = "Umapstoplot", label = "Sample to plot",choices = list() ),
                             conditionalPanel(
                               condition = "input.AreaAnalysis != 'No'",
                               radioButtons(inputId="RetainSubset2", label="Subset?", choices= list("Area","Inverse", "All"), selected = "All")
                             ),
                             
                             sliderInput(inputId = "sizeumapindsamples",label = "Size",min = 0,max = 5,value = 2.4,step = 0.1),
                             actionButton("ploteaumapindsamples", "Plot")
                           ) # Closing fluidPage of CITE UMAP of samples plotted individually
          ), # Closing conditionalPanel of CITE UMAP of samples plotted individually
          
          # CITE independent Tab
          menuItem("CITE independent", tabName = "CITEindep", icon = icon("code-branch")),
          conditionalPanel("input.sidebarmenu === 'CITEindep'",
          fluidPage(
            # File Input 1
            selectInput(inputId = "highlight", label = "Category to highlight",choices = list() ),
            selectInput(inputId = "color", label = "Color",choices = list("darkred", "red", "blue","darkblue", "black"),selected = "darkred" ),
            sliderInput(inputId = "size",label = "Size",min = 0,max = 5,value = 2.4,step = 0.1),
            radioButtons("SelectiononTop", "Selection on Top:",
                         c("Yes" = "yes",
                           "No" = "no")),
            conditionalPanel(
              condition = "input.AreaAnalysis != 'No'",
              radioButtons(inputId="RetainSubset3", label="Subset?", choices= list("Area","Inverse", "All"), selected = "All")
            ),
            actionButton("ploteaCITEIndep", "Plot")
                           ) # Closing fluidPage of CITE independent Tab  
          ), # Closing conditionalPanel of CITE independent Tab
          
          # CITEtoCODEXspatial extrapolation tab
          menuItem("CITEtoCODEXspatial extrapolation", tabName = "spatial1", icon = icon("code-branch")),
          conditionalPanel("input.sidebarmenu === 'spatial1'",
                           fluidPage(
                             textInput(inputId="IdentificationThreshold", label="Probability of Identification (only for RNA)", value= 0),
                             # File Input 1
                             selectInput(inputId = "cluster1", label = "Classification category",choices = list() ),
                             sliderInput(inputId = "size2",label = "Size",min = 0,max = 5,value = 1.3,step = 0.1),
                             conditionalPanel(
                               condition = "input.AreaAnalysis != 'No'",
                               radioButtons(inputId="RetainSubset4", label="Subset?", choices= list("Area","Inverse", "All"), selected = "All")
                             ),
                             actionButton("ploteaExtrapolationTab", "Plot")
                           ) # Closing fluidPage of CITEtoCODEXspatial extrapolation tab  
          ), # Closing conditionalPanel of CITEtoCODEXspatial extrapolation tab
          
          # CITEtoCODEXspatial independent tab
          menuItem("CITEtoCODEXspatial independent", tabName = "spatial2", icon = icon("code-branch")),
          conditionalPanel("input.sidebarmenu === 'spatial2'",
                           fluidPage(
                             
                             selectInput(inputId = "cluster2", label = "Classification category",choices = list("predicted.celltype.l2","Higher_Hierarchy_grouping") ),
                             selectInput(inputId = "clusterincolor", label = "Colored category(ies)",choices = list(), multiple = TRUE ),
                             sliderInput(inputId = "sizespatial2",label = "Size",min = 0,max = 5,value = 2.2,step = 0.1),
                             conditionalPanel(
                               condition = "input.AreaAnalysis != 'No'",
                               radioButtons(inputId="RetainSubset5", label="Subset?", choices= list("Area","Inverse", "All"), selected = "All")
                             ),
                             actionButton("plotea", "Plot")
                           ) # Closing fluidPage of CITEtoCODEXspatial independent tab  
          ), # Closing conditionalPanel of CITEtoCODEXspatial independent tab
          
          # DEs table and Umaps Analysis
          menuItem("DE analysis 1", tabName = "DEsTabandUMaps", icon = icon("code-branch")) ,
          conditionalPanel("input.sidebarmenu === 'DEsTabandUMaps'",
                           fluidPage(
                             
                             selectInput(inputId = "colbyDE", label = "Classification category",choices = list() ),
                             
                             selectizeInput(inputId="selectDE", label="DE analysis", choices = list(), options = list(maxItems = 2)),
                             
                             radioButtons("datatypeDE", "Data type:",
                                          c("RNA" = "cite_mRNA_norm",
                                            "Protein" = "cite_protein"),selected = "RNA"),
                             actionButton("calculaDE", "Calculate"),
                             
                             selectizeInput(inputId="NonDEfeatures", label="Add non DE features", choices = NULL, multiple = TRUE),
                             
                             selectizeInput(inputId="featureUMAPDE", label="Select up to 9 features", choices = list(), options = list(maxItems = 9), selected=NULL),
                             
                             sliderInput(inputId = "sizeDE",label = "Size",min = 0,max = 5,value = 1.2,step = 0.1),
                             
                             conditionalPanel(
                               condition = "input.datatypeDE != 'Protein'",
                               sliderInput(inputId="minco", label="Minimum cutoff (only for RNA)", min = 0,max=1,value=0,step=0.1) ),
                             
                             
                             selectInput(inputId = "colorUMAPDE", label = "High positive Color",choices = list("darkred", "red","dodgerblue2", "blue","darkblue", "black", "gray95","gray75"),selected = "darkred" ),
                             selectInput(inputId = "colorUMAPDE2", label = "Low Color",choices = list("darkred", "red","dodgerblue2", "blue","darkblue", "black", "gray95","gray75"),selected = "gray95" ),
                             selectInput(inputId = "colorUMAPDE3", label = "High negative Color",choices = list("darkred", "red","dodgerblue2", "blue","darkblue", "black", "gray95","gray75"),selected = "dodgerblue2" ),
                             actionButton("ploteaUMAPDE", "Plot")
                             
                           ) # Closing fluidPage of DEs table and Umaps Analysis
          ), # Closing conditionalPanel of DEs table and Umaps Analysis
          
          # DE Heatmap, Ridge and Violin plots
          menuItem("DE analysis 2", tabName = "DEHeatmapRidgeandViolin", icon = icon("code-branch")) ,
          conditionalPanel("input.sidebarmenu === 'DEHeatmapRidgeandViolin'",
                           fluidPage(
                             
                             radioButtons(
                               inputId="allorsomefeatures",
                               label="Select all or some features for Heatmap",
                               choices=list(
                                 "All",
                                 "Manual Selection"
                               ),
                               selected="All"),
                             
                             conditionalPanel(
                               condition = "input.allorsomefeatures != 'All'",
                               selectizeInput(inputId="selectfeaturesforheatmap", label="Select Features for Heatmap", choices = list() , multiple = TRUE ) ),
                             
                             radioButtons(
                               inputId="allorsomeGroups",
                               label="Select all or some categories for Heatmap",
                               choices=list(
                                 "All",
                                 "Manual Selection"
                               ),
                               selected="All"),
                             
                             conditionalPanel(
                               condition = "input.allorsomeGroups != 'All'",
                               selectizeInput(inputId="selectgroupsforheatmap", label="Select categories for Heatmap", choices = list(), multiple = TRUE )
                             ),
                             selectizeInput(inputId="selectgroupsorder", label="Sort categories as desired", choices = list(), multiple = TRUE ),
                             numericInput(inputId="cellnumber",label = "Number of Cells to plot", value = 10,min = 10,max=300,step = 1),
                             numericInput(inputId="lowerbreak",label = "Lower Break", value = -5,min = -100,max=100,step = 1),
                             numericInput(inputId="upperbreak",label = "Upper Break", value = 5,min = -100,max=100,step = 1),
                             
                             numericInput(inputId="heatmapPlotsize",label = "Heatmap size (px)", value = 1000,min = 100,max=3000,step = 1),
                             numericInput(inputId="heatmapLegendsize",label = "Legend size", value = 0.8,min = 0.1,max=3,step = 0.1),
                             numericInput(inputId="heatmapFeaturesize",label = "Feature size", value = 0.8,min = 0.1,max=3,step = 0.1),
                             
                             
                             
                             radioButtons("heatmapcolor", "Heatmap color",
                                          c("Purple/Yellow" = "purpleyellow",
                                            "Blue/Red" = "bluered")),
                             
                             actionButton("ploteaHeatmap", "Plot"),
                             
                             selectizeInput(inputId="selectfeaturesforridgeandviolin", label="Select Features for Ridge/Violin", choices = list() , options = list(maxItems = 12)),
                             selectizeInput(inputId="selectgroupsforridgeandviolin", label="Select categories to plot", choices = list(), multiple = TRUE ),
                             actionButton("plotearidgeandviolin", "Plot Ridge/Violin")
                             
                             
                           ) # Closing fluidPage of DE Heatmap, Ridge and Violin plots
          ), # Closing conditionalPanel of DE Heatmap, Ridge and Violin plots
          
          # Expression CODEX spatial Analysis
          menuItem("Spatial Expression", tabName = "CODEXspatialexpression", icon = icon("code-branch")) ,
          conditionalPanel("input.sidebarmenu === 'CODEXspatialexpression'",
                           fluidPage(
                             selectInput(inputId = "clusterSpatialExpression", label = "Classification category",choices = list("predicted.id","predicted.celltype.l2","Higher_Hierarchy_grouping") ),
                             textInput(inputId="IdentificationThreshold2", label="Probability of Identification (only for RNA)", value= 0),
                             selectInput(inputId = "clusterincolorSpatialExpression", label = "Colored category(ies)",choices = list(), multiple = TRUE ),
                             radioButtons("datatypeCODEXspatial", "Data type:",
                                          c("RNA" = "RNA",
                                            "Protein" = "protein"),selected = "Protein"),
                             selectizeInput(inputId="selectfeaturesCODEXspatial", label="Select one or two features", choices = NULL, options = list(maxItems = 2)),
                             sliderInput(inputId = "sizefeatCODEXspatial",label = "Size",min = 0,max = 5,value = 0.8,step = 0.1),
                             radioButtons(inputId="maximize_differences", label="Maximize differences", 
                                            choices = c("Yes" = "yes",
                                                        "No" = "no"),selected = "No"),
                             selectInput(inputId = "colorCODEXspatial", label = "High positive Color",choices = list("darkred", "red","dodgerblue2", "blue","darkblue", "black", "gray95","gray75"),selected = "darkred" ),
                             selectInput(inputId = "colorCODEXspatial2", label = "Low Color",choices = list("darkred", "red","dodgerblue2", "blue","darkblue", "black", "gray95","gray75"),selected = "gray95" ),
                             selectInput(inputId = "colorCODEXspatial3", label = "High negative Color",choices = list("darkred", "red","dodgerblue2", "blue","darkblue", "black", "gray95","gray75"),selected = "dodgerblue2" ),
                             conditionalPanel(
                               condition = "input.datatypeCODEXspatial != 'Protein'",
                               sliderInput(inputId="mincoCODEXspatial", label="Minimum cutoff (only for RNA)", min = 0,max=1,value=0.60,step=0.01) ),
                             conditionalPanel(
                               condition = "input.AreaAnalysis != 'No'",
                               radioButtons(inputId="RetainSubset6", label="Subset?", choices= list("Area","Inverse", "All"), selected = "All"),
                               radioButtons(inputId="Context", label="Isolated or in space?", choices= list("Isolated","In_Space"), selected = "In_Space")
                             ),
                             conditionalPanel(
                               condition = "input.clusterincolorSpatialExpression.length == 2 & input.selectfeaturesCODEXspatial.length == 1 & input.RetainSubset6 !='All'",
                               radioButtons(inputId="TGOFanalysis", label="TwoGroups_OneFeature Analysis?", choices= list("Yes","No"), selected = "No")
                             ),
                             conditionalPanel(
                               condition = "input.TGOFanalysis == 'Yes'",
                               radioButtons(inputId="OnlyShowGroups", label="Show only the groups?", choices= list("Yes","No"), selected = "No")
                             ),                             
                             actionButton("ploteaCODEXspatial", "Plot")
                             
                           ) # Closing fluidPage of Expression CODEX spatial Analysis
          ), # Closing conditionalPanel of Expression CODEX spatial Analysis
          
          # Expression CODEX UMAP Analysis
          menuItem("CODEX UMAP Clusters", tabName = "CODEXUMAPclusters", icon = icon("code-branch")) ,
          conditionalPanel("input.sidebarmenu === 'CODEXUMAPclusters'",
                           fluidPage(
                             selectizeInput(inputId="selectclusteringtypeUMAPclusters", label="Select clustering type", choices = list(), options = list(maxItems = 1),selected = "orig.ident" ),
                             textInput(inputId="IdentificationThreshold3", label="Probability of Identification (only for RNA-related features)", value= 0),
                             sliderInput(inputId = "sizefeatCODEXUMAPclusters",label = "Size",min = 0,max = 5,value = 0.8,step = 0.1),
                             radioButtons("SelectiononTopCODEXUMAPAnalysis", "Selection on Top:",
                                          c("Yes" = "yes",
                                            "No" = "no")),
                             
                             actionButton("ploteaCODEXUMAPclusters", "Plot")

                           ) # Closing fluidPage of Expression CODEX UMAP clusters Analysis
          ), # Closing conditionalPanel of Expression CODEX UMAP clusters Analysis

          # CODEX independent Tab
          menuItem("CODEX independent", tabName = "CODEXindep", icon = icon("code-branch")),
          conditionalPanel("input.sidebarmenu === 'CODEXindep'",
                           fluidPage(
                             # File Input 1
                             selectInput(inputId = "highlightCODEXindep", label = "Category to highlight",choices = list() ),
                             selectInput(inputId = "colorCODEXindep", label = "Color",choices = list("darkred", "red", "blue","darkblue", "black"),selected = "darkred" ),
                             textInput(inputId="IdentificationThreshold4", label="Probability of Identification (only for RNA-related features)", value= 0),
                             sliderInput(inputId = "sizeCODEX",label = "Size",min = 0,max = 5,value = 1.2,step = 0.1),
                             radioButtons("SelectiononTopCODEX", "Selection on Top:",
                                          c("Yes" = "yes",
                                            "No" = "no")),
                             actionButton("ploteaCODEXindep", "Plot")
                           ) # Closing fluidPage of CODEX independent Tab  
          ), # Closing conditionalPanel of CODEX independent Tab
          
          
          # CODEX Adjacency analysis by feature
          menuItem("Adjacency analysis by feature", tabName = "AdjacencybyFeature", icon = icon("code-branch")) ,
          conditionalPanel("input.sidebarmenu === 'AdjacencybyFeature'",
                           fluidPage(
                             radioButtons("datatypeAdjacencybyFeature", "Data type:",
                                          c("RNA" = "RNA",
                                            "Protein" = "protein")),
                             conditionalPanel(
                               condition = "input.AreaAnalysis != 'No'",
                               radioButtons(inputId="RetainSubset9", label="Subset?", choices= list("Area","Inverse","Both", "All Data"), selected = "All Data"),
                               
                             ),
                             
                             selectizeInput(inputId="selectfeaturesAdjacencybyFeature", label="Select features", choices = NULL, options = list(maxItems = 200))  ,
                             selectInput(inputId = "color1AdjacencybyFeature", label = "Low color",choices = list("darkred", "red","dodgerblue2", "blue","darkblue", "black","white","yellow", "gray95","gray75"),selected = "white" ),
                             
                             selectInput(inputId = "color2AdjacencybyFeature", label = "High color",choices = list("darkred", "red","dodgerblue2", "blue","darkblue", "black","white","yellow", "gray95","gray75"),selected = "darkred" ),
                             
                             actionButton("ploteaAdjacencybyFeature", "Plot")
                             
                           ) # Closing fluidPage of CODEX Adjacency analysis by feature
          ), # Closing conditionalPanel of CODEX Adjacency analysis by feature
          
          # CODEX Adjacency analysis by cluster
          menuItem("Adjacency analysis by cluster", tabName = "AdjacencybyCluster", icon = icon("code-branch")) ,
          conditionalPanel("input.sidebarmenu === 'AdjacencybyCluster'",
                           fluidPage(
                             
                               selectizeInput(inputId="selectclusteringtype", label="Select clustering type", choices = list(), options = list(maxItems = 1),selected = "orig.ident" ),
                               selectInput(inputId = "clusterincolorbyCluster", label = "Select groups",choices = list(), multiple = TRUE ),
                               
                               conditionalPanel(
                                 condition = "input.AreaAnalysis != 'No'",
                                 radioButtons(inputId="RetainSubset10", label="Subset?", choices= list("Area","Inverse","Both", "All Data"), selected = "All Data"),
                               ),
                               
                               selectInput(inputId = "color1AdjacencybyCluster", label = "Low color",choices = list("darkred", "red","dodgerblue2", "blue","darkblue", "black","white","yellow", "gray95","gray75"),selected = "white" ),
                             
                             selectInput(inputId = "color2AdjacencybyCluster", label = "High color",choices = list("darkred", "red","dodgerblue2", "blue","darkblue", "black","white","yellow", "gray95","gray75"),selected = "darkred" ),
                             
                             actionButton("ploteaAdjacencybyCluster", "Plot")
                             
                           ) # Closing fluidPage of CODEX Adjacency analysis by cluster
          ) # Closing conditionalPanel of CODEX Adjacency analysis by cluster
          
         )# Closing sidebarMenu
    ), #Closing dashboardSidebar
    
    ### BODY
    dashboardBody(
      tabItems(
        tabItem(tabName = "fileupload",
                fluidRow(plotOutput ("p1", height = 1000), style =  "width: 95%"
                ), #Closing the fluid Row of CITE UMAP analysis
                br(),
                fluidRow(DT::DTOutput ("tproportions"), style =  "font-size: 95%; width: 95%"
                ) #Closing the fluid Row of the table of proportions of classification categories
        ), # Closing tabItem of CITE UMAP analysis  
        
        tabItem(tabName = "UMAPIndSamples",
                fluidRow(plotOutput ("p5", height = 1000), style =  "width: 95%"
                ) #Closing the fluid Row of CITE UMAP of samples plotted individually
        ),  # Closing tabItem of CITE UMAP of samples plotted individually
        
        tabItem(tabName = "CITEindep",
                fluidRow(plotOutput ("p2", height = 1000), style =  "width: 95%"
                ) #Closing the fluid Row of CITE independent
        ), # Closing tabItem of CITE independent       
        
        tabItem(tabName = "spatial1",
                fluidRow(plotOutput ("p3", height = 1000), style =  "width: 95%"
                         #fluidRow(textOutput ("p3")
                ),  #Closing the fluid Row of CITEtoCODEXspatial extrapolation
                br(),
                fluidRow(DT::DTOutput ("tproportionsCODEX"), style =  "font-size: 95%; width: 95%"
                ) #Closing the fluid Row of the table of proportions of classification categories of CODEX cells
        ),  # Closing tabItem of CITEtoCODEXspatial extrapolation
        
        tabItem(tabName = "spatial2",
                fluidRow(plotOutput ("p4", height = 1000), style =  "width: 95%"
                         #fluidRow(textOutput ("p3")
                ) #Closing the fluid Row of CITEtoCODEXspatial independent
        ),  # Closing tabItem of CITEtoCODEXspatial independent
        
        tabItem(tabName = "DEsTabandUMaps",
                fluidRow(DT::DTOutput ("tDE"), style =  "font-size: 95%; width: 95%"
                ), #Closing the fluid Row of the table of DEsTabandUMaps
                fluidRow(plotOutput ("p6", height = 1000), style =  "width: 95%"
                         #fluidRow(textOutput ("p3")
                )#Closing the fluid Row of the UMAP plots of DEsTabandUMaps
        ), # Closing tabItem of DEsTabandUMaps
        
        tabItem(tabName = "DEHeatmapRidgeandViolin",
                fluidRow(DT::DTOutput ("tDE2"), style =  "font-size: 95%; width: 95%"
                ), #Closing the fluid Row of the table of DEHeatmapRidgeandViolin
                fluidRow(plotOutput ("p7",height = "1500px"), style =  "width: 95%"
                ), #Closing the fluid Row of the heatmap plot of DEHeatmapRidgeandViolin
                br(),
                fluidRow(plotOutput ("p8", height = 1000), style =  "width: 95%"
                ), #Closing the fluid Row of the ridge plot 1 of DEHeatmapRidgeandViolin
                br(),
                fluidRow(plotOutput ("p9", height = 1000), style =  "width: 95%"
                ), #Closing the fluid Row of the violin plot 2 of DEHeatmapRidgeandViolin
                br(),
                fluidRow(plotOutput ("p10", height = 1000), style =  "width: 95%"
                ), #Closing the fluid Row of the violin plot 3 of DEHeatmapRidgeandViolin
                br(),
                fluidRow(plotOutput ("p11", height = 1000), style =  "width: 95%"
                ), #Closing the fluid Row of the violin plot 1 of DEHeatmapRidgeandViolin
                br(),
                fluidRow(plotOutput ("p12", height = 1000), style =  "width: 95%"
                ), #Closing the fluid Row of the violin plot 2 of DEHeatmapRidgeandViolin
                br(),
                fluidRow(plotOutput ("p13", height = 1000), style =  "width: 95%"
                ) #Closing the fluid Row of the violin plot 3 of DEHeatmapRidgeandViolin
        ),  # Clossing tabItem of DEHeatmapRidgeandViolin
        
        tabItem(tabName = "CODEXspatialexpression",
                
                fluidRow(plotOutput ("p14", height = 1000), style =  "width: 95%"
                         
                )#Closing the fluid Row of the Spatial plots of CODEXspatialexpression
        ),  # Clossing tabItem of CODEXspatialexpression
        
        tabItem(tabName = "CODEXUMAPclusters",
                fluidRow(plotOutput ("p19", height = 1000), style =  "width: 95%"
                ), #Closing the fluid Row of the CODEXUMAPclusters
                br(),
                fluidRow(DT::DTOutput ("tproportionsCODEX2"), style =  "font-size: 95%; width: 95%"
                ) #Closing the fluid Row of the table of proportions of classification categories of CODEX cells
        ), # Closing tabItem of the CODEXUMAPclusters
        
        tabItem(tabName = "CODEXindep",
                fluidRow(plotOutput ("p20", height = 1000), style =  "width: 95%"
                ) #Closing the fluid Row of CITE independent
        ), # Closing tabItem of CITE independent     
        
        
        
        tabItem(tabName = "AdjacencybyFeature",
                        fluidRow(plotOutput ("p15", height = 1000), style =  "width: 95%"
                        ), #Closing the fluid Row of the adjacency heatmap by feature
                        br(),
                        fluidRow(plotOutput ("p16", height = 1000), style =  "width: 95%"
                        ), #Closing the fluid Row of the Pearson correlation
                        br(),
                        fluidRow(DT::DTOutput ("AdjacencyTable"), style =  "font-size: 95%; width: 95%" 
                        ), #Closing the fluid Row of the adjacency table by feature
                        br(),
                        fluidRow(plotOutput ("p17", height = 1000), style =  "width: 95%"
                        ), #Closing the fluid Row of the adjacency heatmap of an Area by feature
                        br(),
                        fluidRow(plotOutput ("p18", height = 1000), style =  "width: 95%"
                        ), #Closing the fluid Row of the Pearson correlation of an Area
                        br(),
                        fluidRow(DT::DTOutput ("AdjacencyTableArea"), style =  "font-size: 95%; width: 95%" 
                        ), #Closing the fluid Row of the adjacency table of an Area by feature
                        br(),
                        fluidRow(plotOutput ("p23", height = 1000), style =  "width: 95%"
                        ), #Closing the fluid Row of the adjacency heatmap of an Area by feature
                        br(),
                        fluidRow(plotOutput ("p24", height = 1000), style =  "width: 95%"
                        ), #Closing the fluid Row of the Pearson correlation of an Area
                        br(),
                        fluidRow(DT::DTOutput ("AdjacencyTableANTIArea"), style =  "font-size: 95%; width: 95%" 
                        ) #Closing the fluid Row of the adjacency table of an Area by feature
        ), # Closing tabItem of AdjacencybyFeature
        
        tabItem(tabName = "AdjacencybyCluster",
                        fluidRow(plotOutput ("p21", height = 1000), style =  "width: 95%"
                        ), #Closing the fluid Row of the adjacency heatmap by feature
                        br(),
                        fluidRow(plotOutput ("p22", height = 1000), style =  "width: 95%"
                        ), #Closing the fluid Row of the Pearson correlation
                        br(),
                        fluidRow(DT::DTOutput ("AdjacencyTableCluster"), style =  "font-size: 95%; width: 95%"
                        ), #Closing the fluid Row of the adjacency table by cluster
                        br(),
                        fluidRow(plotOutput ("p25", height = 1000), style =  "width: 95%"
                        ), #Closing the fluid Row of the adjacency heatmap of an Area by feature
                        br(),
                        fluidRow(plotOutput ("p26", height = 1000), style =  "width: 95%"
                        ), #Closing the fluid Row of the Pearson correlation of an Area
                        br(),
                        fluidRow(DT::DTOutput ("AdjacencyTableClusterArea"), style =  "font-size: 95%; width: 95%" 
                        ), #Closing the fluid Row of the adjacency table of an Area by feature
                        br(),
                        fluidRow(plotOutput ("p27", height = 1000), style =  "width: 95%"
                        ), #Closing the fluid Row of the adjacency heatmap of an Area by feature
                        br(),
                        fluidRow(plotOutput ("p28", height = 1000), style =  "width: 95%"
                        ), #Closing the fluid Row of the Pearson correlation of an Area
                        br(),
                        fluidRow(DT::DTOutput ("AdjacencyTableClusterANTIArea"), style =  "font-size: 95%; width: 95%" 
                        ) #Closing the fluid Row of the adjacency table of an Area by feature
                
        ) # Closing tabItem of AdjacencybyCluster       
        
      ) # Closing TabItems
    ) # Closing dashboardBody
  ) # Closing dashboardPage
) # Closing shinyUI


# END
