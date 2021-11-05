#############################################
########     InterroGATOR_ui.R     ##########
#############################################
# Author: Leandro Balzano-Nogueira
# Diabetes Institute, University of Florida (Gainesville)
# Last update: October/19/2021

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
library(DT)


#############################################
## ui.R ##
options(shiny.maxRequestSize=120000*1024^2)

#eltitle <- tags$a("InteGrATOR",tags$img(src="Gator.jpg", height="40", width= "50"))
upload_data <- tags$a(span("InterroGATOR", style = "color: white; font-size: 26px"),tags$img(src="Gator.jpg", height="40", width= "50"))
# dashboardHeader(title = span("Basic ", 
#                              span("dashboard", 
#                                   style = "color: red; font-size: 28px")))
shinyUI(
  dashboardPage(skin = "blue", #skin color
   
    ### HEADER
    dashboardHeader(title=upload_data), ## Header
    #dashboardHeader(title = span("InteGrATOR", style = "color: white; font-size: 28px")),
    
    ### Side Bar
    dashboardSidebar(
      sidebarMenu(id = "sidebarmenu", ### create a menu in the sidebar
          # CITE UMAP Analysis
          menuItem("File Upload", tabName = "fileupload", icon = icon("code-branch")) ,
          conditionalPanel("input.sidebarmenu === 'fileupload'",
          fluidPage(
              # File Input 1
              
              fileInput("file1", label = h4("Insert Seurat rds file")),
              
                                   
              
              fileInput("file2", label = h4("Insert STvEA RData object")),
              
              selectInput(inputId = "colby", label = "Classification category",choices = list() ),             
              sliderInput(inputId = "size1",label = "Size",min = 0,max = 5,value = 1.2,step = 0.1)
              
              ) # Closing fluidPage of input  
            ), # Closing conditionalPanel of Input
          
          # CITE UMAP of samples plotted individually
          menuItem("UMAP individual samples", tabName = "UMAPIndSamples", icon = icon("code-branch")),
          conditionalPanel("input.sidebarmenu === 'UMAPIndSamples'",
                           fluidPage(
                             # File Input 1
                             selectInput(inputId = "Umapstoplot", label = "Sample to plot",choices = list() ),
                             selectInput(inputId = "clusterUMAPIndSamples", label = "Classification category",choices = list("predicted.celltype.l2","Higher_Hierarchy_grouping") ),
                             sliderInput(inputId = "sizeumapindsamples",label = "Size",min = 0,max = 5,value = 1.2,step = 0.1),
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
            sliderInput(inputId = "size",label = "Size",min = 0,max = 5,value = 1.2,step = 0.1),
            radioButtons("SelectiononTop", "Selection on Top:",
                         c("Yes" = "yes",
                           "No" = "no")),
            actionButton("ploteaCITEIndep", "Plot")
                           ) # Closing fluidPage of CITE independent Tab  
          ), # Closing conditionalPanel of CITE independent Tab
          
          # CITEtoCODEXspatial extrapolation tab
          menuItem("CITEtoCODEXspatial extrapolation", tabName = "spatial1", icon = icon("code-branch")),
          conditionalPanel("input.sidebarmenu === 'spatial1'",
                           fluidPage(
                             # File Input 1
                             selectInput(inputId = "cluster1", label = "Classification category",choices = list() ),
                             sliderInput(inputId = "size2",label = "Size",min = 0,max = 5,value = 1.2,step = 0.1)
                           ) # Closing fluidPage of CITEtoCODEXspatial extrapolation tab  
          ), # Closing conditionalPanel of CITEtoCODEXspatial extrapolation tab
          
          # CITEtoCODEXspatial independent tab
          menuItem("CITEtoCODEXspatial independent", tabName = "spatial2", icon = icon("code-branch")),
          conditionalPanel("input.sidebarmenu === 'spatial2'",
                           fluidPage(
                             
                             selectInput(inputId = "cluster2", label = "Classification category",choices = list("predicted.celltype.l2","Higher_Hierarchy_grouping") ),
                             selectInput(inputId = "clusterincolor", label = "Colored category(ies)",choices = list(), multiple = TRUE ),
                             sliderInput(inputId = "sizespatial2",label = "Size",min = 0,max = 5,value = 1.2,step = 0.1),
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
                                            "Protein" = "cite_protein")),
                             actionButton("calculaDE", "Calculate"),
                             
                             selectizeInput(inputId="NonDEfeatures", label="Add non DE features", choices = NULL, multiple = TRUE),
                             
                             selectizeInput(inputId="featureUMAPDE", label="Select up to 9 features", choices = list(), options = list(maxItems = 9), selected=NULL),
                             sliderInput(inputId = "sizeDE",label = "Size",min = 0,max = 5,value = 1.2,step = 0.1),
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
                             radioButtons("heatmapcolor", "Heatmap color",
                                          c("Purple/Yellow" = "purpleyellow",
                                            "Blue/Red" = "bluered")),
                             
                             actionButton("ploteaHeatmap", "Plot"),
                             
                             selectizeInput(inputId="selectfeaturesforridgeandviolin", label="Select Features for Ridge/Violin", choices = list() , options = list(maxItems = 12)),
                             selectizeInput(inputId="selectgroupsforridgeandviolin", label="Select categories to plot", choices = list(), multiple = TRUE ),
                             actionButton("plotearidgeandviolin", "Plot Ridge/Violin")
                             
                             
                           ) # Closing fluidPage of DE Heatmap, Ridge and Violin plots
          ) # Closing conditionalPanel of DE Heatmap, Ridge and Violin plots
          
          
         )# Closing sidebarMenu
    ), #Closing dashboardSidebar
    
    ### BODY
    dashboardBody(
      tabItems(
        tabItem(tabName = "fileupload",
                fluidRow(plotOutput ("p1")
                ), #Closing the fluid Row of CITE UMAP analysis
                br(),
                fluidRow(DT::DTOutput ("tproportions"), style =  "font-size: 95%; width: 95%"
                ) #Closing the fluid Row of the table of proportions of classification categories
        ), # Closing tabItem of CITE UMAP analysis  
        
        tabItem(tabName = "UMAPIndSamples",
                fluidRow(plotOutput ("p5")
                ) #Closing the fluid Row of CITE UMAP of samples plotted individually
        ),  # Closing tabItem of CITE UMAP of samples plotted individually
        
        tabItem(tabName = "CITEindep",
                fluidRow(plotOutput ("p2")
                ) #Closing the fluid Row of CITE independent
        ), # Closing tabItem of CITE independent       
        
        tabItem(tabName = "spatial1",
                fluidRow(plotOutput ("p3")
                         #fluidRow(textOutput ("p3")
                ) #Closing the fluid Row of CITEtoCODEXspatial extrapolation
        ),  # Closing tabItem of CITEtoCODEXspatial extrapolation
        
        tabItem(tabName = "spatial2",
                fluidRow(plotOutput ("p4")
                         #fluidRow(textOutput ("p3")
                ) #Closing the fluid Row of CITEtoCODEXspatial independent
        ),  # Closing tabItem of CITEtoCODEXspatial independent
        
        tabItem(tabName = "DEsTabandUMaps",
                fluidRow(DT::DTOutput ("tDE"), style =  "font-size: 95%; width: 95%"
                ), #Closing the fluid Row of the table of DEsTabandUMaps
                fluidRow(plotOutput ("p6")
                         #fluidRow(textOutput ("p3")
                )#Closing the fluid Row of the UMAP plots of DEsTabandUMaps
        ), # Closing tabItem of DEsTabandUMaps
        
        tabItem(tabName = "DEHeatmapRidgeandViolin",
                fluidRow(DT::DTOutput ("tDE2"), style =  "font-size: 95%; width: 95%"
                ), #Closing the fluid Row of the table of DEHeatmapRidgeandViolin
                fluidRow(plotOutput ("p7")
                ), #Closing the fluid Row of the heatmap plot of DEHeatmapRidgeandViolin
                fluidRow(plotOutput ("p8")
                ), #Closing the fluid Row of the ridge plot 1 of DEHeatmapRidgeandViolin
                fluidRow(plotOutput ("p9")
                ), #Closing the fluid Row of the violin plot 2 of DEHeatmapRidgeandViolin
                fluidRow(plotOutput ("p10")
                ), #Closing the fluid Row of the violin plot 3 of DEHeatmapRidgeandViolin
                fluidRow(plotOutput ("p11")
                ), #Closing the fluid Row of the violin plot 1 of DEHeatmapRidgeandViolin
                fluidRow(plotOutput ("p12")
                ), #Closing the fluid Row of the violin plot 2 of DEHeatmapRidgeandViolin
                fluidRow(plotOutput ("p13")
                ) #Closing the fluid Row of the violin plot 3 of DEHeatmapRidgeandViolin
        )  # Clossing tabItem of DEHeatmapRidgeandViolin
        
        
      ) # Closing TabItems
    ) # Closing dashboardBody
  ) # Closing dashboardPage
) # Closing shinyUI


# END
