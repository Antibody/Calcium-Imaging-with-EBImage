library(shiny)
library("EBImage")
library(pheatmap)
library("magrittr")
library("tibble")
library("ggplot2")
#library("genefilter")
library("GGally")
library(geomtextpath)
library(reshape2)
library(pracma)
#library(dplyr)
library('shinythemes')
library('shinycssloaders')
library(shinybusy)
library(BiocManager)
options(repos = BiocManager::repositories())
knitr::opts_chunk$set(echo = TRUE)

ui <- fluidPage(theme = shinytheme("paper"),
                titlePanel(h6("Calcium Imaging v.0.1.4")),
                
                #++++++++++++Busy indicator +++++++++++++++
                add_busy_spinner(spin = "radar", margins = c(400, 1000), color ="darkblue"),
                #++++++++++++CSS for notifications +++++++++++++++
                tags$head(
                  tags$style(
                    HTML(".shiny-notification {
             position:fixed;
             top: calc(1%);
             left: calc(10%);
             }
             "
                    )
                  )
                ),
                
                sidebarLayout(
                  sidebarPanel(
                    fluidRow( 
                      fileInput(inputId = 'files', 
                                label = 'Select an Image',
                                multiple = TRUE,
                                accept=c('image/tif', 'image/jpeg')),
                      
                      
                      wellPanel(h5("Image mask adjustments"),
                                
                                # A parameter for EBImage's makeBrush
                                numericInput(
                                  inputId="sigma",
                                  label="Image blurring (sigma)",
                                  value=7,
                                  min = 1,
                                  max = 111,
                                  step = 2
                                ),
                                
                                
                                numericInput(
                                  inputId="size",
                                  label="Size (px)",
                                  value=101,
                                  min = 1,
                                  max = 1011,
                                  step = 2
                                ),
                                
                                sliderInput("offset", h6("offset"), 0, 1, 0.07),
                                
                                numericInput(
                                  inputId="rad",
                                  label="masking radius",
                                  value=9,
                                  min = 1,
                                  max = 111,
                                  step = 2
                                ),
                                
                      )
                      
                    ),
                  ),
                  mainPanel(
                    
                    
                    tabsetPanel(
                      
                      tabPanel("Mask", 
                               br(),
                               plotOutput("rasterMask"), 
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               actionButton("addLabels", "Show Mask and Labels"), 
                               br(),
                               
                               checkboxInput(inputId = "rmvCell",
                                             label = "Remove cell(s) from the calculations?",
                                             value = FALSE),
                               conditionalPanel(
                                 condition = "input.rmvCell == true",
                                 tags$div(
                                   HTML(paste("Add cell numbers separated by the coma. Start with the ", tags$span(style="color:blue", "highest cell number e.g. 41, 23, 1"), sep = ""))
                                 ),
                                 textInput("removeCell", "",  value = "1000")
                               ),
                               
                      ),
                      
                      tabPanel("Plot",
                               
                               
                               
                               
                               br(),
                               
                               fluidRow( column(4, wellPanel(h5("Plot and data options"),
                               checkboxInput(inputId = "detrend",
                                             label = "Detrend my data",
                                             value = FALSE),
                               
                               conditionalPanel(
                                 condition = "input.detrend == true",
                               checkboxInput(inputId = "break_point",
                                             label = "Does data contain break points",
                                             value = FALSE),
                               
                               conditionalPanel(
                                 condition = "input.break_point == true",
                                 numericInput("bp", "Break point value (it should be higher than your peaks base) :", value=10, min = 1, max = 10000),
                               
                              
                               
                               )),
                               
                               br(),
                               actionButton("buildPlot", "Build plot"), 
                               
                               )),
                               
                               column(4, wellPanel(h5(""),
                                                             
                                                               
                                                               checkboxInput(inputId = "subtract",
                                                                       label = "Subtract background",
                                                                       value = FALSE),            
                                                   
                                                   
                                                               checkboxInput(inputId = "select",
                                                                             label = "Select individual trace",
                                                                             value = FALSE),
                                                               conditionalPanel(
                                                                 condition = "input.select == true",
                                                                 numericInput("trace", "Trace number:", value=1, min = 1, max = 1000),
                                                                    checkboxInput(inputId = "deltaF",
                                                                               label = "Calculate deltaF/F(0)",
                                                                               value = FALSE),
                                                               
                                                             ),
                                                             
                               ))),
                               
                               br(),
                               plotOutput("plot", height = "100%"),
                               
                               br(),
                               downloadButton("downloadPlotData", "Download data"),
                               
                               
                      ),
                      tabPanel("Correlations", 
                               br(),
                               br(),
                               h4(tags$b("Correlation Table", style = "color:darkblue")),
                               sliderInput("adjPval", h6("Adjusted p-value threshold"), 0, 0.05, 0.01),
                               sliderInput("corr", h6("Correlation coefficient threshold"), 0, 1, 0.8),
                               actionButton("buildTable", "Calculate correlation table"),
                               downloadButton("downloadTable", "Download"),
                               verbatimTextOutput('table'),
                               
                               h4(tags$b("Build heatmap", style = "color:darkblue")),
                               actionButton("buildHeat", "Build"),
                               br(),
                               plotOutput("heat", height = "100%"),
                               
                               
                               
                      ),
                      
                      
                      
                      
                      tabPanel("Interactive image browser", displayOutput("widget")) 
                      
                    )
                  )
                )
)

server <- function(input, output, session) {
  
  output$files <- renderTable(input$files)
  
  files <- reactive({
    files <- input$files
    files$datapath <- gsub("\\\\", "/", files$datapath)
    files
  })
  
  #+++++++++++++++++++++++++ This function extracts numbers from text +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
  numbersFromText <- function(text) {        # this function extracts numbers from text e.g. for the removal of unwanted cells
    text <- gsub(" ", "", text)
    split <- strsplit(text, ",", fixed = FALSE)[[1]]
    as.numeric(split)
  }
  
  
  #+++++++++++++++++++++++ send message about subtract function ++++++++++++++++++++++++++#
  
  observeEvent(input$subtract , {
    if (input$subtract==TRUE)  {
      showNotification("Background is calculated from all of the area that is considered as non-cells. It also adds background values as the last line in the downlodable report table", duration = 10, type = "message")
    }                                      
  })
  
  
  # #+++++++++++++++++++++++ send message negative min value in the trace that is used to calculate dF/F ++++++++++++++++++++++++++#
  # observeEvent(input$buildPlot , {
  #   if (minDataID < 0)  {
  #     showNotification("For the calculation of dF/F a negative minimum value was used", duration = 10, type = "message")
  #   }                                      
  # })
  
  
  
  #+++++++++++++++++++++++ unchecks subtract box if individual traces are chosen and vice-versa +++++++++++++++++++#
  observeEvent(input$deltaF == TRUE ,
               if (input$deltaF == TRUE)  {
                 { #updatedValue = !input$nonParametricTests
                   updateCheckboxInput(session, "subtract",  value=FALSE)
                   updateCheckboxInput(session, "deltaF",  value=TRUE)
                 }
               }
  )
  
  
  observeEvent(input$subtract == TRUE ,
               if (input$subtract == TRUE)  {
                 { #updatedValue = !input$nonParametricTests
                   updateCheckboxInput(session, "deltaF",  value=FALSE)
                   updateCheckboxInput(session, "subtract",  value=TRUE)
                 }
               }
  )
  
  ##########################################################################################################
  ################# Build Plot #############################################################################
  
  
  observeEvent(input$buildPlot, {
    output$plot <- renderPlot({
      img3 = readImage(files = input$files$datapath)
      img3_F1 = readImage((files = input$files$datapath)[1]) # first frame only
      
      img3 <- resize(img3, 512, 512)
      img3_F1 <- resize(img3_F1, 512, 512)
      
      
      cellsSmooth = Image(dim = dim(img3_F1))
      
      sigma <- input$sigma
      size <- input$size
      
      cellsSmooth = filter2(
        img3_F1,
        filter = makeBrush(size = size, shape = "gaussian",
                           sigma=sigma)
      )
      
      disc = makeBrush(size, "disc")
      disc = disc / sum(disc)
      offset = input$offset
      rad <- input$rad
      rmCell <- numbersFromText(input$removeCell)
      
      nucThresh = (cellsSmooth - filter2( cellsSmooth, disc ) > offset)
      
      nucOpened = EBImage::opening(nucThresh, kern = makeBrush(rad, shape = "disc"))
      
      nucSeed = bwlabel(nucOpened)
      
      nucMask = cellsSmooth - filter2(cellsSmooth, disc) > 0
      nucMask = fillHull(nucMask)
      nuclei = propagate(cellsSmooth, nucSeed, mask = nucMask)
      #EBImage::display(nuclei,method = "raster")
      
      
      nucSegOnNuc  = paintObjects(nuclei, tgt = toRGB(img3_F1), col = "#ffff00")
      #EBImage::display(nucSegOnNuc,method = "raster")
      
      
      fr1 = computeFeatures(nuclei,     img3_F1, xname = "Frame1",  refnames = "c1") # this is used to determine how many ROI were detected in the first frame
      
      
      ########################################################################################
      ################# background correction ################################################
      if (input$subtract == T) {
        
        
        data <- data.frame(col1 = rep(NA, dim(fr1)[1]))
        
        
        for(i in 1:dim(img3)[3]) {                             # Head of for-loop
          
          
          new_col <- computeFeatures(nuclei,     img3[,,i], xname = "",   
                                     refnames = "fr_")                      # Creating new variable
          data[ , i] <- new_col[,12]                     # Adding new variable to data
          colnames(data)[i] <- paste0("", i)    # Renaming new variable
        }
        
        whiteImg <- matrix(1, dim(img3_F1), dim(img3_F1)) # create an array of 1s, which will be rendered white
        
        bkg <- whiteImg - nuclei #subtract cells masks from white image to get a mask o a background
        
        #display(bkg, "raster")
        
        data_bkg <- data.frame(col1 = rep(NA, 1))
        
        
        
        for(i in 1:dim(img3)[3]) {                             # A for-loop to create a vector with background intensities timeline
          
          
          new_col <- computeFeatures(bkg,     img3[,,i], xname = "",   
                                     refnames = "fr_")                      # Creating new column for each cell
          data_bkg[ , i] <- new_col[,12]                     # Adding new columns with next timeframe data to a DF
          # colnames(data_bkg)[i] <- paste0("", i)    # Adding a first letter to column names (currently empty)
        }
        
        data_subtr <- sweep(as.matrix(data), MARGIN = 2, as.matrix(data_bkg)) # use SWEEP function to subtract a vector (has to be converted 
        # to a matrix to work with sweep) with background from data DF (also as a matrix) row by row (this is why I use "MARGIN=2") 
        
        data <- as.data.frame(data_subtr)
      }
      
      ###################### If background subtraction is NOT checked ###############################################
      
      if (input$subtract == F) {
        data <- data.frame(col1 = rep(NA, dim(fr1)[1]))       # Creating a new placeholder dataframe with number of rows calulated by: computeFeatures(nuclei,     img3_F1, xname = "Frame1",  refnames = "c1")
        
        
        for(i in 1:dim(img3)[3]) {                             # Head of for-loop
          
          # new dataframe with ROI calculations over entire image stack
          new_col <- computeFeatures(nuclei,     img3[,,i], xname = "P",   # xname for naming the first name of various parameters columns
                                     refnames = "fr_")                      # Creating new variable
          data[ , i] <- new_col[,12]                     # Adding column with intensities [12] to a new dataframe
          colnames(data)[i] <- paste0("", i)    # Renaming new variable
          data
        }
      }
      
      ###################### remove unwanted cells ###############################################
      
      
      if (input$rmvCell == T) {
        
        data <- data[-rmCell, ]
        data
      }
      
   

     
      
      
      ##++++++++++++++++++++++++++++ Detrending ++++++++++++++++++++++++++++++++++++++++###
      if (input$detrend == TRUE) {
       
       
        ###################### break points ############################################### 
        if (input$break_point == F && input$detrend == TRUE) {
          # bp <- input$bp
          # break.points <- seq(from=bp,to=dim(data)[2], by=bp)

          data.t <- t(data) # matrix transposition
          colnames(data.t) <- rownames(data) #restoration of column names
          data.t.df <- as.data.frame(data.t) # conversion to the list formate

          data.t.df.dt <- data.t.df

          #detrending without breaking points###
          for(i in 1:dim(data)[1]){                                       # run each Sample in a loop
            data.t.df.dt[,i] <- detrend(data.t.df[,i], tt = 'linear', bp=c())       # detrend the data
                                   }

          data.dt <- t(data.t.df.dt)
          data.dt.df <- as.data.frame(data.dt)
          dataID <- data.dt.df

                                                              }

      else if (input$break_point == T && input$detrend == TRUE) {
          #bp <- input$bp
          break.points <- seq(from=input$bp,to=dim(data)[2], by=input$bp)
          
          data.t <- t(data) # matrix transposition
          colnames(data.t) <- rownames(data) #restoration of column names
          data.t.df <- as.data.frame(data.t) # conversion to the list formate
          
          data.t.df.dt <- data.t.df 
          
          #detrending with breaking points### 
          for(i in 1:dim(data)[1]){                                       # run each Sample in a loop
            data.t.df.dt[,i] <- detrend(data.t.df[,i], bp=break.points)       # detrend the data
                                  }
          
          data.dt <- t(data.t.df.dt)
          data.dt.df <- as.data.frame(data.dt)
          dataID <- data.dt.df
          dataID
                                 }
        
        
              
              #++++++++++++++++++++ select a single trace ++++++++++++++#
              if (input$select == T) {
                dataID <- (dataID[input$trace, ])
                #++++++++++++++++++++ calculate deltaF/F based on trace minimal value++++++++++++++#
                if (input$deltaF == T) {
                  #dataID <- (dataID[input$trace, ])
                  minDataID <- min(dataID)
                  if (minDataID > 0) {
                    dataID <- (dataID - minDataID)/minDataID
                  }
                  else if (minDataID < 0) {
                    
                    dataID = (dataID - minDataID)/minDataID
                    dataID = -dataID
                    
                  }
                }
                dataID
              }
              
              
              dataID$id = 1:nrow(dataID)
              
              
              
              dataMelt = melt(dataID, variable.name = "Frame", value.name = "Intensity", id.vars = "id")
                                         }
      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++###

      else if (input$detrend == FALSE) {
        
        
        
        # dat31 <- c()
        # for(i in 1:dim(data)[2]){
        #   tmp <- data[,i]
        #   dat31 <- cbind(dat31,tmp)
        # }
        # dat31 <- as.data.frame(dat31)
        # 
        # colnames(dat31) <- colnames(data) # this line is for the returning original column names
        dataID <- data
        
        
        #++++++++++++++++++++ select a single trace ++++++++++++++#
        if (input$select == T) {
          dataID <- (dataID[input$trace, ])
          #++++++++++++++++++++ calculate deltaF/F based on trace minimal value++++++++++++++#
          if (input$deltaF == T) {
            #dataID <- (dataID[input$trace, ])
            minDataID <- min(dataID)
            
              dataID <- (dataID - minDataID)/minDataID
            
                                           }
          dataID
                                 }
        
        dataID$id = 1:nrow(dataID)
      
        dataMelt = melt(dataID, variable.name = "Frame", value.name = "Intensity", id.vars = "id")
      }
    #####################################################################################################
      #++++++++++++++++++++++++++++ build graph  either for deltaF or for Intensity ++++++++++++++++++++#  
      if (input$deltaF == T) {
        dataMelt = melt(dataID, variable.name = "Frame", value.name = "deltaF_over_F", id.vars = "id")
        
        p <-      ggplot(dataMelt, aes(x=Frame, y=deltaF_over_F, group=id, color=id)) +
          geom_textline(aes(label = id), hjust = 1) +
          theme_bw() +
          theme(legend.position = "none")
        return(p)
      }
      #++++++++++++++++++++++++++++ build graph  for Intensity ++++++++++++++++++++#  
      if (input$deltaF == F) {
      p <-      ggplot(dataMelt, aes(x=Frame, y=Intensity, group=id, color=id)) +
        geom_textline(aes(label = id), hjust = 1) +
        theme_bw() +
        theme(legend.position = "none")
      return(p)
                              }
      
      
    } , height = 800, width = 1200)
  })
  ############################################################################ 
  ############# preparing a table for downloading plot data ################
  tablePlotData <- reactive({
    # #+++++++++++++++++++++++++++ add busy indicator for download+++++++++++++++++++++#
    # show_modal_gif( src = "https://jeroen.github.io/images/banana.gif" ) 
    
    img3 = readImage(files = input$files$datapath)
    img3_F1 = readImage((files = input$files$datapath)[1]) # first frame only
    
    img3 <- resize(img3, 512, 512)
    img3_F1 <- resize(img3_F1, 512, 512)
    
    
    cellsSmooth = Image(dim = dim(img3_F1))
    
    sigma <- input$sigma
    size <- input$size
    
    cellsSmooth = filter2(
      img3_F1,
      filter = makeBrush(size = size, shape = "gaussian",
                         sigma=sigma)
    )
    
    disc = makeBrush(size, "disc")
    disc = disc / sum(disc)
    offset = input$offset
    rad <- input$rad
    rmCell <- numbersFromText(input$removeCell)
    
    nucThresh = (cellsSmooth - filter2( cellsSmooth, disc ) > offset)
    
    nucOpened = EBImage::opening(nucThresh, kern = makeBrush(rad, shape = "disc"))
    
    nucSeed = bwlabel(nucOpened)
    
    nucMask = cellsSmooth - filter2(cellsSmooth, disc) > 0
    nucMask = fillHull(nucMask)
    nuclei = propagate(cellsSmooth, nucSeed, mask = nucMask)
    #EBImage::display(nuclei,method = "raster")
    
    
    nucSegOnNuc  = paintObjects(nuclei, tgt = toRGB(img3_F1), col = "#ffff00")
    #EBImage::display(nucSegOnNuc,method = "raster")
    
    
    fr1 = computeFeatures(nuclei,     img3_F1, xname = "Frame_",  refnames = "c1") # this is used to determine how many ROI were detected in the first frame
    
    
    ########################################################################################
    ################# background correction ################################################
    if (input$subtract == T) {
      
      
      data <- data.frame(col1 = rep(NA, dim(fr1)[1]))
      
      
      for(i in 1:dim(img3)[3]) {                             # Head of for-loop
        
        
        new_col <- computeFeatures(nuclei,     img3[,,i], xname = "",   
                                   refnames = "fr_")                      # Creating new variable
        data[ , i] <- new_col[,12]                     # Adding new variable to data
        colnames(data)[i] <- paste0("Frame_", i)    # Renaming new variable
      }
      
      whiteImg <- matrix(1, dim(img3_F1), dim(img3_F1)) # create an array of 1s, which will be rendered white
      
      bkg <- whiteImg - nuclei #subtract cells masks from white image to get a mask of a background
      
      data_bkg <- data.frame(col1 = rep(NA, 1))
      
      
      
      for(i in 1:dim(img3)[3]) {                             # A for-loop to create a vector with background intensities timeline
        
        
        new_col <- computeFeatures(bkg,     img3[,,i], xname = "",   
                                   refnames = "fr_")                      # Creating new variable
        data_bkg[ , i] <- new_col[,12]                     # Adding new columns with next timeframe data to a data
        colnames(data_bkg)[i] <- paste0("Frame_", i)    # Renaming new variable
      }
      
      data_subtr <- sweep(as.matrix(data), MARGIN = 2, as.matrix(data_bkg)) # SWEEP function subtracts a vector (has to be converted 
      # to a matrix to work with sweep) with background from data DF (also as a matrix) row by row (this is why I use "MARGIN=2") 
      
      data <- as.data.frame(data_subtr) #convert matrix back into a DF
      
    }
    
    ###################### If background subtraction is NOT checked ###############################################
    
    if (input$subtract == F) {
      data <- data.frame(col1 = rep(NA, dim(fr1)[1]))       # Creating a new placeholder dataframe with number of rows calulated by: computeFeatures(nuclei,     img3_F1, xname = "Frame1",  refnames = "c1")
      
      
      for(i in 1:dim(img3)[3]) {                             # Head of for-loop
        
        # new dataframe with ROI calculations over entire image stack
        new_col <- computeFeatures(nuclei,     img3[,,i], xname = "P",   # xname for naming the first name of various parameters columns
                                   refnames = "fr_")                      # Creating new variable
        data[ , i] <- new_col[,12]                     # Adding column with intensities [12] to a new dataframe
        colnames(data)[i] <- paste0("Frame_", i)    # Renaming new variable
      }
      data
    }
    
    ###################### remove unwanted cells ###############################################
    
    
    if (input$rmvCell == T) {
      
      data <- data[-rmCell, ]
      data
    }
    
    
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++###
    
    if (input$detrend == FALSE && input$subtract == T) {
      
      dataID <- data
      #++++++++++++++++++++ select a single trace ++++++++++++++#
      if (input$select == T) {
        dataID <- (dataID[input$trace, ])
        #++++++++++++++++++++ calculate deltaF/F based on trace minimal value++++++++++++++#
        if (input$deltaF == T) {
          #dataID <- (dataID[input$trace, ])
          minDataID <- min(dataID)
          dataID <- (dataID - minDataID)/minDataID
        }
        
      }
      
      dataID$Cell_number = 1:nrow(dataID)
      
      ###### add last row to the DF with background
      
      data_bkg$bkg <- "background"
      dataID[dim(dataID)[1]+1, ] <- data_bkg[1, ]
  
      
      dataID
      
    }
    
    else if (input$detrend == FALSE && input$subtract == F) {
      
      dataID <- data
      
      #++++++++++++++++++++ select a single trace ++++++++++++++#
      if (input$select == T) {
        dataID <- (dataID[input$trace, ])
        #++++++++++++++++++++ calculate deltaF/F based on trace minimal value++++++++++++++#
        if (input$deltaF == T) {
          #dataID <- (dataID[input$trace, ])
          minDataID <- min(dataID)
          dataID <- (dataID - minDataID)/minDataID
        }
      }
      dataID$Cell_number = 1:nrow(dataID)
      
        dataID
      
    }
    
    else if (input$detrend == TRUE) {
      
      
      
      ###################### break points ############################################### 
      if (input$break_point == F && input$detrend == TRUE) {
        # bp <- input$bp
        # break.points <- seq(from=bp,to=dim(data)[2], by=bp)
        
        data.t <- t(data) # matrix transposition
        colnames(data.t) <- rownames(data) #restoration of column names
        data.t.df <- as.data.frame(data.t) # conversion to the list formate
        
        data.t.df.dt <- data.t.df
        
        #detrending without breaking points###
        for(i in 1:dim(data)[1]){                                       # run each Sample in a loop
          data.t.df.dt[,i] <- detrend(data.t.df[,i], tt = 'linear', bp=c())       # detrend the data
        }
        
        data.dt <- t(data.t.df.dt)
        data.dt.df <- as.data.frame(data.dt)
        dataID <- data.dt.df
        
      }
      
      else if (input$break_point == T && input$detrend == TRUE) {
        #bp <- input$bp
        break.points <- seq(from=input$bp,to=dim(data)[2], by=input$bp)
        
        data.t <- t(data) # matrix transposition
        colnames(data.t) <- rownames(data) #restoration of column names
        data.t.df <- as.data.frame(data.t) # conversion to the list formate
        
        data.t.df.dt <- data.t.df 
        
        #detrending with breaking points### 
        for(i in 1:dim(data)[1]){                                       # run each Sample in a loop
          data.t.df.dt[,i] <- detrend(data.t.df[,i], bp=break.points)       # detrend the data
        }
        
        data.dt <- t(data.t.df.dt)
        data.dt.df <- as.data.frame(data.dt)
        dataID <- data.dt.df
        dataID
      }
      
      #++++++++++++++++++++ select a single trace ++++++++++++++#
      if (input$select == T) {
        dataID <- (dataID[input$trace, ])
        #++++++++++++++++++++ calculate deltaF/F based on trace minimal value++++++++++++++#
        if (input$deltaF == T) {
          #dataID <- (dataID[input$trace, ])
          minDataID <- min(dataID)
          if (minDataID > 0) {
            dataID <- (dataID - minDataID)/minDataID
          }
          else if (minDataID < 0) {
            
            dataID = (dataID - minDataID)/minDataID
            dataID = -dataID
            
          }
        }
      }
      
      dataID$Cell_number = 1:nrow(dataID)
      
      ############### recalculate background for download
      whiteImg <- matrix(1, dim(img3_F1), dim(img3_F1)) # create an array of 1s, which will be rendered white
      
      bkg <- whiteImg - nuclei #subtract cells masks from white image to get a mask of a background
      
      
      
      data_bkg <- data.frame(col1 = rep(NA, 1))
      
      
      
      for(i in 1:dim(img3)[3]) {                             # A for-loop to create a vector with background intensities timeline
        
        
        new_col <- computeFeatures(bkg,     img3[,,i], xname = "",   
                                   refnames = "fr_")                      # Creating new variable
        data_bkg[ , i] <- new_col[,12]                     # Adding new columns with next timeframe data to a data
        colnames(data_bkg)[i] <- paste0("Frame_", i)    # Renaming new variable
      }
      
      ####### add last row to the DF with background
      data_bkg$bkg <- "background"
      dataID[dim(dataID)[1]+1, ] <- data_bkg[1, ]
      
   
      dataID
    }
    # #+++++++++++++++++++++++++++ remove busy indicator after the download+++++++++++++++++++++#
    # remove_modal_gif() 
    
  })
  
  #####################################################################################
  ############################## build correlation HEATMAP #############################
  observeEvent(input$buildHeat, {
    output$heat <- renderPlot({
      img3 = readImage(files = input$files$datapath)
      img3_F1 = readImage((files = input$files$datapath)[1]) # first frame only
      img3 <- resize(img3, 512, 512)
      img3_F1 <- resize(img3_F1, 512, 512)
      
      cellsSmooth = Image(dim = dim(img3_F1))
      
      sigma <- input$sigma
      size <- input$size
      
      size <- input$size
      cellsSmooth = filter2(
        img3_F1,
        filter = makeBrush(size = size, shape = "gaussian",
                           sigma=sigma)
      )
      
      disc = makeBrush(size, "disc")
      disc = disc / sum(disc)
      offset = input$offset
      rad <- input$rad
      rmCell <- numbersFromText(input$removeCell)
      
      nucThresh = (cellsSmooth - filter2( cellsSmooth, disc ) > offset)
      
      nucOpened = EBImage::opening(nucThresh, kern = makeBrush(rad, shape = "disc"))
      
      nucSeed = bwlabel(nucOpened)
      
      nucMask = cellsSmooth - filter2(cellsSmooth, disc) > 0
      nucMask = fillHull(nucMask)
      nuclei = propagate(cellsSmooth, nucSeed, mask = nucMask)
      #EBImage::display(nuclei,method = "raster")
      
      
      nucSegOnNuc  = paintObjects(nuclei, tgt = toRGB(img3_F1), col = "#ffff00")
      
      
      
      fr1 = computeFeatures(nuclei,     img3_F1, xname = "Frame1",  refnames = "c1") # this is used to determine how many ROI were detected in the first frame
      data <- data.frame(col1 = rep(NA, dim(fr1)[1]))
      
      for(i in 1:dim(img3)[3]) {                             # Head of for-loop
        
        
        new_col <- computeFeatures(nuclei,     img3[,,i], xname = "Frame_",
                                   refnames = "fr_")                      # Creating new variable
        data[ , i] <- new_col[,12]                     # Adding new variable to data
        colnames(data)[i] <- paste0("", i)    # Renaming new variable
      }
      
      
      
      data.c <- cbind(Cell = c(1:dim(fr1)[1]), data)
      
      if (input$rmvCell == T) {
        
        data.c <- data.c[-rmCell, ]
        data.c
      }
 
      
      dat1.t <- t(data.c)
      dat1.t.NoC <- dat1.t[-1,]
      colnames(dat1.t.NoC) <- c(1:dim(dat1.t.NoC)[2])
      data.t.c <- cbind(Time = c(1:dim(dat1.t.NoC)[1]), dat1.t.NoC)
      rownames(data.t.c) <- NULL
      dat1 <- data.t.c
      
      
      
      
      cl <- c("Time_Frame",paste0("Cell.0",1:9),paste0("Cell.",10:(dim(dat1)[2]-1)))
      colnames(dat1) <- cl
      dat1 <- as.data.frame.array(dat1)
      
      dat1.dt <- dat1                                                    #new dataframe for detrended results 
      
      for(i in 1:dim(dat1)[2]){                                          # run each Sample in a loop
        dat1.dt[,i] <- detrend(dat1[,i], tt = 'linear', bp=c())       # detrend the data using linear model
      }; dat1.dt[,1] <- dat1[,1] 
      

      #################################################################
      ### Find out the range of the data and set the minimum value to zero
      drange <- data.frame(row.names = cl, "min" = apply(dat1.dt,2,min), "max" = apply(dat1.dt,2,max))
      min.vals <- matrix(rep(drange$min, dim(dat1)[1]),  ncol=dim(dat1)[2], byrow=T)
      dat1.dt.z <- dat1.dt-min.vals
      drange.z <- data.frame(row.names = cl, "min" = apply(dat1.dt.z,2,min), "max" = apply(dat1.dt.z,2,max))
      #################################################################

      #################################################################
      ####          Build "pretty" heatmap            ####
      # cut the heatmap into segments (first identify by eye how many true sub-clusters you have
      cuts <- 3
      pheatmap(cor(dat1.dt.z[,-1]), margins=c(10,10), cutree_rows=cuts, cutree_cols=cuts)
      
      # dataRN <- data 
      # rownames(dataRN) <- c(paste0("Neuron_", 1:dim(fr1)[1])) # naming the rows in a new dataframe
      # 
      # tDataRN <- t(dataRN) # transpose matrix
      # 
      # p <- pheatmap(cor(tDataRN, use = "pairwise.complete.obs"), cutree_rows=3, cutree_cols=3, treeheight_row=0, treeheight_col=0)                     
      # return(p)
      
      
      
    } , height = 800, width = 800)
  })
  
  imagesUploaded <- reactive({
    
    readImage(files = input$files$datapath)
  })
  
  
  imgUploadedFrame1 <- reactive({
    
    img3_F1 <- readImage((files = input$files$datapath)[1])
    
    img3_F1 <- resize(img3_F1, 512, 512)
    img3_F1
  })
  
  
  imgUploadedMask <- reactive({
    
    img3 = readImage(files = input$files$datapath)
    img3_F1 = readImage((files = input$files$datapath)[1]) # first frame only
    img3 <- resize(img3, 512, 512)
    img3_F1 <- resize(img3_F1, 512, 512)
    cellsSmooth = Image(dim = dim(img3_F1))
    
    sigma <- input$sigma
    size <- input$size
    
    cellsSmooth = filter2(
      img3_F1,
      filter = makeBrush(size = size, shape = "gaussian",
                         sigma=sigma)
    )
    
    disc = makeBrush(size, "disc")
    disc = disc / sum(disc)
    offset = input$offset
    rad <- input$rad
    
    nucThresh = (cellsSmooth - filter2(cellsSmooth, disc ) > offset)
    
    nucOpened = EBImage::opening(nucThresh, kern = makeBrush(rad, shape = "disc"))
    
    nucSeed = bwlabel(nucOpened)
    
    nucMask = cellsSmooth - filter2(cellsSmooth, disc) > 0
    nucMask = fillHull(nucMask)
    nuclei = propagate(cellsSmooth, nucSeed, mask = nucMask)
    #EBImage::display(nuclei,method = "raster")
    
    
    nucSegOnNuc  = paintObjects(nuclei, tgt = toRGB(img3_F1), col = "#ffff00")
    
    #fts = computeFeatures.moment(nuclei)
    
    f <- nucSegOnNuc
    f
  })
  
  tableInput <- reactive({
    
    #+++++++++++++++++++++++++++ add busy indicator for download+++++++++++++++++++++#
    show_modal_gif( src = "https://jeroen.github.io/images/banana.gif" ) 
    
    img3 = readImage(files = input$files$datapath)
    img3_F1 = readImage((files = input$files$datapath)[1]) # first frame only
    img3 <- resize(img3, 512, 512)
    img3_F1 <- resize(img3_F1, 512, 512)
    cellsSmooth = Image(dim = dim(img3_F1))
    
    sigma <- input$sigma
    size <- input$size
    
    cellsSmooth = filter2(
      img3_F1,
      filter = makeBrush(size = size, shape = "gaussian",
                         sigma=sigma)
    )
    
    disc = makeBrush(size, "disc")
    disc = disc / sum(disc)
    offset = input$offset
    rad <- input$rad
    rmCell <- numbersFromText(input$removeCell)
    
    nucThresh = (cellsSmooth - filter2(cellsSmooth, disc ) > offset)
    
    nucOpened = EBImage::opening(nucThresh, kern = makeBrush(rad, shape = "disc"))
    
    nucSeed = bwlabel(nucOpened)
    
    nucMask = cellsSmooth - filter2(cellsSmooth, disc) > 0
    nucMask = fillHull(nucMask)
    nuclei = propagate(cellsSmooth, nucSeed, mask = nucMask)
    
    
    
    #nucSegOnNuc  = paintObjects(nuclei, tgt = toRGB(img3_F1), col = "#ffff00")
    
    #fts = computeFeatures.moment(nuclei)
    
    
    fr1 = computeFeatures(nuclei,     img3_F1, xname = "Frame1",  refnames = "c1") # this is used to determine how many ROI were detected in the first frame
    data <- data.frame(col1 = rep(NA, dim(fr1)[1]))
    
    for(i in 1:dim(img3)[3]) {                             # Head of for-loop
      
      
      new_col <- computeFeatures(nuclei,     img3[,,i], xname = "Frame_",
                                 refnames = "fr_")                      # Creating new variable
      data[ , i] <- new_col[,12]                     # Adding new variable to data
      colnames(data)[i] <- paste0("", i)    # Renaming new variable
    }
    
    
    data.c <- cbind(Cell = c(1:dim(fr1)[1]), data)
    
    
    if (input$rmvCell == T) {
      
      data.c <- data.c[-rmCell, ]
      data.c
    }
    
    dat1.t <- t(data.c)
    dat1.t.NoC <- dat1.t[-1,]
    colnames(dat1.t.NoC) <- c(1:dim(dat1.t.NoC)[2])
    data.t.c <- cbind(Time = c(1:dim(dat1.t.NoC)[1]), dat1.t.NoC)
    rownames(data.t.c) <- NULL
    dat1 <- data.t.c
    
    cl <- c("Time_Frame",paste0("Cell.0",1:9),paste0("Cell.",10:(dim(dat1)[2]-1)))
    colnames(dat1) <- cl
    dat1 <- as.data.frame.array(dat1)
    
    dat1.dt <- dat1                                                    #new dataframe for detrended ersults
    
    for(i in 1:dim(dat1)[2]){                                          # run each Sample in a loop
      dat1.dt[,i] <- detrend(dat1[,i], tt = 'linear', bp = c())       # detrend the data using linear model
    }; dat1.dt[,1] <- dat1[,1]
    
    
    ### Find out the range of the data and set the minimum value to zero
    drange <- data.frame(row.names = cl, "min" = apply(dat1.dt,2,min), "max" = apply(dat1.dt,2,max))
    min.vals <- matrix(rep(drange$min, dim(dat1)[1]),  ncol=dim(dat1)[2], byrow=T)
    dat1.dt.z <- dat1.dt-min.vals
    drange.z <- data.frame(row.names = cl, "min" = apply(dat1.dt.z,2,min), "max" = apply(dat1.dt.z,2,max))
    #################################################################
   
    ####################################################################
    ################ correlation table summary generation ##
    
    pairwise.combn <- t(combn(2:dim(dat1)[2],2))
    
    p.val <- c()
    rq <- c()
    cors <- c()
    for(i in 1:dim(pairwise.combn)[1]){
      lm1 <- summary(lm(dat1.dt.z[,pairwise.combn[i,1]] ~ dat1.dt.z[,pairwise.combn[i,2]] ))
      p.val <- c(p.val,lm1$coefficients[2,4])
      rq <- c(rq,lm1$r.squared)
      cors <- c(cors,cor(dat1.dt.z[,pairwise.combn[i,1]] , dat1.dt.z[,pairwise.combn[i,2]] ))
    }
    
    comb <- data.frame("X" = colnames(dat1)[pairwise.combn[,1]],
                       "Y" = colnames(dat1)[pairwise.combn[,2]],
                       "correlation" = cors,
                       "R-squared" = rq,
                       "P.value" = p.val)
    comb$Adj.P.Value <- p.adjust(comb$P.value)
    
    cat("Number of correlating pairs with adj. p-value<", input$adjPval, "& correlation>", input$corr, "::",
        length(comb[comb$Adj.P.Value<input$adjPval&comb$correlation>input$corr, 1]), "\n",
        "Percentage of correlating pairs, of all possible pairs:", 
        length(comb[comb$Adj.P.Value<input$adjPval&comb$correlation>input$corr, 1])*100/length(comb$X), "%", "\n")
    
    
    cat("Correlation table filtered by adj. p-value<", input$adjPval, "& correlation>", input$corr, ":", "\n")
    print(comb[comb$Adj.P.Value<input$adjPval&comb$correlation>input$corr,]  )
    
    #comb[comb$Adj.P.Value<0.01&comb$correlation>0.8,]
    
  })
  
  observeEvent(input$buildTable, {
    output$table = renderPrint({                      # create output for displaying a correlation table summary in UI
      tableInput <- reactive({
        img3 = readImage(files = input$files$datapath)
        img3_F1 = readImage((files = input$files$datapath)[1]) # first frame only
        img3 <- resize(img3, 512, 512)
        img3_F1 <- resize(img3_F1, 512, 512)
        cellsSmooth = Image(dim = dim(img3_F1))
        
        sigma <- input$sigma
        size <- input$size
        
        cellsSmooth = filter2(
          img3_F1,
          filter = makeBrush(size = size, shape = "gaussian",
                             sigma=sigma)
        )
        
        disc = makeBrush(size, "disc")
        disc = disc / sum(disc)
        offset = input$offset
        rad <- input$rad
        rmCell <- numbersFromText(input$removeCell)
        
        nucThresh = (cellsSmooth - filter2(cellsSmooth, disc ) > offset)
        
        nucOpened = EBImage::opening(nucThresh, kern = makeBrush(rad, shape = "disc"))
        
        nucSeed = bwlabel(nucOpened)
        
        nucMask = cellsSmooth - filter2(cellsSmooth, disc) > 0
        nucMask = fillHull(nucMask)
        nuclei = propagate(cellsSmooth, nucSeed, mask = nucMask)
        
        
        
        #nucSegOnNuc  = paintObjects(nuclei, tgt = toRGB(img3_F1), col = "#ffff00")
        
        #fts = computeFeatures.moment(nuclei)
        
        
        fr1 = computeFeatures(nuclei,     img3_F1, xname = "Frame1",  refnames = "c1") # this is used to determine how many ROI were detected in the first frame
        data <- data.frame(col1 = rep(NA, dim(fr1)[1]))
        
        for(i in 1:dim(img3)[3]) {                             # Head of for-loop
          
          
          new_col <- computeFeatures(nuclei,     img3[,,i], xname = "Frame_",
                                     refnames = "fr_")                      # Creating new variable
          data[ , i] <- new_col[,12]                     # Adding new variable to data
          colnames(data)[i] <- paste0("", i)    # Renaming new variable
        }
        
        
        data.c <- cbind(Cell = c(1:dim(fr1)[1]), data)
        
        
        if (input$rmvCell == T) {
          
          data.c <- data.c[-rmCell, ]
          data.c
        }
        
        dat1.t <- t(data.c)
        dat1.t.NoC <- dat1.t[-1,]
        colnames(dat1.t.NoC) <- c(1:dim(dat1.t.NoC)[2])
        data.t.c <- cbind(Time = c(1:dim(dat1.t.NoC)[1]), dat1.t.NoC)
        rownames(data.t.c) <- NULL
        dat1 <- data.t.c
        
        cl <- c("Time_Frame",paste0("Cell.0",1:9),paste0("Cell.",10:(dim(dat1)[2]-1)))
        colnames(dat1) <- cl
        dat1 <- as.data.frame.array(dat1)
        
        dat1.dt <- dat1                                                    #new dataframe for detrended ersults
        
        for(i in 1:dim(dat1)[2]){                                          # run each Sample in a loop
          dat1.dt[,i] <- detrend(dat1[,i], tt = 'linear', bp = c())       # detrend the data using linear model
        }; dat1.dt[,1] <- dat1[,1]
        
        
        ### Find out the range of the data and set the minimum value to zero
        drange <- data.frame(row.names = cl, "min" = apply(dat1.dt,2,min), "max" = apply(dat1.dt,2,max))
        min.vals <- matrix(rep(drange$min, dim(dat1)[1]),  ncol=dim(dat1)[2], byrow=T)
        dat1.dt.z <- dat1.dt-min.vals
        drange.z <- data.frame(row.names = cl, "min" = apply(dat1.dt.z,2,min), "max" = apply(dat1.dt.z,2,max))
        #################################################################
        
        
        
        
        ####################################################################
        ################ correlation table summary generation ##
        
        pairwise.combn <- t(combn(2:dim(dat1)[2],2))
        
        p.val <- c()
        rq <- c()
        cors <- c()
        for(i in 1:dim(pairwise.combn)[1]){
          lm1 <- summary(lm(dat1.dt.z[,pairwise.combn[i,1]] ~ dat1.dt.z[,pairwise.combn[i,2]] ))
          p.val <- c(p.val,lm1$coefficients[2,4])
          rq <- c(rq,lm1$r.squared)
          cors <- c(cors,cor(dat1.dt.z[,pairwise.combn[i,1]] , dat1.dt.z[,pairwise.combn[i,2]] ))
        }
        
        comb <- data.frame("X" = colnames(dat1)[pairwise.combn[,1]],
                           "Y" = colnames(dat1)[pairwise.combn[,2]],
                           "correlation" = cors,
                           "R-squared" = rq,
                           "P.value" = p.val)
        comb$Adj.P.Value <- p.adjust(comb$P.value)
        
        cat("Number of correlating pairs with adj. p-value<", input$adjPval, "& correlation>", input$corr, "::",
            length(comb[comb$Adj.P.Value<input$adjPval&comb$correlation>input$corr, 1]), "\n",
            "Percentage of correlating pairs, of all possible pairs:", 
            length(comb[comb$Adj.P.Value<input$adjPval&comb$correlation>input$corr, 1])*100/length(comb$X), "%", "\n")
        
        
        cat("Correlation table filtered by adj. p-value<", input$adjPval, "& correlation>", input$corr, ":", "\n")
        print(comb[comb$Adj.P.Value<input$adjPval&comb$correlation>input$corr,]  )
        
        #comb[comb$Adj.P.Value<0.01&comb$correlation>0.8,]
        
      })
      
      tableInput()
    })
    
   
  })
  
  
  
  output$downloadPlotData <- downloadHandler(
    
    
    filename = 'Calcium_traces.csv',
    
    content = function(file) {
      
      #+++++++++++++++++++++++++++ add busy indicator for download+++++++++++++++++++++#
      show_modal_gif( src = "https://jeroen.github.io/images/banana.gif" )
      
      write.csv(tablePlotData(), file, row.names = FALSE)
      
      #+++++++++++++++++++++++++++ remove busy indicator after the download+++++++++++++++++++++#
      remove_modal_gif() 
      #dev.off()
    }              )
  
  
  output$downloadTable <- downloadHandler(
    
    filename = 'Correlaition_Table.csv',
    
    content = function(file) {
      #+++++++++++++++++++++++++++ add busy indicator for download+++++++++++++++++++++#
      show_modal_gif( src = "https://jeroen.github.io/images/banana.gif" )
      write.csv(tableInput(), file, row.names = FALSE)
      
      #+++++++++++++++++++++++++++ remove busy indicator after the download+++++++++++++++++++++#
      remove_modal_gif() 
      #dev.off()
    }              
    
    
    
    
    )
  
  output$widget <- renderDisplay({
    display(imagesUploaded(), method = 'browser')
  })
  
  output$rasterUploaded1 <- renderPlot({
    plot(imgUploadedFrame1())
  })
  
  observeEvent (input$addLabels==TRUE, {
    
    
    if (input$addLabels==TRUE) {
      output$rasterMask <- renderPlot({
        plot(imgUploadedMask())
        img3 = readImage(files = input$files$datapath)
        img3_F1 = readImage((files = input$files$datapath)[1]) # first frame only
        img3 <- resize(img3, 512, 512)
        img3_F1 <- resize(img3_F1, 512, 512)
        cellsSmooth = Image(dim = dim(img3_F1))
        
        sigma <- input$sigma
        size <- input$size
        
        cellsSmooth = filter2(
          img3_F1,
          filter = makeBrush(size = size, shape = "gaussian",
                             sigma=sigma)
        )
        
        disc = makeBrush(size, "disc")
        disc = disc / sum(disc)
        offset = input$offset
        rad <- input$rad
        rmCell <- numbersFromText(input$removeCell)
        
        nucThresh = (cellsSmooth - filter2(cellsSmooth, disc ) > offset)
        
        nucOpened = EBImage::opening(nucThresh, kern = makeBrush(rad, shape = "disc"))
        
        nucSeed = bwlabel(nucOpened)
        
        nucMask = cellsSmooth - filter2(cellsSmooth, disc) > 0
        nucMask = fillHull(nucMask)
        nuclei = propagate(cellsSmooth, nucSeed, mask = nucMask)
        #EBImage::display(nuclei,method = "raster")
        
        nucSegOnNuc  = paintObjects(nuclei, tgt = toRGB(img3_F1), col = "#ffff00")
        
        #EBImage::display(nucSegOnNuc,method = "raster")
        
        fts = computeFeatures.moment(nuclei)
        
        #if(is.na(rmCell)) return(NULL) 
        if (input$rmvCell == T) {
          
          fts <- fts[-rmCell, ]
          fts
        }
        
        text(fts[,"m.cx"], fts[,"m.cy"], labels=seq_len(nrow(fts)), col="red", cex=1)
        
        
      }, height = 512, width = 512)
      
      
    }
    
    else if (input$addLabels==FALSE) {
      
      output$rasterMask <- renderPlot({
        req(input$files$datapath)
        plot(imgUploadedFrame1())
        
        
      }, height = 512, width = 512)
    }
  })
  
  output$raster <- renderPlot({
    plot(imgUploadedFrame1())
  })
  
} 
shinyApp(ui = ui, server = server)
