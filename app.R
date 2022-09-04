library(shiny)
library("EBImage")
library(pheatmap)
library("magrittr")
library("tibble")
library("ggplot2")
library("genefilter")
library("GGally")
library(geomtextpath)
library(reshape2)
library(pracma)
#library(dplyr)
library('shinythemes')
library('shinycssloaders')
library(BiocManager)
options(repos = BiocManager::repositories())
knitr::opts_chunk$set(echo = TRUE)

ui <- fluidPage(theme = shinytheme("paper"),
                titlePanel(h6("Calcium Imaging v.0.0.1")),
                
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
        
        
        
    
      )
    ),
    mainPanel(
      
      h4(tags$b(span("Plot", style = "color:darkblue"))),
      tabsetPanel(
        tabPanel("Plot",
                 br(),
                 actionButton("buildPlot", "Build"), 
                 checkboxInput(inputId = "detrend",
                               label = "Detrend my data",
                               value = FALSE),
                 br(),
                 plotOutput("plot")
                 
                 
                 ),
        tabPanel("Heatmap", 
                 br(),
                 actionButton("buildHeat", "Heatmap"),
                 br(),
                 plotOutput("heat"),
                 ),
        tabPanel("Mask", 
                 br(),
                 actionButton("addLabels", "Show labels"),
                 br(),
                 plotOutput("rasterMask"), 
                 ),
        
        
        tabPanel("Interactive image browser", displayOutput("widget"))
      
    )
  )
)
)

server <- function(input, output) {
  
  output$files <- renderTable(input$files)
  
  files <- reactive({
    files <- input$files
    files$datapath <- gsub("\\\\", "/", files$datapath)
    files
    
    
  })
  
  
  output$images <- renderUI({
    if(is.null(input$files)) return(NULL)
    image_output_list <- 
      lapply(1:nrow(files()),
             function(i)
             {
               imagename = paste0("image", i)
               imageOutput(imagename)
             })
    
    do.call(tagList, image_output_list)
  })
  
  observe({
    if(is.null(input$files)) return(NULL)
    for (i in 1:nrow(files()))
    {
      print(i)
      local({
        my_i <- i
        imagename = paste0("image", my_i)
        print(imagename)
        output[[imagename]] <- 
          renderImage({
            list(src = files()$datapath[my_i],
                 alt = "Image failed to render")
          }, deleteFile = FALSE)
        #img3 = readImage(files()$datapath[my_i])
       #
        
      })
      
      
    }
    
    
  })
  
   
    
  
   
  
  observeEvent(input$buildPlot, {
            output$plot <- renderPlot({
              img3 = readImage(files = input$files$datapath)
              
              # w = makeBrush(size = 101, shape = "gaussian", sigma = 7)
              # tibble(w = w[(nrow(w)+1)/2, ])
              # #%>%
              #   #ggplot(aes(y = w, x = seq(along = w))) + geom_point()
              # 
              # nucSmooth = filter2(getFrame(img3, 1), w)
              
              cellsSmooth = Image(dim = dim(img3))
              
              sigma <- rep (7, dim(img3)[3])
              for(i in seq_along(sigma)) {
                cellsSmooth[,,i] = filter2(
                  img3[,,i],
                  filter = makeBrush(size = 101, shape = "gaussian",
                                     sigma = sigma[i])
                )
              }
              
              disc = makeBrush(101, "disc")
              disc = disc / sum(disc)
              offset = 0.07
              nucThresh = (cellsSmooth[,,1] - filter2( cellsSmooth[,,1], disc ) > offset)
              
              nucOpened = EBImage::opening(nucThresh, kern = makeBrush(3, shape = "disc"))
              
              nucSeed = bwlabel(nucOpened)
              
              nucMask = cellsSmooth[,, 1] - filter2(cellsSmooth[ , , 1], disc) > 0
              nucMask = fillHull(nucMask)
              nuclei = propagate(cellsSmooth[,,1], nucSeed, mask = nucMask)
            
              
              nucSegOnNuc  = paintObjects(nuclei, tgt = toRGB(img3[, , 1]), col = "#ffff00")
              
            
              

              fr1 = computeFeatures(nuclei,     img3[,,1], xname = "Frame1",  refnames = "c1") # this is used to determine how many ROI were detected in the first frame
              data <- data.frame(col1 = rep(NA, dim(fr1)[1]))

              for(i in 1:dim(img3)[3]) {                             # Head of for-loop


                new_col <- computeFeatures(nuclei,     img3[,,i], xname = "Frame_",
                                           refnames = "fr_")                      # Creating new variable
                data[ , i] <- new_col[,12]                     # Adding new variable to data
                colnames(data)[i] <- paste0("", i)    # Renaming new variable
              }
              
              ##++++++++++++++++++++++++++++ Detrending ++++++++++++++++++++++++++++++++++++++++###
              if (input$detrend == TRUE) {
              
              #cn <- colnames(data) # for later use in the returning original column names
              
              #library(pracma)
              
              # br <- 4                                                        # You can try different break points
              # break.points <- seq(from=br,to=dim(data)[1], by=br)
              # data.dt <- data                                                 # create a duplicate data frame
              
                dat31 <- c()
                for(i in 1:dim(data)[2]){
                  tmp <- detrend(data[,i], tt = 'linear', bp = c()) # fits a "moving" linear regression to the data and subracts the "time" component from the total intensities
                  dat31 <- cbind(dat31,tmp)
                }
                dat31 <- as.data.frame(dat31)
                
                cn <- colnames(data) # for later use in the returning original column names
                colnames(dat31) <- cn
                dataID <- dat31
              
              # dataID <- cbind(cellID = c(1:dim(fr1)[1]), dataID)
              # 
              # dataID <- dataID[,-1]
              dataID$id = 1:nrow(dataID)
              dataMelt = melt(dataID, variable.name = "Frame", value.name = "Intensity", id.vars = "id")
                                            }
              #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++###
              
              else if (input$detrend == FALSE) {
              dataID <- data

              dataID <- cbind(cellID = c(1:dim(fr1)[1]), dataID)

              dataID <- dataID[,-1]
              dataID$id = 1:nrow(dataID)

              #require(reshape2)
              dataMelt = melt(dataID, variable.name = "Frame", value.name = "Intensity", id.vars = "id")
              }
            
              
              
              
                p <-      ggplot(dataMelt, aes(x=Frame, y=Intensity, group=id, color=id)) +
                geom_textline(aes(label = id), hjust = 1) +
                theme_bw() +
                theme(legend.position = "none")
                return(p)
              
               
      
    } , height = 800, width = 1200)
  })
  
  observeEvent(input$buildHeat, {
    output$heat <- renderPlot({
      img3 = readImage(files = input$files$datapath)
      
      # w = makeBrush(size = 101, shape = "gaussian", sigma = 7)
      # tibble(w = w[(nrow(w)+1)/2, ])
      # %>%
      # ggplot(aes(y = w, x = seq(along = w))) + geom_point()
      
      #nucSmooth = filter2(getFrame(img3, 1), w)
      
      cellsSmooth = Image(dim = dim(img3))
      
      sigma <- rep (7, dim(img3)[3])
      for(i in seq_along(sigma)) {
        cellsSmooth[,,i] = filter2(
          img3[,,i],
          filter = makeBrush(size = 101, shape = "gaussian",
                             sigma = sigma[i])
        )
      }
      
      disc = makeBrush(101, "disc")
      disc = disc / sum(disc)
      offset = 0.07
      nucThresh = (cellsSmooth[,,1] - filter2( cellsSmooth[,,1], disc ) > offset)
      
      nucOpened = EBImage::opening(nucThresh, kern = makeBrush(3, shape = "disc"))
      
      nucSeed = bwlabel(nucOpened)
      
      nucMask = cellsSmooth[,, 1] - filter2(cellsSmooth[ , , 1], disc) > 0
      nucMask = fillHull(nucMask)
      nuclei = propagate(cellsSmooth[,,1], nucSeed, mask = nucMask)
      #EBImage::display(nuclei,method = "raster")
      
     
      
      nucSegOnNuc  = paintObjects(nuclei, tgt = toRGB(img3[, , 1]), col = "#ffff00")
      
      
      
      
      fr1 = computeFeatures(nuclei,     img3[,,1], xname = "Frame1",  refnames = "c1") # this is used to determine how many ROI were detected in the first frame
      data <- data.frame(col1 = rep(NA, dim(fr1)[1]))
      
      for(i in 1:dim(img3)[3]) {                             # Head of for-loop
        
        
        new_col <- computeFeatures(nuclei,     img3[,,i], xname = "Frame_",
                                   refnames = "fr_")                      # Creating new variable
        data[ , i] <- new_col[,12]                     # Adding new variable to data
        colnames(data)[i] <- paste0("", i)    # Renaming new variable
      }
      
      
      
      dataRN <- data 
      rownames(dataRN) <- c(paste0("Neuron_", 1:dim(fr1)[1])) # naming the rows in a new dataframe
      
      tDataRN <- t(dataRN) # transpose matrix
      
      p <- pheatmap(cor(tDataRN, use = "pairwise.complete.obs"), cutree_rows=3, cutree_cols=3, treeheight_row=0, treeheight_col=0)                     
      return(p)
      
      
      
    } , height = 800, width = 1200)
  })
  
  imgUploadedFrame1 <- reactive({
    
    readImage(files = input$files$datapath)
  })
  
  
  imgUploadedMask <- reactive({
    
    img3 = readImage(files = input$files$datapath)
    
   
    
       
    cellsSmooth = Image(dim = dim(img3))
    
    sigma <- rep (7, dim(img3)[3])
    for(i in seq_along(sigma)) {
      cellsSmooth[,,i] = filter2(
        img3[,,i],
        filter = makeBrush(size = 101, shape = "gaussian",
                           sigma = sigma[i])
      )
    }
    
    disc = makeBrush(101, "disc")
    disc = disc / sum(disc)
    offset = 0.07
    nucThresh = (cellsSmooth[,,1] - filter2( cellsSmooth[,,1], disc ) > offset)
    
    nucOpened = EBImage::opening(nucThresh, kern = makeBrush(3, shape = "disc"))
    
    nucSeed = bwlabel(nucOpened)
    
    nucMask = cellsSmooth[,, 1] - filter2(cellsSmooth[ , , 1], disc) > 0
    nucMask = fillHull(nucMask)
    nuclei = propagate(cellsSmooth[,,1], nucSeed, mask = nucMask)
  
    
    nucSegOnNuc  = paintObjects(nuclei, tgt = toRGB(img3[, , 1]), col = "#ffff00")
    
   
    
    f <- nucSegOnNuc
    f
  })
  
  
  
  output$widget <- renderDisplay({
    display(imgUploadedFrame1())
  })
  
  output$rasterUploaded1 <- renderPlot({
    plot(imgUploadedFrame1())
  })
  
  observe ({
  if (input$addLabels==TRUE) {
  output$rasterMask <- renderPlot({
    plot(imgUploadedMask())
    fts = computeFeatures.moment(nuclei)
   
     text(fts[,"m.cx"], fts[,"m.cy"], labels=seq_len(nrow(fts)), col="white", cex=.8)
                              
  })
  }
  
  else if (input$addLabels==FALSE) {
    output$rasterMask <- renderPlot({
      plot(imgUploadedFrame1())
      
      
    })
  }
  })
    
  output$raster <- renderPlot({
    plot(imgUploadedFrame1())
  })
  
} 
shinyApp(ui = ui, server = server)
