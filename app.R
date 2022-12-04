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
library(BiocManager)
options(repos = BiocManager::repositories())
knitr::opts_chunk$set(echo = TRUE)

ui <- fluidPage(theme = shinytheme("paper"),
                titlePanel(h6("Calcium Imaging v.0.1.2")),
                
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
                        value=3,
                        min = 1,
                        max = 111,
                        step = 2
                      ),
                      
                      
                      numericInput(
                        inputId="size",
                        label="Size (px)",
                        value=51,
                        min = 1,
                        max = 1011,
                        step = 2
                      ),
                      
                      sliderInput("offset", h6("offset"), 0, 1, 0.05),
                      
                      numericInput(
                        inputId="rad",
                        label="masking radius",
                        value=3,
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
                               actionButton("buildPlot", "Build"), 
                               checkboxInput(inputId = "detrend",
                                             label = "Detrend my data",
                                             value = FALSE),
                               br(),
                               plotOutput("plot", height = "100%")
                               
                               
                      ),
                      tabPanel("Correlations", 
                               br(),
                               h4(tags$b("Build heatmap", style = "color:darkblue")),
                               actionButton("buildHeat", "Build"),
                               br(),
                               plotOutput("heat", height = "100%"),
                               
                               br(),
                               h4(tags$b("Correlation Table", style = "color:darkblue")),
                               sliderInput("adjPval", h6("Adjusted p-value threshold"), 0, 0.05, 0.01),
                               sliderInput("corr", h6("Correlation coefficient threshold"), 0, 1, 0.8),
                               actionButton("buildTable", "Calculate correlation table"),
                               downloadButton("downloadTable", "Download"),
                               verbatimTextOutput('table'),
                               
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
       
      })
      
    }
    
  })
  
  numbersFromText <- function(text) {        # this functipon extracts numbers from text
    text <- gsub(" ", "", text)
    split <- strsplit(text, ",", fixed = FALSE)[[1]]
    as.numeric(split)
  }
  
  # exists_as_number <- function(item) {
  #   !is.null(item) && !is.na(item)
  # }
  # 
  
  
  
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
      
      
      data <- data.frame(col1 = rep(NA, dim(fr1)[1]))
      
      
      for(i in 1:dim(img3)[3]) {                             # Head of for-loop
        
        
        new_col <- computeFeatures(nuclei,     img3[,,i], xname = "Frame_",
                                   refnames = "fr_")                      # Creating new variable
        data[ , i] <- new_col[,12]                     # Adding new variable to data
        colnames(data)[i] <- paste0("", i)    # Renaming new variable
      }
      
      
      if (input$rmvCell == T) {

        data <- data[-rmCell, ]
        data
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
        
        cn <- colnames(data) # for use in the next line for the returning original column names
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
        
      
        
        dat31 <- c()
        for(i in 1:dim(data)[2]){
          tmp <- data[,i] # fits a "moving" linear regression to the data and subracts the "time" component from the total intensities
          dat31 <- cbind(dat31,tmp)
        }
        dat31 <- as.data.frame(dat31)
        
        cn <- colnames(data) # for use in the next line for the returning original column names
        colnames(dat31) <- cn
        dataID <- dat31
        
        # dataID <- cbind(cellID = c(1:dim(fr1)[1]), dataID)
        # 
        # dataID <- dataID[,-1]
        dataID$id = 1:nrow(dataID)
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
      
      dat1.dt <- dat1                                                    #new dataframe for detrended ersults 
      
      for(i in 1:dim(dat1)[2]){                                          # run each Sample in a loop
        dat1.dt[,i] <- detrend(dat1[,i], tt = 'linear', bp = c())       # detrend the data using linear model
      }; dat1.dt[,1] <- dat1[,1] 
      
      
      #################################################################
      ### Find out the range of the data and set the minimum value to zero
      drange <- data.frame(row.names = cl, "min" = apply(dat1.dt,2,min), "max" = apply(dat1.dt,2,max))
      min.vals <- matrix(rep(drange$min, dim(dat1)[1]),  ncol=dim(dat1)[2], byrow=T)
      dat1.dt.z <- dat1.dt-min.vals
      drange.z <- data.frame(row.names = cl, "min" = apply(dat1.dt.z,2,min), "max" = apply(dat1.dt.z,2,max))
      #################################################################
      
      
      #################################################################
      ####           detect peaks & smooth out the data            ####
      #
      # smoothing is needed because the raw data is very "spiky", this makes peak detection hard
      # peak is defined as high intensity: normalized intensity > 2x sd &
      # peak has >2 successive measurements increase before peak maximum &
      # peak has >2 successive measurements decrease after peak maximum
      smoothing.parameter <- 0.02                                     # Smooting parameter for lowess (larger values make smoothing rougher)
      successive.points <- 3                                          # how many successive points need to increase or decrease
      
      dat1.peaks <- list()                                            # Peaks are stores in a list
      dat1.dt.z.s <- dat1.dt.z                                        # The smoothed values go here
      
      for(i in 2:ncol(dat1.dt)){
        tmp.data <- dat1.dt.z[,i]
        noiselevel <- 2*sd(tmp.data)
        tmp.data <- lowess(tmp.data, f=smoothing.parameter)$y
        dat1.dt.z.s[,i] <- tmp.data
        tmp <- findpeaks(tmp.data, nups=successive.points, ndowns=successive.points, minpeakheight=noiselevel )
        if(is.null(tmp[,1])) { tmp <- matrix(NA,ncol=4) }
        else { tmp <- tmp }
        tmp <- data.frame("Intensity" = tmp[,1], "Peak.max"=tmp[,2], "Peak.start"=tmp[,4], "Peak.end"=tmp[,4], "Noiselevel" = noiselevel)
        dat1.peaks[[i-1]] <- tmp
      }; names(dat1.peaks) <- colnames(dat1.dt)[-1]
      
      
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
    
   
    
    nucSegOnNuc  = paintObjects(nuclei, tgt = toRGB(img3_F1), col = "#ffff00")
    
    fts = computeFeatures.moment(nuclei)
  
    
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
    
    
    #################################################################
    ####           detect peaks & smooth out the data            ####
    #
    # smoothing is needed because the raw data is very "spiky", this makes peak detection hard
    # peak is defined as high intensity: normalized intensity > 2x sd &
    # peak has >2 successive measurements increase before peak maximum &
    # peak has >2 successive measurements decrease after peak maximum
    smoothing.parameter <- 0.02                                     # Smooting parameter for lowess (larger values make smoothing rougher)
    successive.points <- 3                                          # how many successive points need to increase or decrease
    
    dat1.peaks <- list()                                            # Peaks are stores in a list
    dat1.dt.z.s <- dat1.dt.z                                        # The smoothed values go here
    
    for(i in 2:ncol(dat1.dt)){
      tmp.data <- dat1.dt.z[,i]
      noiselevel <- 2*sd(tmp.data)
      tmp.data <- lowess(tmp.data, f=smoothing.parameter)$y
      dat1.dt.z.s[,i] <- tmp.data
      tmp <- findpeaks(tmp.data, nups=successive.points, ndowns=successive.points, minpeakheight=noiselevel )
      if(is.null(tmp[,1])) { tmp <- matrix(NA,ncol=4) }
      else { tmp <- tmp }
      tmp <- data.frame("Intensity" = tmp[,1], "Peak.max"=tmp[,2], "Peak.start"=tmp[,4], "Peak.end"=tmp[,4], "Noiselevel" = noiselevel)
      dat1.peaks[[i-1]] <- tmp
    }; names(dat1.peaks) <- colnames(dat1.dt)[-1]
    
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
      tableInput()
      
    })
    
  })
  
  
  output$downloadTable <- downloadHandler(
    
    
    filename = 'Correlaition_Table.csv',
    
    content = function(file) {
      write.csv(tableInput(), file, row.names = FALSE)
      #dev.off()
    }              )
  
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
