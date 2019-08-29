library(shinydashboard)
library(NMF)
library(ggplot2)
library(fields)
library(png)
library(shinyBS)
nmf.options(grid.patch=TRUE)
#source("/Users/eqiao14/Desktop/Visualization/helpers.R")

#Helper Functions
adjust_lower = function (input_table, threshold) {
  
  for ( i in 1:ncol(input_table)) {
    for( j in 1:nrow(input_table)) {
      
      if(input_table[j,i] < threshold) {
        input_table[j,i] = threshold
        
      }
    }
  }
  return(input_table)
}

adjust_upper = function (input_table, threshold) {
  
  
  for ( i in 1:ncol(input_table)) {
    for( j in 1:nrow(input_table)) {
      
      if(input_table[j,i] > threshold) {
        input_table[j,i] = threshold
        
      }
    }
  }
  
  return(input_table)
}

get_genes = function(input_table) {
  
  genes = input_table[,1]
  genes = as.character((genes))
  
  return(genes)
}

format = function(input_table) {
  
  input_table[,1] = NULL
  
  return(input_table)
  
}

normalize = function(input_table1, annotations) {
  # Tie in the annotaitons with RNA-seq values
  input_table = cbind(annotations, input_table1)
  HK_Pos = c()
  #Total number of columns
  index = 1
  for (i in 1:(ncol(input_table) - 2)) {
    HK_Pos[index] = i + 2
    index = index + 1
  }
  
  attach(input_table)
  #generates the HK table
  newtable = subset(input_table, input_table[,1] == "HK", select=HK_Pos)
  detach(input_table)
  
  #generate a list of the normalized HK means
  meanvalues = c()
  index = 1
  upperindex = 1
  list = c()
  
  #convert newtable values to numeric
  for (i in 1:ncol(newtable)) {
    for (j in 1:nrow(newtable)) {
      list[index] = newtable[j,i]
      index = index + 1
    }
    meanvalues[upperindex] = mean(list)
    upperindex = upperindex + 1
    list = c()
    index = 1
  }
  
  #Generate new table
  input_table = format(input_table)
  input_table[,1] = NULL
  for ( i in 1:ncol(input_table)) {
    for (j in 1:nrow(input_table)){
      input_table[j,i] = input_table[j,i] - meanvalues[i]
    }
  }
  
  return(input_table)
}

annotation_numbers = function(anno) {
  
  colors = c()
  
  if(anno == 2) {
    colored = c("blue", "green")
  }
  else if(anno == 3) {
    colored = c("blue", "green", "red")
  }
  else if(anno == 4) {
    colored = c("blue", "green", "red", "yellow")
  }
  else if(anno == 5) {
    colored = c("blue", "green", "red", "yellow", "plum")
  }
  else if(anno == 6) {
    colored = c("blue", "green", "red", "yellow", "plum", "gray0")
  }
  else if (anno == 7) {
    colored = c("blue", "green", "red", "yellow", "plum", "gray0",
                "aquamarine")
  }
  else if (anno == 8) {
    colored = c("blue", "green", "red", "yellow", "plum", "gray0",
                "aquamarine","white")
  }
  else if (anno == 9) {
    colored = c("blue", "green", "red", "yellow", "plum", "gray0",
                "aquamarine","white", "orange")
  }
  else if (anno == 10) {
    colored = c("blue", "green", "red", "yellow", "plum", "gray0",
                "aquamarine","white", "orange", "lightslateblue")
  }
  else if (anno == 11) {
    colored = c("blue", "green", "red", "yellow", "plum", "gray0",
                "aquamarine","white", "orange", "lightslateblue", "burlywood")
  }
  else if (anno == 12) {
    colored = c("blue", "green", "red", "yellow", "plum", "gray0",
                "aquamarine","white", "orange", "lightslateblue", "burlywood", 
                "darkgrey")
  }
  else if (anno == 13) {
    colored = c("blue", "green", "red", "yellow", "plum", "gray0",
                "aquamarine","white", "orange", "lightslateblue", "burlywood", 
                "darkgrey","deepskyblue", "darkolivegreen","darkolivegreen")
  }
  else if (anno == 14) {
    colored = c("blue", "green", "red", "yellow", "plum", "gray0",
                "aquamarine","white", "orange", "lightslateblue", "burlywood", 
                "darkgrey","deepskyblue", "darkolivegreen")
  }
  else if (anno == 15) {
    colored = c("blue", "green", "red", "yellow", "plum", "gray0",
                "aquamarine","white", "orange", "lightslateblue", "burlywood", 
                "darkgrey", "deepskyblue", "darkolivegreen", "darkolivegreen1")
  }
  else {
    colored = NULL
  }
  
  return (list(colored))
  
}

#Gets the annotation information from the submitted files 
get_annotate = function(input_table) {
  
  input_table = data.frame(input_table)
  
  return(input_table)
  
}

median_center = function(input_table) {
  num_col = ncol(input_table)
  num_row = nrow(input_table)
  
  #List of median values 
  median_values = c()
  
  for (i in 1:num_row){
    rowlist = c()
    #generate list for each column to find the median
    for (j in 1:num_col){
      rowlist[j] = input_table[i,j]
    }
    median_values[i] = median(rowlist)
  }
  #create new table
  for(i in 1:num_row){
    for(j in 1:num_col) {
      input_table[i,j] = input_table[i,j] - median_values[i]
    }
  }
  
  return(input_table)
}

#####################################################################
#Start of Server.R
#####################################################################

function(input, output, session) {
  
  HK_Normal = reactive({
    if (is.null(input$file1))
      return(NULL)
    table = input$file1
    
    
  })

  anno_1 = reactive({
    
    table = input$annotate_1
    
    if (is.null(table))
      return(NULL)
    
    table$datapath
    
    
  })
  
  anno_4 = reactive({
    
    table = input$annotate_4
    
    if (is.null(table))
      return(NULL)
    
    table$datapath
    
    
  })
  
  anno_list_1 = reactive({
    if(is.null(input$annotate_1)){
      return(NULL)
    }
    table = input$annotate_1
    anno_1 = read.csv(table$datapath)
    Anno =  get_annotate(anno_1)
    Anno
    
  })
  
  anno_list_4 = reactive({
    if(is.null(input$annotate_4)){
      return(NULL)
    }
    table = input$annotate_4
    anno_4 = read.csv(table$datapath)
    Anno =  get_annotate(anno_4)
    Anno
    
  })
  
  data1 = reactive({
    
    table = input$file1
    if (is.null(table))
      return(NULL)
    
    table$gene = NULL
    
    table$datapath
    
  })
  
  gene_list = reactive({
    if (is.null(input$file1))
      return(NULL)
    table = input$file1
    gene = read.csv(table$datapath)
    gene = get_genes(gene)
    gene
  })
  
  myplot = function(){
    
    if (is.null(data1()))
      return(NULL)
    
    #Annotation variables to be set
    colann = NULL
    colrow = NULL
    colored = c()
    num_col = 0
    num_row = 0
    
    #Read in and format the csv Table
    gene_table = read.csv(data1())
    gene_table = format(gene_table)
    
    #If normalization desired, then generate a new table
    #to use with the normalized values
    if(input$Normalize) {
      
      #Check if Transcript Annotations Present for normalization
      if(is.null(input$annotate_4)){
        createAlert(session, "no_HK", "exampleAlert", title = "Error",
                    content = "Annotation File Missing", append = FALSE)
        return(NULL)
      }
      
      #reset gene_table to fit the table and normalize
      gene_table = read.csv(data1())
      xy = read.csv(anno_4())
      gene_table = normalize(gene_table, xy)
    }
    
    if(input$Normalize_med){
      gene_table = median_center(gene_table)
    }
    
    if(input$Lower & input$Lower_T) {
      gene_table = adjust_lower(gene_table, input$Lower_T)
    }
    
    if(input$Upper) {
      gene_table = adjust_upper(gene_table, input$Upper_T)
    }
    
    #Check if Col Annotations are present 
    if(!is.null(input$annotate_1)){
      
      #Generate list of annotations 
      colann = anno_list_1()
      num_col = ncol(read.csv(anno_1()))
      
      for ( i in 1:num_col) {
        colored[i] = annotation_numbers(nlevels(colann[,i]))
      }
    }
    
    if(!is.null(input$annotate_4)){
      
      #Still use ncol
      num_row = ncol(read.csv(anno_4()))
      colrow = anno_list_4()
      
      for (i in 1:num_row) {
        colored[num_col + i] = annotation_numbers(nlevels(colrow[,i]))
      }
    }
    
    aheatmap(gene_table, labRow = gene_list(), 
             annCol = colann, annRow = colrow, annColors = colored)
    
  }
  
  output$plot_b = renderPlot({

  myplot()
  
  })
  
  output$zoom <- renderUI({
    plotOutput("plot", width = paste0(input$Zoom), height = input$Zoom)
  })
  
  output$plot <- renderPlot({
    myplot()
  })
  

  final_gene_table = reactive({
    
    #Read in and format the csv Table
    gene_table = read.csv(data1())
    gene_table = format(gene_table)
    
    #If normalization desired, then generate a new table
    #to use with the normalized values
    if(input$Normalize) {
      
      #Check if Transcript Annotations Present for normalization
      if(is.null(input$annotate_4)){
        createAlert(session, "no_HK", "exampleAlert", title = "Error",
                    content = "Annotation File Missing", append = FALSE)
        return(NULL)
      }
      
      #reset gene_table to fit the table and normalize
      gene_table = read.csv(data1())
      xy = read.csv(anno_4())
      gene_table = normalize(gene_table, xy)
    }
    
    if(input$Normalize_med){
      gene_table = median_center(gene_table)
    }
    
    if(input$Lower & input$Lower_T) {
      gene_table = adjust_lower(gene_table, input$Lower_T)
    }
    
    if(input$Upper) {
      gene_table = adjust_upper(gene_table, input$Upper_T)
    }
    
  gene_table = cbind(gene_list(),gene_table)
  gene_table
    
  })
  
  output$downloadData = downloadHandler(
    filename = function(){
      paste(input$text, ".pdf", sep="")
    },
    content = function(file) {
      pdf(file)
      myplot()
      dev.off()
    }
  )
  outputOptions(output, "downloadData", suspendWhenHidden=FALSE)
  
  output$downloadData_2 = downloadHandler(
    filename = function(){
      paste(input$csv, ".csv", sep="")
    },
    content = function(file) {
      write.csv(final_gene_table(), file)
    }
  )
  
  

  outputOptions(output, "downloadData_2", suspendWhenHidden=FALSE)
  
  output$zip = downloadHandler(
    filename = function() {
      paste("Sample_Data", "zip", sep=".")
    },
    
    content = function(file) {
      file.copy("Sample_Data.zip", file)
    }
  )
  
  outputOptions(output, "zip", suspendWhenHidden=FALSE)
}

##################################################
# Old Code, Might be Useful later
##################################################
  
  

  #image = reactive({
  #  outfile = tempfile(fileext='.pdf')
  # if(is.null(outfile)){
  #    return(NULL)
  #  }
  #  pdf(outfile)
  #  myplot()
  #  dev.off()
  #  list(src = outfile, width = 500, height = 500,
  #       alt = "This is alternate text")
    
  #})
  
  #output$zoom = renderImage({
  #  outfile = tempfile(fileext='.pdf')
  #  if(is.null(outfile)){
  #     return(NULL)
  #    }
  #    pdf(outfile)
  #   myplot()
  #    dev.off()
  #  list(src = outfile, width = input$Zoom, height = input$Zoom,
  #         alt = "This is alternate text")
  #  })
  
  

  
  #output$zoom = renderImage({
    
      # Get width and height of image output
      #width  <- session$clientData$output_image1_width
      #height <- session$clientData$output_image1_height
      #npixels <- width * height
    
      # Fill the pixels for R, G, B
      #m <- matrix(1, nrow = height, ncol = width)
      
      #im(m)

      # Convert the vector to an array with 3 planes
      #img <- array(c(m, m, m), dim = c(height, width, 3))
      
      # Write it to a temporary file
      #outfile <- tempfile(fileext = ".png")
      #writePNG(img, target = outfile)
      
      # Return a list containing information about the image
      #list(
      #  src = outfile,
      ##  contentType = "image/png",
      #  width = width,
      #  height = height,
      #  alt = "This is alternate text"
      #)
  #  })
  
#ranges2 <- reactiveValues(x = NULL, y = NULL)

#observe({
#  brush <- input$plot_b
#  if (!is.null(brush)) {
#    ranges2$x <- c(brush$xmin, brush$xmax)
#    ranges2$y <- c(brush$ymin, brush$ymax)
    
#  } else {
#    ranges2$x <- NULL
#    ranges2$y <- NULL
#  }
#})



