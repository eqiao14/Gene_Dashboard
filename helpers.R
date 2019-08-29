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




