library(readxl)
#' Read in data from the li-6800
#'
#' This function accepts a licor data file location (xlsx or no-extension (tsv)), and returns a tibble of all the data without the chaff
#' @param location the path to the data file
#' @param purgeComments Removes comments from the file recommended to leave this TRUE - will still work with FALSE but there will be issues
#' @param makeConstCol turns S and Oxygen constants into a column. Currently doesn't work with .xlsx files
#' @param makeCommentsCol Turns comments into their own column
#' @param xlsxIndex which sheet to read on an xlsx.
#' @name licorData
#' @export

licorData <- function(location, purgeComments = T, makeConstCol = F, makeCommentsCol=T,xlsxIndex=1){
  excel <- regexpr(".xlsx$",location)>=0
  if(excel){
    suppressMessages(data <- readxl::read_excel(path = location,sheet = xlsxIndex,col_names = F))
    maxCols <- ncol(data)
    data <- tibble::as_tibble(data)
    makeConstCol <- F #makeconstcol currently does not function with excel imports
    
  } else {
    maxCols <- max(utils::count.fields(file=location,sep = "\t",quote = "")) #This logic is needed to get around the long, narrow header at the top of the TSV
    data <- utils::read.table(location,sep="\t",quote="",dec=".",stringsAsFactors = FALSE,blank.lines.skip = FALSE,skipNul = FALSE,comment.char="",row.names = NULL,col.names = 1:maxCols)
    if(is.na(data[1,length(data)])){  #Terminal \t at the end of every line
      data <- data[,-length(data)] #remove blank column
      maxCols <- maxCols -1
    }
    data <- tibble::as_tibble(data)
  }
  
  if(makeConstCol){ #Store in where and what changes (sysconsts)
    oxylocs <- dfPullConst(df=data,label = "SysConst:Oxygen")
    slocs <- dfPullConst(df=data,label = "Const:S")
  }
  #cat("end of sysconst check ",Sys.time(),"\n")
  counter <- 0 #track how many lines of text we delete for adjustments to the sysconst columns later
  tester <- F
  
  
  if((maxCols-10)<0)
    temp <- maxCols
  else
    temp <- maxCols - 10
  
  rnum <- dplyr::first(which(!(is.na(data[,temp]) | data[,temp]=="")))
  if(!is.na(rnum)){
    data <- data[-c(1:rnum),] #also deletes the category row
    counter <- counter+rnum
  } else{
    stop("There is no data here.")
  }
  
  colnames(data)<-as.character(unlist(data[1,]))
  data<- data[-1,]
  data<- data[-1,] #delete header and unit lines
  counter <- counter+2
  data <- tibble::set_tidy_names(data,quiet = T) #there are like 5 columns all called "time." this renames them to time..2,time..3 etc
  
  if(makeConstCol){
    slocs$row <- slocs$row - counter
    oxylocs$row <- oxylocs$row - counter
    if(nrow(oxylocs)>0){
      data <- dfIncorporate(df = data, newDat = oxylocs, label = "Oxygen")
    }
    if(nrow(slocs)>0){
      data <- dfIncorporate(df = data, newDat = slocs, label = "Leaf Area")
    }
  }
  
  if(purgeComments & !excel){
    #remove comments from the dataframe if requested (default TRUE)
    counter <- 1
    comlocs <- c()
    coms <- c()
    comtimes<- c()
    #I bet we can speed this up by vectorizing
    comlocs <- which(data[,maxCols]=="")-1 # the -1 tells you where to put the resulting comments later
    if(makeCommentsCol & length(comlocs)!=0){
      comtimes <- data[comlocs+1,1]
      coms <- data[comlocs+1,2]
    }
    
    if(length(comlocs)!=0){
      data <- data[-c(comlocs+1),]
    }
    
    if(length(comlocs)!=0 & makeCommentsCol){
      commentdf <- data.frame("hhmmss" = comtimes, "Comments" = coms,stringsAsFactors = F)
      colnames(commentdf) <- c("hhmmss", "Comments")
      colnames(data)[grep("hhmmss",colnames(data))[1]]<- "hhmmss" #rename first instance of hhmmss to just hhmmss for sorting
      data <- tibble::add_column(data,"Comments"=NA,.before=2)
      counter <- 1
      while(counter<= length(comlocs)){
        #add the comments back in as a new row, with NA everywhere except in the comment column and the time column
        #because they are added back in as new rows entirely, they do not have to continue to increment
        data <- tibble::add_row(data,Comments = commentdf[counter,2],hhmmss=commentdf[counter,1],.after = comlocs[counter])
        counter <- counter+1
      }
    }
  }
  
  
  
  if(makeCommentsCol & excel){
    colnames(data)[grep("hhmmss",colnames(data))[1]]<- "hhmmss" #rename first instance of hhmmss for sorting
    data2 <- suppressMessages(readxl::read_excel(path = location,sheet = 2,col_names = F)) #in the xlsx comments are stored on page 2
    commentlocs <- grep(pattern = "^[\\d]{2}",unlist(data2[,1]),perl=T) #only the comments have got numbers at the front of them on page 2
    comments <- unlist(unname(data2[,2][commentlocs,]))
    data <- tibble::add_column(data,"Comments"=NA,.before=2)
    comdf <- data.frame("Comments" = comments,"hhmmss" = unlist(data2[,1])[commentlocs],stringsAsFactors = F)
    data <- dplyr::bind_rows(data,comdf)
    data <- data[order(as.numeric(strptime(data$hhmmss,format="%H:%M:%S"))),] #strptime converts the hhmmss numerics into a sortable time
  }
  
  if("ETR" %in% colnames(data)){
    data$ETR[is.na(data$ETR)] <- ""
  }
  
  suppressMessages(data <- readr::type_convert(data))
  
  return(data)
}

#' Find where constants are changed, and gather those changes into a data frame to return.
#' 
#' @param df The dataframe to search through
#' @param label What to search the dataframe for
#' @param limit_search If TRUE, will limit search to only the first few columns (recommended for large datasets)
#' @name dfPullConst

dfPullConst <- function(df, label,limit_search=T){
  if(limit_search==TRUE){
    constlocs <- tibble::as_tibble(which(arr.ind=T,df[,1:4]==label))
  }else{
    constlocs <- tibble::as_tibble(which(arr.ind=T,df==label))
  }
  if(nrow(constlocs)>0){
    constlocs[,2]<- constlocs[,2]+1
    vals <- c()
    for(i in 1:nrow(constlocs)){
      id <- constlocs[i,]
      vals <- rbind(vals,df[id$row,id$col])
    }
    names(vals) <- "value"
    constlocs <- dplyr::bind_cols(constlocs,vals)
  }
  return(constlocs)
}

#' Incorporate a pulled constant into the dataframe
#'
#' @param df The dataframe to incorporate into
#' @param newDat The new dataframe to add in, with "value"s and "row"s
#' @param label The name of the new column
#' @name dfIncorporate

dfIncorporate <- function(df, newDat, label){
  #Create a vector of the initialization value
  dat <- rep(newDat$value[1],nrow(df))
  #If the pulled constant never changes, we do not need to incorporate the changes.
  if(length(newDat$value)>1){
    #loop through all the new values of the constant.
    for(i in 2:length(newDat$value)){
      #each time, set the value of everything beneath - so, iteratively change the "bottom" part of the vector
      dat[newDat$row[i]:length(dat)] <- newDat$value[i]
      #the other approach, setting the value for only newDat$row[i]:newDat$row[i+1],
      #works if you leave the loop before the last one and execute code similar to above.
      #but this is pretty parsimonious
    }
  }
  return(
    #label is a variable, must unquote it. := is necessary for unquoted "assignments"
    tibble::add_column(df, !!(label):= dat)
  )
}

dat_6_14_1 <- licorData(location = "LC_2023/LiCOR/2023-06-014_logdata_6_14_1_CH.txt",makeConstCol = TRUE)
dat_6_14_2 <- licorData(location = "LC_2023/LiCOR/2023-06-014_logdata_6_14_2_CH.txt",makeConstCol = TRUE)
dat_6_15 <- licorData(location = "LC_2023/LiCOR/2023-06-15-1343_logdata_ch_6_15", makeConstCol = TRUE)
dat_6_15_2 <- licorData(location = "LC_2023/LiCOR/2023-06-15-2054_logdata_ch_6_15_2", makeConstCol = TRUE)
dat_6_16 <- licorData(location = "LC_2023/LiCOR/2023-06-16-0823_logdata_ch_6_16", makeConstCol = TRUE)
