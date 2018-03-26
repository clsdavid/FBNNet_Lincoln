#'Split a compressed gene regulate function into symbols
#'
#'@param expression A regulate function like GeneA&GeneB&!GeneC
#'@param outputType A output type with two options 1 the splition included "!", others not included "!"
#'@param lowerCase Is lower case required?
#'@return The return is a vector of the input regulate function, a example is "GeneA" "&" "GeneB" "&" "!" "GeneC"
#'@examples
#' splitExpression("GeneA&GeneB&!GeneC")
#' splitExpression("GeneA&GeneB&!GeneC",1)
#' splitExpression("GeneA&GeneB&!GeneC",0)
#'@export
splitExpression <- function(expression, outputType=1,lowerCase=FALSE)
{
  tryCatch(
    {
      if(identical(expression,1)|identical(expression,0))
      {
        return (c(expression))
      }

      if (lowerCase)
      {
        expression <- tolower(expression)
      }
      if(outputType==1)
      {
        # add extra whitespace to brackets and operators
        expression <- gsub("(\\(|\\)|\\[|\\]|\\||\\&|\\!|\\=|\\.\\.|[a-zA-Z0-9_]*)"," \\1 ", expression)
      }else
      {
        # add extra whitespace to brackets and operators
        expression <- gsub("(\\(|\\)|\\[|\\]|\\||\\&|\\!|\\=|\\.\\.|[a-zA-Z0-9_!]*)"," \\1 ", expression)
      }


      # strip multiple whitespace characters
      expression <- gsub("[ ]+", " ", expression)
      expression <- gsub("^[ ]+", "", expression)
      expression <- gsub("[ ]+$", "", expression)

      # split up at whitespace positions
      res <- strsplit(expression, " ", fixed=TRUE)[[1]]
    },
    error = function(e)
    {
      sprintf("Error converting to FBN Network \"%s\": %s", expression, e$message)
      stop(sprintf("Error converting to FBN Network \"%s\": %s", expression, e$message))

    })
  return(res)
}

#'Remove duplidate timeseries sample data
#'
#'@param timeseriescube A list of time series data
#'@return The return a list of time series data
#'@examples
#' mat1<-matrix(c("1","2","3","4","5","6","7","8","9"),3,3)
#' mat2<-matrix(c("1","2","3","4","5","6","7","8","9"),3,3)
#' listtest<-list(mat1,mat2)
#' FBNDataReduction(listtest)
#'@export
FBNDataReduction<-function(timeseriescube)
{
  shrinkTimeseriesCube<-function(timeseriescube)
  {
    duplicateIndexes <- duplicated(timeseriescube)
    deduplicatedCube <- timeseriescube[!duplicateIndexes]
    return(deduplicatedCube)
  }
  return(shrinkTimeseriesCube(timeseriescube))
}

#'@export
removeduplicateCol<-function(mat)
{
  res<-c()
  cols<-dim(mat)[2]
  index<-1L
  res<-lapply(1:cols,function(index,mat){
    if(index>1)
    {
      cond<-abs(mat[,index-1]-mat[,index])
      if(!all(cond%in%0))
      {
        res<-mat[,index]
      }
    }else
    {
      res<-mat[,index]
    }
    return(res)
  },mat)

  return(do.call(cbind,res[!length(res)==0]))
}

#'@export
dissolve <- function(x){
  operator <- function(x,name=NULL){
    if(is.list(x)){
      for(i in seq(x)){
        operator(x[[i]],names(x)[[i]])
      }
    }else{
      combi[[length(combi)+1]] <<- x
      names(combi)[[length(combi)]]<<-name
    }
  }
  combi=list()
  operator(x)
  return(combi)
}

cosine<-function(x,y)
{
  return(crossprod(x, y)/sqrt(crossprod(x) * crossprod(y)))
}

#'@export
convertMatrixIntoListByRow<-function(mat)
{
  res<- lapply(seq_len(nrow(mat)), function(i) mat[i,])
  names(res)<-c(rownames(mat))
  return(res)
}

#'@export
similarityBetweenMatrix<-function(timeseries1,timeseries2,index)
{
  if(!identical(dim(timeseries1),dim(timeseries2)))
  {
    stop("The two matrixes must have the same dimensions")
  }

  differ<-abs(timeseries1-timeseries2)

  #correlation<-cor(c(timeseries1),c(timeseries2))
  zerosum<-length(differ[differ==0])
  correlation<-zerosum/length(differ)

  if(is.na(correlation) | is.null(correlation))
  {
    correlation<-0
  }

  if(correlation<=0.2)
    return (c("veryunlikely",correlation,index))

  if(correlation>0.2&correlation<=0.4)
    return ((c("unlikely",correlation,index)))

  if(correlation>0.4&correlation<=0.6)
    return ((c("likely",correlation,index)))

  if(correlation>0.6&correlation<=0.8)
    return ((c("similar",correlation,index)))

  if(correlation>0.8)
    return ((c("verysimilar",correlation,index)))
}

#'Check the similarity between time series
#'
#'@param originalTimeseriesCube The original data set that contains samples and each sample contains genes and time points
#'@param reconstructedTimeSeriesCube The reconstructed data set that contains samples and each sample contains genes and time points
#'@return similarity report
#'@examples
#' coming later
#'@export
checkSimilarity<-function(originalTimeseriesCube,reconstructedTimeSeriesCube)
{
  res<-list()
  if(!identical(length(originalTimeseriesCube),length(reconstructedTimeSeriesCube)))
  {
    stop("The length of each cube must be identical")
  }

  for(i in seq_along(originalTimeseriesCube))
  {
    if(!identical(dim(originalTimeseriesCube[[i]]),dim(reconstructedTimeSeriesCube[[i]])))
    {
      stop("The dimension of each cube must be identical")
    }
    res[[i]]<-similarityBetweenMatrix(originalTimeseriesCube[[i]],reconstructedTimeSeriesCube[[i]],i)
  }

  return(res)
}
getSubsetBasedOnSimilarity<-function(timeseries,similarityreport)
{
  res<-list()
  for(i in seq_along(similarityreport))
  {
    entry<-similarityreport[[i]]
    index<-as.numeric(entry[[3]])
    res[[i]]<-timeseries[[index]]
  }

  return(res)
}

#'Generate similarity report
#'
#'@param similarityreport The raw similarity report which was created by the function checkSimilarity
#'@return An organized similarity report
#'@examples
#' coming later
#'@export
generateSimilarReport<-function(similarityreport)
{
  cond1<-sapply(similarityreport,function(entry)as.numeric(entry[[2]]) >=0.9)
  cond2<-sapply(similarityreport,function(entry)as.numeric(entry[[2]])>=0.8 & as.numeric(entry[[2]])<0.9)
  cond3<-sapply(similarityreport,function(entry)as.numeric(entry[[2]])>=0.7 & as.numeric(entry[[2]])<0.8)
  cond4<-sapply(similarityreport,function(entry)as.numeric(entry[[2]])>=0.6 & as.numeric(entry[[2]])<0.7)
  cond5<-sapply(similarityreport,function(entry)as.numeric(entry[[2]])>=0.5 & as.numeric(entry[[2]])<0.6)
  cond6<-sapply(similarityreport,function(entry)as.numeric(entry[[2]])>=0.4 & as.numeric(entry[[2]])<0.5)
  cond7<-sapply(similarityreport,function(entry)as.numeric(entry[[2]])>=0.3 & as.numeric(entry[[2]])<0.4)
  cond8<-sapply(similarityreport,function(entry)as.numeric(entry[[2]])>=0.2 & as.numeric(entry[[2]])<0.3)
  cond9<-sapply(similarityreport,function(entry)as.numeric(entry[[2]])>=0.1 & as.numeric(entry[[2]])<0.2)
  cond10<-sapply(similarityreport,function(entry)as.numeric(entry[[2]])<0.1)

  res<-list()
  res[[1]]<-similarityreport[cond1]
  res[[2]]<-similarityreport[cond2]
  res[[3]]<-similarityreport[cond3]
  res[[4]]<-similarityreport[cond4]
  res[[5]]<-similarityreport[cond5]
  res[[6]]<-similarityreport[cond6]
  res[[7]]<-similarityreport[cond7]
  res[[8]]<-similarityreport[cond8]
  res[[9]]<-similarityreport[cond9]
  res[[10]]<-similarityreport[cond10]
  names(res)[[1]]<-"A"
  names(res)[[2]]<-"B"
  names(res)[[3]]<-"C"
  names(res)[[4]]<-"D"
  names(res)[[5]]<-"E"
  names(res)[[6]]<-"F"
  names(res)[[7]]<-"G"
  names(res)[[8]]<-"H"
  names(res)[[9]]<-"I"
  names(res)[[10]]<-"J"
  return(res)
}

#'Remove some time points out from the target time series
#'
#'@param timeseriesCube A data set that contains samples and each sample contains genes and time points
#'@param reducedTimePoints How many time points reduced between two observation time points
#'@return reduced time series cube
#'@examples
#' coming later
#'@export
timeseriesReduction<-function(timeseriesCube,reducedTimePoints)
{
  tryCatch(
    {
      #get FBN matrix
      res<-lapply(timeseriesCube,function(series){
        subres<-series
        colHeaders<-colnames(series)
        rowHeaders<-rownames(series)

        timpoints<-length(colHeaders)
        if(timpoints%%2!=0)
        {
          timpoints<-timpoints-1
        }

        #remove the beginning and end of time points
        timpoints<-timpoints-2
        interval<-timpoints/reducedTimePoints

        #if(timpoints%%reducedTimePoints!=0 | (length(colHeaders)%%2==0 & interval%%2!=0))
        #{
        #stop(sprintf("The time sereis can not be reduced by %s equally",reducedTimePoints))
        #}

        if(interval%%2!=0)
        {
          interval
        }

        i<-1
        start<-i+1

        while(i <=(reducedTimePoints-1))
        {
          end<-start+interval-1
          subres<-subres[,-(start:end)]
          i<-i+1
          start<-i+1
          #colnames(subres[,i])<-colHeaders[[i]]
        }
        return(subres)
      })
      return(res)
    },
    error = function(e)
    {
      stop(sprintf("Error executing Time sereis reduction: %s", e$message))
    })

}

timeseriesReductionBasedOnVector<-function(timeseriesCube,vector)
{
  tryCatch(
    {
      #get FBN matrix
      res<-lapply(timeseriesCube,function(series){
        return(series[,vector])
      })
      return(res)
    },
    error = function(e)
    {
      stop(sprintf("Error executing Time sereis reduction: %s", e$message))
    })

}
extractTimeseriesDataByColumn<-function(timeseriesmatrix,numofCulumns)
{
  if(!is.matrix(timeseriesmatrix))
  {
    stop("The input data is not a matrix")
  }
  columnames<-colnames(timeseriesmatrix)

  if(length(columnames)%%numofCulumns!=0)
  {
    stop("The data cannot be split equallly")
  }

  numofsplits<-length(columnames)/numofCulumns
  res<-list()
  i<-1
  j<-1
  while(i<=numofsplits)
  {
    res[[i]]<-timeseriesmatrix[,j:(i*numofCulumns)]
    j<-1+(i*numofCulumns)
    i<-i+1
  }
  return(res)
}

##configulation
#parameters for Cube
#' @export
p_value_FBNCube<-0.05

#' @export
thredsod_FBNCube<-0.5 #dependes on the level of confidence

#' @export
Mutual_toleranceRate<--1

#' @export
minMutual_range<-0.5 #1 for normal time series maximum is one

#' @export
enableclusterpruning<-TRUE # FALSE for short time series?
#minSupport_FBN<-0.001#support and confidence framework has limitations
#parameters for miner
#defualt

#' @export
error_toleranceRate<-0

#' @export
maxFBNRulesOfActivator<-5

#' @export
maxFBNRulesOfInhibitor<-5

#' @export
configHasLoaded<-FALSE
#' @export
readConfigFile<-function()
{
  if(configHasLoaded==FALSE)
  {
    library(yaml)
    config = yaml.load_file("config/config.yml")
    p_value_FBNCube<-as.numeric(config$cube$p_value_FBNCube)
    thredsod_FBNCube<-as.numeric(config$cube$thredsod_FBNCube)
    Mutual_toleranceRate<-as.numeric(config$cube$Mutual_toleranceRate)
    minMutual_range<-as.numeric(config$cube$minMutual_range)
    enableclusterpruning<-as.numeric(config$cube$enableclusterpruning)
    error_toleranceRate<-as.numeric(config$network$error_toleranceRate)
    maxFBNRulesOfActivator<-as.numeric(config$network$maxFBNRulesOfActivator)
    maxFBNRulesOfInhibitor<-as.numeric(config$network$maxFBNRulesOfInhibitor)

    options("p_value_FBNCube"=p_value_FBNCube)  #using getOption("mypkg-myval") to retreive the option
    options("thredsod_FBNCube"=thredsod_FBNCube)
    options("Mutual_toleranceRate"=Mutual_toleranceRate)
    options("minMutual_range"=minMutual_range) #1 for normal time series maximum is one
    options("enableclusterpruning"=enableclusterpruning) # FALSE for short time series?
    options("error_toleranceRate"=error_toleranceRate)
    options("maxFBNRulesOfActivator"=maxFBNRulesOfActivator)
    options("maxFBNRulesOfInhibitor"=maxFBNRulesOfInhibitor)

    configHasLoaded<-TRUE
  }
  #return (config)
}

#'Retrieving genome sequence data using SeqinR. The code is a sample from the book "a little book of r for bioinformatic"
#'
#'@param accession NBCI accession number
#'@return genome sequence
#'@examples
#' coming later
#'@export
getncbiseq <- function(accession)
{
  require("seqinr") # this function requires the SeqinR R package
  # first find which ACNUC database the accession is stored in:
  dbs <- c("genbank","refseq","refseqViruses","bacterial")
  numdbs <- length(dbs)
  for (i in 1:numdbs)
  {
    db <- dbs[i]
    choosebank(db)
    # check if the sequence is in ACNUC database ’db’:
    resquery <- try(query(".tmpquery", paste("AC=", accession)), silent = TRUE)
    if (!(inherits(resquery, "try-error")))
    {
      queryname <- "query2"
      thequery <- paste("AC=", accession, sep="")
      query('queryname','thequery')
      # see if a sequence was retrieved:
      seq <- getSequence(query2$req[[1]])
      closebank()
      return(seq)
    }
    closebank()
  }
  print(paste("ERROR: accession",accession,"was not found"))
}

#'@export
outputFasta<-function(names,sequenceData)
{
  newfile<-paste(names, ".fasta", sep="")
  write.fasta(names=names, sequences=dengueseq, file.out=newfile)
  print(paste("Output sequence to",newfile))
}

#'@export
readFasta<-function(fastaFileName)
{
  require("seqinr")
  fastadata <- read.fasta(file = fastaFileName)
  return(fastadata)
}

#'@export
slidingwindowplot <- function(windowsize, inputseq)
{
  starts <- seq(1, length(inputseq)-windowsize, by = windowsize)
  n <- length(starts) # Find the length of the vector "starts"
  chunkGCs <- numeric(n) # Make a vector of the same length as vector "starts", but just containing
  for (i in 1:n) {
    chunk <- inputseq[starts[i]:(starts[i]+windowsize-1)]
    chunkGC <- GC(chunk)
    print(chunkGC)
    chunkGCs[i] <- chunkGC
  }
  plot(starts,chunkGCs,type="b",xlab="Nucleotide start position",ylab="GC content")
}

#'@export
rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}

#'@export
rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}


