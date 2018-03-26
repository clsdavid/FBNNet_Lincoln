## try http:// if https:// URLs are not supported, use follow to install required bioconductor packages
## source("https://bioconductor.org/biocLite.R")
## biocLite(c("affy", "limma"))
## help(package="limma")
## installed.packages() %>% .[, c("Package", "LibPath")] to find all libary path
## test "C:\\Users\\chenl\\Dropbox\\FBNNet\\ChildhoodLeukeamiaDataFile\\GSE2677_RAW"

#'A benchmark type function to test a complete process of FBN model
#'Step 1, read Affy data
#' @export
getRelatedAffyRawData<-function(cellDirectory)
{
  # 1 Read in probe level data
  # The affy package will automatically download the appropriate array annotation
  # when you require it. However, if you wish you may download and install the
  # cdf environment you need from http://www.bioconductor.org/packages/release/data/annotation/
  # manually. If there is no cdf environment currently built for your particular chip and you
  # have access to the CDF file then you may use the makecdfenv package to create one
  # yourself. To make the cdf packaes, Microsoft Windows users will need to use the tools
  # described here: http://cran.r-project.org/bin/windows/rw-FAQ.html

  #FIRST solution
  #mycdf <- read.cdffile(cdfFileName)
  #source("http://www.bioconductor.org/biocLite.R")
  #biocLite("affy")
  #biocLite(cdfname)
  require("affy")
  #library(affy)
  #biocLite("HG-U133_Plus_2cdf")
  listCellFiles <- dir(path=cellDirectory, pattern="*\\.CEL", full.names=TRUE)
  affydata <- ReadAffy(filenames=listCellFiles)
  # raw expression data
  expdata <- exprs(affydata)
  samp <- sampleNames(affydata)
  probes <- featureNames(affydata)

  res<-list()
  res[[1]]<-affydata
  names(res)[[1]]<-"RawData"

  res[[2]]<-expdata
  names(res)[[2]]<-"ExpressionData"

  res[[3]]<-samp
  names(res)[[3]]<-"SampleName"

  res[[4]]<-probes
  names(res)[[4]]<-"Probes"

  #cells<-read.affybatch()
  ##alternatively sample from https://www.bioconductor.org/help/workflows/arrays/
  #require(affy) # Affymetrix pre-processing
  #require(limma)#two-color pre-processing; differential expression
  ## import "phenotype" data, describing the experimental design
  #phenoData <- read.AnnotatedDataFrame(system.file("extdata", "pdata.txt",package="arrays"))

  ## RMA normalization
  #celfiles <- system.file("extdata", package="arrays")

  #eset <- justRMA(phenoData=phenoData,celfile.path=celfiles)

  ## differential expression
  #combn <- factor(paste(pData(phenoData)[,1],pData(phenoData)[,2], sep = "_"))
  #fit <- lmFit(eset, design)  # fit each probeset to model
  #efit <- eBayes(fit)        # empirical Bayes adjustment
  #res<-topTable(efit, coef=2)
  return(res)
}

#'A benchmark type function to test a complete process of FBN model
#'Step 2, normalizing data
#' @export
normalizeTimesereisRawData<-function(rawData)
{
  #source("http://www.bioconductor.org/biocLite.R")
  #biocLite("affyPLM")
  require(affyPLM)
  # Normalizing Data
  #
  # The Affy package has implementations of a number of normalization methods
  # for single-channel arrays. This includes (among others):
  #- mas5() - Affymetrix's own normalization program
  #- rma() - 'Robust Multi-Chip' average
  #- gcrma() - A bias-corrected RMA
  # GCRMA is good but takes significantly longer than RMA, so RMA is the
  # most commonly used
  #cat("test")
  #rawData <- as.numeric(as.character(rawData[[1]]))
  #nvals <- rma(rawData)
  nvals <- gcrma(rawData$RawData)
  # normalised expression data
  ned <- exprs(nvals)
  nsamp <- sampleNames(nvals)
  nprobes <- featureNames(nvals)

  res<-list()
  res[[1]]<-nvals
  names(res)[[1]]<-"NormalizedData"

  res[[2]]<-ned
  names(res)[[2]]<-"NormalizedExpData"

  res[[3]]<-nsamp
  names(res)[[3]]<-"SampleName"

  res[[4]]<-nprobes
  names(res)[[4]]<-"FeatureNames"
  return(res)
}

## function(x) paste(x[1],x[2],x[3],sep="-")
#' @export
convertIntoSampleTimeSeries<-function(normalizedData,func=function(x) paste(x[1],x[2],x[3],sep="-"), splitor=as.character("-"))
{

  mat<-normalizedData
  #Obtain the last part of each column names
  groups <- sapply(strsplit(x = colnames(mat), split = splitor), func)
  print(groups)
  #Go through each unique column name and extract the corresponding columns
  res<-lapply(unique(groups), function(x) mat[,which(groups == x)])
  names(res)<-unique(groups)

  #need to convert column name into numbers
  #str_extract_all(names(df), '[0-9]+')
  return(res)
}

#' @export
reorderSampleTimeSeries<-function(convertedTimeSeries,func=function(x) x[length(x)],splitor=as.character("-"))
{
  res<-lapply(1:length(convertedTimeSeries),function(index)
  {
    result<-convertedTimeSeries[[index]]
    ocolnames<-colnames(result)
    ocolnames<-gsub("\\..*","",ocolnames)
    ocolnames <- sapply(strsplit(x = ocolnames, split = splitor),func)
    ocolnames <- as.integer(str_extract_all(ocolnames, '[0-9]+'))
    colnames(result)<-ocolnames
    return(result[, order(as.numeric(colnames(result)))])
  })
  names(res)<-names(convertedTimeSeries)
  return(res)
}
#'A benchmark type function to test a complete process of FBN model
#'Step 3, identify significantly expressed genes, which are strongly related with the samples / study purposes
#' @export
identifyDifferentiallyExpressedGenes<-function(orderSampleTimeSeries)
{
  #The expression values from RMA are log2 transformed, so to calculate
  #the log ratio you simply subtract one from the other.
  #log2(x) - log2(y) = log2(x/y)
  #Note here the the log ratio gives you the fold change directly. A log
  #ratio of 1 = 2-fold up regulated and a log ratio of -1 = 2-fold down
  #regulated (when comparing x vs y).
  tempdata<-lapply(orderSampleTimeSeries,function(input){
     cols<-colnames(input)
     if(length(cols)<2)
     {
       stop("The input must have columns more than one")
     }

     res<-round(log2(input[,2]/input[,1]),5)

     for(i in seq_along(cols))
     {
       if(i>2)
       {
         res<-cbind(res,round(log2(input[,i]/input[,i-1]),5))
       }
     }
     return(res)
      #res1<-cbind(res)
      #res2<-cbind(round(log2(input[,2]/input[,1]),5),round(log2(input[,3]/input[,2]),5))
      #return(res2)
    })
  #tempdata<-lapply(orderSampleTimeSeries,function(input)cbind(round(log2(input[,2]/input[,1]),4),round(log2(input[,3]/input[,2]),4)))
  diffExpressed<-lapply(tempdata,function(subdata) subdata[which(!apply(subdata,1,FUN = function(x){
    rowvalue<-as.numeric(x)
    #print(paste(rowvalue," length=",length(dim(rowvalue)),sep=""))
    #if(!is.vector(rowvalue))
    #{
      #print(is.vector(rowvalue))
      #return(TRUE)
    #}

    all(rowvalue > -1 & rowvalue<1)
    })),]) #find differently expressed genes
  alldifferExpressednames<-unique(unlist(sapply(diffExpressed,function(x)rownames(x)))) #find the maximum common gene set
  probesetMapping<-mapToGeneNameWithhgu133plus2Db(alldifferExpressednames)
  filteredDataKnown<-lapply(orderSampleTimeSeries,function(subdata)subdata[probesetMapping$MappedProbeset[,1],])
  filteredDataUnKnown<-lapply(orderSampleTimeSeries,function(subdata)subdata[probesetMapping$UnMappedProbeset[,1],])

  #update known probset names to gene names
  filteredDataKnown<-lapply(filteredDataKnown,function(mtx){
    geneNames<-probesetMapping$MappedProbeset[which(probesetMapping$MappedProbeset$PROBEID%in%rownames(mtx)),]$SYMBOL
    rownames(mtx)<-geneNames
    return(mtx)
  })

  #find all common genes in the result, the result should be
  for(i in seq_along(filteredDataKnown))
  {
    subset<-filteredDataKnown[[i]]
    if(i>1)
    {
      subset<-subset[rownames(subset)%in%commonNames,]
      commonNames<-rownames(subset)
    }else
    {
      commonNames<-rownames(subset)
    }
  }
  #allCommonGeneNames<-unique(unlist(sapply(filteredDataKnown,function(x)rownames(x)))) #find the maximum common gene set
  filteredDataKnown<-lapply(filteredDataKnown,function(subdata)subdata[commonNames,])
  sampleNames<-names(filteredDataKnown)
  allDiffData<-lapply(sampleNames,function(sampleName)rbind(filteredDataKnown[[sampleName]],filteredDataUnKnown[[sampleName]]))
  names(allDiffData)<-sampleNames
  filteredData<-list()
  filteredData[[1]]<-filteredDataKnown
  filteredData[[2]]<-filteredDataUnKnown
  filteredData[[3]]<-allDiffData
  filteredData[[4]]<-probesetMapping
  names(filteredData)[[1]]<-"KnownDiffExpressedTimeSeriesData"
  names(filteredData)[[2]]<-"UnKnownDiffExpressedTimeSeriesData"
  names(filteredData)[[3]]<-"AllDiffExpressedTimeSeriesData"
  names(filteredData)[[4]]<-"ProbesetMapping"
  return(filteredData)
}

#' @export
mapPropersetsWithhgu133plus2Db<-function(orderSampleTimeSeries)
{
  commonNames<-c()
  for(i in seq_along(orderSampleTimeSeries))
  {
    subset<-orderSampleTimeSeries[[i]]
    if(i>1)
    {
      subset<-subset[rownames(subset)%in%commonNames,]
      commonNames<-rownames(subset)
    }else
    {
      commonNames<-rownames(subset)
    }
  }
  alldifferExpressednames<-commonNames
  #alldifferExpressednames<-unique(unlist(sapply(orderSampleTimeSeries,function(x)rownames(x)))) #find the maximum common gene set
  probesetMapping<-mapToGeneNameWithhgu133plus2Db(alldifferExpressednames)
  filteredDataKnown<-lapply(orderSampleTimeSeries,function(subdata)subdata[probesetMapping$MappedProbeset[,1],])
  filteredDataUnKnown<-lapply(orderSampleTimeSeries,function(subdata)subdata[probesetMapping$UnMappedProbeset[,1],])

  #update known probset names to gene names
  filteredDataKnown<-lapply(filteredDataKnown,function(mtx){
    geneNames<-probesetMapping$MappedProbeset[which(probesetMapping$MappedProbeset$PROBEID%in%rownames(mtx)),]$SYMBOL
    rownames(mtx)<-geneNames
    return(mtx)
  })


  #remove duplicate gene names
  filteredDataKnown<-lapply(filteredDataKnown,function(mtx){
    return(mtx[!duplicated(mtx),])
  })

  #find all common genes in the result, the result should be
  for(i in seq_along(filteredDataKnown))
  {
    subset<-filteredDataKnown[[i]]
    if(i>1)
    {
      subset<-subset[rownames(subset)%in%commonNames,]
      commonNames<-rownames(subset)
    }else
    {
      commonNames<-rownames(subset)
    }
  }
  #allCommonGeneNames<-unique(unlist(sapply(filteredDataKnown,function(x)rownames(x)))) #find the maximum common gene set
  filteredDataKnown<-lapply(filteredDataKnown,function(subdata)subdata[commonNames,])

  filteredData<-list()
  filteredData[[1]]<-filteredDataKnown
  filteredData[[2]]<-filteredDataUnKnown
  filteredData[[3]]<-probesetMapping
  names(filteredData)[[1]]<-"KnownGenesExp"
  names(filteredData)[[2]]<-"UnKnownGenesExp"
  names(filteredData)[[3]]<-"ProbesetMapping"
  return(filteredData)
}

#data mapping
#internal method
#' @export
mapToGeneNameWithhgu133plus2Db<-function(probeset)
{
  #if the following package hasen't been installed, install them
  #source("http://bioconductor.org/biocLite.R")
  #biocLite("annotate")
  #biocLite("hgu133plus2.db")
  require("AnnotationDbi")
  require("hgu133plus2.db")
  #find all types
  keyTypes<-keytypes(hgu133plus2.db)

  res<-list()
  getGeneNames<-select(hgu133plus2.db, probeset, c("SYMBOL","ENTREZID", "GENENAME"))
  mapped<-getGeneNames[!duplicated(getGeneNames[,"ENTREZID"]),] #de-duplicated
  mapped<-mapped[!is.na(mapped[,"ENTREZID"]),]
  unmapped<-getGeneNames[is.na(getGeneNames[,"ENTREZID"]),] #get all unknown probeset
  res[[1]]<-mapped
  res[[2]]<-unmapped
  names(res)[[1]]<-"MappedProbeset"
  names(res)[[2]]<-"UnMappedProbeset"
  return(res)
  #find all types
  #keytypes(hgu133plus2.db)
}

#This method is used to construct time series data for Leukeamia data from GSE2677
#Experimental method
#'@export
constructTestGSE2677Data<-function(files)
{
  #files<-c("C:\\Users\\chenl\\Dropbox\\FBNNet\\ChildhoodLeukeamiaDataFile\\GSE2677_RAW","C:\\Users\\chenl\\Dropbox\\FBNNet\\Genome Data\\GSE13670_RAW\\GSE13670_RAW","C:\\Users\\chenl\\Dropbox\\FBNNet\\Genome Data\\GSE20489_RAW\\GSE20489_RAW","C:\\Users\\chenl\\Dropbox\\FBNNet\\Genome Data\\GSE42088_RAW\\GSE42088_RAW"，"C:\\Users\\chenl\\Dropbox\\FBNNet\\Genome Data\\GSE54992_RAW\\GSE54992_RAW"，"C:\\Users\\chenl\\Dropbox\\FBNNet\\Genome Data\\GSE57194_RAW\\GSE57194_RAW","C:\\Users\\chenl\\Dropbox\\FBNNet\\Genome Data\\GSE62533_RAW\\GSE62533_RAW")

  if(!is_vector(files))
  {
    stop("The parameter files must be a vector of string, i.e. a list of data files")
  }

  combinedTimeseries<-list()
  for(filename in files)
  {
    print("Start to process file:")
    print(filename)
    rawdata<-getRelatedAffyRawData(filename)
    nrawdata<-normalizeTimesereisRawData(rawdata)$NormalizedExpData
    stimeseries<-convertIntoSampleTimeSeries(nrawdata)
    sortedtimeseries<-reorderSampleTimeSeries(stimeseries)
    combinedTimeseries[[length(combinedTimeseries)+1]]<-sortedtimeseries
    print("Finish")
  }

  res<-dissolve(combinedTimeseries)
  return(res)
  #diffgenes<-identifyDifferentiallyExpressedGenes(sortedtimeseries)
  #timeseries<-discreteTimeSeries(diffgenes$KnownGenesExp,method="average")
  #getClusteredTimeseries<-dividedDiscreteDataintosmallgroups(diffgenes$KnownGenesExp,timeseries)
}
