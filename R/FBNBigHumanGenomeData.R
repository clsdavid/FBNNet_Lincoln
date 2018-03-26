
#'@export
downloadGenomeData<-function()
{
  ## try http:// if https:// URLs are not supported
  source("https://bioconductor.org/biocLite.R")
  biocLite("GEOquery")

  #use gse2553 <- getGEO('GSE2677',GSEMatrix=TRUE) to download GSE2677 data for example.
  #To detect time of day dependent gene expression in human epidermis suction blister samples from 20 healthy subjects were obtained at three different time points throughout the day. RNA from 20 subjects were used to perform whole genome microarray analysis. Microarrays from 19 subjects showed sufficient quality to perform analysis for differential gene expression. We detected significant differential expression levels for several canonical clock genes such as Per1, Per2, Per3, Bmal1 and Rev-Erb_alpha throughout the day. In total we identified 294 genes that showed significant circadian gene expression including several transcription factors and rate limiting enzymes. To our knowledge this is the first study to investigate genome wide circadian gene expression in human epidermis.


  #1)
  #gathering data from gse35635
  #https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE35635
  gse35635 <- getGEO("GSE35635", GSEMatrix = TRUE)
  eset_gse35635<-gse35635[[1]]
  assayData_gse35635<-assayDataElement(eset_gse35635,'exprs')
  for(i in seq(1:20)){colnames(assayData_gse35635)[[i]]<-paste("T-",colnames(assayData_gse35635)[[i]],"-",i,"-","0", sep="",collapse="")}
  for(i in seq(21:40)){colnames(assayData_gse35635)[[i]]<-paste("T-",colnames(assayData_gse35635)[[i]],"-",i,"-","5", sep="",collapse="")}
  for(i in seq(41:60)){colnames(assayData_gse35635)[[i]]<-paste("T-",colnames(assayData_gse35635)[[i]],"-",i,"-","10", sep="",collapse="")}

  #https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6480
  #search key word GPL570[Accession] time course
  #1)
  #Gene expression data from S. aureus-exposed macrophages
  #GSE13670
  #GPL570 	[HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array

  #*must be GPL570
  #2)
  #	GPL570: [HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array
  #DataSet Record GDS6177:
  #use GSE20489
  #gse <- getGEO("GSE20489") 22.2MB
  #	Acute alcohol consumption effect on whole blood (control group): time course


  #3)
  #GPL570 	[HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array
  #GSE62533
  #Inhibitor of apoptosis proteins as promising therapeutic targets in chronic lymphocytic leukemia

  #4)
  #GPL570 	[HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array
  #GSE57194
  #In Vitro Transformation of Primary Human CD34+ Cells by AML Fusion Oncogenes: Early Gene Expression Profiling Reveals Possible Drug Target in AML

  #5)
  # GPL570 	[HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array
  #Expression data from Leishmania major infected human dendritic cells
  #GSE42088

  #6）
  #GSE41300
}

#'@export
loadExperimentData<-function()
{
  #to remove all variables we can use the following command
  # rm(list=ls(all=TRUE))

  print("load synchronous training series data -> synchronoustrainingseries. \n\n")
  load("synchronoustrainingseries.RDATA",envir = parent.frame())

  print("load synchronous cellcycle cube -> synchronous_cellcycle_cube. \n\n")
  load("synchronous_cellcycle_cube.RDATA",envir = parent.frame())

  print("load asynchronous training series data -> asynchronoustrainingseries. \n\n")
  load("asynchronoustrainingseries.RDATA",envir = parent.frame())

  print("load asynchronous cellcycle cube -> asynchronous_cellcycle_cube. \n\n")
  load("asynchronous_cellcycle_cube.RDATA",envir = parent.frame())

  print("load leukeamia data -> leukeamia. \n\n")
  load("leukeamia.RDATA",envir = parent.frame())

  print("load cellcycle genes -> cellcyclegenes \n\n")
  load("cellcyclegenes.RDATA")

  genes<<-rownames(synchronoustrainingseries[[1]])
  synfbnnetwork<<-mineFBNNetwork(synchronous_cellcycle_cube,genes)
  asynfbnnetwork<<-mineFBNNetwork(asynchronous_cellcycle_cube,genes)
  startstates<<-lapply(synchronoustrainingseries,function(x)x[,1])
  attractor<<-searchForAttractors(synfbnnetwork,startstates,genes)
}

