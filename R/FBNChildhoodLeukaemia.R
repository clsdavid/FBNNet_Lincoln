## try http:// if https:// URLs are not supported, use follow to install required bioconductor packages
## source("https://bioconductor.org/biocLite.R")
## biocLite(c("affy", "limma"))
## help(package="limma")
## installed.packages() %>% .[, c("Package", "LibPath")] to find all libary path
## test "C:\\Users\\chenl\\Dropbox\\FBNNet\\ChildhoodLeukeamiaDataFile\\GSE2677_RAW"


#'@export
leukaemiaExperiment<-function(cellDirectory)
{
  rawdata<-getRelatedAffyRawData(cellDirectory)
  nrawdata<-normalizeTimesereisRawData(rawdata)$NormalizedExpData
  stimeseries<-convertIntoSampleTimeSeries(nrawdata)
  sortedtimeseries<-reorderSampleTimeSeries(stimeseries)
  diffgenes<-identifyDifferentiallyExpressedGenes(sortedtimeseries)
  timeseries<-discreteTimeSeries(diffgenes$KnownGenesExp,method="average")

  getClusteredTimeseries<-dividedDiscreteDataintosmallgroups(diffgenes$KnownGenesExp,timeseries)

  #todo:
  #build all cubes for all clusters
  cubeLeukaemia<-constructFBNCubeAndNetworkInClusters(getClusteredTimeseries,5,1,TRUE)
  #merge all networks? seperate by clusters?
  networks<-mergeClusterNetworks(cubeLeukaemia)
  #remove nodes that donot have any connections
  #draw dynamic networks
  #
  return(timeseries)
}
