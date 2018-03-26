
#' @export
constructFBNCubeAndNetworkInClusters<-function(clusteredTimeseriesData,maxK=5,temporal=1,useParallel=FALSE)
{
  print("Starting.....")
  clusteredFBNCube<-lapply(clusteredTimeseriesData,function(timeseriesData,maxK=5,temporal=1,useParallel=FALSE){
    res<-list()
    #get row values i.e. genes
    print("Creating sub cube")
    genes<-rownames(timeseriesData[[1]])
    subcube<-constructFBNCube(genes,genes,timeseriesData,maxK,temporal,useParallel)
    res[[1]]<-subcube
    names(res)[[1]]<-"Cube"
    res[[2]]<-mineFBNNetwork(subcube,genes)
    names(res)[[2]]<-"Network"
    #print("Created sub cube")
    return(res)
  },maxK,temporal,useParallel)
  class(clusteredFBNCube)<-"ClusteredFBNCube"
  print("save the cube.")
  save(clusteredFBNCube, file='clusteredFBNCube.Rdata')
  print("End.....")
  return(clusteredFBNCube)
}

#'Create an Orchard cube
#'
#'@param targetgenes A vector of genes that will be treated as target genes
#'@param allgenes All genes that are available for building up the cube
#'@param timeseriesCube A list of samples in which a sample is a matrix that contains gene states where genes in rows and time points in columms
#'@param maxK The maximum level the cube can dig in
#'@param temporal A value that used to be 1 indicates the previous steps the current one can depend on
#'@param useParallel If it is TRUE, the constructing will run it in parallel, otherwise in a singl thread
#'@return A Orcahrd cube that contains all precomputed measures
#'@examples
#' mat1<-matrix(c("1","0","0","1","0","0","0","1","1"),3,3, dimnames=list(c("gene1","gene2","gene3"),c("1","2","3")))
#' mat2<-matrix(c("1","1","0","1","0","1","1","1","0"),3,3, dimnames=list(c("gene1","gene2","gene3"),c("1","2","3")))
#' listtest<-list(mat1,mat2)
#' constructFBNCube(c("gene1","gene2"),c("gene1","gene2","gene3"),listtest,4,1,FALSE)
#' @export
constructFBNCube<-function(targetgenes,allgenes,timeseriesCube,maxK=5,temporal=1,useParallel=FALSE)
{
  #construct gene tree by timeseries * samples * timepoints(columns)
  internalloopByWhole<-function(i,targetgenes,allgenes,clusters,timeseriesCube,maxK,temporal,getCurrentState,getpreviousState){
    getStatisticInternal<-function(result)
    {
      subAvgErrorT<-unlist(lapply(result,function(rule)
        return (as.numeric(rule$Statistic$avgErrorP))
      ))

      subAvgErrorF<-unlist(lapply(result,function(rule)
        return (as.numeric(rule$Statistic$avgErrorN))
      ))

      subAvgSignalT<-unlist(lapply(result,function(rule)
        return (as.numeric(rule$Statistic$avgSignalP))
      ))

      subAvgSignalF<-unlist(lapply(result,function(rule)
        return (as.numeric(rule$Statistic$avgSignalN))
      ))

      minAvgErrorT<-min(subAvgErrorT)
      maxAvgErrorT<-max(subAvgErrorT)
      minAvgErrorF<-min(subAvgErrorF)
      maxAvgErrorF<-max(subAvgErrorF)

      avgErrorT<-sum(subAvgErrorT)/length(subAvgErrorT)
      avgErrorF<-sum(subAvgErrorF)/length(subAvgErrorF)

      stdErrorT<-sd(subAvgErrorT)
      stdErrorF<-sd(subAvgErrorF)

      avgSignalT<-sum(subAvgSignalT)/length(subAvgSignalT)
      avgSignalF<-sum(subAvgSignalF)/length(subAvgSignalF)

      stdSignalT<-sd(subAvgSignalT)
      stdSignalF<-sd(subAvgSignalF)

      res<-list()
      index<-length(res)+1
      res[[index]]<-avgErrorT
      names(res)[[index]]<-"avgErrorP"
      index<-index+1

      res[[index]]<-avgErrorF
      names(res)[[index]]<-"avgErrorN"
      index<-index+1

      res[[index]]<-minAvgErrorT
      names(res)[[index]]<-"minAvgErrorP"
      index<-index+1

      res[[index]]<-maxAvgErrorT
      names(res)[[index]]<-"maxAvgErrorP"
      index<-index+1

      res[[index]]<-minAvgErrorF
      names(res)[[index]]<-"minAvgErrorN"
      index<-index+1

      res[[index]]<-maxAvgErrorF
      names(res)[[index]]<-"maxAvgErrorN"
      index<-index+1

      res[[index]]<-avgSignalT
      names(res)[[index]]<-"avgSignalP"
      index<-index+1

      res[[index]]<-avgSignalF
      names(res)[[index]]<-"avgSignalN"
      index<-index+1

      res[[index]]<-stdErrorT
      names(res)[[index]]<-"stdErrorP"
      index<-index+1

      res[[index]]<-stdErrorF
      names(res)[[index]]<-"stdErrorN"
      index<-index+1

      res[[index]]<-stdSignalT
      names(res)[[index]]<-"stdSignalP"
      index<-index+1

      res[[index]]<-stdSignalF
      names(res)[[index]]<-"stdSignalN"
      index<-index+1

      return(res)
    }
    try({
      res<-list()
      res[[i]]<-list()
      names(res)[[i]]<-targetgenes[[i]]

      gene<-targetgenes[[i]]

      if(length(clusters[[gene]])==0)
      {
        return(NULL)
      }
      if(enableclusterpruning)
      {
        allgenes<-names(clusters[[gene]])
      }

      allRelatedTargetGeneMeasures<-clusteringOnTargetGene(gene,getCurrentState,getpreviousState,timeseriesCube,allgenes,temporal,FALSE)

      if(is.null(allRelatedTargetGeneMeasures)|length(allRelatedTargetGeneMeasures)==0)
      {
        return(NULL)
      }
      allRelatedConditionalGenesMeasures<-clusters
      response<-buildProbabilityTreeOnTargetGene(gene,getCurrentState,getpreviousState,timeseriesCube,allgenes,NULL,NULL,maxK,temporal,allRelatedTargetGeneMeasures,allRelatedConditionalGenesMeasures,FALSE)

      cond1<-sapply(response,function(entry)!is.null(entry))
      response<-(response[cond1][unlist(lapply(response[cond1],length)!=0)])

      if(is.null(response))
      {
        return(NULL)
      }
      res[[i]][[1]]<-response
      names(res[[i]])[[1]]<-"SubGenes"
      res[[i]][[2]]<-getStatisticInternal(response)
      names(res[[i]])[[2]]<-"Statistic"

      return(res)
    })}

  time1<-as.numeric(Sys.time())
  #genes<-fbnnetwork$genes
  myCube<-timeseriesCube

  cols<-colnames(timeseriesCube[[1]])

  res<-list()
  reducedCube<-FBNDataReduction(myCube)

  completeState<-extractGeneStateFromTimeSeriesCube(timeseriesCube)
  getCurrentState<-completeState
  getpreviousState<-completeState

  #construct data
  getCurrentState<-getCurrentState[,-1]#remove the first one
  getpreviousState<-getpreviousState[,-length(colnames(getpreviousState))] #remove the last one

  #clustering
  if(useParallel)
  {
    clusters<-doParallelWork(function(iIndex,allgenes,timeseriesCube,temporal,geneState){
      allRelatedConditionalGenesMeasures<-list()
      res<-clusteringOnTargetGene(allgenes[[iIndex]],geneState,geneState,timeseriesCube,allgenes,temporal,FALSE)
      if(!is.null(res)|length(res)>0)
      {
        allRelatedConditionalGenesMeasures[[1]]<-res
        names(allRelatedConditionalGenesMeasures)[[1]]<-allgenes[[iIndex]]
        return(allRelatedConditionalGenesMeasures)
      }else
      {
        return(NULL)
      }

    },allgenes,reducedCube,temporal,completeState)

  }else
  {
    clusters<-doNonParallelWork(function(iIndex,allgenes,timeseriesCube,temporal,geneState){
      allRelatedConditionalGenesMeasures<-list()
      res<-clusteringOnTargetGene(allgenes[[iIndex]],geneState,geneState,timeseriesCube,allgenes,temporal,FALSE)
      if(!is.null(res)|length(res)>0)
      {
        allRelatedConditionalGenesMeasures[[1]]<-res
        names(allRelatedConditionalGenesMeasures)[[1]]<-allgenes[[iIndex]]
        return(allRelatedConditionalGenesMeasures)
      }else
      {
        return(NULL)
      }
    },allgenes,reducedCube,temporal,completeState)
  }

  if(useParallel)
  {
    res<-doParallelWork(internalloopByWhole,targetgenes,allgenes,clusters,reducedCube,maxK,temporal,getCurrentState,getpreviousState)

  }else
  {
    res<-doNonParallelWork(internalloopByWhole,targetgenes,allgenes,clusters,reducedCube,maxK,temporal,getCurrentState,getpreviousState)
  }


  if(!is.null(res)&length(res)>0)
  {
    cond1<-sapply(res,function(entry)!is.null(entry))
    res<-(res[cond1][unlist(lapply(res[cond1],length)!=0)])
    class(res)<-c("FBNCube")
  }

  if(useParallel)
  {
    closeAllConnections()
  }
  time2<-as.numeric(Sys.time())
  print(paste("Total cost ", time2-time1, " seconds to construct gene cube", sep='',collapse = ''))
  return(res)
}
