#source http://www.sthda.com/english/articles/25-cluster-analysis-in-r-practical-guide/111-types-of-clustering-methods-overview-and-quick-start-r-code/
#install.packages("factoextra")
#install.packages("cluster")
#install.packages("magrittr")
require("cluster")
require("factoextra")
require("magrittr")

##kmeansParameters<-list(type="kmeans",numOfClusters=10,nstart=100,iter.max=1000)
##class(kmeansParameters)="ModelParameter"
##hierarchicalParameters<-list(type="hierarchical",distmethod="euclidean",hclustmethod="ward.D2")
##class(hierarchicalParameters)="ModelParameter"
##nbclustParameters<-list(type="nbclust",distmethod="euclidean",min.nc=2,max.nc=10,nbmethod="complete",nbindex="all")
##class(nbclustParameters)="ModelParameter"
##fuzzyParameters<-list(type="fanny",metric = "euclidean", stand = FALSE)
##class(fuzzyParameters)="ModelParameter"
#' @export
clusterTimeSeries <- function(timeSeriesData, method=c("kmeans","hierarchical","diana","fanny","nbclust"),modelParameters)
{
  if (is.matrix(timeSeriesData))
    clusterData <- timeSeriesData
  else
    stop("The type of timeSeriesData must be matrix")

  #switch between the different methods
  switch(match.arg(method),
         kmeans={
           set.seed(123)
           if (!inherits(modelParameters,"ModelParameter"))
             stop("modelParameters must be inherited from ModelParameter")

           if(modelParameters$type!="kmeans")
             stop("The type of the ModelParameter cannot be applied to kmeans")

           numOfCluster<-modelParameters$numOfClusters
           if(is.na(numOfCluster) | is.null(numOfCluster))
             stop("The parameter numOfCluster is missing")

           nstart<-modelParameters$nstart
           iter.max<-modelParameters$iter.max

           thisCluster<-clusterData %>%
             na.omit() %>%  # Remove missing values (NA)
             scale() %>% # standardize variables
             kmeans(numOfCluster,nstart=nstart,iter.max=iter.max)

           #visualize the cluster
           fviz_cluster(thisCluster, clusterData,
                        ellipse.type = "convex",
                        palette = "jco",
                        ggtheme = theme_minimal())

           return(list(type="kmeans",cluster=thisCluster))
         },
         hierarchical={ #It does not require to pre-specify the number of clusters to be generated.
           if (!inherits(modelParameters,"ModelParameter"))
             stop("modelParameters must be inherited from ModelParameter")

           if(modelParameters$type!="hierarchical")
             stop("The type of the ModelParameter cannot be applied to hierarchical")

           distmethod<-modelParameters$distmethod
           hclustmethod<-modelParameters$hclustmethod
           numOfCluster<-modelParameters$numOfClusters
           if(is.na(numOfCluster) | is.null(numOfCluster))
             stop("The parameter numOfCluster is missing")

           thisCluster<- clusterData %>%
             na.omit() %>%          # Remove missing values (NA)
             scale() %>%                    # Scale the data # standardize variables
             dist(method = distmethod) %>% # Compute dissimilarity matrix
             hclust(method = hclustmethod)     # Compute hierachical clustering

           clusterCut <- cutree(thisCluster, numOfCluster)

           # Visualize using factoextra
           # Cut in 4 groups and color by groups
           library(factoextra)
           fviz_dend(thisCluster, k = numOfCluster, # Cut in four groups
                     cex = 0.5, # label size
                     k_colors = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
                     color_labels_by_k = TRUE, # color labels by groups
                     rect = TRUE # Add rectangle around groups
           )

           return(list(type=match.arg(method),cluster=thisCluster,clusterCut=clusterCut))

         },
         diana={
           if (!inherits(modelParameters,"ModelParameter"))
             stop("modelParameters must be inherited from ModelParameter")

           if(modelParameters$type!="diana")
             stop("The type of the ModelParameter cannot be applied to diana")

           numOfCluster<-modelParameters$numOfClusters
           if(is.na(numOfCluster) | is.null(numOfCluster))
             stop("The parameter numOfCluster is missing")

           # Compute diana()
           library(cluster)
           thisCluster<- clusterData %>%
             na.omit() %>%          # Remove missing values (NA)
             scale() %>%                    # Scale the data # standardize variables
             diana(stand = TRUE)

           clusterCut <- cutree(thisCluster, numOfCluster)

           # Plot the dendrogram
           library(factoextra)
           fviz_dend(thisCluster, cex = 0.5,
                     k = numOfCluster, # Cut in four groups
                     palette = "jco" # Color palette
           )

           return(list(type=match.arg(method),cluster=thisCluster,clusterCut=clusterCut))
         },
         fanny={
           #fuzzy cluster
           library(cluster)
           if (!inherits(modelParameters,"ModelParameter"))
             stop("modelParameters must be inherited from ModelParameter")

           if(modelParameters$type!="fanny")
             stop("The type of the ModelParameter cannot be applied to fanny")

           metric<-modelParameters$metric
           numOfCluster<-modelParameters$numOfClusters
           stand <-modelParameters$stand #True or false
           thisCluster<- clusterData %>%
             na.omit() %>%          # Remove missing values (NA)
             scale() %>%                    # Scale the data # standardize variables
             fanny(k=numOfCluster,metric = metric, stand = stand)  # Compute fuzzy clustering with k = 2

           library(factoextra)
           fviz_cluster(thisCluster, ellipse.type = "norm", repel = TRUE,
                        palette = "jco", ggtheme = theme_minimal(),
                        legend = "right")
           return(list(type=match.arg(method),cluster=thisCluster))

         },
         nbclust={ #Determining the optimal number of clusters
           if (!inherits(modelParameters,"ModelParameter"))
             stop("modelParameters must be inherited from ModelParameter")

           if(modelParameters$type!="nbclust")
             stop("The type of the ModelParameter cannot be applied to nbclust")

           require("NbClust")
           distmethod<-modelParameters$distmethod #distance method like euclidean
           min.nc<-modelParameters$min.nc
           max.nc<-modelParameters$max.nc
           nbmethod<-modelParameters$nbmethod
           nbindex<-modelParameters$nbindex
           thisCluster <- clusterData %>%
             na.omit() %>%          # Remove missing values (NA)
             scale() %>% # standardize variables
             NbClust(distance = distmethod,
                     min.nc = min.nc, max.nc = max.nc,
                     method = nbmethod, index =nbindex)
           library(factoextra)
           fviz_nbclust(thisCluster, ggtheme = theme_minimal())

           return(list(type=match.arg(method),cluster=thisCluster))

         },
         mclust={
           #model based cluster

         },
         stop("'method' must be one of \"kmeans\",\"edgeDetector\",\"scanStatistic\"")
  )
}

#' @export
fuzzyTreeCluster<-function(requireddata, minElements,maxElements,useParallel=FALSE)
{
  if (is.matrix(requireddata))
    fullData <- requireddata
  else if(is.list(requireddata))
    # in case of list, paste all matrices before clustering
  {
    fullData<-do.call(cbind,requireddata)
  }else
  {
    stop("The type of time seriesdata is not supported")
  }

  rowValues<-rownames(fullData)

  fuzzyParameters<-list(type="fanny",metric = "euclidean", stand = FALSE,numOfClusters=2)
  class(fuzzyParameters)="ModelParameter"

  clusterobject<-clusterTimeSeries(fullData,method="fanny",fuzzyParameters)

  clusters<-clusterobject$cluster$clustering

  fac<-factor(clusters)
  dat <- data.frame(cluster=clusters)
  groups<-split.data.frame(dat,dat$cluster)

  groupsdata<-list()
  #ToDo: calculate group data
  #
  groupsdata<-lapply(groups,function(g,fullData){
    gnames<-rownames(g)
    return(fullData[rownames(fullData) %in% gnames,])
  },fullData)
  #load configulation
  readConfigFile()

  internalloop<-function(i,groupsdata,minElements,maxElements,fuzzyParameters){


      subgroupData<-groupsdata[[i]]
      res<-list()
      res[[i]]<-list()
      rowValues<-rownames(subgroupData)
      lenOfsub<-length(rowValues)

      if(lenOfsub<=maxElements)
      {
        res[[i]]<-subgroupData
        names(res)[[i]]<-"leaf"
        return(res)
      }

      tryCatch(
      {
        clusterobject<-clusterTimeSeries(subgroupData,method="fanny",fuzzyParameters)
        clusters<-clusterobject$cluster$clustering


        fac<-factor(clusters)
        dat <- data.frame(cluster=clusters)
        groups<-split.data.frame(dat,dat$cluster)

        #ToDo: calculate group data
        #
        groupsdata<-lapply(groups,function(g,subgroupData){
          gnames<-rownames(g)
          return(subgroupData[rownames(subgroupData) %in% gnames,])
        },subgroupData)

        res[[i]]<-doNonParallelWork(internalloop,groupsdata,minElements,maxElements,fuzzyParameters)
        names(res)[[i]]<-"sub"
        return(res)
      },
      error = function(e)
      {
        emas<-e
        stop(sprintf("Error executing FBN Cluster: %s", e$message))
      })
  }

  time1<-as.numeric(Sys.time())

  res<-list()

  if(length(rowValues)==0)
  {
    return(list())
  }
  if(useParallel)
  {
    res<-doParallelWork(internalloop,groupsdata,minElements,maxElements,fuzzyParameters)

  }else
  {
    res<-doNonParallelWork(internalloop,groupsdata,minElements,maxElements,fuzzyParameters)
  }

  cond1<-sapply(res,function(entry)!is.null(entry))
  res<-(res[cond1][unlist(lapply(res[cond1],length)!=0)])



  if(useParallel)
  {
    closeAllConnections()
  }

  time2<-as.numeric(Sys.time())
  print(paste("Total cost ", time2-time1, " seconds to run FuzzyTree algorithm ", sep='',collapse = ''))
  return(res)
}
