#' @export
discreteTimeSeries <- function(timeSeriesData, method=c("average"))
{
  res<-timeSeriesData

  if (is.matrix(timeSeriesData))
    fullData <- timeSeriesData
  else if(is.list(timeSeriesData))
    # in case of list, paste all matrices before clustering
  {
    fullData<-do.call(cbind,timeSeriesData)
  }else
  {
    stop("The type of time seriesdata is not supported")
  }

  #use average as a thread to discrete time series
  switch(match.arg(method),
         average={

            fullData<-abs(fullData)+1 #add 1 to avoid 0
            rowmean<-t(rowMeans(fullData))
            colmean<-t(rowmean)
           if (is.matrix(res))
           {
             numofcol<-ncol(fullData)
             matrixMean<-rep.col(colmean,numofcol)
             colNames<-colnames(fullData)
             rowNames<-rownames(fullData)
             binarizedTimeSeries <- sapply(as.data.frame(fullData-(fullData%%matrixMean)), function(vr)as.numeric(vr>0))
             rownames(binarizedTimeSeries)<-rowNames
             colnames(binarizedTimeSeries)<-colNames
           }
           else if(is.list(res))
             # in case of list, paste all matrices before clustering
             # all matrixs in the list must have the same dimensions
           {
             binarizedTimeSeries <- lapply(res,function(m){
                   numofcol<-ncol(m)
                   matrixMean<-rep.col(colmean,numofcol)
                   fullData<-abs(m)+1 #add 1 to avoid 0
                   colNames<-colnames(m)
                   rowNames<-rownames(m)
                   sub_res<-sapply(as.data.frame(fullData-(fullData%%matrixMean)), function(vr)as.numeric(vr>0))
                   rownames(sub_res)<-rowNames
                   colnames(sub_res)<-colNames
                   return(sub_res)
                 }
               )
           }
        },
        stop("'method' must be one of \"average\",\"kmeans\"")
  )
}

#' @export
dividedintosmallgroups<-function(timeSeriesData)
{
  if (is.matrix(timeSeriesData))
    fullData <- timeSeriesData
  else if(is.list(timeSeriesData))
    # in case of list, paste all matrices before clustering
  {
    fullData<-do.call(cbind,timeSeriesData)
  }else
  {
    stop("The type of time seriesdata is not supported")
  }
  fuzzgroups<-fuzzyTreeCluster(fullData,10,30)
  return(fuzzgroups)
}

#' @export
dividedDiscreteDataintosmallgroups<-function(originalTimeSeriesData,discretedTimeSeriesdata)
{
  findAllLeaves<-function(treegoups)
  {
    res1<-list()
    len<-length(treegoups)

    nm <- names(treegoups)[1]
    if(nm=="leaf")
    {
      res1[length(res1)+1]<-treegoups[1]
    }
    else
    {
      newlen<-length(res1)+1
      res1[[newlen]]<-list()
      res1[[newlen]]<-findAllLeaves(treegoups[[1]])
    }

    if(len==2)
    {
      nm <- names(treegoups)[2]
      if(nm=="leaf")
      {
        res1[length(res1)+1]<-treegoups[2]
      }
      else
      {
        newlen2<-length(res1)+1
        res1[[newlen2]]<-list()
        res1[[newlen2]]<-dissolve(findAllLeaves(treegoups[[2]]))
      }
    }

    return(res1)
  }

  fuzzygroups<-dividedintosmallgroups(originalTimeSeriesData)

  groups<-dissolve(findAllLeaves(fuzzygroups))

  res<-lapply(groups,function(subgroup,discretedTimeSeriesdata){
      subnames<-rownames(subgroup)
      res<-lapply(discretedTimeSeriesdata,function(mtx,subnames)
        {
          res1<-mtx[rownames(mtx)%in%subnames,]
          return(res1)
        },subnames)
      return(res)
  },discretedTimeSeriesdata)
  names(res)<-c(1:length(res))
  class(res)<-"ClusteredTimeseriesData"
  return(res)
}

