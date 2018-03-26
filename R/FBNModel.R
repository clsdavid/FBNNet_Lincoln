

#'This method is used to reconstruct time series data
#'
#'@param fbnCube A pre constructed Orchard cube
#'@param fbnnetwork An object of FBNNetwork
#'@param timeseriesCube An original time series data that contains samples, genes and time points. It can be as list of samples contain genes and the initial time point.
#'@param interval An interval is the time gap between two time points
#'@return A list object that contains reconstructed time series and FBN network
#'@examples
#' mat1<-matrix(c("1","0","0","1","0","0","0","1","1"),3,3, dimnames=list(c("gene1","gene2","gene3"),c("1","2","3")))
#' mat2<-matrix(c("1","1","0","1","0","1","1","1","0"),3,3, dimnames=list(c("gene1","gene2","gene3"),c("1","2","3")))
#' listtestOriginal<-list(mat1,mat2)
#' cube<-constructFBNCube(c("gene1","gene2"),c("gene1","gene2","gene3"),listtestOriginal,4,1,FALSE)
#' network<-mineFBNNetwork(cube,c("gene1","gene2"))
#' mat1<-matrix(c("1","1","0"),3,1, dimnames=list(c("gene1","gene2","gene3"),c("1")))
#' mat2<-matrix(c("1","1","1"),3,1, dimnames=list(c("gene1","gene2","gene3"),c("1")))
#' listtestInitial<-list(mat1,mat2)
#' result<-reconstructTimeseries(cube,network,listtestInitial,1,1,3)
#' @export
reconstructTimeseries<-function(fbnCube,fbnnetwork,timeseriesCube,interval=1,maxTimepoints=100)
{

  if(!is.numeric(interval) | interval<=0L | !all.equal(interval, as.integer(interval)))
  {
    stop("interval must be integer")
  }

  if (!(inherits(fbnnetwork,"FundamentalBooleanNetwork")))
    stop("Network must be inherited from FundamentalBooleanNetwork")

  genes<-fbnnetwork$genes
  res<-list()
  #code
  tryCatch(
    {
        reconstructed<-lapply(timeseriesCube,function(timeseries){
          return (reconstructMissingTimepoints(interval,timeseries,fbnnetwork,genes))
        })


        res[[1]]<-reconstructed
        names(res)[[1]]<-"reconstructed"
        res[[2]]<-fbnCube
        names(res)[[2]]<-"FBN"
    },
    error = function(e)
    {
      emas<-e
      stop(sprintf("Error executing FBN model: %s", e$message))
    })

  return (res)
}

#'A random selection function denoted as P[]
#'
#'@param probability A probability of an event
#'@return TRUE or FALSE indicating whether or not the event has been selected
#'@examples
#' coming later
#' @export
randomSelection<-function(probability)
{
  tryCatch(
    {
      #index<-1L
      #while(index<=10L)
      #{
      set.seed(1)
      res<-sample(c(FALSE,TRUE), size=1, replace=TRUE, prob=c(1-probability,probability))
      #index<-index+1L
      #}
    },
    error = function(e)
    {
      pissue<-probability
      mes<-e$message
      stop(sprintf("Error randomSelection: %s, probability=%s", e$message,probability))
    })

  return(res)
}

reconstructMissingTimepoints<-function(interval,targetSeriesMatrix,fbnNetwork,genes)
{

  if(!is.matrix(targetSeriesMatrix))
  {
    stop("The target series must be a type of matrix i.e. a table has two dimensions")
  }

  tryCatch(
    {
      numrow<-dim(targetSeriesMatrix)[[1]]
      rowNames<-rownames(targetSeriesMatrix)
      originalColNames<-colnames(targetSeriesMatrix)
      colNames<-as.numeric(originalColNames)
      if(all(is.na(colNames)))
      {
        colnames(targetSeriesMatrix)<-1:length(colNames)
        colNames<-as.numeric(colnames(targetSeriesMatrix))
      }
      i<-1L
      recCol<-c(colNames[[i]])
      while(i<length(colNames))
      {
        j<-colNames[[i]]
        while(j<colNames[[i+1]])
        {
          j<-j+interval
          recCol<-c(recCol,j)
        }
        i<-i+1
      }
      mat<-matrix(0,nrow=numrow,ncol=length(recCol),byrow=FALSE,dimnames<-list(rowNames,recCol))
      eColNames<-colnames(mat)
      #get initial state
      mat[,1]<-targetSeriesMatrix[,1]
      k<-2
      premat<-mat[,k-1]
      decayIndex<-c()
      timestepTrack<-list()
      while(k<=length(eColNames))
      {
        nextState<-getFBMSuccessor(fbnNetwork,premat,rowNames,type="synchronous",decayIndex,timestepTrack)
        mat[,k]<-nextState$nextState
        premat<-mat[,k]
        decayIndex<-nextState$decayIndex
        timestepTrack<-nextState$timestepTrack
        k<-k+1
      }

      return(mat)
    },
    error = function(e)
    {
      mes<-e$message
      stop(sprintf("Error executing FBN model: %s", e$message))
    })
}

#'This method is used to calculate the next state
#'
#'@param fbnnetwork An object of FBNNetwork
#'@param currentState A vector of current gene state
#'@param genes a list of genes which index order must match with the current state
#'@param type A type of Boolean network update schema choosen from synchronous, asynchronous and time step based
#'@param timedecay An value indicates the period of time when to degrade an activated gene if no activators presented. It is usually one time step
#'@param timestepTrack A list of processed time steps for each FBM function. This parameter is used with the type of timestep
#'@return A list object that contains reconstructed time series and FBN network
#'@examples
#' require(BoolNet)
#' require(FBNNet)
#' network<-loadNetwork("testthat/others/cellcycle.txt")
#' trainingseries<-FBNDataReduction(generateTimeSeries(network,2000,43))
#' cube<-constructFBNCube(network$genes,network$genes,trainingseries,4,1,TRUE)
#' NETWORK2<-mineFBNNetwork(cube,network$genes)
#' state<-c("0","1","1","0","1","1","1","1","0","0")
#' names(state)<-c("CycD","Rb","E2F","CycE","CycA","p27","Cdc20","Cdh1","UbcH10","CycB")
#' getFBMSuccessor(NETWORK2,state,names(state))
#' @export
getFBMSuccessor<-function(fbnNetwork,currentState,genes,type=c("synchronous","asynchronous"),decayIndex=c(),timestepTrack=list())
{
  internalFun<-function(gene,interactions,geneState,currentState,timedecay,decayIndex=1,timestepTrack=list())
  {
    genefunctions<-interactions[[gene]]
    ini<-currentState[[gene]]==1L
    decay<-timedecay
    if(decay<1)decay<-1L
    if(length(genefunctions)>0)
    {
      #find all activators' probabilities
      condOfActivation <- sapply(genefunctions, function(activator) activator$type == 1L)
      funcOfActivators<-genefunctions[condOfActivation]
      #find all inhibitors' probabilities
      condOfInhibitors<- sapply(genefunctions, function(inhibitor) inhibitor$type == 0L)
      funcOfInhibitors<-genefunctions[condOfInhibitors]
    }else
    {
      funcOfActivators<-list()
      funcOfInhibitors<-list()
    }

    prFA<-FALSE
    probabilityFA<-0.0
    if(length(funcOfActivators)>0)
    {
      for(fa in seq_along(funcOfActivators))
      {
        fbnName<-names(funcOfActivators)[[fa]]

        if(is.null(timestepTrack[[fbnName]]))
        {
          newindex<-length(timestepTrack)+1L
          timestepTrack[[newindex]]<-1L
          names(timestepTrack)[[newindex]]<-fbnName
        }

        if(as.numeric(timestepTrack[[fbnName]])>=as.numeric(funcOfActivators[[fa]]$timestep))
        {
          pregeneInput<-dissolve(lapply(funcOfActivators[[fa]]$input,function(geneindex){
            res<-list()
            res[[1]]<-currentState[[geneindex]]
            names(res)[[1]]<-genes[[geneindex]]
            return (res)
          }))

          probability<-getProbabilityFromCubeBasedOnInput(gene,1,funcOfActivators[[fa]]$expression,funcOfActivators[[fa]]$probability,pregeneInput)

          if(!length(probability)==0)
          {
            prFA<-prFA | randomSelection(probability)
          }
          timestepTrack[[fbnName]]<-1L
        }else
        {
          timestepTrack[[fbnName]]<-as.numeric(timestepTrack[[fbnName]])+1L
        }
      }
    }

    prFD<-FALSE
    probabilityFD<-0.0
    if(length(funcOfInhibitors)>0)
    {
      for(fd in seq_along(funcOfInhibitors))
      {
        fbnName<-names(funcOfInhibitors)[[fd]]
        if(is.null(timestepTrack[[fbnName]]))
        {
          newindex<-length(timestepTrack)+1L
          timestepTrack[[newindex]]<-1L
          names(timestepTrack)[[newindex]]<-fbnName
        }

        if(as.numeric(timestepTrack[[fbnName]])>=as.numeric(funcOfInhibitors[[fd]]$timestep))
        {
          pregeneInput<-dissolve(lapply(funcOfInhibitors[[fd]]$input,function(geneindex){
            res<-list()
            res[[1]]<-currentState[[geneindex]]
            names(res)[[1]]<-genes[[geneindex]]
            return (res)
          }))

          probability<-getProbabilityFromCubeBasedOnInput(gene,0,funcOfInhibitors[[fd]]$expression,funcOfInhibitors[[fd]]$probability,pregeneInput)

          if(!length(probability)==0)
          {
            prFD<-prFD| randomSelection(probability)
          }
          timestepTrack[[fbnName]]<-1L
        }else
        {
          timestepTrack[[fbnName]]<-as.numeric(timestepTrack[[fbnName]])+1L
        }
      }
    }

    #nothing happen then deactivate when the decay time is due
    if(decay>0)
    {
      if(prFA | prFD)
      {
        decayIndex<-1L
      }else
      {
        #check decay
        if(as.numeric(decayIndex)>=as.numeric(decay))
        {
          ini<-FALSE
          decayIndex<-1L
        }else
        {
          decayIndex<-as.numeric(decayIndex)+1L
        }
      }
    }

    result<-(ini|prFA)&(!prFD)
    res<-list()
    res[[1]]<-result
    names(res)[[1]]<-"nextState"
    res[[2]]<-decayIndex
    names(res)[[2]]<-"decayIndex"
    res[[3]]<-timestepTrack
    names(res)[[3]]<-"timestepTrack"
    #result<-(prFA)&(!prFD)
    #nextState<-c(nextState,as.numeric(result))
    return(res)
  }
  interactions <-fbnNetwork$interactions
  names(currentState)<-genes

  #decayIndex<-0L
  nextState<-c()
  fixedgeneIndex<-which(fbnNetwork$fixed!=-1)
  fixedgenes<-fbnNetwork$genes[fixedgeneIndex]
  type<-match.arg(type,c("synchronous","asynchronous"))

  if(length(decayIndex)==0)
  {
    decayIndex<-rep(1L,length(genes))
    names(decayIndex)<-genes
  }
  if(type=="asynchronous")
  {
    nonfixedgenes<-genes[!genes%in%fixedgenes]
    j<-0
    while(j <=10)
    {
      gene<-sample(nonfixedgenes, 1)
      j<-j+1
    }
    timedecay<-fbnNetwork$timedecay[[gene]]
    ini<-currentState[[gene]]==1L
    result<-internalFun(gene,interactions,ini,currentState,timedecay,decayIndex[[gene]],timestepTrack)
    decayIndex[[gene]]<-result[[2]]
    timestepTrack<-result[[3]]
    nextState<-currentState
    #update next state for the random choosen gene
    nextState[[gene]]<-as.numeric(result[[1]])
  }else
  {
    for(j in seq_along(genes))
    {
      gene<-genes[[j]]
      ini<-currentState[[gene]]==1L
      timedecay<-fbnNetwork$timedecay[[gene]]
      #if gene is fixed then get the current gene state
      if(gene%in%fixedgenes)
      {
        nextState<-c(nextState,as.numeric(ini))
        next
      }
      result<-internalFun(gene,interactions,ini,currentState,timedecay,decayIndex[[gene]],timestepTrack)
      decayIndex[[gene]]<-result[[2]]
      timestepTrack<-result[[3]]
      nextState<-c(nextState,as.numeric(result[[1]]))
    }
  }
  res<-list()
  res[[1]]<-nextState
  names(res)[[1]]<-"nextState"
  res[[2]]<-decayIndex
  names(res)[[2]]<-"decayIndex"
  res[[3]]<-timestepTrack
  names(res)[[3]]<-"timestepTrack"
  return(res)
}
