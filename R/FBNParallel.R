#The follow functions should be used internally and none of them should be exposed to outside user.
#'Executing a list of processes in parallel
#'
#'@param parallelFuc A process that will be run in parallel
#'@param listitems The main list for each process
#'@param ... other parameters
#'@return No return
#'@examples
#'coming later
#' @export
doParallelWork <- function(parallelFuc,listitems,...)
{
  if(!is.function(match.fun(parallelFuc)))
  {
    stop("The parameter is not a function")
  }

  tryCatch({
    require(parallel)
    require(foreach)
    require(doParallel)
    cl <- makeCluster(detectCores()[1] - 1)
    registerDoParallel(cl)
    #* Do the work .....**
    res<-foreach(k = 1:length(listitems),.combine=c,.packages = c("foreach","FBNNet"))%dopar% {
      try({
        result<-parallelFuc(k,listitems,...)
        #reclaim memory
        gc()
        return(result)
      })
    }
    stopCluster(cl)
    #closeAllConnections()
    cond1<-sapply(res,function(entry)!is.null(entry))
    #remove the unwant outer list
    #return(unlist(res[!is.null(res)],recursive = FALSE))
    return(res[cond1][unlist(lapply(res[cond1],length)!=0)])
  },
  error = function(e)
  {
    emas<-e
    stop(sprintf("Error doParallelWork: %s", e$message))
  })

}

#'Executing a list of processes not in parallel
#'
#'@param parallelFuc A process that will be run in parallel
#'@param listitems The main list for each process
#'@param ... other parameters
#'@return No return
#'@examples
#'coming later
#' @export
doNonParallelWork <- function(parallelFuc,listitems,...)
{
  #foreach environment, the option stringsAsFactors is required to be set to FALSE
  options(stringsAsFactors = FALSE)
  if(!is.function(match.fun(parallelFuc)))
  {
    stop("The parameter is not a function")
  }
  tryCatch({
    require(parallel)
    require(foreach)
    require(doParallel)
    res<-(foreach(k = 1:length(listitems),.combine=c)%do% {
      result<-parallelFuc(k,listitems,...)
      #reclaim memory
      #gc()
      return(result)
    })
    cond1<-sapply(res,function(entry)!is.null(entry))
    #remove the unwant outer list
    return(res[cond1][unlist(lapply(res[cond1],length)!=0)])
  },
  error = function(e)
  {
    emas<-e
    stop(sprintf("Error doNonParallelWork: %s", e$message))
  })
}

#'Executing a list of processes in parallel with decreasing list
#'
#'@param parallelFuc A process that will be run in parallel
#'@param listitems The main list for each process
#'@param unprocessedListitems The remain list items that haven't been processed yet
#'@param ... other parameters
#'@return No return
#'@examples
#' coming later
#' @export
doNonParallelWorkDecrease <- function(parallelFuc,listitems,unprocessedListitems,...)
{
  #foreach environment, the option stringsAsFactors is required to be set to FALSE
  options(stringsAsFactors = FALSE)
  if(!is.function(match.fun(parallelFuc)))
  {
    stop("The parameter is not a function")
  }
  tryCatch({
    require(parallel)
    require(foreach)
    require(doParallel)
    res<-(foreach(k = 1:length(listitems),.combine=c)%do% {

      reindex<-which(listitems[k]%in%unprocessedListitems)
      unprocessedListitems2<-unprocessedListitems[-reindex]

      result<-parallelFuc(k,listitems,unprocessedListitems2,...)

      if(!is.null(result))
      {
        unprocessedListitems<-unprocessedListitems2
      }
      #reclaim memory
      #gc()
      return(result)
    })
    cond1<-sapply(res,function(entry)!is.null(entry))
    #remove the unwant outer list
    return(res[cond1][unlist(lapply(res[cond1],length)!=0)])
  },
  error = function(e)
  {
    emas<-e
    stop(sprintf("Error doNonParallelWork: %s", e$message))
  })
}
unregister <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
