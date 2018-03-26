
#state should be a vector with names corresponding to genes
#genesOn is a vector of genes that are marked as on
#genesOff is a vector of genes that are marked as off
#the indexes of genes must be corresponding to the state
#startstates<-lapply(trainingseries,function(x)x[,1])
#' @export
searchForAttractors<-function(fbnnetwork,startStates=list(),genes,timedecay=1, type=c("synchronous","asynchronous"), genesOn=c(),genesOff=c(),maxSearch=1000)
{
  if (!(inherits(fbnnetwork,"FundamentalBooleanNetwork")))
    stop("Network must be inherited from FundamentalBooleanNetwork")

  type<-match.arg(type,c("synchronous","asynchronous"))

  #random generated test data
  if(length(startStates)==0)
  {
    startStates[[length(startStates)+1]]<-sample(0:1, length(genes), replace=T,prob=c(0.5,0.5))
  }

  if(length(genesOn)>0 & is.vector(genesOn))
  {
    genesOnIndexes<-which(genes%in%genesOn)
    lapply(startStates,function(state)state[genesOnIndexes]<-1)
  }

  if(length(genesOff)>0 & is.vector(genesOff))
  {
    genesOffIndexes<-which(genes%in%genesOff)
    lapply(startStates,function(state)state[genesOffIndexes]<-0)
  }
  resultList<-list()
  basinStates<-list()
  stateAssighed<-list()
  tempbasinStates<-list()
  tryCatch(
    {
        for(s in seq_along(startStates))
        {
          iniState<-startStates[[s]]
          searchedStates<-list()
          names(iniState)<-genes

          if(length(stateAssighed)>0)
          {
            if(Position(function(x) identical(x, iniState), stateAssighed, nomatch = 0)>0)
            {
              next
            }
          }

          found<-FALSE
          maxS<-1
          searchedStates[[length(searchedStates)+1]]<-iniState
          tempbasinStates[[length(tempbasinStates)+1]]<-iniState
          decayIndex<-c()
          timestepTrack<-list()
          while(!found & maxS<=maxSearch)
          {
            result<-getFBMSuccessor(fbnnetwork,iniState,genes,type=type,decayIndex,timestepTrack)
            nextState<-result$nextState
            decayIndex<-result$decayIndex
            timestepTrack<-result$timestepTrack
            names(nextState)<-genes
            #check if an attractor is found
            if(length(stateAssighed)>0)
            {
              if(Position(function(x) identical(x, nextState), stateAssighed, nomatch = 0)>0)
              {
                break
              }
            }
            #correction
            if(length(basinStates)>0)
            {
              foundI<-0
              lapply(seq_along(basinStates),function(statesIndex){
                states<-basinStates[[statesIndex]]
                if(Position(function(x) identical(x, nextState), states, nomatch = 0)>0)
                {
                  foundI<-statesIndex
                }
              })
              if(length(tempbasinStates)>0 & foundI>0)
              {
                basinStates[[foundI]]<-unique(dissolve(list(dissolve(basinStates[[foundI]]),tempbasinStates)))
                break
              }
            }

            # find vector in list of vectors
            searchedIndex<-Position(function(x) identical(x, nextState), searchedStates, nomatch = 0)

            if(searchedIndex>0)
            {
              found<-TRUE
              #check if the attractor has been found
              attractors<-lapply(resultList,function(result)result[[1]])
              if(!Position(function(x) identical(x, nextState), attractors, nomatch = 0)>0)
              {
                currentlength<-length(resultList)+1
                resultList[[currentlength]]<-list()
                #resultList[[currentlength]][[1]]<-nextState

                resultList[[currentlength]]<-dissolve(list(dissolve(searchedStates[searchedIndex:length(searchedStates)]),nextState))
                #resultList[[currentlength]][[3]]<-basinStates
                basinStates[[currentlength]]<-list()
                temp<-dissolve(tempbasinStates)
                basinStates[[currentlength]]<-temp[!temp%in%resultList[[currentlength]]]
                names(basinStates)[[currentlength]]<-currentlength
                stateAssighed<-unique(dissolve(list(stateAssighed,resultList[[currentlength]])))
              }
            }else
            {
              found<-FALSE
              searchedStates[[length(searchedStates)+1]]<-nextState
              tempbasinStates[[length(tempbasinStates)+1]]<-nextState
            }
            iniState<-nextState

            maxS<-maxS+1
          }
        }
      res<-list()
      res[[1]]<-resultList
      names(res)[[1]]<-"Attractors"
      res[[2]]<-genes
      names(res)[[2]]<-"Genes"
      res[[3]]<-basinStates
      names(res)[[3]]<-"BasinOfAttractor"
      class(res)<-c("FBMAttractors")
      return(res)
    },
    error = function(e)
    {
      mes<-e$message
      stop(sprintf("Error executing FBN model: %s", e$message))
    })
}



networkFixUpdate<-function (network, fixIndices, values)
{
  if (!(inherits(fbnnetwork,"FundamentalBooleanNetwork")))
    stop("Network must be inherited from FundamentalBooleanNetwork")

  if (length(fixIndices) != length(values) && length(values) !=
      1)
    stop("fixIndices and values must have the same number of elements, or values must have 1 element!")

  if (any(is.na(network$fixed[fixIndices])))
    stop("fixIndices contains invalid indices!")

  if (any(values != 0 & values != 1 & values != -1))
    stop("Please supply only 0, 1, or -1 in values!")

  network$fixed[fixIndices] <- as.integer(values)

  return(network)
}
