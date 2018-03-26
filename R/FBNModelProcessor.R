
#'A benchmark type function to test a complete process of FBN model
#'
#'@param booleanNetwork A Boolean Network object
#'@param shorttrainingtimeSeries A training timeseries matrix data set
#'@param shorttestseries A test timeseries matrix data set
#'@param originaltestseries A original timeseries matrix data set
#'@param interval The internal gat between time points
#'@param decay The decay time step for protein decay
#'@param temporal A temporal setting with default =1
#'@param convertToFBN if it is a need to convert the boolean network to FBN Network. If the boolean network object is not a FBNNetwork object, then it must be set to TRUE
#'@param existingFBNCube An existing FBN cube object
#'@return A list of object
#'@examples
#' coming later
#' @export
FBNBenchmark<-function(booleanNetwork,shorttrainingtimeSeries,shorttestseries,originaltestseries,interval=1,temporal=1,convertToFBN=TRUE,existingFBNCube=NULL)
{
  #validate network types
  if (!(inherits(booleanNetwork,"BooleanNetworkCollection")))
  {
    stop("booleanNetwork must be inherited from BooleanNetwork")
  }

  if(!is.numeric(interval) | interval<=0L | !all.equal(interval, as.integer(interval)))
  {
    stop("interval must be integer")
  }

  if(convertToFBN)
  {
    network<-convertToFBNNetwork(booleanNetwork)
  }else
  {
    network<-booleanNetwork
  }

  if(is.null(existingFBNCube))
  {
    fbncube<-constructFBNCube(network$genes,network$genes,shorttrainingtimeSeries,3,temporal,FALSE)
  }else
  {
    fbncube<-existingFBNCube
  }

  result<-reconstructTimeseries(fbncube,network,shorttestseries,interval)

  #we don't need to include all observed time points except for the first one
  includedObservedTimepoint<-FALSE
  if(includedObservedTimepoint)
  {
    similar<-checkSimilarity(originaltestseries,result$reconstructed)
  }else
  {

    verifiedfiles<-lapply(1:length(result$reconstructed),function(index)
      {
        ocolnames<-colnames(originaltestseries[[index]])
        if(all(!is.numeric(ocolnames)))
        {
          ocolnames<-1:length(ocolnames)
        }
        rcolnames<-colnames(result$reconstructed[[index]])
        return(result$reconstructed[[index]][,which(rcolnames%in%ocolnames)])
        #which(eColNames==colNames[[i]]
     })
    similar<-checkSimilarity(originaltestseries,verifiedfiles)
  }

  cond1<-sapply(similar,function(entry)entry[[1]]=="similar")
  cond2<-sapply(similar,function(entry)entry[[1]]=="likely")
  cond3<-sapply(similar,function(entry)entry[[1]]=="verysimilar")
  cond4<-sapply(similar,function(entry)entry[[1]]=="unlikely")
  cond5<-sapply(similar,function(entry)entry[[1]]=="veryunlikely")

  res<-list()
  res[[1]]<-result
  res[[2]]<-similar


  res[[3]]<-similar[cond1]
  res[[4]]<-similar[cond2]
  res[[5]]<-similar[cond3]
  res[[6]]<-similar[cond4]
  res[[7]]<-similar[cond5]

  names(res)[[1]]<-"FBNResult"
  names(res)[[2]]<-"SimilarityReport"
  names(res)[[3]]<-"Similar"
  names(res)[[4]]<-"Likely"
  names(res)[[5]]<-"verysimilar"
  names(res)[[6]]<-"unlikely"
  names(res)[[7]]<-"veryunlikely"

  #benchmark result
  pmcond<-sapply(similar,function(entry)as.numeric(entry[[2]])==1)
  mmcond<-sapply(similar,function(entry)as.numeric(entry[[2]])<1)

  pm<-similar[pmcond]
  mm<-similar[mmcond]

  mmvalues<-unlist(lapply(mm,function(cond)as.numeric(cond[[2]])))
  avgMMvalues<-ifelse(length(mm)>0,sum(mmvalues)/length(mm),0)

  errorRate<-  ((1-avgMMvalues)*length(mm))/(length(pm)+length(mm))
  perfectMatchedRate<-length(pm)/(length(pm)+length(mm))
  missMatchedRate<-length(mm)/(length(pm)+length(mm))
  accurateRate<-(length(pm)+avgMMvalues*length(mm))/(length(pm)+length(mm))

  res[[8]]<-errorRate
  res[[9]]<-accurateRate
  res[[10]]<-missMatchedRate
  res[[11]]<-perfectMatchedRate
  names(res)[[8]]<-"ErrorRate"
  names(res)[[9]]<-"AccurateRate"
  names(res)[[10]]<-"MissMatchedRate"
  names(res)[[11]]<-"PrfectMatchedRate"
  return (res)
}

MiningFromShortTimeSeriesData<-function(trainingseries,testseries,validateseries,genes,haltrate=0,reconstructTimeInterval=1,maxTryRate=1,maxk=4,multithread=FALSE,timeporal=1,decay=1)
{

  cube<-constructFBNCube(genes,genes,trainingseries,maxk,timeporal,multithread)

  testErrorTry<-0

  bestResult<-NULL
  bestNetwork<-NULL
  bestsetting<-""
  tries<-1
  if(maxTryRate<testErrorTry)
    maxTryRate<-testErrorTry

  while(testErrorTry<=maxTryRate)
  {
    error_toleranceRate<<-testErrorTry #update global values
    print(paste("Execute ",tries, " times with error_toleranceRate=",error_toleranceRate, " and confidence_toleranceRate=",confidence_toleranceRate, sep="", collapse = ""))
    extractedNetwork<-mineFBNNetwork(cube,genes)
    resultfile<-FBNBenchmark(extractedNetwork,trainingseries,testseries,validateseries,reconstructTimeInterval,decay,timeporal,FALSE,cube)

    if(!is.null(bestResult))
    {
      if(resultfile$ErrorRate < bestResult$ErrorRate)
      {
        bestResult<-resultfile
        bestNetwork<-extractedNetwork
        bestsetting<-error_toleranceRate
      }
    }else
    {
      bestResult<-resultfile
      bestNetwork<-extractedNetwork
      bestsetting<-error_toleranceRate
    }
    print(paste("The best error rate is ",bestResult$ErrorRate, " And compare error rate is ",resultfile$ErrorRate , sep="", collapse = ""))
    if(bestResult$ErrorRate==0.0 | bestResult$ErrorRate<haltrate)
      break;
    tries<-tries+1
    #testConfidence<-testConfidence+0.1
    #}
    #testConfidence<-1 #reset
    testErrorTry<-testErrorTry+0.05
  }


  if(bestResult$ErrorRate<=0.0 | bestResult$ErrorRate<=haltrate | maxTryRate==testErrorTry)
  {
    cube<-constructFBNCube(genes,genes,bestResult$FBNResult$reconstructed,maxk,timeporal,multithread)
    extractedNetwork<-mineFBNNetwork(cube,genes)
    resultfile<-FBNBenchmark(extractedNetwork,trainingseries,testseries,testseries,reconstructTimeInterval,decay,timeporal,FALSE,cube)
    bestResponse<-list()
    bestResponse[[1]]<-bestNetwork
    bestResponse[[2]]<-bestResult
    bestResponse[[3]]<-bestsetting
    names(bestResponse)[[1]]<-"BestNetwork"
    names(bestResponse)[[2]]<-"BestResult"
    names(bestResponse)[[3]]<-"BestSetting"
  }else
  {
    bestResponse<-MiningFromShortTimeSeriesData(bestResult$FBNResult$reconstructed,testseries,validateseries,genes,haltrate,reconstructTimeInterval,maxTryRate,maxk,multithread,timeporal,decay)
  }

  return(bestResponse)
}

#' coming later
#' @export
MiningFromShortTimeSeriesDataUnsupervisionTwoStages<-function(trainingseries,testseries,genes,haltrate=0,reconstructTimeInterval=1,maxTryRate=1,maxk=4,multithread=FALSE,timeporal=1,decay=1)
{

  cube<-constructFBNCube(genes,genes,trainingseries,maxk,timeporal,multithread)

  testErrorTry<-0
  testminmutual_range<-0.5
  maxmutual_range<-1

  maxRules<-5
  maxFBNRulesOfActivator<<-maxRules
  maxFBNRulesOfInhibitor<<-maxRules
  bestResult<-NULL
  bestNetwork<-NULL
  bestsetting<-""
  tries<-1
  if(maxTryRate<testErrorTry)
    maxTryRate<-testErrorTry
  #Stage 1
  foundbest<-FALSE
  while(testErrorTry<=maxTryRate)
  {
    error_toleranceRate<<-testErrorTry #update global values
    while(testminmutual_range<=maxmutual_range)
    {
      minMutual_range<<-testminmutual_range

      print(paste("Execute ",tries, " times with error_toleranceRate=",error_toleranceRate, " and minimum minMutual_range=",minMutual_range, sep="", collapse = ""))
      extractedNetwork<-mineFBNNetwork(cube,genes)
      #print(extractedNetwork)
      resultfile<-FBNBenchmark(extractedNetwork,NULL,testseries,testseries,reconstructTimeInterval,decay,timeporal,FALSE,cube)

      if(!is.null(bestResult))
      {
        if(resultfile$ErrorRate < bestResult$ErrorRate)
        {
          print("The current best network")
          print(extractedNetwork)
          bestResult<-resultfile
          bestNetwork<-extractedNetwork
          bestsetting<-error_toleranceRate
        }
      }else
      {
        bestResult<-resultfile
        bestNetwork<-extractedNetwork
        bestsetting<-error_toleranceRate
      }
      print(paste("The best error rate is ",bestResult$ErrorRate, " And compare error rate is ",resultfile$ErrorRate , sep="", collapse = ""))
      if(bestResult$ErrorRate==0.0 | bestResult$ErrorRate<haltrate)
      {
        foundbest<-TRUE
        break;
      }
      tries<-tries+1
      testminmutual_range<-testminmutual_range+0.1
    }

    if(foundbest)
    {
      break;
    }
    testminmutual_range<-0.5
    #testConfidence<-testConfidence+0.1
    #}
    #testConfidence<-1 #reset
    testErrorTry<-testErrorTry+0.1
  }


  if(bestResult$ErrorRate<=0.0 &  error_toleranceRate==0 | bestResult$ErrorRate<=haltrate | maxTryRate==testErrorTry)
  {
    #output result
    cube<-constructFBNCube(genes,genes,bestResult$FBNResult$reconstructed,maxk,timeporal,multithread)
    extractedNetwork<-mineFBNNetwork(cube,genes)
    resultfile<-FBNBenchmark(extractedNetwork,NULL,testseries,bestResult$FBNResult$reconstructed,reconstructTimeInterval,decay,timeporal,FALSE,cube,FALSE)
    bestsetting<-c("error_toleranceRate"=error_toleranceRate,"maxFBNRulesOfActivator"=maxFBNRulesOfActivator,"maxFBNRulesOfInhibitor"=maxFBNRulesOfInhibitor)
    bestResponse<-list()
    bestResponse[[1]]<-extractedNetwork
    bestResponse[[2]]<-resultfile
    bestResponse[[3]]<-bestsetting
    names(bestResponse)[[1]]<-"BestNetwork"
    names(bestResponse)[[2]]<-"BestResult"
    names(bestResponse)[[3]]<-"BestSetting"
    print("Best Network generated:")
    print(bestResponse[[3]])
  }else
  {
    #Stage
    bestResponse<-MiningFromShortTimeSeriesDataUnsupervisionTwoStages(bestResult$FBNResult$reconstructed,testseries,genes,haltrate,reconstructTimeInterval,maxTryRate,maxk,multithread,timeporal,decay)
  }
  error_toleranceRate<<-0
  maxMutual_range<<-1
  minMutual_range<<-1
  return(bestResponse)
}

MiningFromShortTimeSeriesDataUnsupervisionThreeStages<-function(trainingseries,testseries,genes,haltrate=0,reconstructTimeInterval=1,maxTryRate=1,maxk=4,multithread=FALSE,timeporal=1,decay=1)
{

  cube<-constructFBNCube(genes,genes,trainingseries,maxk,timeporal,multithread)

  testErrorTry<-0
  maxRules<-5
  maxFBNRulesOfActivator<<-maxRules
  maxFBNRulesOfInhibitor<<-maxRules
  bestResult<-NULL
  bestNetwork<-NULL
  bestsetting<-""
  tries<-1
  if(maxTryRate<testErrorTry)
    maxTryRate<-testErrorTry
  #Stage 1
  while(testErrorTry<=maxTryRate)
  {
    error_toleranceRate<<-testErrorTry #update global values
    print(paste("Execute ",tries, " times with error_toleranceRate=",error_toleranceRate, " and confidence_toleranceRate=",confidence_toleranceRate, sep="", collapse = ""))
    extractedNetwork<-mineFBNNetwork(cube,genes)

    resultfile<-FBNBenchmark(extractedNetwork,NULL,testseries,testseries,reconstructTimeInterval,decay,timeporal,FALSE,cube)

    if(!is.null(bestResult))
    {
      if(resultfile$ErrorRate < bestResult$ErrorRate)
      {
        bestResult<-resultfile
        bestNetwork<-extractedNetwork
        bestsetting<-error_toleranceRate
      }
    }else
    {
      bestResult<-resultfile
      bestNetwork<-extractedNetwork
      bestsetting<-error_toleranceRate
    }
    print(paste("The best error rate is ",bestResult$ErrorRate, " And compare error rate is ",resultfile$ErrorRate , sep="", collapse = ""))
    if(bestResult$ErrorRate==0.0 | bestResult$ErrorRate<haltrate)
      break;
    tries<-tries+1
    #testConfidence<-testConfidence+0.1
    #}
    #testConfidence<-1 #reset
    testErrorTry<-testErrorTry+0.1
  }

  bestresultgroups<-list()
  maxTestTime<-10

  #Stage 2 and 3
  #State 3
  print("Stage 3 Testing")
  i<-1
  while(i<=maxTestTime)
  {
    print(paste("test=",i,sep="",collapse = ""))
    bestResponse<-MiningFromShortTimeSeriesDataUnsupervisionTwoStages(bestResult$FBNResult$reconstructed,testseries,genes,haltrate,reconstructTimeInterval,maxTryRate,maxk,multithread,timeporal,decay)
    bestResult<-bestResponse$BestResult
    bestresultgroups<-list(bestresultgroups,bestResponse$BestResult$FBNResult$reconstructed)
    i<-i+1
  }

  bestresultgroups<-shrinkTimeseriesCube(dissolve(bestresultgroups))
  cube<-constructFBNCube(genes,genes,bestresultgroups,maxk,timeporal,multithread)
  extractedNetwork<-mineFBNNetwork(cube,genes)
  resultfile<-FBNBenchmark(extractedNetwork,NULL,testseries,bestResult$FBNResult$reconstructed,reconstructTimeInterval,decay,timeporal,FALSE,cube)

  #final
  bestsetting<-c("error_toleranceRate"=error_toleranceRate,"maxFBNRulesOfActivator"=maxFBNRulesOfActivator,"maxFBNRulesOfInhibitor"=maxFBNRulesOfInhibitor)
  bestResponse<-list()
  bestResponse[[1]]<-extractedNetwork
  bestResponse[[2]]<-resultfile
  bestResponse[[3]]<-bestsetting
  names(bestResponse)[[1]]<-"BestNetwork"
  names(bestResponse)[[2]]<-"BestResult"
  names(bestResponse)[[3]]<-"BestSetting"
  print("Best Network generated:")
  print(bestResponse[[3]])


  return(bestResponse)
}

MiningMinimumRulesFromShortTimeSeriesDataUnsupervision<-function(trainingseries,testseries,genes,haltrate=0,reconstructTimeInterval=1,maxTryRate=1,maxk=4,multithread=FALSE,timeporal=1,decay=1)
{

  cube<-constructFBNCube(genes,genes,trainingseries,maxk,timeporal,multithread)

  testErrorTry<-0
  maxRules<-5
  maxFBNRulesOfActivator<<-maxRules
  maxFBNRulesOfInhibitor<<-maxRules
  bestResult<-NULL
  bestNetwork<-NULL
  bestsetting<-""
  tries<-1
  if(maxTryRate<testErrorTry)
    maxTryRate<-testErrorTry

  while(testErrorTry<=maxTryRate)
  {
    error_toleranceRate<<-testErrorTry #update global values
    print(paste("Execute ",tries, " times with error_toleranceRate=",error_toleranceRate, " and confidence_toleranceRate=",confidence_toleranceRate, sep="", collapse = ""))
    extractedNetwork<-mineFBNNetwork(cube,genes)

    resultfile<-FBNBenchmark(extractedNetwork,trainingseries,testseries,testseries,reconstructTimeInterval,decay,timeporal,FALSE,cube)

    if(!is.null(bestResult))
    {
      if(resultfile$ErrorRate < bestResult$ErrorRate)
      {
        bestResult<-resultfile
        bestNetwork<-extractedNetwork
        bestsetting<-error_toleranceRate
      }
    }else
    {
      bestResult<-resultfile
      bestNetwork<-extractedNetwork
      bestsetting<-error_toleranceRate
    }
    print(paste("The best error rate is ",bestResult$ErrorRate, " And compare error rate is ",resultfile$ErrorRate , sep="", collapse = ""))
    if(bestResult$ErrorRate==0.0 | bestResult$ErrorRate<haltrate)
      break;
    tries<-tries+1
    #testConfidence<-testConfidence+0.1
    #}
    #testConfidence<-1 #reset
    testErrorTry<-testErrorTry+0.1
  }


  if(bestResult$ErrorRate<=0.0 &  error_toleranceRate==0 | bestResult$ErrorRate<=haltrate | maxTryRate==testErrorTry)
  {

    cube<-constructFBNCube(genes,genes,bestResult$FBNResult$reconstructed,maxk,timeporal,multithread)
    extractedNetwork<-mineFBNNetwork(cube,genes)
    bestResult<-FBNBenchmark(extractedNetwork,NULL,testseries,bestResult$FBNResult$reconstructed,reconstructTimeInterval,decay,timeporal,FALSE,cube)

    print("Verify error")
    print(bestResult$ErrorRate)

    #find the best number of rules
    startActivator<-1
    startInhibitor<-1
    tries<-1
    while(startActivator<=maxRules)
    {
      needToStop<-FALSE
      maxFBNRulesOfActivator<<-startActivator
      while(startInhibitor<=maxRules)
      {
        maxFBNRulesOfInhibitor<<-startInhibitor
        print(paste("Execute ",tries, " times with maxFBNRulesOfActivator=",maxFBNRulesOfActivator, " and maxFBNRulesOfInhibitor=",maxFBNRulesOfInhibitor, sep="", collapse = ""))

        extractedNetwork<-mineFBNNetwork(cube,genes)
        resultfile<-FBNBenchmark(extractedNetwork,NULL,testseries,bestResult$FBNResult$reconstructed,reconstructTimeInterval,decay,timeporal,FALSE,cube)
        print(paste("The error rate is ",resultfile$ErrorRate, sep="", collapse = ""))
        if(resultfile$ErrorRate<=0.0 &  error_toleranceRate==0 | resultfile$ErrorRate<=haltrate)
        {
          print(bestResult$ErrorRate)
          print(extractedNetwork)
          bestResult<-resultfile
          bestNetwork<-extractedNetwork
          bestsetting<-c("error_toleranceRate"=error_toleranceRate,"maxFBNRulesOfActivator"=maxFBNRulesOfActivator,"maxFBNRulesOfInhibitor"=maxFBNRulesOfInhibitor)
          needToStop<-TRUE
          break
        }
        tries<-tries+1
        startInhibitor<-startInhibitor+1
      }
      startInhibitor<-1
      if(needToStop)
        break;

      startActivator<-startActivator+1
    }
    startActivator<-1

    bestResponse<-list()
    bestResponse[[1]]<-bestNetwork
    bestResponse[[2]]<-bestResult
    bestResponse[[3]]<-bestsetting
    names(bestResponse)[[1]]<-"BestNetwork"
    names(bestResponse)[[2]]<-"BestResult"
    names(bestResponse)[[3]]<-"BestSetting"
    print("Best Network generated:")
    print(bestResponse[[3]])
  }else
  {
    bestResponse<-MiningMinimumRulesFromShortTimeSeriesDataUnsupervision(bestResult$FBNResult$reconstructed,testseries,genes,haltrate,reconstructTimeInterval,maxTryRate,maxk,multithread,timeporal,decay)
  }

  return(bestResponse)
}
