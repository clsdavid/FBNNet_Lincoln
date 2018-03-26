mineFBNNetworkCore<-function(fbnGeneCube,genes,useParallel=FALSE)
{
  #load configulation
  readConfigFile()

  internalloop<-function(i,genes,fbnGeneCube){
    recursiveMiningFBNFunction<-function(fbnGeneCubeStem,targetgene,currentsubgene,mainStatistic,groupbyGene)
    {
      if(is.null(groupbyGene))
      {
        groupbyGene<-currentsubgene
      }
      chiSQP<-fbnGeneCubeStem[[currentsubgene]][["Statistic"]][["chiSQ"]]
      chiSQ_c<-fbnGeneCubeStem[[currentsubgene]][["Statistic"]][["chiSQ_c"]]
      fisherTest<-fbnGeneCubeStem[[currentsubgene]][["Statistic"]][["fisherTest"]]
      fisherTest_c<-fbnGeneCubeStem[[currentsubgene]][["Statistic"]][["fisherTest_c"]]

      Noise_P<-as.numeric(fbnGeneCubeStem[[currentsubgene]][["Statistic"]]["Noise_P"])
      Noise_N<-as.numeric(fbnGeneCubeStem[[currentsubgene]][["Statistic"]]["Noise_N"])

      avgErrorP<-as.numeric(fbnGeneCubeStem[[currentsubgene]][["Statistic"]]["avgErrorP"])
      avgErrorN<-as.numeric(fbnGeneCubeStem[[currentsubgene]][["Statistic"]]["avgErrorN"])

      avgErrorPsd<-as.numeric(fbnGeneCubeStem[[currentsubgene]][["Statistic"]]["stdErrorP"])
      avgErrorNsd<-as.numeric(fbnGeneCubeStem[[currentsubgene]][["Statistic"]]["stdErrorN"])

      pvalue<-min(chiSQP,fisherTest)
      pvalue_c<-min(chiSQ_c,fisherTest_c)

      parentStatic<-fbnGeneCubeStem[[currentsubgene]][["Statistic"]]

      isNegativeCorrelated<-fbnGeneCubeStem[[currentsubgene]][["Statistic"]][["isNegativeCorrelated"]]
      isPossitiveCorrelated<-fbnGeneCubeStem[[currentsubgene]][["Statistic"]][["isPossitiveCorrelated"]]

      input<-sort(fbnGeneCubeStem[[currentsubgene]][["Statistic"]][["Input"]])

      isWeakLink<-fbnGeneCubeStem[[currentsubgene]][["Statistic"]][["IsWeakLink"]]

      #identity<-paste(input,"_",sep="",collapse = '')
      res<-list()
      resT<-list()
      resF<-list()
      resTMerged<-list()
      resFMerged<-list()
      se_T<-0
      se_F<-0
      pickTcons<-0
      pickFcons<-0
      pickTvalue<-0
      pickFvalue<-0
      pickTTransitionP<-0
      pickFTransitionP<-0
      pickTallconfidence<-0
      pickFallconfidence<-0
      #defauilt time step is 1 (synchronous update schema)
      timestepT<-1
      timestepF<-1

      pickT<-fbnGeneCubeStem[[currentsubgene]][["ActivatorAndInhibitor"]][["Activator"]]
      pickF<-fbnGeneCubeStem[[currentsubgene]][["ActivatorAndInhibitor"]][["Inhibitor"]]


      #if(targetgene=="Gene2" & pickT[1]=="Gene1&!Gene4&Gene5")
      #{
        #print(pickT)
        #print(pvalue)
      #}

      #if(pvalue<chi_P_FBN)#is significant enough
      #{
        if(!is.null(pickT))
        {
          pickTvalue<-as.numeric(pickT["Confidence"]) #avgSignal_T + sdSignal_T#
          pickTallconfidence<-as.numeric(pickT["all_confidence"])
          pickTmaxconfidence<-as.numeric(pickT["max_confidence"])
          pickTTransitionP<-as.numeric(pickT["TransitionP"])
          pickTMutualInfo<-abs(as.numeric(pickT["MutualInfo"]))
          pickTcosine<-as.numeric(pickT["cosine"])
          pickTsupport<-as.numeric(pickT["support"])
          pickTcons<-as.numeric(pickT["rate"])#*******
          se_T<-as.numeric(pickT["Noise"])
          identityT<-pickT["Identity"]
          error_Entropy_P<-pickT["Entropy"]
          pickTType<-pickT["type"]
          preSupportT<-pickT["preSupport"]
          timestepT<-as.numeric(pickT["timestep"])
        }

        if(!is.null(pickF))
        {
          pickFvalue<-as.numeric(pickF["Confidence"]) #avgSignal_F + sdSignal_F#
          pickFallconfidence<-as.numeric(pickF["all_confidence"])
          pickFmaxconfidence<-as.numeric(pickF["max_confidence"])
          pickFTransitionP<-as.numeric(pickF["TransitionP"])
          pickFMutualInfo<-abs(as.numeric(pickF["MutualInfo"]))
          pickFcosine<-as.numeric(pickF["cosine"])
          pickFsupport<-as.numeric(pickF["support"])
          pickFcons<-as.numeric(pickF["rate"]) #******
          se_F<-as.numeric(pickF["Noise"])
          identityF<-pickF["Identity"]
          error_Entropy_N<-pickF["Entropy"]
          pickFType<-pickF["type"]
          preSupportF<-pickF["preSupport"]
          timestepF<-as.numeric(pickF["timestep"])
        }
        errorActivator<-se_T
        errorInhibitor<-se_F

        #pickTcosine
        subapplied<-FALSE
        if(!is.null(fbnGeneCubeStem[[currentsubgene]]$SubGenesT))
        {
          subapplied<-TRUE
          fbnGeneCubeStemT<-fbnGeneCubeStem[[currentsubgene]]$SubGenesT
          nextgenesT<-names(fbnGeneCubeStemT)
          indexT<-1
          for(j in seq_along(nextgenesT))
          {
            nextgene<-nextgenesT[[j]]
            resTU<-dissolve(recursiveMiningFBNFunction(fbnGeneCubeStemT,targetgene,nextgene,mainStatistic,groupbyGene))
            if(length(resTU)>0)
            {
              condT<-sapply(resTU,function(entry)!is.null(entry))
              resT[[indexT]]<-resTU[unlist(lapply(resTU[condT],length)!=0)]
              indexT<-indexT+1
            }
          }
          if(length(resT)>0)
          {
            resT<-dissolve(resT)
          }
        }
        if(!is.null(fbnGeneCubeStem[[currentsubgene]]$SubGenesF))
        {
          subapplied<-TRUE
          fbnGeneCubeStemF<-fbnGeneCubeStem[[currentsubgene]]$SubGenesF
          nextgenesF<-names(fbnGeneCubeStemF)
          indexF<-1
          for(j in seq_along(nextgenesF))
          {
            nextgene<-nextgenesF[[j]]
            resFU<-dissolve(recursiveMiningFBNFunction(fbnGeneCubeStemF,targetgene,nextgene,mainStatistic,groupbyGene))
            if(length(resFU)>0)
            {
              condF<-sapply(resFU,function(entry)!is.null(entry))
              resF[[indexF]]<-resFU[unlist(lapply(resFU[condF],length)!=0)]
              indexF<-indexF+1
            }
          }
          if(length(resF)>0)
          {
            resF<-dissolve(resF)
          }
        }

        if(pvalue<p_value_FBNCube & !is.null(pickT) & pickTcons>=1  & pickTMutualInfo>=minMutual_range & pickTMutualInfo<=1 & pickTvalue>thredsod_FBNCube)
        {
          #if(targetgene=="Gene2" & pickT[1]=="Gene1&!Gene4&Gene5")
          #{
          #print(pickT)
          #print(pvalue)
          #}

          getT<-pickT
        }else
        {
          getT<-NULL
        }

        #if(!is.null(pickF) & pickFcons>=1)  #original work
        if(pvalue<p_value_FBNCube& !is.null(pickF) & pickFcons>=1 & pickFMutualInfo>=minMutual_range & pickTMutualInfo<=1 & pickFvalue>thredsod_FBNCube)
        {
          getF<-pickF
        }else
        {
          getF<-NULL
        }
        newindex<-length(res)+1


        if(!is.null(getT))
        {

          res[[newindex]]<-c("targets"=targetgene,"factor"=getT[1],"type"=1,"identity"=identityT, "error"=round(errorActivator,5),"P"=round(pickTvalue,4),"support"=pickTsupport,"timestep"=timestepT,"input"=paste(input,collapse = ","),"numOfInput"=length(input), "causeRate"=round(pickTcons,5),"GroupBy"=groupbyGene, "tRate"=pickTTransitionP,"all_confidence"=pickTallconfidence,"MutualInfo"=pickTMutualInfo,"avgErrorP"=avgErrorP,"Entropy_Error"=error_Entropy_P,"dimensionType"=pickTType)

          newindex<-length(res)+1

        }

        if(!is.null(getF))
        {
          #identityF<-paste(identity,0,sep="",collapse = '')
          res[[newindex]]<-c("targets"=targetgene,"factor"=getF[1],"type"=0,"identity"=identityF, "error"=round(errorInhibitor,5),"P"=round(pickFvalue,4),"support"=pickFsupport,"timestep"=timestepF,"input"=paste(input,collapse = ","),"numOfInput"=length(input), "causeRate"=round(pickFcons,5),"GroupBy"=groupbyGene, "tRate"=pickFTransitionP,"all_confidence"=pickFallconfidence,"MutualInfo"=pickFMutualInfo,"avgErrorN"=avgErrorN,"Entropy_Error"=error_Entropy_N,"dimensionType"=pickFType)

          newindex<-length(res)+1

        }
      #}


      if(length(resT)>0 & length(resF))
      {
        weedResT<-list()
        weedResF<-list()

        for(i in seq_along(resT))
        {
          for(j in seq_along(resF))
          {
            if(resT[[i]][["input"]]==resF[[j]][["input"]] & resT[[i]][["dimensionType.type"]]==resF[[j]][["dimensionType.type"]] & as.numeric(resT[[i]][["MutualInfo"]])==as.numeric(resF[[j]][["MutualInfo"]]))
            {
              weedResT[length(weedResT)+1]<-resT[i]
              weedResF[length(weedResF)+1]<-resF[j]
            }
          }
        }
        weedResT<-unique(weedResT)
        weedResF<-unique(weedResF)
        resT<-resT[!resT%in%weedResT]
        resF<-resF[!resF%in%weedResF]
      }
      if(length(resT)>0)
      {
        res<-list(res,resT)
      }

      if(length(resF)>0)
      {
        res<-list(res,resF)
      }

      res<-dissolve(res)

      if(length(res)>0)
      {
        cond1<-sapply(res,function(entry)!is.null(entry))
        res<-(res[cond1][unlist(lapply(res[cond1],length)!=0)])
      }


      return(res)
    }

    removeDuplicates<-function(factors)
    {
      cond1<-sapply(factors,function(x)x[[4]])
      #getAllDuplicate<-factors[duplicated(cond1)]
      return(factors[!duplicated(cond1)])
    }

    try({
      targetGene<-genes[[i]]
      res<-list()
      res[[i]]<-list()
      names(res)[[i]]<-targetGene
      currentStem<-fbnGeneCube[[targetGene]]$SubGenes
      Statistic<-fbnGeneCube[[targetGene]]$Statistic

      nextgenes<-names(currentStem)

      for(k in seq_along(nextgenes))
      {
        currentGene<-nextgenes[[k]]
        subres<-dissolve(recursiveMiningFBNFunction(currentStem,targetGene,currentGene,Statistic,NULL))
        if(length(subres)>0)
        {
          cond1<-sapply(subres,function(entry)!is.null(entry))
          subres<-(subres[cond1][unlist(lapply(subres[cond1],length)!=0)])
          if(!is.null(subres))
          {
            if(length(subres)>0)
            {
              index<-length(res[[i]])+1
              #get result and remove duplicates
              resultsub<-dissolve(subres)
              res[[i]][[index]]<-resultsub
              names(res[[i]])[[index]]<-currentGene
            }
          }
        }
      }

      #pre pruning
      maxNumOfRule<-50

      preresponse<-removeDuplicates(dissolve(res[[i]]))

      if(length(preresponse)==0)
      {
        return(list())
      }
      cond2<-sapply(preresponse,function(entry)!is.null(entry))
      preresponse<-(preresponse[cond2][unlist(lapply(preresponse[cond2],length)!=0)])

      activators<-lapply(preresponse,function(rule)
      {
        if(is.null(rule))return(NULL)

        if(as.numeric(rule["type"])==1)
        {
          #if(rule["targets"]=="Gene2" & rule["factor.factor"]=="Gene1&!Gene4&Gene5")
          #{
            #print(rule)
          #}
          return (rule)
        }
      }
      )

      inhibitors<-lapply(preresponse,function(rule)
      {
        if(is.null(rule))return(NULL)
        if(as.numeric(rule["type"])==0)
        {
          return (rule)
        }
      }
      )

      cond2<-sapply(activators,function(entry)!is.null(entry))
      activators<-(activators[cond2][unlist(lapply(activators[cond2],length)!=0)])

      cond2<-sapply(inhibitors,function(entry)!is.null(entry))
      inhibitors<-(inhibitors[cond2][unlist(lapply(inhibitors[cond2],length)!=0)])


      finalRes<-list()
      finalRes[[1]]<-activators
      finalRes[[2]]<-inhibitors
      finalRes<-dissolve(finalRes)

      res[[i]]<-finalRes
      #res[[i]]<-dissolve(res[[i]])
      return(res)
    })}

  time1<-as.numeric(Sys.time())

  res<-list()

  if(length(genes)==0)
  {
    return(list())
  }
  if(useParallel)
  {
    res<-doParallelWork(internalloop,genes,fbnGeneCube)

  }else
  {
    res<-doNonParallelWork(internalloop,genes,fbnGeneCube)
  }

  cond1<-sapply(res,function(entry)!is.null(entry))
  res<-(res[cond1][unlist(lapply(res[cond1],length)!=0)])



  if(useParallel)
  {
    closeAllConnections()
  }

  time2<-as.numeric(Sys.time())
  print(paste("Total cost ", time2-time1, " seconds to mine Fundamental Boolean Functions ", sep='',collapse = ''))
  return(res)
}

mineFBNNetworkStage2<-function(res)
{

  if(is.null(res)|length(res)==0)
  {
    return(list())
  }



  filteredres<-lapply(res,function(entry){

    allerrors<-unlist(lapply(entry,function(rule)
      {
        return (as.numeric(rule["error"]))
      }
    ))


    if(is.null(allerrors))
      allerrors<-c(0)


    allMutualInfo<-unlist(lapply(entry,function(rule)
        return (as.numeric(rule["MutualInfo"]))
    ))



    if(is.null(allMutualInfo))
      allMutualInfo<-c(0)


    #error
    tolerancesofErrorsSd<-sd(allerrors)
    tolerancesofErrorsSd<-ifelse(is.na(tolerancesofErrorsSd),0,tolerancesofErrorsSd)

    toleranceOfErrors<-round(min(allerrors) + error_toleranceRate*tolerancesofErrorsSd,5)#toleranceRate*(max(errorsActivator)-min(errorsActivator))

    #confidence

    toleranceOfMutualInfoSD<-sd(allMutualInfo)

    toleranceOfMutualInfoSD<-ifelse(is.na(toleranceOfMutualInfoSD),0,toleranceOfMutualInfoSD)

    toleranceOfMutualInfo<-round(sum(allMutualInfo)/length(allMutualInfo)+Mutual_toleranceRate * toleranceOfMutualInfoSD,5)


    if(toleranceOfMutualInfo>1)toleranceOfMutualInfo<-1


    res<-list()
    for(e in seq_along(entry))
    {
      #if(entry[[e]]["targets"]=="Gene2" & entry[[e]]["factor.factor"]=="Gene1&!Gene4&Gene5")
      #{
        #print(entry[[e]])
      #}

      if(as.numeric(entry[[e]]["MutualInfo"])>=toleranceOfMutualInfo & as.numeric(entry[[e]]["error"])<=toleranceOfErrors & as.numeric(entry[[e]]["type"])==1)
      {
        res[[length(res)+1]]<-c(entry[[e]],"tolerance_error"=toleranceOfErrors)
      }

      if(as.numeric(entry[[e]]["MutualInfo"])>=toleranceOfMutualInfo & as.numeric(entry[[e]]["error"])<=toleranceOfErrors & as.numeric(entry[[e]]["type"])==0)
      {
        res[[length(res)+1]]<-c(entry[[e]],"tolerance_error"=toleranceOfErrors)
      }
    }

    if(length(res)==0)
    {
      return(res)
    }

    activators<-lapply(res,function(rule)
    {
      if(is.null(rule))return(NULL)
      if(as.numeric(rule["type"])==1)
      {
        return (rule)
      }
    }
    )

    inhibitors<-lapply(res,function(rule)
    {
      if(is.null(rule))return(NULL)
      if(as.numeric(rule["type"])==0)
      {
        return (rule)
      }
    }
    )

    cond1<-sapply(activators,function(entry)!is.null(entry))
    activators<-(activators[cond1][unlist(lapply(activators[cond1],length)!=0)])

    cond1<-sapply(inhibitors,function(entry)!is.null(entry))
    inhibitors<-(inhibitors[cond1][unlist(lapply(inhibitors[cond1],length)!=0)])

    activators<-activators[order(-as.numeric(sapply(activators,"[[","MutualInfo")),as.numeric(sapply(activators,"[[","error")),-as.numeric(sapply(activators,"[[","support")), as.numeric(sapply(activators,"[[","numOfInput")))]
    if(length(activators)>maxFBNRulesOfActivator)
    {
      activators<-activators[1:maxFBNRulesOfActivator]
    }
    #order by less input
    activators<-activators[order(as.numeric(sapply(activators,"[[","error")),as.numeric(sapply(activators,"[[","numOfInput")))]

    inhibitors<-inhibitors[order(-as.numeric(sapply(inhibitors,"[[","MutualInfo")),as.numeric(sapply(inhibitors,"[[","error")),-as.numeric(sapply(inhibitors,"[[","support")),as.numeric(sapply(inhibitors,"[[","numOfInput")))]
    if(length(inhibitors)>maxFBNRulesOfInhibitor)
    {
      inhibitors<-inhibitors[1:maxFBNRulesOfInhibitor]
    }
    #order by less input
    inhibitors<-inhibitors[order(as.numeric(sapply(inhibitors,"[[","error")),as.numeric(sapply(inhibitors,"[[","numOfInput")))]

    finalRes<-list()
    finalRes[[1]]<-activators
    finalRes[[2]]<-inhibitors
    finalRes<-dissolve(finalRes)
    #futher remove based on maximum activators and maximum inhibitors
    return(finalRes)
  })


  cond1<-sapply(filteredres,function(entry)!is.null(entry))
  filteredres<-(filteredres[cond1][unlist(lapply(filteredres[cond1],length)!=0)])

  finalFilteredlist<-list()

  targetgenes<-names(filteredres)
  for(i in seq_along(targetgenes))
  {
    target<-targetgenes[i]
    finalFilteredlist[[i]]<-list()
    names(finalFilteredlist)[[i]]<-target

    ruleset<-filteredres[[target]]
    for(j in seq_along(ruleset))
    {
      rule<-ruleset[[j]]
      for(k in seq_along(ruleset))
      {
        rule2<-ruleset[[k]]
        if(as.numeric(rule[["numOfInput"]])<as.numeric(rule2[["numOfInput"]]) & as.numeric(rule[["type"]])==as.numeric(rule2[["type"]]) & all(splitExpression(rule[["input"]],2)%in%splitExpression(rule2[["input"]],2)==TRUE))
        {
          finalFilteredlist[[i]][[length(finalFilteredlist[[i]])+1]]<-rule2
        }
      }
    }

    filteredres[[target]]<-filteredres[[target]][!filteredres[[target]]%in%finalFilteredlist[[i]]]
  }

  class(filteredres)<-c("FBNTrueCubeMiner")
  return(filteredres)
}

#'Mine FBN Networks from an Orchard cube
#'
#'@param fbnGeneCube A pre constructed Orchard cube
#'@param genes The target genes in the output
#'@param useParallel An option turns on parallel
#'@return A object of FBN network
#'@examples
#' mat1<-matrix(c("1","0","0","1","0","0","0","1","1"),3,3, dimnames=list(c("gene1","gene2","gene3"),c("1","2","3")))
#' mat2<-matrix(c("1","1","0","1","0","1","1","1","0"),3,3, dimnames=list(c("gene1","gene2","gene3"),c("1","2","3")))
#' listtest<-list(mat1,mat2)
#' cube<-constructFBNCube(c("gene1","gene2"),c("gene1","gene2","gene3"),listtest,4,1,FALSE)
#' network<-mineFBNNetwork(cube,c("gene1","gene2"))
#' @export
mineFBNNetwork<-function(fbnGeneCube,genes,useParallel=FALSE)
{
  res<-mineFBNNetworkCore(fbnGeneCube,genes,useParallel)
  res<-mineFBNNetworkStage2(res)
  finalresult<-convertMinedResultToFBNNetwork(res,genes)
  if(length(finalresult)>0)
  {
    cond1<-sapply(finalresult,function(entry)!is.null(entry))
    finalresult<-(finalresult[cond1][unlist(lapply(finalresult[cond1],length)!=0)])
  }
  class(finalresult)<-c("FundamentalBooleanNetwork","BooleanNetworkCollection")
  return(finalresult)
}
