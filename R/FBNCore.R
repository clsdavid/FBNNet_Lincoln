#'An internal method to extract gene state from timeseries cube
#'
#'@param timeSeriesCube A list object contains samples in which a sample contains genes in row and time points in column
#'@return A combined matrix
#'@examples
#' mat1<-matrix(c("1","0","0","1","0","0","0","1","1"),3,3, dimnames=list(c("gene1","gene2","gene3"),c("1","2","3")))
#' mat2<-matrix(c("1","1","0","1","0","1","1","1","0"),3,3, dimnames=list(c("gene1","gene2","gene3"),c("1","2","3")))
#' listtest<-list(mat1,mat2)
#' extractGeneStateFromTimeSeriesCube(listtest)
#' @export
extractGeneStateFromTimeSeriesCube<-function(timeSeriesCube)
{
  genes<-rownames(timeSeriesCube[[1]])

  newseries<-lapply(1:length(timeSeriesCube),function(entryindex,ngenes){
    entry<-timeSeriesCube[[entryindex]]
    if(entryindex<length(timeSeriesCube))
    {
      entry2<-cbind(entry,rep(9,ngenes))
      entry<-entry2
    }
    return(entry)
  },length(genes))
  combined<-do.call(cbind,newseries)

  m <- mapply(combined, FUN=as.numeric)
  combined <- matrix(data=m, ncol=dim(combined)[[2]],nrow=dim(combined)[[1]], dimnames = list(rownames(combined),colnames(combined)))
  #can't remove duplicate as it will impact on the final result
  return(combined)
}

extractGenesFromTimeSeriesCube<-function(timeSeriesCube,targetgenes,geneState=NULL)
{
  rownames<-rownames(timeSeriesCube[[1]])

  if(!all(targetgenes%in%rownames))
  {
    stop("All or some part of the target genes are not founded in the timeseries cube")
  }

  filter<-rownames[rownames%in%targetgenes]
  if(is.null(geneState))
  {
    thisgenestate<-extractGeneStateFromTimeSeriesCube(timeSeriesCube)
  }else
  {
    thisgenestate<-geneState
  }
  numOfColumns<-dim(thisgenestate)[2]
  extracted<-thisgenestate[rownames(thisgenestate)%in%filter]

  mat<-matrix(extracted,nrow=length(filter),ncol=numOfColumns,byrow=FALSE,dimnames<-list(filter,c(1:numOfColumns)))
  return(mat)

}

#'A core method to calculate all measures for a target gene with a condition gene
#'
#'@param timeSeriesCube A list object contains samples in which a sample contains genes in row and time points in column
#'@param fixedgenestate prefixed gene state
#'@param targetgene The tage gene's name
#'@param newgene The new added conditional gene
#'@param geneStates The current combined gene state
#'@param preGeneStates The previous combined gene state
#'@param temporal The temporal setting and the default value is 1
#'@param temporalGeneStates The combined gene state for temporal
#'@return A list of measures
#'@examples
#'genesInput<-c("CycD","p27","CycE","E2F")
#'testseries<-list()
#'testseries[[1]]<-matrix(c(1,1,1,1,0,1,0,1,1,1,1,1,0,0,0,1),nrow=4,ncol=4,byrow=FALSE,dimnames<-list(genesInput,c("1","2","3","4")))
#'testseries[[2]]<-matrix(c(1,1,0,1,0,1,0,1,0,1,0,1,1,0,0,1),nrow=4,ncol=4,byrow=FALSE,dimnames<-list(genesInput,c("1","2","3","4")))
#'testseries[[3]]<-matrix(c(1,0,0,1,1,1,0,1,1,1,0,1,1,1,0,1),nrow=4,ncol=4,byrow=FALSE,dimnames<-list(genesInput,c("1","2","3","4")))
#'
#'getCurrentState<-extractGeneStateFromTimeSeriesCube(testseries)
#'getCurrentState<-getCurrentState[,-1]#remove the first one
#'
#'getpreviousState<-extractGeneStateFromTimeSeriesCube(testseries)
#'getpreviousState<-getpreviousState[,-length(colnames(getpreviousState))] #remove the last one
#'probability<-getGenePrababilities(testseries,NULL,"CycD","p27",getCurrentState,getpreviousState,1)
#' @export
getGenePrababilities<-function(timeSeriesCube,fixedgenestate,targetgene,newgene,geneStates=NULL,preGeneStates=NULL,temporal=1,temporalGeneStates=NULL)
{
  #P(B|A)=P(A & B)/P(A)=(frq(A&B)/n)/(frq(A)/n)=frq(A&B)/frq(A)=P(A&D)/(P(A&D)+P(!A & D)) Bayers rule
  #order of targetgenes is very important
  rownames<-rownames(timeSeriesCube[[1]])

  if(is.null(fixedgenestate))
  {
    fixedgenestate<-list()
    conditionalgenes<-c(newgene)
  }else
  {
    conditionalgenes<-names(fixedgenestate)

    if(!all(conditionalgenes%in%rownames))
    {
      stop("All or some part of the target genes are not founded in the timeseries cube")
    }

    if(any(c(newgene)%in%conditionalgenes))
    {
      conditionalgenes<-conditionalgenes[conditionalgenes!=newgene]
      fixedgenestate<-fixedgenestate[names(fixedgenestate)!=newgene]
    }else
    {
      conditionalgenes<-c(conditionalgenes,newgene)
    }
  }

  #initialize
  TTgenestatesCond<-fixedgenestate
  FFgenestatesCond<-fixedgenestate

  newindex<-length(fixedgenestate)+1
  #set target gene
  TTgenestatesCond[newindex]<-1
  names(TTgenestatesCond)[newindex]<-newgene

  FFgenestatesCond[newindex]<-0
  names(FFgenestatesCond)[newindex]<-newgene

  #get states in order
  filter<-rownames[rownames%in%conditionalgenes]
  stateTCond<-unlist(lapply(filter,function(x)TTgenestatesCond[[x]]))
  stateFCond<-unlist(lapply(filter,function(x)FFgenestatesCond[[x]]))

  finalStateTTCound<-c(stateTCond,1)
  finalStateTFCound<-c(stateFCond,1)
  finalStateFTCound<-c(stateTCond,0)
  finalStateFFCound<-c(stateFCond,0)

  #Test with temporal Boolean functions
  if(is.null(temporalGeneStates))
  {
    getAllTemporalStates<-generateTemporalGeneStates(timeSeriesCube,targetgene,conditionalgenes,temporal,geneStates,preGeneStates)
  }else
  {
    getAllTemporalStates<-temporalGeneStates
  }

  totalTCond <-0
  totalFCond <-0
  totalTTarget<-0
  totalFTarget<-0
  totalTT<-0
  totalTF<-0
  totalFT<-0
  totalFF<-0
  countOfDifferentTemporal<-length(getAllTemporalStates)
  resultGroup<-list()
  for(i in seq_along(getAllTemporalStates))
  {
    #condition
    extracedCond<-getAllTemporalStates[[i]]$ConditionState
    extracedTarget<-getAllTemporalStates[[i]]$TargetState

    concatenatedMatrix<-rbind(extracedCond,extracedTarget)

    midTTCond<-abs(finalStateTTCound-concatenatedMatrix)
    midTFCond<-abs(finalStateTFCound-concatenatedMatrix)
    midFTCond<-abs(finalStateFTCound-concatenatedMatrix)
    midFFCond<-abs(finalStateFFCound-concatenatedMatrix)


    resTTCond<-which(colSums(midTTCond)==0)
    resTFCond<-which(colSums(midTFCond)==0)
    resFTCond<-which(colSums(midFTCond)==0)
    resFFCond<-which(colSums(midFFCond)==0)

    lenTT<-length(resTTCond)
    lenTF<-length(resTFCond)
    lenFT<-length(resFTCond)
    lenFF<-length(resFFCond)

    totalTT<-totalTT+lenTT
    totalTF<-totalTF + lenTF
    totalFT<-totalFT + lenFT
    totalFF<-totalFF +lenFF

    #calculate total target. This should calculate once and need to revise
    if(length(stateTCond)>1)
    {
      concatenatedMatrix<-concatenatedMatrix[(length(rownames(concatenatedMatrix))-1):length(rownames(concatenatedMatrix)),]
      midTT<-abs(c(1,1)-concatenatedMatrix)
      midTF<-abs(c(0,1)-concatenatedMatrix)
      midFT<-abs(c(1,0)-concatenatedMatrix)
      midFF<-abs(c(0,0)-concatenatedMatrix)

      resTT<-which(colSums(midTT)==0)
      resTF<-which(colSums(midTF)==0)
      resFT<-which(colSums(midFT)==0)
      resFF<-which(colSums(midFF)==0)

      totalTTarget<-length(resTT)+length(resTF)
      totalFTarget<-length(resFT)+length(resFF)
    }else
    {
      totalTTarget<-totalTT+totalTF
      totalFTarget<-totalFT+lenFF
    }

    totalTCond<-totalTT+totalFT
    totalFCond<-totalTF+totalFF
    totalsamples<-totalTTarget+totalFTarget
    conditionP_TT<-totalTCond/totalsamples
    conditionP_FF<-totalFCond/totalsamples
    targetP_TT<-totalTTarget/totalsamples
    targetP_FF<-totalFTarget/totalsamples

    #final p(B if A)=p(B and A)/p(A), A=conditions, B is the target
    conditionTT<-ifelse(totalTCond>0,totalTT/totalTCond,0)
    conditionFT<-ifelse(totalTCond>0,totalFT/totalTCond,0)
    conditionTF<-ifelse(totalFCond>0,totalTF/totalFCond,0)
    conditionFF<-ifelse(totalFCond>0,totalFF/totalFCond,0)

    #final p(A if B)=p(B and A)/p(B), A=conditions, B is the target
    targetTT<-ifelse(totalTTarget>0,totalTT/totalTTarget,0)
    targetFT<-ifelse(totalFTarget>0,totalFT/totalFTarget,0)
    targetTF<-ifelse(totalTTarget>0,totalTF/totalTTarget,0)
    targetFF<-ifelse(totalFTarget>0,totalFF/totalFTarget,0)

    variable1<-c(totalTT+1,totalTF+1)
    variable2<-c(totalFT+1,totalFF+1)
    pTable<-as.table(cbind(variable1,variable2))
    chiSQ<-chisq.test(pTable,correct=FALSE,simulate.p.value = TRUE)

    fisher<-fisher.test(pTable,simulate.p.value = TRUE)

    noise_T<-min(conditionTT,conditionFT)
    noise_F<-min(conditionTF,conditionFF)

    signal_T<-max(conditionTT,conditionFT)
    signal_F<-max(conditionTF,conditionFF)

    #from the book Data mining concepts and techniques
    isNegativeCorrelated<-ifelse(totalsamples==0,0,((totalTF/totalsamples)*(totalFT/totalsamples))>((totalTT/totalsamples)*(totalFF/totalsamples))) #mainly inhibiting the target gene
    isPossitiveCorrelated<-ifelse(totalsamples==0,0,((totalTF/totalsamples)*(totalFT/totalsamples))<((totalTT/totalsamples)*(totalFF/totalsamples)))#mainly activating the target gene

    #error transformation
    errorActivator<-0
    errorInhibitor<-0

    if(max(conditionTT,conditionTF)==conditionTT)
    {
      errorActivator<-noise_T #no problem
    }else
    {
      errorActivator<-noise_F #no problem
    }

    if(max(conditionFT,conditionFF)==conditionFT)
    {
      errorInhibitor<-noise_T #no problem
    }else
    {
      errorInhibitor<-noise_F #no problem
    }

    noise_T<-errorActivator
    noise_F<-errorInhibitor

    if(any(c(conditionTT,conditionTF,conditionFT,conditionFF)>1))
    {
      stop("probability should be less or equal to one")
    }

    #Shannon entropy (Shannon & Weaver, 1963) and REVEAL (Liang,1998)
    #X is conditional, Y is target
    p_x1<-totalTCond/totalsamples
    p_x2<-totalFCond/totalsamples
    HX<--sum(p_x1*log2(p_x1)+p_x2*log2(p_x2))

    p_y1<-totalTTarget/totalsamples
    p_y2<-totalFTarget/totalsamples
    HY<--sum(p_y1*log2(p_y1)+p_y2*log2(p_y2))

    HX<-ifelse(is.na(HX),0,HX)
    HY<-ifelse(is.na(HY),0,HY)
    #H(X|Y)
    #HX_Y<--sum(conditionTT*log2(conditionTT)+conditionFT*log2(conditionFT)+conditionTF*log2(conditionTF)+conditionFF*log2(conditionFF))
    #HY_X<--sum(targetTT*log2(targetTT)+targetFT*log2(targetFT)+targetTF*log2(targetTF)+targetFF*log2(targetFF))
    #X is conditional inputs and Y is the target gene
    HXT_YT<- -conditionTT*log(conditionTT)
    HXF_YT<- -conditionTF*log(conditionTF)
    HXT_YF<- -conditionFT*log(conditionFT)
    HXF_YF<- -conditionFF*log(conditionFF)

    HXT_YT<-ifelse(is.na(HXT_YT),0,HXT_YT)
    HXF_YT<-ifelse(is.na(HXF_YT),0,HXF_YT)
    HXT_YF<-ifelse(is.na(HXT_YF),0,HXT_YF)
    HXF_YF<-ifelse(is.na(HXF_YF),0,HXF_YF)

    #Mutual Information
    MXT_YT <- HX - HXT_YT
    MXF_YT <- HX - HXF_YT
    MXT_YF <- HX - HXT_YF
    MXF_YF <- HX - HXF_YF

    dXT_YT<-abs(round(MXT_YT,5)/round(HX,5))
    dXF_YT<-abs(round(MXF_YT,5)/round(HX,5))
    dXT_YF<-abs(round(MXT_YF,5)/round(HX,5))
    dXF_YF<-abs(round(MXF_YF,5)/round(HX,5))

    #at moment use totalsamples receive the best result with 0.00013
    supportTT<-ifelse(totalsamples>0,totalTT/totalsamples,0)
    supportFT<-ifelse(totalsamples>0,totalFT/totalsamples,0)
    supportTF<-ifelse(totalsamples>0,totalTF/totalsamples,0)
    supportFF<-ifelse(totalsamples>0,totalFF/totalsamples,0)

    result<-list("TT"=conditionTT,"TF"=conditionTF,"FT"=conditionFT,"FF"=conditionFF,
                "Noise_P"=noise_T,"Noise_N"=noise_F,"Signal_P"=signal_T,"Signal_N"=signal_F,
                "chiSQ"=chiSQ$p.value, "fisherTest"=fisher$p.value,"conditionTCount"=totalTCond,"conditionFCount"=totalFCond,"conditionT"=conditionP_TT,"conditionF"=conditionP_FF,
                "targetT"=targetP_TT,"targetF"=targetP_FF,"targetTCount"=totalTTarget,"targetFCount"=totalFTarget,
                "lenTT"=totalTT,"lenTF"=totalTF,"lenFT"=totalFT,"lenFF"=totalFF,
                "targetTT"=targetTT,"targetTF"=targetTF,"targetFT"=targetFT,"targetFF"=targetTT,
                "isNegativeCorrelated"=isNegativeCorrelated,"isPossitiveCorrelated"=isPossitiveCorrelated,
                "dXT_YT"=dXT_YT,"dXF_YT"=dXF_YT,"dXT_YF"=dXT_YF,"dXF_YF"=dXF_YF, "HInput"=HX,"HOutput"=HY,
                "supportTT"=supportTT,"supportFT"=supportFT,"supportTF"=supportTF,"supportFF"=supportFF,"totalvalues"=totalsamples,
                "HXT_YT"=HXT_YT,"HXF_YT"=HXF_YT,"HXT_YF"=HXT_YF,"HXF_YF"=HXF_YF,
                "MXT_YT"=MXT_YT,"MXF_YT"=MXF_YT,"MXT_YF"=MXT_YF,"MXF_YF"=MXF_YF,"timestep"=i)

    resultGroup[[length(resultGroup)+1]]<-result
  }

 #subresults<-lapply(resultGroup,function(result)max((result[["TT"]]+result[["supportTT"]]),(result[["TF"]]+result[["supportTF"]]),(result[["FT"]]+result[["supportFT"]]),(result[["FF"]]+result[["supportFF"]])))
  subresults<-lapply(resultGroup,function(result)max(result[["TT"]],result[["TF"]],result[["FT"]],result[["FF"]]))
 if(length(subresults)<=1)
 {
   return(resultGroup[[1]])
 }else
 {
   bestindexs<-which(subresults==max(unlist(subresults)))
   bestgroup<-resultGroup[bestindexs]
   subresults<-lapply(bestgroup,function(result)max(result[["supportTT"]],result[["supportTF"]],result[["supportFT"]],result[["supportFF"]]))
   if(length(subresults)<=1)
   {
     return(bestgroup[[1]])
   }else
   {
     bestindex<-which(subresults==max(unlist(subresults)))
     return(bestgroup[[bestindex]])
   }
 }

}


generateTemporalGeneStates<-function(timeSeriesCube,targetgene,conditionalgenes,temporal=1,geneStates=NULL,preGeneStates=NULL)
{
  genestates<-list()
  #Test 1

  if(is.null(geneStates))
  {
    getCurrentState<-extractGeneStateFromTimeSeriesCube(timeSeriesCube)
  }else
  {
    getCurrentState<-geneStates
  }

  if(is.null(preGeneStates))
  {
    getpreviousState<-extractGeneStateFromTimeSeriesCube(timeSeriesCube)
  }else
  {
    getpreviousState<-preGeneStates
  }

  index<-temporal

  stateindex<-1
  while(index>0)
  {

    extracedCondition<-extractGenesFromTimeSeriesCube(timeSeriesCube,conditionalgenes,getpreviousState)
    extracedTargetGeneState<-extractGenesFromTimeSeriesCube(timeSeriesCube,c(targetgene),getCurrentState)

    genestates[[stateindex]]<-pairlist("ConditionState"=extracedCondition,"TargetState"=extracedTargetGeneState)
    index<-index-1
    stateindex<-stateindex+1

    getCurrentState<-getCurrentState[,-1]#remove the first one
    lastindex<-length(colnames(getpreviousState))
    getpreviousState<-getpreviousState[,-lastindex] #remove the last one
  }

  return(genestates)
}

#'An utility function to verify whether or not the input expression and the required gene input are marched
#'
#'@param geneState Pre gene input state
#'@param expression The expression of the target regulatory function
#'@return TRUE or FALSE
#'@examples
#' coming later
#' @export
isSatisfied<-function(geneState,expression)
{
  res<-TRUE
  #inistate<-geneState[[1]][[1]]
  if(identical(expression,"1"))
  {
    #over expressed?
    return(TRUE)
  }

  if(identical(expression,"0"))
  {
    #over inhibited
    return(FALSE)
  }

  for(i in seq_along(geneState))
  {
    index<-which(expression==names(geneState)[[i]])[1]
    if(length(index)==0 | length(index)>1)
    {
      res<-FALSE
    }

    if(index>1)
    {
      #if contain negation
      if(identical(expression[index-1],"!"))
      {
        res<-res& !(as.numeric(geneState[[i]])==1L)
      }else
      {
        if(!identical(expression[index],"!")&!identical(expression[index],"&")&!identical(expression[index],"(")&!identical(expression[index],")"))
        {
          res<-res& (as.numeric(geneState[[i]])==1L)
        }
      }
    }
    else
    {
      if(!identical(expression[index],"!")&!identical(expression[index],"&")&!identical(expression[index],"(")&!identical(expression[index],")"))
      {
        res<-res& (as.numeric(geneState[[i]]) ==1L)
      }
    }
  }
  return (res)
}

getGeneInputState<-function(geneInput,expression)
{
  res<-list()
  #inistate<-geneState[[1]][[1]]
  if(identical(expression,"1"))
  {
    res[[1]]<-1
    names(res)[[1]]<-geneInput[[1]]
    #over expressed?
    return(res)
  }

  if(identical(expression,"0"))
  {
    #over inhibited
    res[[1]]<-0
    names(res)[[1]]<-geneInput[[1]]
    #over expressed?
    return(res)
  }

  for(i in seq_along(geneInput))
  {
    gene<-geneInput[[i]]
    index<-which(expression==gene)[1]
    if(length(index)==0 | length(index)>1)
    {
      res[[i]]<-0
    }

    if(index>1)
    {
      #if contain negation
      if(identical(expression[index-1],"!"))
      {
        res[[i]]<-0
      }else
      {
        if(!identical(expression[index],"!")&!identical(expression[index],"&")&!identical(expression[index],"(")&!identical(expression[index],")"))
        {
          res[[i]]<-1
        }
      }
    }
    else
    {
      if(!identical(expression[index],"!")&!identical(expression[index],"&")&!identical(expression[index],"(")&!identical(expression[index],")"))
      {
        res[[i]]<-1
      }
    }
    names(res)[[i]]<-gene
  }
  return (res)
}

'do we need this'
getProbabilityFromGeneCube<-function(FBNGenecube,targetGene,geneInputState,type)
{
  getInputGenes<-names(geneInputState)

  getInitStem<-FBNGenecube[[targetGene]]$SubGenes
  getInitStatis<-FBNGenecube[[targetGene]]$Statistic


  #need to conside the different combinations
  res<-c(0)
  probability<-0
  for(k in seq_along(getInputGenes))
  {
    getfirstgene<-getInputGenes[k]
    getsubgenes<-getInputGenes[getInputGenes!=getfirstgene]
    getStem<-getInitStem[[getfirstgene]]
    if(length(getsubgenes)>0)
    {
      if(geneInputState[[getfirstgene]]==1)
      {
        getStem<-getStem$SubGenesT
      }else
      {
        getStem<-getStem$SubGenesF
      }

      for(i in seq_along(getsubgenes))
      {
        subgene<-getsubgenes[[i]]
        getsubStem<-getStem[[subgene]]
        if(is.null(getsubStem))
        {
          break;
        }

        if(geneInputState[[subgene]]==1)
        {
          if(type==1)
          {
            probability<-getsubStem$Statistic$TargetTrue$True[2] #TT
          }else
          {
            probability<-getsubStem$Statistic$TargetFalse$True[2] #FT
          }
        }else
        {
          if(type==1)
          {
            probability<-getsubStem$Statistic$TargetTrue$False[2] #TF
          }else
          {
            probability<-getsubStem$Statistic$TargetFalse$False[2] #FF
          }
        }

        if(i<length(getsubgenes))
        {
          if(geneInputState[[subgene]]==1)
          {
            getStem<-getsubStem$SubGenesT
          }else
          {
            getStem<-getsubStem$SubGenesF
          }
          if(is.null(getStem))
            break;
        }
      }
    }else
    {
      if(geneInputState[[getfirstgene]]==1)
      {
        if(type==1)
        {
          probability<-getStem$Statistic$TargetTrue$True[2]#TT
        }else
        {
          probability<-getStem$Statistic$TargetFalse$True[2]#FT
        }
      }else
      {
        if(type==1)
        {
          probability<-getStem$Statistic$TargetTrue$False[2]#TF
        }else
        {
          probability<-getStem$Statistic$TargetFalse$False[2]#FF
        }
      }
    }
    if(is.na(probability) | is.null(probability))
      probability<-0

    res<-c(res,probability)
  }


  return(max(res))
}
