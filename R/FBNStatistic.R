
#'Find all genes that are related to the target gene
#'
#'@param targetGene The target gene
#'@param currentState current combined gene state i.e. the gene state at t
#'@param previousState  previous combined gene state i.e. the gene state at t-1
#'@param timeseriesCube a list of samples in which a sample contains gene and time points
#'@param genes a list of genes that the related genes re mined from
#'@param matchedgenes a list of pre matched genes
#'@param matchedexpression previous matched expression string
#'@param maxK The maximum depth the cube can go to. The default value is 5
#'@param temporal 1: d means temporal = t-1 ... t-d
#'@param allRelatedTargetGeneMeasures all measures that are related to the target genes.
#'@param allRelatedConditionalGenesMeasures all measures that are related to the conditional genes
#'@param An option to turn on parallel
#'@return A list of measures
#'@examples
#' mat1<-matrix(c("1","0","0","1","0","0","0","1","1"),3,3, dimnames=list(c("gene1","gene2","gene3"),c("1","2","3")))
#' mat2<-matrix(c("1","1","0","1","0","1","1","1","0"),3,3, dimnames=list(c("gene1","gene2","gene3"),c("1","2","3")))
#' listtest<-list(mat1,mat2)
#' completeState<-extractGeneStateFromTimeSeriesCube(listtest)
#' getCurrentState<-completeState
#' getpreviousState<-completeState
#' getCurrentState<-getCurrentState[,-1]#remove the first one
#' getpreviousState<-getpreviousState[,-length(colnames(getpreviousState))] #remove the last one
#' clusteringOnTargetGene("gene1",getCurrentState,getpreviousState,listtest,c("gene1","gene2","gene3"),1,FALSE)
#' @export
buildProbabilityTreeOnTargetGene<-function(targetGene,currentState,previousState,timeseries,genes,matchedgenes,matchedexpression,maxK=5,temporal=1,allRelatedTargetGeneMeasures,allRelatedConditionalGenesMeasures,useParallel=FALSE)
{
  #tree mining stop when 1 or 0

  res<-list()
  readConfigFile()
  internalloop<-function(g,pgenes,unprocessedGenes,ptimeseriesCube,ptargetGene,currentState,previousState,pmatchedgenes,pmatchedexpression,pmaxK,pPrePrababilities,temporal,allRelatedTargetGeneMeasures,allRelatedConditionalGenesMeasures){
    try({
      res<-list()
      gene<-pgenes[g]
      timestep<-1 #default as synchronous
      res[[g]]<-list()
      names(res)[[g]]<-gene
      expression<-c()
      expressionT<-c()
      expressionF<-c()
      preprocessed<-c()
      isFirst<-FALSE

      if(!is.null(pmatchedgenes))
      {
        newMatchedGenesT<-pmatchedgenes
        newMatchedGenesF<-pmatchedgenes
        preprocessed<-names(pmatchedgenes)
      }else
      {
        newMatchedGenesT<-list()
        newMatchedGenesF<-list()
        isFirst<-TRUE
      }

      #update for next level
      if(!is.null(pmatchedexpression))
      {
        expression<-splitExpression(pmatchedexpression)
        expressionT<-c(expression,"&",gene)
        expressionF<-c(expression,"&","!",gene)

        if(is.null(newMatchedGenesT[[gene]]))
        {
          index<-length(newMatchedGenesT)+1
          newMatchedGenesT[[index]]<-1
          names(newMatchedGenesT)[[index]]<-gene
        }

        if(is.null(newMatchedGenesF[[gene]]))
        {
          index<-length(newMatchedGenesF)+1
          newMatchedGenesF[[index]]<-0
          names(newMatchedGenesF)[[index]]<-gene
        }

      }else
      {
        newMatchedGenesT<-list(gene=1)
        names(newMatchedGenesT)[[1]]<-gene
        newMatchedGenesF<-list(gene=0)
        names(newMatchedGenesF)[[1]]<-gene
        expressionT<-c(gene)
        expressionF<-c("!",gene)
      }

      expT<-paste(expressionT,sep='',collapse = '')
      expF<-paste(expressionF,sep='',collapse = '')

      inputgenes<-sort(names(newMatchedGenesT))

      #for debug
      testprint<-paste("Target=",ptargetGene," current=",gene, " processed=",paste(names(pmatchedgenes),sep=",",collapse = ","),sep="",collapse = "")
      #print(pmaxK)
      #print(testprint)


      #if(ptargetGene=="Gene2" & gene=="Gene1" & is.null(pmatchedgenes))
      #if(testprint=="Target=Gene2 current=Gene1 processed=Gene1")
      #if(ptargetGene=="Gene2" & gene=="Gene1" & is.null(pmatchedgenes))
      if(testprint=="Target=Gene2 current=Gene5 processed=Gene1,Gene4")
      {
        print(testprint)
      }
      #if(testprint=="Target=E2F current=p27 processed=Rb")
      #{
      #print("Stop")
      #}

      #if(testprint=="Target=E2F current=Rb processed=p27")
      #{
      #print("Stop")
      #}
      #debug end

      if(any(c(gene) %in% preprocessed))
      {
        return(NULL)
      }

      if(!is.null(gene) & !is.na(gene))
      {
        if(!is.null(pmatchedgenes))
        {
          probabilityOfFourCombines<-getGenePrababilities(ptimeseriesCube,pmatchedgenes,ptargetGene,gene,currentState,previousState,temporal)
        }else
        {
          probabilityOfFourCombines<-allRelatedTargetGeneMeasures[[gene]]
        }


        chiSQ<-probabilityOfFourCombines$chiSQ
        fisherTest<-probabilityOfFourCombines$fisherTest
        timestep<-probabilityOfFourCombines$timestep
        #P(conditions(t)|target(t+1))
        TT<-as.numeric(probabilityOfFourCombines$TT)
        FT<-as.numeric(probabilityOfFourCombines$FT)
        TF<-as.numeric(probabilityOfFourCombines$TF)
        FF<-as.numeric(probabilityOfFourCombines$FF)

        supportTT<-as.numeric(probabilityOfFourCombines$supportTT)
        supportFT<-as.numeric(probabilityOfFourCombines$supportFT)
        supportTF<-as.numeric(probabilityOfFourCombines$supportTF)
        supportFF<-as.numeric(probabilityOfFourCombines$supportFF)

        noise_P<-as.numeric(probabilityOfFourCombines$Noise_P)  #min(TT,FT)
        noise_N<-as.numeric(probabilityOfFourCombines$Noise_N)  #min(TF,FF)

        signal_P<-as.numeric(probabilityOfFourCombines$Signal_P)#max(TT,FT)
        signal_N<-as.numeric(probabilityOfFourCombines$Signal_N)#max(TT,FT)

        conditionTCount<-as.numeric(probabilityOfFourCombines$conditionTCount)
        conditionFCount<-as.numeric(probabilityOfFourCombines$conditionFCount)

        #P(target(t+1))
        targetT<-as.numeric(probabilityOfFourCombines$targetT)
        targetF<-as.numeric(probabilityOfFourCombines$targetF)
        conditionT<-as.numeric(probabilityOfFourCombines$conditionT)
        conditionF<-as.numeric(probabilityOfFourCombines$conditionF)

        targetTT<-as.numeric(probabilityOfFourCombines$targetTT)
        targetFT<-as.numeric(probabilityOfFourCombines$targetFT)
        targetTF<-as.numeric(probabilityOfFourCombines$targetTF)
        targetFF<-as.numeric(probabilityOfFourCombines$targetFF)


        targetTCount<-probabilityOfFourCombines$targetTCount
        targetFCount<-probabilityOfFourCombines$targetFCount

        dXT_YT<-as.numeric(probabilityOfFourCombines$dXT_YT)
        dXF_YT<-as.numeric(probabilityOfFourCombines$dXF_YT)
        dXT_YF<-as.numeric(probabilityOfFourCombines$dXT_YF)
        dXF_YF<-as.numeric(probabilityOfFourCombines$dXF_YF)

        lenTT =as.numeric(probabilityOfFourCombines$lenTT)
        lenFT =as.numeric(probabilityOfFourCombines$lenFT)
        lenTF =as.numeric(probabilityOfFourCombines$lenTF)
        lenFF =as.numeric(probabilityOfFourCombines$lenFF)

        isNegativeCorrelated<-probabilityOfFourCombines$isNegativeCorrelated
        isPossitiveCorrelated<-probabilityOfFourCombines$isPossitiveCorrelated

        #first level prunning

        if(is.null(isNegativeCorrelated)|is.null(isPossitiveCorrelated))
        {
          return(NULL)
        }
        if(!isNegativeCorrelated & !isPossitiveCorrelated)
        {
          return(NULL)
        }

        #pvalue to reject null hypothese
        pvalue<-min(chiSQ,fisherTest)

        #we cannot do this because it will upset short time series data
        #if(pvalue>p_value_FBNCube)
        #{
        #return(NULL)
        #}

        #The counter envicende to valify which cause which
        counterevidenceOfFourCombines<-getGenePrababilities(ptimeseriesCube,pmatchedgenes,gene,ptargetGene,currentState,previousState,temporal)

        noise_P_c<-as.numeric(counterevidenceOfFourCombines$Noise_P)  #min(TT,FT)
        noise_N_c<-as.numeric(counterevidenceOfFourCombines$Noise_N)  #min(TF,FF)


        #P(target(t)|conditions(t+1))
        TT_c<-as.numeric(counterevidenceOfFourCombines$TT)
        FT_c<-as.numeric(counterevidenceOfFourCombines$FT)
        TF_c<-as.numeric(counterevidenceOfFourCombines$TF)
        FF_c<-as.numeric(counterevidenceOfFourCombines$FF)

        supportTT_c<-as.numeric(counterevidenceOfFourCombines$supportTT)
        supportFT_c<-as.numeric(counterevidenceOfFourCombines$supportFT)
        supportTF_c<-as.numeric(counterevidenceOfFourCombines$supportTF)
        supportFF_c<-as.numeric(counterevidenceOfFourCombines$supportFF)

        chiSQ_c<-counterevidenceOfFourCombines$chiSQ
        fisherTest_c<-counterevidenceOfFourCombines$fisherTest

        conditionTCount_c<-counterevidenceOfFourCombines$conditionTCount
        conditionFCount_c<-counterevidenceOfFourCombines$conditionFCount
        targetTCount_c<-counterevidenceOfFourCombines$targetTCount
        targetFCount_c<-counterevidenceOfFourCombines$targetFCount
        targetT_c<-counterevidenceOfFourCombines$targetT
        targetF_c<-counterevidenceOfFourCombines$targetF
        lenTT_c =counterevidenceOfFourCombines$lenTT
        lenFT_c =counterevidenceOfFourCombines$lenFT
        lenTF_c =counterevidenceOfFourCombines$lenTF
        lenFF_c =counterevidenceOfFourCombines$lenFF

        maxTT<-round(max(TT,TT_c),5)
        minTT<-round(min(TT,TT_c),5)
        maxTF<-round(max(TF,TF_c),5)
        minTF<-round(min(TF,TF_c),5)
        maxFT<-round(max(FT,FT_c),5)
        minFT<-round(min(FT,FT_c),5)
        maxFF<-round(max(FF,FF_c),5)
        minFF<-round(min(FF,FF_c),5)

        if(!is.null(pPrePrababilities))
        {
          finalTT<-as.numeric(pPrePrababilities$transitionT)*TT
          finalTF<-as.numeric(pPrePrababilities$transitionT)*TF
          finalFT<-as.numeric(pPrePrababilities$transitionF)*FT
          finalFF<-as.numeric(pPrePrababilities$transitionF)*FF

          preSE_T<-as.numeric(pPrePrababilities$Noise_P)
          preSE_F<-as.numeric(pPrePrababilities$Noise_N)
          parentGene<-pPrePrababilities$ParentGene
          presupportP<-as.numeric(pPrePrababilities$supportP)
          presupportN<-as.numeric(pPrePrababilities$supportN)

        }else
        {
          finalTT<-TT * targetT
          finalTF<-TF * targetT
          finalFT<-FT * targetF
          finalFF<-FF * targetF

          preSE_T<-noise_P
          preSE_F<-noise_N
          parentGene<-NULL

          presupportP<-NULL
          presupportN<-NULL

        }

        #Test

        #rateTT<-round((TT*supportTT)/(TT_c*supportTT_c),2)#very important and effective measures, >=1 means the cause relationship is from the condition to target
        rateTT<-round(TT/TT_c,2)
        if(is.infinite(rateTT))
        {
          rateTT<-ifelse(TT==1,100,TT)
        }

        if(is.na(rateTT))
        {
          rateTT<-0
        }


        #rateTF<-round((TF*supportTF)/(TF_c*supportTF_c),2)#very important and effective measures
        rateTF<-round(TF/TF_c,2)
        if(is.infinite(rateTF))
        {
          rateTF<-ifelse(TF==1,100,TF)
        }

        if(is.na(rateTF))
        {
          rateTF<-0
        }

        #rateFT<-round((FT*supportFT)/(FT_c*supportFT_c),2)#very important and effective measures
        rateFT<-round(FT/FT_c,2)
        if(is.infinite(rateFT))
        {
          rateFT<-ifelse(FT==1,100,FT)
        }

        if(is.na(rateFT))
        {
          rateFT<-0
        }


        #rateFF<-round((FF*supportFF)/(FF_c*supportFF_c),2)#very important and effective measures
        rateFF<-round(FF/FF_c,2)
        if(is.infinite(rateFF))
        {
          rateFF<-ifelse(FF==1,100,FF)
        }

        if(is.na(rateFF))
        {
          rateFF<-0
        }

        #transition probabilities from the parent
        prePrababilitiesT=list("transitionT"=finalTT,"transitionF"=finalFT, "supportP"=supportTT,"supportN"= supportFT, "Noise_P"=noise_P,"Noise_N"=noise_N,"ParentGene"=gene)
        prePrababilitiesF=list("transitionT"=finalTF,"transitionF"=finalFF, "supportP"=supportTF,"supportN"= supportFF, "Noise_P"=noise_P,"Noise_N"=noise_N,"ParentGene"=gene)


        cosineTT<-sqrt(TT * TT_c)
        cosineTF<-sqrt(TF * TF_c)
        cosineFT<-sqrt(FT * FT_c)
        cosineFF<-sqrt(FF * FF_c)

        #calculate the scores
        #scoreP<-5
        #scoreN<-5

        Noise_P_ZToY<-as.numeric(preSE_T)
        Noise_N_ZToY<-as.numeric(preSE_F)
        HNoiseP<-0
        HNoiseN<-0

        stopForNextLevel<- FALSE #!isCausedByConditions

        #ARACNE test mutual information remove if I(X,Y)<= min(I(X,Z),I(Z,Y))-tolerance threshold
        if(!is.null(parentGene))
        {
          #P(A&B)>P(A)P(B)
          if(!is.null(allRelatedConditionalGenesMeasures[[parentGene]][[gene]]))
          {
            XToZ<-allRelatedConditionalGenesMeasures[[parentGene]][[gene]]
          }else
          {
            XToZ<-getGenePrababilities(ptimeseriesCube,NULL,parentGene,gene,currentState,previousState,temporal)
          }
          Noise_P_XToZ<-as.numeric(XToZ$Noise_P)
          Noise_N_XToZ<-as.numeric(XToZ$Noise_N)

          XToY<-allRelatedTargetGeneMeasures[[gene]]
          Noise_P_XToY<-as.numeric(XToY$Noise_P)
          Noise_N_XToY<-as.numeric(XToY$Noise_N)


          #Information theory https://en.wikipedia.org/wiki/Entropy
          #information entropy
          sp1<-noise_P * log(noise_P)
          sp2<-Noise_P_ZToY * log(Noise_P_ZToY)
          sp3<-Noise_P_XToZ * log(Noise_P_XToZ)
          sp4<-Noise_P_XToY * log(Noise_P_XToY)

          sn1<-noise_N * log(noise_N)
          sn2<-Noise_N_ZToY * log(Noise_N_ZToY)
          sn3<-Noise_N_XToZ * log(Noise_N_XToZ)
          sn4<-Noise_N_XToY * log(Noise_N_XToY)

          HNoiseP<--1 * sum(sp1,sp2,sp3,sp4)#ientropy_P)
          HNoiseN<--1 * sum(sn1,sn2,sn3,sn4)#ientropy_N)

          HNoiseP<-ifelse(is.na(HNoiseP),0,HNoiseP)
          HNoiseN<-ifelse(is.na(HNoiseN),0,HNoiseN)

        }

        if(max(TT,TF)==TT)
        {
          pickT<-TT
          pickTRate<-rateTT
          pickTMutual<-dXT_YT
        }else
        {
          pickT<-TF
          pickTRate<-rateTF
          pickTMutual<-dXF_YT
        }

        if(max(FT,FF)==FT)
        {
          pickF<-FT
          pickFRate<-rateFT
          pickFMutual<-dXT_YF
        }else
        {
          pickF<-FF
          pickFRate<-rateFF
          pickFCosine<-cosineFF
          pickFMutual<-dXF_YF
        }

        #Third level prunning
        #(scoreP==0 & scoreN==0)
        #{
        #return(NULL)
        #}

        #Forth level prunning
        #very important
        if(HNoiseP>0 & HNoiseN>0)
        {
          stopForNextLevel<-TRUE
        }

        subindex<-length(res[[g]])+1


        errorTdata<-c(noise_P)
        errorFdata<-c(noise_N)
        signalTdata<-c(signal_P)
        signalFdata<-c(signal_N)

        subresultT<-NULL
        subresultF<-NULL
        if(!stopForNextLevel)#primary selection
        {
          if(pmaxK>1)
          {
            pmaxK<-pmaxK-1
            exlcudedSubgenes<-c(names(newMatchedGenesT))
            nextGenes<-unprocessedGenes[!unprocessedGenes%in%exlcudedSubgenes]
            if(length(nextGenes)>0)
            {
              if(pmaxK>length(nextGenes))
              {
                pmaxK<-length(nextGenes)
              }

              res[[g]][[subindex]]<-list()
              if(isFirst)
              {
                subresultT<-doNonParallelWorkDecrease(internalloop,nextGenes,nextGenes,ptimeseriesCube,ptargetGene,currentState,previousState,newMatchedGenesT,expT,pmaxK,prePrababilitiesT,temporal,allRelatedTargetGeneMeasures,allRelatedConditionalGenesMeasures)
              }else{
                subresultT<-doNonParallelWork(internalloop,nextGenes,nextGenes,ptimeseriesCube,ptargetGene,currentState,previousState,newMatchedGenesT,expT,pmaxK,prePrababilitiesT,temporal,allRelatedTargetGeneMeasures,allRelatedConditionalGenesMeasures)
              }

              if(!is.null(subresultT)&length(subresultT)>0)
              {
                subAvgErrorT<-unlist(lapply(subresultT,function(rule)
                  return (as.numeric(rule$Statistic$avgErrorP))
                ))

                subAvgErrorF<-unlist(lapply(subresultT,function(rule)
                  return (as.numeric(rule$Statistic$avgErrorN))
                ))

                subAvgSignalT<-unlist(lapply(subresultT,function(rule)
                  return (as.numeric(rule$Statistic$avgSignalP))
                ))

                subAvgSignalF<-unlist(lapply(subresultT,function(rule)
                  return (as.numeric(rule$Statistic$avgSignalN))
                ))

                errorTdata<-c(errorTdata,subAvgErrorT)
                errorFdata<-c(errorFdata,subAvgErrorF)
                signalTdata<-c(signalTdata,subAvgSignalT)
                signalFdata<-c(signalFdata,subAvgSignalF)

              }
            }
            exlcudedSubgenes<-c(names(newMatchedGenesF))
            nextGenes<-unprocessedGenes[!unprocessedGenes%in%exlcudedSubgenes]
            if(length(nextGenes)>0)
            {
              if(pmaxK>length(nextGenes))
              {
                pmaxK<-length(nextGenes)
              }

              res[[g]][[subindex]]<-list()
              if(isFirst)
              {
                subresultF<-doNonParallelWorkDecrease(internalloop,nextGenes,nextGenes,ptimeseriesCube,ptargetGene,currentState,previousState,newMatchedGenesF,expF,pmaxK,prePrababilitiesF,temporal,allRelatedTargetGeneMeasures,allRelatedConditionalGenesMeasures)
              }else
              {
                subresultF<-doNonParallelWork(internalloop,nextGenes,nextGenes,ptimeseriesCube,ptargetGene,currentState,previousState,newMatchedGenesF,expF,pmaxK,prePrababilitiesF,temporal,allRelatedTargetGeneMeasures,allRelatedConditionalGenesMeasures)
              }


              if(!is.null(subresultF)&length(subresultF)>0)
              {
                subAvgErrorT<-unlist(lapply(subresultF,function(rule)
                  return (as.numeric(rule$Statistic$avgErrorP))
                ))

                subAvgErrorF<-unlist(lapply(subresultF,function(rule)
                  return (as.numeric(rule$Statistic$avgErrorN))
                ))

                subAvgSignalT<-unlist(lapply(subresultF,function(rule)
                  return (as.numeric(rule$Statistic$avgSignalP))
                ))

                subAvgSignalF<-unlist(lapply(subresultF,function(rule)
                  return (as.numeric(rule$Statistic$avgSignalN))
                ))


                errorTdata<-c(errorTdata,subAvgErrorT)
                errorFdata<-c(errorFdata,subAvgErrorF)
                signalTdata<-c(signalTdata,subAvgSignalT)
                signalFdata<-c(signalFdata,subAvgSignalF)
              }
            }
          }
        }

        avgErrorT<-sum(errorTdata)/length(errorTdata)
        avgErrorF<-sum(errorFdata)/length(errorFdata)

        avgSignalT<-sum(signalTdata)/length(signalTdata)
        avgSignalF<-sum(signalFdata)/length(signalFdata)

        stdErrorT<-sd(errorTdata)
        stdErrorF<-sd(errorFdata)

        stdSignalT<-sd(signalTdata)
        stdSignalF<-sd(signalFdata)

        stdErrorT<-ifelse(is.na(stdErrorT),0,stdErrorT)
        stdErrorF<-ifelse(is.na(stdErrorF),0,stdErrorF)
        stdSignalT<-ifelse(is.na(stdSignalT),0,stdSignalT)
        stdSignalF<-ifelse(is.na(stdSignalF),0,stdSignalF)


        #TTvalue<-max(ifelse(avgSignalT+stdSignalT>1,1,avgSignalT+stdSignalT),TT)

        #TFvalue<-max(ifelse(avgSignalT+stdSignalT>1,1,avgSignalT+stdSignalT),TF)

        #FTvalue<-max(ifelse(avgSignalF+stdSignalF>1,1,avgSignalF+stdSignalF),FT)

        #FFvalue<-max(ifelse(avgSignalF+stdSignalF>1,1,avgSignalF+stdSignalF),FF)

        noiseP_value<-round(min(ifelse(avgErrorT+stdErrorT>1,1,avgErrorT+stdErrorT),noise_P),5)
        noiseN_value<-round(min(ifelse(avgErrorF+stdErrorF>1,1,avgErrorF+stdErrorF),noise_N),5)

        #TTvalue<-TT#ifelse(testTT>0.1,TT,0)

        #TFvalue<-TF#ifelse(testTF>0.1,TF,0)

        #FTvalue<-FT#ifelse(testFT>0.1,FT,0)

        #FFvalue<-FF#ifelse(testFF>0.1,FF,0)

        #noiseP_value<-noise_P
        #noiseN_value<-noise_N


        dontGoNextLevelForT<-FALSE
        dontGoNextLevelForF<-FALSE


        #Activator and Inhibitor
        res[[g]][[subindex]]<-list()
        names(res[[g]])[[subindex]]<-"ActivatorAndInhibitor"
        if(pickT==TT)
        {
          pattern<-sapply(inputgenes,function(gene)paste(gene,"$",newMatchedGenesT[[gene]],sep="",collapse = ""))
          identity<-paste(pattern,sep="_",collapse = "_")
          identity<-paste(identity,"_Activator",sep="",collapse = "")

          res[[g]][[subindex]][[1]]<-c("factor"=expT,"Confidence"=round(TT,5),"ConfidenceCounter"=round(TT_c,5),"all_confidence"=minTT,"max_confidence"=maxTT, "TransitionP"=round(finalTT,6),"cosine"=round(cosineTT,3),"support"=round(supportTT,5),"rate"=rateTT, "Noise"=round(noiseP_value,5),"Entropy"=round(HNoiseP,5),"Identity"=identity,"MutualInfo"=round(dXT_YT,5),"type"="TT","timestep"=timestep)
          names(res[[g]][[subindex]])[[1]]<-"Activator"
        }else
        {
          pattern<-sapply(inputgenes,function(gene)paste(gene,"$",newMatchedGenesF[[gene]],sep="",collapse = ""))
          identity<-paste(pattern,sep="_",collapse = "_")
          identity<-paste(identity,"_Activator",sep="",collapse = "")

          res[[g]][[subindex]][[1]]<-c("factor"=expF,"Confidence"=round(TF,5),"ConfidenceCounter"=round(TF_c,5),"all_confidence"=minTF,"max_confidence"=maxTF,"TransitionP"=round(finalTF,6), "cosine"=round(cosineTF,3),"support"=round(supportTF,5),"rate"=rateTF, "Noise"=round(noiseP_value,5),"Entropy"=round(HNoiseP,5),"Identity"=identity,"MutualInfo"=round(dXF_YT,5),"type"="TF","timestep"=timestep)
          names(res[[g]][[subindex]])[[1]]<-"Activator"
        }

        if(pickF==FT)
        {
          pattern<-sapply(inputgenes,function(gene)paste(gene,"$",newMatchedGenesT[[gene]],sep="",collapse = ""))
          identity<-paste(pattern,sep="_",collapse = "_")
          identity<-paste(identity,"_Inhibitor",sep="",collapse = "")

          res[[g]][[subindex]][[2]]<-c("factor"=expT,"Confidence"=round(FT,5),"ConfidenceCounter"=round(FT_c,5),"all_confidence"=minFT,"max_confidence"=maxFT,"TransitionP"=round(finalFT,6),"cosine"=round(cosineFT,3),"support"=round(supportFT,5),"rate"=rateFT, "Noise"=round(noiseN_value,5),"Entropy"=round(HNoiseN,5),"Identity"=identity,"MutualInfo"=round(dXT_YF,5),"type"="FT","timestep"=timestep)
          names(res[[g]][[subindex]])[[2]]<-"Inhibitor"
        }else
        {
          pattern<-sapply(inputgenes,function(gene)paste(gene,"$",newMatchedGenesF[[gene]],sep="",collapse = ""))
          identity<-paste(pattern,sep="_",collapse = "_")
          identity<-paste(identity,"_Inhibitor",sep="",collapse = "")

          res[[g]][[subindex]][[2]]<-c("factor"=expF,"Confidence"=round(FF,5),"ConfidenceCounter"=round(FF_c,5),"all_confidence"=minFF,"max_confidence"=maxFF,"TransitionP"=round(finalFF,6),"cosine"=round(cosineFF,3),"support"=round(supportFF,5),"rate"=rateFF, "Noise"=round(noiseN_value,5),"Entropy"=round(HNoiseN,5),"Identity"=identity,"MutualInfo"=round(dXF_YF,5),"type"="FF","timestep"=timestep)
          names(res[[g]][[subindex]])[[2]]<-"Inhibitor"
        }
        subindex<-subindex+1

        #add sub genes
        if(!is.null(subresultT))
        {
          res[[g]][[subindex]]<-subresultT
          names(res[[g]])[[subindex]]<-"SubGenesT"
          subindex<-subindex+1
        }


        if(!is.null(subresultF))
        {
          res[[g]][[subindex]]<-subresultF
          names(res[[g]])[[subindex]]<-"SubGenesF"
          subindex<-subindex+1
        }

        #statistic methods
        res[[g]][[subindex]]<-list()
        names(res[[g]])[[subindex]]<-"Statistic"

        statindex<-1



        res[[g]][[subindex]][[statindex]]<-list()
        names(res[[g]][[subindex]])[[statindex]]<-"TargetTrue"
        res[[g]][[subindex]][[statindex]][[1]]<-c("factor"=expT,"Confidence"=round(TT,5),"ConfidenceCounter"=round(TT_c,5),"all_confidence"=minTT,"max_confidence"=maxTT, "support"=round(supportTT,5), "supportCounter"=supportTT_c, "rate"=rateTT,"TransitionP"=round(finalTT,5),"type"="TT","MutualInfo"=round(dXT_YT,5),"timestep"=timestep)
        names(res[[g]][[subindex]][[statindex]])[[1]]<-"True"
        res[[g]][[subindex]][[statindex]][[2]]<-c("factor"=expF,"Confidence"=round(TF,5),"ConfidenceCounter"=round(TF_c,5),"all_confidence"=minTF,"max_confidence"=maxTF,"support"=round(supportTF,5), "supportCounter"=supportTF_c,"rate"=rateTF,"TransitionP"=round(finalTF,5),"type"="TF" ,"MutualInfo"=round(dXF_YT,5),"timestep"=timestep)
        names(res[[g]][[subindex]][[statindex]])[[2]]<-"False"
        statindex<-statindex+1


        res[[g]][[subindex]][[statindex]]<-list()
        names(res[[g]][[subindex]])[[statindex]]<-"TargetFalse"
        res[[g]][[subindex]][[statindex]][[1]]<-c("factor"=expT,"Confidence"=round(FT,5),"ConfidenceCounter"=round(FT_c,5),"all_confidence"=minFT,"max_confidence"=maxFT,"support"=round(supportFT,5), "supportCounter"=supportFT_c,"rate"=rateFT,"TransitionP"=round(finalFT,5),"type"="FT" ,"MutualInfo"=round(dXT_YF,5),"timestep"=timestep)
        names(res[[g]][[subindex]][[statindex]])[[1]]<-"True"
        res[[g]][[subindex]][[statindex]][[2]]<-c("factor"=expF,"Confidence"=round(FF,5),"ConfidenceCounter"=round(FF_c,5),"all_confidence"=minFF,"max_confidence"=maxFF,"support"=round(supportFF,5),  "supportCounter"=supportFF_c,"rate"=rateFF,"TransitionP"=round(finalFF,5),"type"="FF" ,"MutualInfo"=round(dXF_YF,5),"timestep"=timestep)
        names(res[[g]][[subindex]][[statindex]])[[2]]<-"False"
        statindex<-statindex+1

        res[[g]][[subindex]][[statindex]]<-inputgenes
        names(res[[g]][[subindex]])[[statindex]]<-"Input"
        statindex<-statindex+1

        res[[g]][[subindex]][[statindex]]<-noise_P
        names(res[[g]][[subindex]])[[statindex]]<-"Noise_P"
        statindex<-statindex+1

        res[[g]][[subindex]][[statindex]]<-noise_N
        names(res[[g]][[subindex]])[[statindex]]<-"Noise_N"
        statindex<-statindex+1

        res[[g]][[subindex]][[statindex]]<-preSE_T
        names(res[[g]][[subindex]])[[statindex]]<-"preNoise_P"
        statindex<-statindex+1

        res[[g]][[subindex]][[statindex]]<-preSE_F
        names(res[[g]][[subindex]])[[statindex]]<-"preNoise_N"
        statindex<-statindex+1

        res[[g]][[subindex]][[statindex]]<-avgErrorT
        names(res[[g]][[subindex]])[[statindex]]<-"avgErrorP"
        statindex<-statindex+1

        res[[g]][[subindex]][[statindex]]<-avgErrorF
        names(res[[g]][[subindex]])[[statindex]]<-"avgErrorN"
        statindex<-statindex+1

        res[[g]][[subindex]][[statindex]]<-avgSignalT
        names(res[[g]][[subindex]])[[statindex]]<-"avgSignalP"
        statindex<-statindex+1

        res[[g]][[subindex]][[statindex]]<-avgSignalF
        names(res[[g]][[subindex]])[[statindex]]<-"avgSignalN"
        statindex<-statindex+1


        res[[g]][[subindex]][[statindex]]<-stdErrorT
        names(res[[g]][[subindex]])[[statindex]]<-"stdErrorP"
        statindex<-statindex+1

        res[[g]][[subindex]][[statindex]]<-stdErrorF
        names(res[[g]][[subindex]])[[statindex]]<-"stdErrorN"
        statindex<-statindex+1

        res[[g]][[subindex]][[statindex]]<-stdSignalT
        names(res[[g]][[subindex]])[[statindex]]<-"stdSignalP"
        statindex<-statindex+1

        res[[g]][[subindex]][[statindex]]<-stdSignalF
        names(res[[g]][[subindex]])[[statindex]]<-"stdSignalN"
        statindex<-statindex+1


        res[[g]][[subindex]][[statindex]]<-isNegativeCorrelated
        names(res[[g]][[subindex]])[[statindex]]<-"isNegativeCorrelated"
        statindex<-statindex+1
        res[[g]][[subindex]][[statindex]]<-isPossitiveCorrelated
        names(res[[g]][[subindex]])[[statindex]]<-"isPossitiveCorrelated"
        statindex<-statindex+1

        res[[g]][[subindex]][[statindex]]<-chiSQ
        names(res[[g]][[subindex]])[[statindex]]<-"chiSQ"
        statindex<-statindex+1

        res[[g]][[subindex]][[statindex]]<-chiSQ_c
        names(res[[g]][[subindex]])[[statindex]]<-"chiSQ_c"
        statindex<-statindex+1

        res[[g]][[subindex]][[statindex]]<-fisherTest
        names(res[[g]][[subindex]])[[statindex]]<-"fisherTest"
        statindex<-statindex+1

        res[[g]][[subindex]][[statindex]]<-fisherTest_c
        names(res[[g]][[subindex]])[[statindex]]<-"fisherTest_c"
        statindex<-statindex+1

        subindex<-length(res[[g]])+1

        res[[g]][[subindex]]<-list()
        names(res[[g]])[[subindex]]<-"Source"
        res[[g]][[subindex]][[1]]<-c("targetT"=round(targetT,5),"targetTCount"=targetTCount,"SupportCount"=lenTT,"ConditionCount"=conditionTCount,"targetT_c"=round(targetT_c,5),"targetTCount_c"=targetTCount_c,"SupportCount_c"=lenTT_c,"ConditionCount_c"=conditionTCount_c)
        names(res[[g]][[subindex]])[[1]]<-"TT"
        res[[g]][[subindex]][[2]]<-c("targetT"=round(targetT,5),"targetTCount"=targetTCount,"SupportCount"=lenTF,"ConditionCount"=conditionFCount,"targetT_c"=round(targetT_c,5),"targetTCount_c"=targetTCount_c,"SupportCount_c"=lenTF_c,"ConditionCount_c"=conditionFCount_c)
        names(res[[g]][[subindex]])[[2]]<-"TF"
        res[[g]][[subindex]][[3]]<-c("targetF"=round(targetF,5),"targetFCount"=targetFCount,"SupportCount"=lenFT,"ConditionCount"=conditionTCount,"targetF_c"=round(targetF_c,5),"targetFCount_c"=targetFCount_c,"SupportCount_c"=lenFT_c,"ConditionCount_c"=conditionTCount_c)
        names(res[[g]][[subindex]])[[3]]<-"FT"
        res[[g]][[subindex]][[4]]<-c("targetF"=round(targetF,5),"targetFCount"=targetFCount,"SupportCount"=lenFF,"ConditionCount"=conditionFCount,"targetF_c"=round(targetF_c,5),"targetFCount_c"=targetFCount_c,"SupportCount_c"=lenFF_c,"ConditionCount_c"=conditionFCount_c)
        names(res[[g]][[subindex]])[[4]]<-"FF"
        subindex<-subindex+1
        return(res)
      }
    })}

  maxK<-ifelse(maxK>length(genes),length(genes),maxK)

  #print(genes)
  time1<-as.numeric(Sys.time())
  res[[1]]<-list()
  if(useParallel)
  {
    res<-doParallelWork(internalloop,genes,genes,timeseries,targetGene,currentState,previousState,matchedgenes,matchedexpression,maxK,NULL,temporal,allRelatedTargetGeneMeasures,allRelatedConditionalGenesMeasures)

  }else
  {
    res<-doNonParallelWorkDecrease(internalloop,genes,genes,timeseries,targetGene,currentState,previousState,matchedgenes,matchedexpression,maxK,NULL,temporal,allRelatedTargetGeneMeasures,allRelatedConditionalGenesMeasures)
  }
  time2<-as.numeric(Sys.time())

  print(paste("Total cost ", time2-time1, " seconds to process genes for the target gene ",targetGene, sep='',collapse = ''))

  return(res)
}

#'Find all genes that are related to the target gene
#'
#'@param targetGene The target gene
#'@param currentState current combined gene state i.e. the gene state at t
#'@param previousState  previous combined gene state i.e. the gene state at t-1
#'@param timeseriesCube a list of samples in which a sample contains gene and time points
#'@param genes a list of genes that the related genes re mined from
#'@param temporal 1: d means temporal = t-1 ... t-d
#'@return A list of genes that are related to the target genes
#'@examples
#' mat1<-matrix(c("1","0","0","1","0","0","0","1","1"),3,3, dimnames=list(c("gene1","gene2","gene3"),c("1","2","3")))
#' mat2<-matrix(c("1","1","0","1","0","1","1","1","0"),3,3, dimnames=list(c("gene1","gene2","gene3"),c("1","2","3")))
#' listtest<-list(mat1,mat2)
#' completeState<-extractGeneStateFromTimeSeriesCube(listtest)
#' getCurrentState<-completeState
#' getpreviousState<-completeState
#' getCurrentState<-getCurrentState[,-1]#remove the first one
#' getpreviousState<-getpreviousState[,-length(colnames(getpreviousState))] #remove the last one
#' clusteringOnTargetGene("gene1",getCurrentState,getpreviousState,listtest,c("gene1","gene2","gene3"),1,FALSE)
#' @export
clusteringOnTargetGene<-function(targetGene,currentState,previousState,timeseriesCube,genes,temporal=1,useParallel=FALSE)
{

  readConfigFile()
  internalloop<-function(g,pgenes,ptimeseriesCube,ptargetGene,currentState,previousState,temporal){
    try({
      res<-list()
      gene<-pgenes[g]


      if(!is.null(gene) & !is.na(gene))
      {
        probabilityOfFourCombines<-getGenePrababilities(ptimeseriesCube,NULL,ptargetGene,gene,currentState,previousState,temporal)

        chiSQ<-probabilityOfFourCombines$chiSQ
        fisherTest<-probabilityOfFourCombines$fisherTest

        #pvalue to reject null hypothese
        pvalue<-min(chiSQ,fisherTest)

        if(pvalue>p_value_FBNCube)
        {
          return(NULL)
        }

        res[[g]]<-probabilityOfFourCombines
        names(res)[[g]]<-gene
        return(res)
      }
    })}

  res<-list()

  time1<-as.numeric(Sys.time())
  res[[1]]<-list()
  if(useParallel)
  {
    res<-doParallelWork(internalloop,genes,timeseriesCube,targetGene,currentState,previousState,temporal)

  }else
  {
    res<-doNonParallelWork(internalloop,genes,timeseriesCube,targetGene,currentState,previousState,temporal)
  }
  return(res)
}

#'An internal method to get the probability of regulatory function from the cube
#'
#'@param FBNCube An object of orchard cube
#'@param targetgene The target gene
#'@param funcType type of function
#'@param FBNExpression expression of function
#'@param preGeneInputs pre input genes' state
#'@return A probablity of the target regulatory function
#'@examples
#' coming later
#' @export
getProbabilityFromCubeBasedOnInput<-function(targetgene,funcType,FBNExpression,FBNProbability,preGeneInputs)
{
  #internal functions
  isInputStateMatchedFBNFunction<-function(inputstate,expression)
  {
    splitedexpression<-splitExpression(expression)

    res<-isSatisfied(preGeneInputs,splitedexpression)
    return(res)
  }

  if(!funcType%in%c(1,0))
  {
    stop("The function type value must be in the range of \"0\" and \"1\"")
  }

  if(isInputStateMatchedFBNFunction(preGeneInputs,FBNExpression))
  {
    probability<-as.numeric(FBNProbability) #getProbabilityFromGeneCube(FBNCube,targetgene,preGeneInputs,funcType)
  }else
  {
    probability<-0
  }

  return(as.numeric(probability))
}
