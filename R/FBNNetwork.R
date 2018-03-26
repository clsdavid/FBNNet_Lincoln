#'Load FBN network data into a FBN network object
#'
#'@param file A file that contains FBN data
#'@param bodySeparator A specific text separator
#'@param lowercaseGenes if all genes are in lower case
#'@return An object of FBN network
#'@examples
#' coming later
#' @export
loadFBNNetwork<-function (file, bodySeparator = ",", lowercaseGenes = FALSE)
{
  #internal function
  matchNames <- function(rule)
  {
    regexpr <- "([_a-zA-Z][_a-zA-Z0-9]*)[,| |\\)|\\||\\&|\\[]"
    rule <- paste(gsub(" ", "", rule, fixed=TRUE)," ",sep="")
    res <- unique(unname(sapply(regmatches(rule,gregexpr(regexpr, rule))[[1]],
                                function(m)
                                {
                                  sapply(regmatches(m,regexec(regexpr,m)),function(x)x[2])
                                })))

    # remove operators
    isOp <- sapply(res, function(x)
    {
      tolower(x) %in% c("all", "any",
                        "sumis", "sumgt", "sumlt",
                        "maj", "timegt", "timelt", "timeis")
    })

    return(res[!isOp])
  }

    func <- readLines(file, -1)
    func <- gsub("#.*", "", trimws(func))
    func <- func[nchar(func) > 0]
    if (length(func) == 0)
        stop("Header expected!")

    header <- func[1]
    header <- tolower(trimws(strsplit(header, bodySeparator)[[1]]))

    if (length(header) < 3 || header[1] != "targets" || !(header[2] %in%
        c("functions", "factors")) || header[3] != "type")
        stop(paste("Invalid header:", func[1]))

    func <- func[-1]
    if (lowercaseGenes)
        func <- tolower(func)
    func <- gsub("[^\\[\\]a-zA-Z0-9_\\|\\&!\\(\\) \t\\-+=.,]+",
        "_", func, perl = TRUE)
    tmp <- unname(lapply(func, function(x) {
        bracketCount <- 0
        lastIdx <- 1
        chars <- strsplit(x, split = "")[[1]]
        res <- c()
        if (length(chars) > 0) {
            for (i in seq_along(chars)) {
                if (chars[i] == "(")
                  bracketCount <- bracketCount + 1
                else if (chars[i] == ")")
                  bracketCount <- bracketCount - 1
                else if (chars[i] == bodySeparator && bracketCount ==
                  0) {
                  res <- c(res, trimws(paste(chars[lastIdx:(i -
                    1)], collapse = "")))
                  lastIdx <- i + 1
                }
            }
            res <- c(res, trimws(paste(chars[lastIdx:length(chars)],
                collapse = "")))
        }
        return(res)
    }))

    targets <- sapply(tmp, function(rule) trimws(rule[1]))
    for (target in targets) {
        if (regexec("^[a-zA-Z_][a-zA-Z0-9_]*$", target)[[1]] ==
            -1)
            stop(paste("Invalid gene name:", target))
    }

    factors <- sapply(tmp, function(rule) trimws(rule[2]))

    types <-sapply(tmp, function(rule) as.numeric(rule[3]))

    temporal <- length(grep("timeis|timelt|timegt|\\[|\\]", factors,
        ignore.case = TRUE) > 0)
    if (temporal && !symbolic) {
        warning("The network contains temporal elements. This requires loading the model with symbolic=TRUE!")
        symbolic <- TRUE
    }

    factors.tmp <- lapply(factors,matchNames)

    genes <- unique(c(targets, unname(unlist(factors.tmp))))


    fixed <- rep(-1, length(genes))
    names(fixed) <- genes
    interactions <- list()


    for (i in seq_along(targets)) {
        target <- targets[i]
        interaction <- generateFBNInteraction(factors[i], genes)
        if (length(interaction$func) == 1) {
           fixed[target] <- interaction$func
        }
        interaction[[length(interaction)+1]]<-types[i]
        names(interaction)[[length(interaction)]]<-"type"
        interactions[[target]][[length(interactions[[target]]) + 1]] <- interaction
    }
    onlyInputs <- setdiff(genes, targets)
    if (length(onlyInputs) > 0) {
        for (gene in onlyInputs) {
            warning(paste("There is no transition function for gene \"",
              gene, "\"! Assuming an input!", sep = ""))

            interactions[[gene]] = list(input = length(interactions) + 1, func = c(0, 1), expression = gene)
        }
    }

    res <- list(interactions = interactions, genes = genes,
        fixed = fixed)

    class(res)<-c("BooleanNetworkCollection")

    #return(convertToFBNNetwork(res))
    return(res)
}

#' @export
convertToBooleanNetworkCollection<-function(network)
{
  #validate network types
  if (!inherits(network,"BooleanNetwork"))
    stop("Network1 must be inherited from BooleanNetwork")

  merges<-mergeInteraction(lapply(network$interactions,convertInteraction),lapply(network$interactions,convertInteraction))
  res <- list(interactions = merges, genes = unique(c(network$genes,network$genes)),
              fixed = unique(c(network$fixed,network$fixed)))

  class(res) <- "BooleanNetworkCollection"
  return(res)
}

#' @export
mergeNetwork<-function(network1,network2)
{
  #validate network types
  if (!(inherits(network1,"BooleanNetwork")||inherits(network1,"BooleanNetworkCollection")))
    stop("Network1 must be inherited from BooleanNetwork or BooleanNetworkCollection")

  if (!(inherits(network2,"BooleanNetwork")||inherits(network2,"BooleanNetworkCollection")))
    stop("Network2 must be inherited from BooleanNetwork or BooleanNetworkCollection")

  merges<-mergeInteraction(lapply(network1$interactions,convertInteraction),lapply(network2$interactions,convertInteraction))

  res <- list(interactions = merges, genes = unique(c(network1$genes,network2$genes)),
              fixed = unique(c(network1$fixed,network2$fixed)))

  class(res) <- "BooleanNetworkCollection"
  return(res)
}

#' @export
convertInteraction<-function(interaction)
{
  res<-list()
  if(mode(interaction)=="list")
  {
    res<-interaction
  }else if(mode(interaction)=="pairlist")
  {
    res<-list(list(input=interaction$input,func=interaction$func,expression=interaction$expression,error=NA))
  }
  return(res)
}

#' @export
mergeInteraction<-function(interactions1,interactions2)
{
  res<-list()
  entry<-list()

  for(name1 in names(interactions1))
  {
    if(!(name1 %in% names(res)))
    {
      res[length(res)+1]<-interactions1[name1]
      names(res)[length(res)]<-name1
    }

    if(name1 %in% names(interactions2))
    {
      duplicateFound<-FALSE
      for(j in length(interactions2[name1]))
      {
        for(i in length(res[name1]))
        {
          if(all(res[name1][[i]]$input==interactions2[name1][[j]]$input))
          {
            duplicateFound<- TRUE
            break;
          }
        }

        if(!duplicateFound)
        {
          res[name1][[length(res[name1])+1]]<-list(input=interactions2[name1][[j]]$input,func=interactions2[name1][[j]]$func,expression=interactions2[name1][[j]]$expression,error=interactions2[name1][[j]]$error)
        }
      }
    }
  }


  for(name2 in names(interactions2))
  {
    if(!(name2 %in% names(res)))
    {
      res[length(res)+1]<-interactions2[name2]
      names(res)[length(res)]<-name2
    }

    if(name2 %in% names(interactions1))
    {
      duplicateFound<-FALSE
      for(j in length(interactions1[name2]))
      {
        for(i in length(res[name2]))
        {
          if(all(res[name2][[i]]$input==interactions1[name2][[j]]$input))
          {
            duplicateFound<- TRUE
            break;
          }
        }

        if(!duplicateFound)
        {
          res[name2][[length(res[name2])+1]]<-list(input=interactions1[name2][[j]]$input,func=interactions1[name2][[j]]$func,expression=interactions1[name2][[j]]$expression,error=interactions1[name2][[j]]$error)
        }
      }
    }
  }
  return(res)
}

#' @export
mergeClusterNetworks<-function(clusteredFBNCube)
{
  #validate network types
  if (!(inherits(clusteredFBNCube,"ClusteredFBNCube")))
    stop("clusteredFBNCube must be inherited from ClusteredFBNCube")

  len<-length(clusteredFBNCube)
  i<-1
  res<-clusteredFBNCube[[1]]$Network
  i<-i+1
  while(i<=len)
  {
    res<-mergeNetwork(res,clusteredFBNCube[[i]]$Network)
    i<-i+1
  }
  return(convertToFBNNetwork(res))
}

#' @export
filterNetworkConnections<-function(networks)
{
  res<-networks
  filterednetworks<-res$interactions[lapply(res$interactions,length)>0]
  res$interactions<-filterednetworks
  return(res)
}
