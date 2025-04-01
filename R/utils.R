#' ·············································································
#' Building Diseasome 
#' R version 4.4.0 (2024-04-24 ucrt)
#' Dec 2024
#' ·············································································
#' utils.R

##····························································· check.packages()
## Check installed packages and load packages
#' @param pkg Character vector contains package names
check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)>0){
    cat(paste0("\nInstalling ", paste(new.pkg,collapse = ", ")))
    if (!requireNamespace("BiocManager", quietly = TRUE)){
      install.packages("BiocManager")
    }
    bioconductor_packages <- BiocManager::available()
    bioconductor_packages <- bioconductor_packages[bioconductor_packages %in% new.pkg]
    if (length(bioconductor_packages) > 0){
      BiocManager::install(bioconductor_packages, dependencies = TRUE,ask=FALSE)
    }
    new.pkg = new.pkg[!(new.pkg %in% bioconductor_packages)]
    if (length(new.pkg)>0){
      install.packages(new.pkg, dependencies = TRUE)
    }
  }
  res <- lapply(pkg,function(x){
    cat(paste0("\nLoading ",x))
    suppressMessages(require(x,character.only = T))
  })
} 
##··············································································

##····································································norm.log()
## Log Normalization
#' @param data Gene expression data.frame
norm.log<-function(data){
  qx <- as.numeric(quantile(data, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  if (LogC) { data[which(data <= 0)] <- NaN
  data <- log2(data) }
  return(data)
}
##··············································································

##····································································summ.var()
## Summary of gene variance
#' @param data Gene expression data.frame
#' @param freqCut cutoff for the ratio of the most common value to the second 
#' most common value
#' @param uniqueCut the cutoff for the percentage of distinct values out of the
#' number of total samples
#' @param minVar threshold of standard deviation to remove genes
summ.var<-function(data,
                   freqCut = 70/30,
                   uniqueCut = 30,
                   minVar = 0.05){ 
  require(caret)
  
  nzv <- data.frame(
    caret::nearZeroVar(
      x=t(data),
      freqCut = freqCut,
      uniqueCut = uniqueCut,
      saveMetrics = TRUE,
      allowParallel = TRUE),
    "sd"=apply(data,1,sd))
  nzv$thrSD<-ifelse(nzv$sd<=minVar,TRUE,FALSE)
  
  return(nzv)
}
##··············································································

##·····························································mergeExpression()
## Merge expression by gene. Internaly used by annotateGenes()
#' @param gene Gene identifier to annotate
#' @param genome Table with gene-probe set connections
#' @param expressionMatrix Gene expression data.frame
#' @param method Method used to merge expression of several probes pointing to 
#' the same gene. "median", "mean", "max", "sum"
mergeExpression = function(gene, 
                           genome,
                           expressionMatrix,
                           method){ 
  
  probe_ids = genome[genome$toGenes==gene,"fromGenes"] ## Select probes for each gene
  if (length(probe_ids)>1){
    method<-ifelse(method %in% c("median","mean","sum","max"),method,"median")
    switch (method,
            median = { res =  colMedians(as.matrix(expressionMatrix[probe_ids,]))},
            mean = { res =  colMeans2(as.matrix(expressionMatrix[probe_ids,]))},
            max = { res<-apply(as.matrix(expressionMatrix[probe_ids,]),2,function(x){
              max(x,na.rm = TRUE)})},
            sum = { res<-apply(as.matrix(expressionMatrix[probe_ids,]),2,function(x){
              sum(x,na.rm = TRUE)})})
  } else{
    res <- as.numeric(expressionMatrix[probe_ids,])  
  }
  return(res)
}
##··············································································

##·······························································annotateGenes()
## Function to annotate genes from a gene-expression matrix 
#' @param data Gene expression data.frame
#' @param toGenes Identifier of db to annotate genes
#' @param fromGenes Identifier of db of probes/genes used in data
#' @param method Method used to merge expression of several probes pointing to 
#' the same gene. "median", "mean", "max", "sum"
#' @param ensembl ensembl object from biomaRt package
annotateGenes<-function(data,
                        toGenes='external_gene_name',
                        fromGenes='ensembl_gene_id',
                        method = "median",
                        ensembl = NULL){
  require("biomaRt")
  require("parallel")
  require("tidyr")
  
  ## Get Gene annotation database
  if(is.null(ensembl)){
    ensembl = useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  }
  genome = getBM(attributes=c(toGenes,fromGenes),mart = ensembl)
  genome <- genome %>% `colnames<-`(c("toGenes","fromGenes")) %>% 
    replace(.=="",NA) %>% drop_na() %>% filter(fromGenes %in% as.character(rownames(data)))
  
  data = data[genome$fromGenes,]
  finalGenes = unique(genome$toGenes)
  nCores<-ifelse(as.character(Sys.info()["sysname"])=="Windows",1,detectCores())
  
  temp = as.data.frame(do.call("rbind", mclapply(finalGenes, mergeExpression,
                                                 genome = genome,
                                                 expressionMatrix = data,
                                                 method = method,
                                                 mc.cores = nCores)))
  rownames(temp) = finalGenes
  colnames(temp) = colnames(data)
  
  return(temp)
}
##··············································································

##······························································· ClusterStats()
## Get cluster stability
#' @param list.i 
#' @param Dist.list 
#' @param combinations
ClusterStats<-function(list.i, ## permutation list
                       Dist.list, ## List of distance matrices
                       combinations){ ## Matrix with parameter combinations
  require("SNFtool")
  require("NbClust")
  cat("Running ClusterStats")
  
  nDatasets<-as.numeric(unlist(list.i))
  Times = 15; 	# Number of Iterations, usually (10~20)
  
  ListaNB<-lapply(1:nrow(combinations),function(i){
    
    ## Parameter and dataset selection
    K<<-combinations[i,1]
    Alpha<<-combinations[i,2]
    sel.samples<-sample(1:length(Dist.list),nDatasets,replace=F)
    
    tmp.list<-list()
    for(l.i in 1:length(Dist.list)){
      data<-Dist.list[[l.i]]
      w<-affinityMatrix(data,K,Alpha)
      tmp.list[[l.i]]<-w
    }
    tmp.list<-tmp.list[sel.samples]
    
    W = SNF(tmp.list, K, Times)
    Nb<-NbClust(data = W, diss = NULL, distance = "euclidean", min.nc = 2, max.nc = 15,
                method = "complete", index = "ch", alphaBeale = 0.1)
    
    return(Nb$All.index)
  })
  
  ## Get summary
  nbclusts<-lapply(ListaNB,function(x){
    res<-as.numeric(names(x[order(x,decreasing=T)][1]))
  })
  nbclusts<-unlist(nbclusts)
  res<-table(nbclusts)
  
  res<-c(as.numeric(res["2"]),as.numeric(res["3"]),as.numeric(res["4"]),as.numeric(res["5"]),
         as.numeric(res["6"]),as.numeric(res["7"]),as.numeric(res["8"]),as.numeric(res["9"]),
         as.numeric(res["10"]),as.numeric(res["11"]),as.numeric(res["12"]),as.numeric(res["13"]),
         as.numeric(res["14"]),as.numeric(res["15"]))
  res[is.na(res)]<-0
  
  res<-data.frame(cbind(rep(nDatasets,14),2:15,res))
  colnames(res)<-c("datasets","clusters","nbindex")
  
  
  return(res)
  
}


##······························································vectSimilarity()
## Return similarity measurements between two vectors
#' @param x numeric or categorical vector
#' @param y numeric or categorical vector
#' @param method one or several metrics of similarity: "euclidean","maximum",
#' "canberra","binary","minkowski","pearson","spearman","kendall","cosine",
#' "jaccard", "cohenk","adjrand"
vectSimilarity<-function(x,y,method="euclidean"){
  
  require("jaccard")
  require("psych")
  require("lsa")
  require("mclust")
  
  ## Build for other similarity metrics
  results<-NULL
  
  if(any(c("euclidean","maximum","canberra","binary","minkowski") %in% method)){
    tmpMethods<-method[method %in% c("euclidean","maximum","canberra","binary","minkowski")]
    res<-unlist(lapply(tmpMethods,function(met){
      dist(rbind(x,y),method = met)}))
    names(res)<-tmpMethods
    results<-c(results,res)
  }
  
  if(any(c("pearson","spearman","kendall") %in% method)){
    tmpMethods<-method[method %in% c("pearson","spearman","kendall")]
    res<-unlist(lapply(tmpMethods,function(met){
      cor.test(x,y,method = met)$estimate}))
    names(res)<-tmpMethods
    results<-c(results,res)
  }
  
  if("cosine" %in% method){
    results<-c(results,cosine(x,y))
    names(results)[length(results)]<-"cosine"
  }
  
  if("jaccard" %in% method){
    results<-c(results,jaccard(x,y))
    names(results)[length(results)]<-"jaccard"
  }
  if("cohenk" %in% method){
    results<-c(results,cohen.kappa(cbind(x,y))$kappa)
    names(results)[length(results)]<-"cohenk"
  }
  if("adjrand" %in% method){
    results<-c(results,adjustedRandIndex(x,y))
    names(results)[length(results)]<-"adjrand"
  }
  
  return(results)
}


##··············································································


##······························································ jaccard_index()
## Get jaccard index for vector of different length
#' @param set1 vector1
#' @param set2 vector2
jaccard_index <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))
  union <- length(union(set1, set2))
  return(intersection / union)
}

# jaccard_index <- function(a, b) {
#   intersection = length(intersect(a, b))
#   union = length(a) + length(b) - intersection
#   return (intersection/union)
# }
##··············································································


##··································································· getNodes()
#' @param listDB database of pathways (nested list)
#' @param simmilarity.threshold threshold of similarity to remove similar
#' pathways 
#' @param ann.db data.frame contains info about pathDB. Must be contain the columns:
#' "TraitID","AnnotationID","term","source" and "coannotation"
#' @param max.length similarity is only measured for genesets with length <= max.length
getNodes<-function(listDB,ann.db,simmilarity.threshold = 0.8,max.length=NULL){
  
  
  if(!is.null(max.length)){
    tmp.db<- Filter(function(x) length(x) <= max.length, listDB)
  }else{
    tmp.db<-listDB
  }
  
  
  ## Get similarity matrix
  similarity_matrix <- matrix(0, length(tmp.db), length(tmp.db),
                              dimnames = list(names(tmp.db),names(tmp.db)))
  
  total_comparisons <- (length(tmp.db) * (length(tmp.db) - 1)) / 2
  pb <- txtProgressBar(min = 0, max = total_comparisons, style = 3)
  progress <- 0
  for (i in 1:(length(tmp.db)-1)) {
    for (j in (i+1):length(tmp.db)) {
      progress <- progress + 1
      setTxtProgressBar(pb, progress)
      similarity_matrix[i, j] <- jaccard_index(tmp.db[[i]], tmp.db[[j]])
      similarity_matrix[j, i] <- similarity_matrix[i, j]
    }
  }
  
  # Similarity matrix to graph using a threshold (0.8)
  adjacency_matrix <- similarity_matrix > simmilarity.threshold
  graph <- graph_from_adjacency_matrix(adjacency_matrix,
                                       mode = "undirected", diag = FALSE)
  
  # Identify nodes
  components <- clusters(graph)
  groups <- split(names(tmp.db), components$membership)
  
  ## Filter out nodes with only one gene
  groups<-groups[sapply(groups, function(x) length(x) > 1)]
  
  ## Filter step
  cat("\nFiltering paths with high Jaccard Index...")
  removePaths<-NULL
  for(g in 1:length(groups)){
    
    x<-groups[[g]]
    tmp<-tmp.db[x]
    tmp <- tmp[order(sapply(tmp, length),decreasing = T)]
    coann<-names(tmp)
    
    for(i in 1:nrow(ann.db)){
      if(ann.db$TraitID[i] %in% coann){
        if(is.na(ann.db[i,"coannotation"])){
          ann.db[i,"coannotation"]<-paste(coann, collapse = ", ")
        }else{
          ann.db[i,"coannotation"]<-paste0(ann.db[i,"coannotation"],
                                           paste(coann, collapse = ", "),collapse = ", ")
        }
      }
    }
    removePaths<-c(removePaths,coann[-1])
  }
  removePaths<-unique(removePaths)
  
  listDB <- listDB[!names(listDB) %in% removePaths]
  
  return(list("db"=listDB,"ann"=ann.db))
}


##··························································· setPackingFilter()
## Evaluate combinations
#' @param pathDB list of all pathways
#' @param mRef output of pathMED::createReference, used for measure gain information
#' @param ann.db data.frame contains info about pathDB. Must be contain the columns:
#' "TraitID","AnnotationID","term","source" and "coannotation"
#' @param coverage.thr (0-1) Percentage of overlap in set packing
#' @param max.combs Maximum number of elements to combine
#' @param gainIndex Minimal difference in percentage (0-100) for gain information
#' for selection of large pathways or the combination of small pathways 
#' @param minCorr Minimal correlation to consider redundant the small and large paths
setPackingFilter<-function(pathDB,
                           mRef,
                           ann.db,
                           coverage.thr=0.8,
                           max.combs =NULL,
                           gainIndex=10,
                           minCorr=0.75){
  
  paths2Remove<-NULL
  pb <- txtProgressBar(min = 1, max = length(pathDB), style = 3)
  for(i in 1:length(pathDB)){
    setTxtProgressBar(pb, i)
    tmp<-pathDB[[i]]
    
    if(length(tmp)>5){
      tmp.list<-pathDB[sapply(pathDB, function(x) length(x) < length(tmp))]
      tmp.list <- tmp.list[sapply(tmp.list, function(x) all(x %in% tmp))]
      
      if(length(tmp.list)>1){
        
        ## The genes of the small traits complete at least 80% of the large trait
        if(length(unique(unlist(tmp.list))) >= length(tmp)*coverage.thr){
          
          #' How many information is gained using the Large or small traits
          #' % of significant pathways across datasets comparing large pathway vs
          #' small combination of pathways
          difs<-unlist(lapply(1:length(mRef$mscores),function(d){
            m<-mRef$mscores[[d]][c(names(pathDB)[i],names(tmp.list)),]
            
            diff<-((sum(abs(m[1,]) > 1.65)/ncol(m))-
                     (sum(apply(m, 2, function(columna) any(abs(columna) > 1.65)))/ncol(m)) ) *100
          }))
          difs<-mean(difs,na.rm=T) 
          #' dif: mean percentage of patients that loss significance if large path 
          #' is removed
          
          
          #' How correlated are large and small traits
          cors<-unlist(lapply(1:length(mRef$mscores),function(d){
            m<-mRef$mscores[[d]][c(names(pathDB)[i],names(tmp.list)),]
            
            mxCor<-max(unlist(lapply(2:nrow(m),function(r){
              cor.tmp<-cor.test(as.numeric(m[1,]),as.numeric(m[r,]))
              return(as.numeric(cor.tmp$estimate))
            })),na.rm = T)
          }))
          cors[is.na(cors)]<-0
          cors[!is.finite(cors)]<-0
          cors<-median(cors,na.rm=T)
          
          print(paste0(i,"|| Corr: ",round(cors,digits = 2),"| Loss: ",round(difs,digits=2)))
          
          ##···························
          if(difs>(-gainIndex) & cors>=minCorr){
            ## Keep Large path and remove small traits:
            #' Use small traits do not improve significantly the information gain
            
            paths2Remove<-c(paths2Remove,names(tmp.list))  
            
            ## Re-annotation
            coann<-c(names(pathDB)[i],names(tmp.list))
            for(anr in 1:nrow(ann.db)){
              
              if(ann.db$TraitID[anr] %in% coann){
                
                if(is.na(ann.db[anr,"coannotation"])){
                  ann.db[anr,"coannotation"]<-paste(coann, collapse = ", ")
                }else{
                  ann.db[anr,"coannotation"]<-paste0(ann.db[anr,"coannotation"],
                                                     paste(coann, collapse = ", "),collapse = ", ")
                }}}
            
            
            ##·····
          }else{
            if(difs<=(-gainIndex)){
              
              best_combinations <- setPacking(vect.list = tmp.list,
                                              vect= tmp,
                                              coverage.thr = coverage.thr,
                                              max.combs = max.combs)
              
              if(!is.null(best_combinations)){ 
                ## Set packing: Remove the Large path if can be formed using the small traits
                #' Measure if the other small traits are gainful to keep or discard
                
                paths2Remove<-c(paths2Remove,names(pathDB)[i]) 
                coann<-names(pathDB)[i]
                
                selected<-best_combinations[1,"paths"]
                selected<-unlist(strsplit(selected, ", "))
                
                
                #' Information gain
                difs<-unlist(lapply(1:length(mRef$mscores),function(d){
                  
                  m<-mRef$mscores[[d]][names(tmp.list),]
                  
                  diff<-((sum(abs(m[selected,]) > 1.65)/ncol(m))-
                           (sum(apply(m, 2, function(columna) any(abs(columna) > 1.65)))/ncol(m)) ) *100
                }))
                difs<-mean(difs,na.rm=T)
                
                
                if(mean(difs,na.rm=T)>(-gainIndex)){
                  coann<-c(coann,names(tmp.list))
                  paths2Remove<-c(paths2Remove,names(tmp.list)[!names(tmp.list) %in% selected])
                }else{
                  coann<-c(coann,selected)
                }
                
                ## Re-annotation
                for(anr in 1:nrow(ann.db)){
                  
                  if(ann.db$TraitID[anr] %in% coann){
                    
                    if(is.na(ann.db[anr,"coannotation"])){
                      ann.db[anr,"coannotation"]<-paste(coann, collapse = ", ")
                    }else{
                      ann.db[anr,"coannotation"]<-paste0(ann.db[anr,"coannotation"],
                                                         paste(coann, collapse = ", "),collapse = ", ")
                    }}}
                
              } ## Else: Not effective set Packing, keep paths without modifications
              
            } ## Else: difs <10 y no cors, keep paths without modifications
            
          } ##···························
          
        } ## The genes of the small traits complete at least 80% of the large trait
      } ## length(tmp.list)<=1
      
    }
    
  } # for loop
  
  
  ## Filter pathDB and save R Object
  pathDB <- pathDB[!names(pathDB) %in% paths2Remove]
  
  return(list("db"=pathDB,"ann"=ann.db))
}

##··············································································


##································································ .setPacking()
## Function to evaluate how combinations of @vect.list elements form @vect
#' @param vect.list list with small genesets (nested list with vector of genes)
#' @param vect vector with genes from a large geneset
#' @param coverage.thr minumun coverage to consider that @vect is formed
#'  using a combination of vect.list elements
#' @param max.combs Maximum number of vect.list elements to combinate 
setPacking <- function(vect.list, vect,coverage.thr=0.8,max.combs=NULL) {
  require("combinat")
  require("dplyr")
  require("parallel")
  require("arrangements")
  
  if(is.null(max.combs)){
    max.combs<-length(vect.list)
  }else{
    if(max.combs>length(vect.list)){
      max.combs<-length(vect.list)
    }
  }
  ## High Computational demanding. Limit of combinations stablished
  while(choose(length(vect.list), max.combs)>50000000){
    max.combs<-max.combs-1
    print("Large number of combinations")
  }
  
    all_combinations<-lapply(1:max.combs,function(i){
      combs <- combinations(1:length(vect.list), k = i,layout = "list")
      #combs <- combn(vect.list, i, simplify = FALSE)
      #print(i)
      
      filtered_combs <- Filter(function(comb) {
        l.tmp<-vect.list[comb]
        combined_genes <- unique(unlist(l.tmp))
        coverage <- length(intersect(combined_genes, vect)) / length(vect)
        coverage >= coverage.thr  
      }, combs)

      return(filtered_combs)
      
      })
    gc()

  all_combinations <- do.call(c, Filter(function(x) length(x) > 0, all_combinations))
  
  all_combinations<-lapply(all_combinations,function(x){
    return(vect.list[x]) })
  
  ## Get Combination stats
  if(length(all_combinations)>0){
    results <- lapply(all_combinations, .evaluate_combination, vect)
    
    m.results<-as.data.frame(do.call("rbind",lapply(1:length(results),function(x){
      r.x<-results[[x]]
      res<-c(r.x$coverage,r.x$extra_genes,
             paste(names(r.x$comb), collapse = ", "))
      
    })))
    colnames(m.results)<-c("coverage","extraGenes","paths")
    m.results$coverage<-as.numeric(m.results$coverage)
    m.results$extraGenes<-as.numeric(m.results$extraGenes)
    
    m.results<-m.results[order(-m.results$coverage,m.results$extraGenes),]
    m.results <- m.results[m.results$coverage >=coverage.thr, ]  # & m.results$extraGenes==0
    
    
    if(nrow(m.results)!=0){
      return(m.results)
    }else{
      return(NULL)
    }
  }else{
    return(NULL)
  }
}
##··············································································

##······················································ .evaluate_combination()
## Evaluate combinations
#' @param comb combination of genesets (vectors of genes)
#' @param vect vector with genes from a large geneset
.evaluate_combination <- function(comb, vect) {
  combined_genes <- unique(unlist(comb))
  covered_genes <- intersect(combined_genes, vect)
  uncovered_genes <- setdiff(combined_genes, vect)
  coverage <- length(covered_genes) / length(vect)
  extra_genes <- length(uncovered_genes)
  return(list(coverage = coverage, extra_genes = extra_genes, comb = comb))
}
##··············································································


##······················································ Diversity indexes
## Evaluate combinations
shannon_entropy <- function(terms) {
  freq_table <- table(terms)  
  probs <- freq_table / sum(freq_table) 
  entropy <- -sum(probs * log2(probs))
  return(entropy)
}

calculate_ZCI <- function(totalTerms, iterTerms) {
  H_total <- shannon_entropy(totalTerms)
  H_iter <- shannon_entropy(iterTerms)
  
  if (H_total == 0) {
    return(0)
  }
  ZCI <- H_iter / H_total
  return(ZCI)
}
##··············································································


##··············································································
find_best_clusters <- function(nbclust_results,metric) {
  library(dplyr)
  
  # Define mode of action of each index
  metric_types <- c("kl"="low",
                    "ch"="high",
                    "hartigan"="low",
                    "cindex"="low",
                    "db"="low",
                    "silhouette"="high",
                    "ratkowsky"="high",
                    "ball"="low",
                    "ptbiserial"="high",
                    "gap"="high",
                    "frey"="change",
                    "mcclain"="low",
                    "dunn"="high",
                    "hubert"="high",
                    "gamma"="high",
                    "gplus"="high",
                    "tau"="high",
                    "sdindex"="low",
                    "dindex"="high",
                    "sdbw"="low")
  
  metric_type<-metric_types[metric]
  values<-nbclust_results$Allindex
  k_values <- nbclust_results$Kclusters
  
  peaks <- c()  
  if (metric_type %in% c("high", "low")) {
    
    for (i in seq_along(values)) {
      if (i == 1) {  # Primer elemento (solo comparar con el siguiente)
        if ((metric_type == "high" && values[i] > values[i + 1]) ||
            (metric_type == "low" && values[i] < values[i + 1])) {
          peaks <- c(peaks, k_values[i])
        }
      } else if (i == length(values)) {  # Último elemento (solo comparar con el anterior)
        if ((metric_type == "high" && values[i] > values[i - 1]) ||
            (metric_type == "low" && values[i] < values[i - 1])) {
          peaks <- c(peaks, k_values[i])
        }
      } else {  # Valores intermedios (comparar con ambos lados)
        if ((metric_type == "high" && values[i] > values[i - 1] && values[i] > values[i + 1]) ||
            (metric_type == "low" && values[i] < values[i - 1] && values[i] < values[i + 1])) {
          peaks <- c(peaks, k_values[i])
        }
      }
    }
    
  } else if (metric_type == "change") {
    diffs <- abs(diff(values))  
    threshold <- quantile(diffs, 0.75, na.rm = TRUE) 
    change_points <- which(diffs >= threshold) + 1 
    
    peaks <- k_values[change_points]
  }
  
  return(peaks)
}
##··············································································


##···································································getSankey()
## Sankey Plot

getSankey<-function(df_data){
  
  require(ggplot2)
  require(ggalluvial)
  require(dplyr)
  require(tidyr)
  
  colores <- colorRampPalette(brewer.pal(12, "Set3"))(length(unique(df_data[,ncol(df_data)])))
  
  df_long <- df_data %>%
    pivot_longer(cols = -ID, names_to = "K", values_to = "Cluster") 
  
  p1<-ggplot(df_long, aes(x = K, stratum = Cluster, alluvium = ID, fill = as.factor(Cluster))) +
    geom_flow(color = "black", size = 0.2) +  
    geom_stratum(alpha = 0.5, color = "black") +  
    geom_text(stat = "stratum", aes(label = Cluster), size = 3, color = "black") +  # ⬅️ Añadir etiquetas
    scale_fill_manual(values = colores) +  
    theme_minimal() +
    labs(title = "Sankey Plot de Clusters", x = "K", fill = "Cluster") +
    theme(legend.position = "none")
  
  return(p1)
}

# ann.db[ann.db$TraitID %in% cp[cp$K25 == 24,"ID"],"term"]
# ann.db[ann.db$TraitID %in% cp[cp$K32 == 21,"ID"],"term"]


##··············································································

##································································· ClusterDot()
## Function to get dotplot of cluster annotations
#' @param paths Dataframe with columns: ClusterID, K and Function
ClusterDot<-function(paths,useTM=FALSE){
  
  require("pheatmap")
  require("ggplot2")
  require("dplyr")
  
  ## Text Mining (optional)
  if(useTM){
    print("Text Mining Not Implemented Yet")
  }else{
    paths$Terms<-paths$Function
  }
  
  ## Get Frecuencies of each term in M with respect to the total
  freqTable<-data.frame("ID"=names(table(paths$Terms)),
                        "n"=as.numeric(table(paths$Terms)))
  
  dat<-do.call("cbind",lapply(unique(paths$K),function(cl){
    tmp<-paths[paths$K==cl,]
    res<-unlist(lapply(1:nrow(freqTable),function(p){
      score<-(nrow(tmp[tmp$Function==freqTable[p,"ID"],]) / freqTable[p,"n"]) *100
      return(score)
    }))
    names(res)<-freqTable$ID
    res<-round(res,digits = 2)
    return(res)
  }))
  colnames(dat)<-paste0("M",unique(paths$K))
  
  ## Get Composition
  freqTable2<-freqTable[freqTable$ID!="TBD" & freqTable$ID!="TBA",]
  
  dat2<-do.call("cbind",lapply(unique(paths$K),function(cl){
    
    tmp<-paths[paths$K==cl,]
    
    res<-unlist(lapply(1:nrow(freqTable2),function(p){
      
      score<-(nrow(tmp[tmp$Function==freqTable2[p,"ID"],]) / nrow(tmp)) *100
      return(score)
    }))
    names(res)<-freqTable2$ID
    res<-round(res,digits = 2)
    return(res)
  }))
  colnames(dat2)<-paste0("M",unique(paths$K))
  
  
  ## Plot conjunto
  paths<-intersect(rownames(dat),rownames(dat2))
  df<-data.frame(reshape2::melt(dat[paths,]),
                 reshape2::melt(dat2[paths,]))
  df<-df[,c(1,2,3,6)]
  colnames(df)<-c("Term","Module","Perc","Comp")
  
  ## Get order
  p.Ord<-pheatmap(dat[rownames(dat)!="TBD" & rownames(dat)!="TBA",],scale="none", show_colnames = T, cluster_cols = F,
                  cluster_rows = T, show_rownames = T, border_color = "#404040",
                  fontsize = 6,
                  breaks=seq(0,100,length.out = 100),
                  color = colorRampPalette(c("white","orange","#ce7457"))(100))
  
  p.Ord<-rownames(dat[rownames(dat)!="TBD" & rownames(dat)!="TBA",])[p.Ord$tree_row$order]
  
  df$Term<-factor(df$Term,levels=p.Ord)
  
  ## Filter data with Perc = 0
  df_filtered <- df %>%
    filter(Perc > 0)
  
  # Paso 2: Hacer el plot con el dataframe filtrado
  p1<- ggplot(df_filtered, aes(x = Module, y = Term, fill = Comp, size = Perc)) +
    geom_point(color = "black", shape = 21, stroke = 0.5) +
    scale_fill_gradient2(mid = "white", high = "#ce7457") + 
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,size = 7)) +
    labs(x = "", y = "", fill = "Composition (%)", size = "Frequency (%)")
  
  return(p1)
  
}

##··············································································

##······························································· explorePaths()
## Function to get dotplot of cluster annotations
#' @param paths Dataframe with columns: ClusterID and K
explorePaths<-function(paths,ann.db){
  
  require(tm)
  require(SnowballC)
  require(wordcloud)
  require(RColorBrewer)
  
  subcl<-unique(paths$K)
  
  res<-lapply(subcl,function(cl){
    
    termList<-paths[paths$K==cl,"ID"]
    coann<-ann.db[ann.db$TraitID %in% termList,"coannotation"]
    coann<-coann[!is.na(coann)]  
    
    if(length(coann)!=0){
      coann<-unlist(strsplit(coann,", "))
      termList<-unique(c(termList,coann))
    } 
    
    termList<-ann.db[ann.db$TraitID %in% termList,"term"]
    termList<-termList[termList!="TBA" & termList!="TBD"]
    
    corpus <- Corpus(VectorSource(termList))
    corpus <- tm_map(corpus, content_transformer(tolower))
    corpus <- tm_map(corpus, removePunctuation)
    corpus <- tm_map(corpus, removeNumbers)
    corpus <- tm_map(corpus, stripWhitespace)
    corpus <- tm_map(corpus, removeWords, stopwords("english"))
    corpus <- tm_map(corpus, stemDocument, language = "english")
    
    dtm <- as.matrix(DocumentTermMatrix(corpus))
    freqs <- colSums(dtm)
    freqs <- sort(freqs, decreasing = TRUE)
    
    freqs<-data.frame("term"=names(freqs),"n"=as.numeric(freqs))
    
    return(list("termList"=termList,"freqs"=freqs)) 
  })
  names(res)<-paste0("K",subcl)
  return(res)
}

##··············································································

##······························································· searchByTerm()
## Retrieve Genesets ID related with a term
searchByTerm<-function(term,ann.db){
  tmp<-ann.db[grepl(tolower(term),tolower(ann.db$term)),]
  
  ids<-unique(c(tmp$TraitID),
              unlist(lapply(strsplit(tmp$coannotation[!is.na(tmp$coannotation)],","),trimws)))
  
  tmp.ann<-ann.db[ann.db$TraitID %in% ids,]
  return(tmp.ann$TraitID)
}

##··············································································
