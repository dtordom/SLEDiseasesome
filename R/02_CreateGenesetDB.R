#' ·············································································
#' Building Diseasome 
#' R version 4.4.0 (2024-04-24 ucrt)
#' Dec 2024
#' ·············································································
#' Dissect and merge pathway databases, and reduce overlapping

## ·············································································
## Set environment ---- 


set.seed(12345)
# memory.limit(size = 1600000000)
setwd("C:/Users/danie/Desktop/WORK/DISEASOME")

source("C:/Users/danie/Desktop/WORK/DISEASOME/Code/utils.R")

check.packages(c("parallel","matrixStats","biomaRt","NOISeq","stringr","dplyr",
                 "tidyr","doParallel","caret","pbapply","BiocParallel","tibble",
                 "pathMED","NbClust","ConsensusClusterPlus","SNFtool","BloodGen3Module","purrr",
                 "igraph","pheatmap","SNFtool","UpSetR","ComplexHeatmap"))


## ·············································································
## Collect all pathway databases ---- 

## Load transcriptome data
DATA<-readRDS("C:/Users/danie/Desktop/WORK/DISEASOME/RData/Datasets.rds")

dbs<-pathMED:::genesetsData[c("tmod","go_bp","go_cc","go_mf","reactome","kegg","wikipathways")]
# tmod includes LI and DC databases

## Add literature-based Genesets
# Add aditional IFN signatures from bibliography
ifnSig<-read.csv("C:/Users/danie/Desktop/WORK/DISEASOME/Data/IFN_Signature.csv",
                 sep=";")[,1:2]
ifnSig <- as.data.frame(ifnSig %>% separate_rows(Gene, sep = " /// "))

knowlge<-list("IFN_I" = as.character(ifnSig[ifnSig$Signature=="IFI_I","Gene"]),
               "IFN_II" = as.character(ifnSig[ifnSig$Signature=="IFI_II","Gene"]),
               "IFN_II_III" = as.character(ifnSig[ifnSig$Signature=="IFI_II_III","Gene"]),
               "IFN_I_II_III" = as.character(ifnSig[ifnSig$Signature=="IFI_I_II_III","Gene"]),
               "T_cell_exhaustion" = c("CTLA4", "IL7R", "LAG3", "PDCD1", "ABCE1"))
 
dbs[["knowlge"]]<-knowlge

## BloodGen3Module
annb3m<-as.data.frame(BloodGen3Module:::Module_listGen3)

b3m<-lapply(unique(annb3m$Module),function(term){
  tmp<-annb3m[annb3m$Module==term,]
  return(as.character(tmp$Gene))
})
names(b3m)<-unique(annb3m$Module)
dbs[["B3M"]]<-b3m

## Single-cell db
xcelldb<-readRDS("C:/Users/danie/Desktop/WORK/DISEASOME/RData/xcellDB.rds")

nonBlood<-c("Astrocytes","Hepatocytes","Adipocytes","aDC","Chondrocytes",
            "Endothelial","Epithelial","Keratinocytes","Melanocytes",
            "Mesangial","MSC","Neurons","Osteoblast","Preadipocytes","muscle")

# Filter non-blood cells
xcelldb<-xcelldb[!sapply(names(xcelldb), function(element, pattern) {
  sapply(pattern, function(patr) grepl(patr, element)) %>% any()
}, nonBlood)]


dbs[["xcell"]]<-xcelldb

## 11 Sources of pathways-gene annotations



## ······································································ Step 1 
## DissectDB: Split and Merge all databases ---- 
#' DATA: nested list with datasets (list with Disease and Healthy matrices)
#' dbs.all: list of pathways and their genes
#' ann.info: dataframe with annotation_id, term and source for each pathway

# Annotation reference object
ann.info<-do.call("rbind",lapply(1:length(dbs),function(i){
  if(names(dbs)[i]=="xcell" | names(dbs)[i]=="knowlge"){
    paths<-names(dbs[[i]])
    tmp<-data.frame("annotation_id"=paths,"term"=paths,
                    "source"= rep(names(dbs)[i],length(paths)))
    return(tmp)
  }
  if(names(dbs)[i]=="B3M"){
    paths<-names(dbs[[i]])
    tmp<-annb3m[annb3m$Module %in% paths,c("Module","Function")]
    tmp<-tmp[!duplicated(tmp),]
    colnames(tmp)<-c("annotation_id","term")
    tmp$source<-rep(names(dbs)[i],nrow(tmp))
    return(tmp)
  }

  if(!names(dbs)[i] %in% c("xcell","B3M","knowlge")){
    paths<-names(dbs[[i]])
    tmp<-pathMED:::ann_info
    tmp<-tmp[tmp$annotation_id %in% paths,]
    tmp$source<-rep(names(dbs)[i],nrow(tmp))
    return(tmp)
  }
}))

## Joint all databases in one
dbs.all<-flatten(dbs)


## DissectDB: split pathways into co-expressed sub-paths
#' Use big datasets (remove datasets with low shared genes across cohorts)
custom.db<-dissectDB(refData = DATA[-c(4,8,12,14)], 
                     geneSets = dbs.all,
                     minPathSize = 8, minSplitSize = 3,maxSplits = NULL,
                     explainedVariance = 70, percSharedGenes = 90)

print(length(custom.db)) ## 47866, 47337

rm(nonBlood,xcelldb,knowlge,ifnSig,b3m,annb3m,dbs)
gc()


save.image("C:/Users/danie/Desktop/WORK/DISEASOME/RData/DB_step1.RData")


## ······································································ Step 2 
## Get scores for all pathways and datasets ---- 

## Get Mscores for all the SLE Datasets
SCORES<-mScores_createReference(refData=DATA,
                                geneSets=custom.db,
                                cores = 12)


## Create a full table to store co-annotation information
ann.db<-as.data.frame(do.call("rbind",lapply(1:length(custom.db),function(p){
  
  nme<-gsub("\\.split.*", "", names(custom.db)[p])
  tmp<-ann.info[ann.info$annotation_id %in% nme,]
  tmp<-c(names(custom.db)[p],as.character(tmp))
  return(tmp)
})))
colnames(ann.db)<-c("TraitID","AnnotationID","term","source")
ann.db$coannotation<-NA


save.image("C:/Users/danie/Desktop/WORK/DISEASOME/RData/DB_step2.RData")


## Check for Outliers Outliers
print(unlist(lapply(SCORES$mscores,max)))


## ······································································ Step 3
## Reducing pathway overlapping using Set theory (Set Packing) ---- 

# Run: 'Set Environment' section
# load("C:/Users/danie/Desktop/WORK/DISEASOME/RData/DB_step2.RData")

## Remove pathways with <= 2 genes
custom.db <- custom.db[!names(custom.db) %in% 
                         names(custom.db[lengths(custom.db) <=2])]
gc()
## 38797


## Run Set packing
time0<-Sys.time()
res.sPack<-setPackingFilter(pathDB=custom.db, mRef=SCORES, ann.db=ann.db,
                           coverage.thr=0.8, max.combs =8, 
                           gainIndex=10, minCorr=0.75)
time1<-Sys.time()
print(as.numeric(difftime(time1,time0,units="min"))) # 243.8337 mins

custom.db<-res.sPack$db
ann.db<-res.sPack$ann

rm(res.sPack,time0,time1)
save.image("C:/Users/danie/Desktop/WORK/DISEASOME/RData/DB_step3.RData")

## 37537

## ······································································ Step 4 
## Reducing pathway overlapping (Jaccard index) ---- 

# Run: 'Set Environment' section
# load("C:/Users/danie/Desktop/WORK/DISEASOME/RData/DB_step3.RData")

## 1. Similarity 0.8 
time0<-Sys.time()
res.JI<-getNodes(listDB = custom.db, ann.db=ann.db,
                 simmilarity.threshold = 0.8, max.length = NULL)
time1<-Sys.time()
print(as.numeric(difftime(time1,time0,units="min"))) ##

custom.db<-res.JI$db
ann.db<-res.JI$ann


rm(res.JI,time0,time1)
save.image("C:/Users/danie/Desktop/WORK/DISEASOME/RData/DB_step4.RData")


## 2. Similarity for Small Genesets (4-5) 
time0<-Sys.time()
res.JI<-getNodes(listDB = custom.db, ann.db=ann.db,
                 simmilarity.threshold = 0.7, max.length = 4)
time1<-Sys.time()
print(as.numeric(difftime(time1,time0,units="min"))) ##

custom.db<-res.JI$db
ann.db<-res.JI$ann


## Remove duplicated values in coannotation
for(i in 1:nrow(ann.db)){
  if(!is.na(ann.db[i,"coannotation"])){
    x<-unique(unlist(strsplit(ann.db[i,"coannotation"],split = ", ")))
    ann.db[i,"coannotation"]<-paste(x, collapse = ", ")
  }
}

rm(list=setdiff(ls(),c("ann.db","custom.db","DATA","SCORES")))

## Save R Object
save.image("C:/Users/danie/Desktop/WORK/DISEASOME/RData/DB_step5.RData")

## 36853

## ······································································ Step 5 
## Tuning parameter for function mScores_filterPaths (pathMED)
#' Get stats based on parameter combination:
#' False Discovery Rate (fdr)
#' Loss of information (Variability across patients)
#' Database diversity (Shannon index)

# load("C:/Users/danie/Desktop/WORK/DISEASOME/RData/DB_step5.RData")

## Get Mscores for all the SLE Datasets
SCORES<-mScores_createReference(refData=DATA,
                                geneSets=custom.db,
                                cores = 16)
gc()

save.image("C:/Users/danie/Desktop/WORK/DISEASOME/RData/DB_step6.RData")


## Step 5.1: Get false discovery rate ················
# load("C:/Users/danie/Desktop/WORK/DISEASOME/RData/DB_step6.RData")

## Create database of random pathways and get scores
all_genes <- unique(unlist(custom.db))
sizes <- c(3,quantile(sapply(custom.db, length),probs = 0.9))
random.db <- lapply(1:1000, function(x) {
  sample(all_genes, sample(sizes[1]:sizes[2], 1))})
names(random.db)<-paste("rp.",1:length(random.db))

SCORES.rd<-mScores_createReference(refData=DATA,
                                   geneSets = random.db,
                                   cores = 16)


## Get fdr for each combination based on mscore significances
#' fdr: number of random paths that achieves significance / total nº of paths
params <- expand.grid(percentage = seq(from = 1, to = 30, by = 1), datasets = 1:15)
params$fdr<-unlist(lapply(1:nrow(params),function(p){
  drp.rp<-mScores_filterPaths(MRef=SCORES.rd,
                       min_datasets=params[p,"datasets"],
                       perc_samples=params[p,"percentage"],
                       Pcutoff=0.05,plotMetrics = FALSE)  
  res<-length(drp.rp)
  return(as.numeric((res/1000)*100))
}))
#params$datasets <- factor(params$datasets, levels = c(1:15))
params$datasets<-as.numeric(params$datasets)

# Heatmap plot
p3<-ggplot(params, aes(x = percentage, y = datasets, fill = fdr)) +
  geom_tile(color="black") +
  labs(title = "Random-based false positives", x = "Percentage of patients",
       y = "Number of datasets",fill = "FDR") +
  scale_fill_gradientn(
    colors = c("lightblue", "gold", "darkorange", "darkorange3", "darkred"),
    values = scales::rescale(c(0, 5, 10, 50, 100)),
    limits = c(0, 100)) +
      theme(legend.text = element_text(size=8),
            axis.text.x=element_text(size=7, color = "black"),
            axis.text.y = element_text(size=7, color = "black"),
            axis.title = element_text(size=8, color = "black", face = "bold"), 
            strip.text = element_text(size=7), 
            strip.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(size=10,hjust=0.5,face="bold"))+
  scale_x_continuous(
    breaks = c(1,5,10,15,20,25,30),labels = paste0(c(1,5,10,15,20,25,30), " %")) +
  scale_y_continuous(breaks = seq(1, 15, by = 1),labels = seq(1, 15, by = 1))
plot(p3)

rm(all_genes,sizes)

save.image("C:/Users/danie/Desktop/WORK/DISEASOME/RData/DB_step7.RData")


## Step 5.2: Calculate the loss of information················

#' Get nº PCAs that explained the 80% of the varianza in patients using all paths
explainedVar<-unlist(lapply(1:length(SCORES$mscore),function(dat){
  pca <- FactoMineR::PCA(t(SCORES$mscores[[dat]]), graph = F)
  pca_eig <- as.data.frame(pca$eig)
  pca_eig <-pca_eig[pca_eig$`cumulative percentage of variance` < 80, ]
  return(nrow(pca_eig) + 1) }))
names(explainedVar)<-names(SCORES$mscore)


pcaVars<-do.call("rbind",lapply(1:nrow(params),function(i){
#' Get nº PCAs that explained the 80% of the varianza in patients using selected
#' paths based on parameter selection

  drp<-mScores_filterPaths(MRef=SCORES,
                           min_datasets=as.numeric(params[i,"datasets"]),
                           perc_samples=params[i,"percentage"],Pcutoff=0.05,
                           plotMetrics = FALSE)
  
  pca.dat<-unlist(lapply(1:length(SCORES$mscore),function(dat){
    pca <- FactoMineR::PCA(t(SCORES$mscores[[dat]][names(drp),]), graph = F)
    pca_eig <- as.data.frame(pca$eig)
    pca_eig <-pca_eig[pca_eig$`cumulative percentage of variance` < 80, ]
    return(nrow(pca_eig) + 1)
  }))
  names(pca.dat)<-names(SCORES$mscores)
  
  res<- (pca.dat / explainedVar ) * 100
  return(res)
}))
names(pcaVars)<-names(SCORES$mscores)


params$lossInfo<-as.numeric(unlist(apply(pcaVars,1,mean)))


p4<-ggplot(params, aes(x = percentage, y = datasets, fill = lossInfo)) +
  geom_tile(color="black") +
  labs(title = "Data Variance Through PCA Reduction", x = "Percentage of patients",
       y = "Number of datasets",fill = "%PCAs") +
  scale_fill_gradientn(
    colors = c("black", "white"),
    values = scales::rescale(c(0, 100)),
    limits = c(0, 100)) +
  theme(legend.text = element_text(size=8),
        axis.text.x=element_text(size=7, color = "black"),
        axis.text.y = element_text(size=7, color = "black"),
        axis.title = element_text(size=8, color = "black", face = "bold"), 
        strip.text = element_text(size=7), 
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=10,hjust=0.5,face="bold"))+
  scale_x_continuous(
    breaks = c(1,5,10,15,20,25,30),labels = paste0(c(1,5,10,15,20,25,30), " %")) +
  scale_y_continuous(breaks = seq(1, 15, by = 1),labels = seq(1, 15, by = 1))
plot(p4)


## Save results
save.image("C:/Users/danie/Desktop/WORK/DISEASOME/RData/DB_step8.RData")


## Step 5.3: Database diversity (Shannon index) ················
library(tidyverse)
library(tidytext)
library(SnowballC)
library(vegan)
library("pathMED")
library(ggpubr)
library(UpSetR)
library(ComplexHeatmap)

## Sparse all database terms
all_terms<-ann.db[ann.db$TraitID %in% names(custom.db),"term"]

all_terms <- all_terms %>% as_tibble() %>% rename(term = value) %>%
  mutate(term = str_to_lower(term)) %>% unnest_tokens(word, term) %>%
  anti_join(stop_words)  

## Words to roots
all_terms <- all_terms %>% mutate(stem = wordStem(word))  


div.res<-do.call("rbind",lapply(1:nrow(params),function(i){
  
  #print(i)
  drp<-mScores_filterPaths(MRef=SCORES,
                           min_datasets=as.numeric(params[i,"datasets"]),
                           perc_samples=params[i,"percentage"],Pcutoff=0.05,
                           plotMetrics = FALSE)
  drp<-names(drp)
  
  func_terms<-ann.db[ann.db$TraitID %in% drp,"term"]
  
  clean_terms <- func_terms %>% as_tibble() %>% rename(term = value) %>%
    mutate(term = str_to_lower(term)) %>% unnest_tokens(word, term) %>%
    anti_join(stop_words)  
  
  ## Words to roots
  clean_terms <- clean_terms %>% mutate(stem = wordStem(word))  
  
  ## Here I can remove rows with specific terms (optional)
  ##
  
  # Shannon
  freq_table <- clean_terms %>% count(stem, sort = TRUE)
  # freq_table
  term_counts <- freq_table$n
  shannon_index <- diversity(term_counts, index = "shannon")

  # ZCI value  
  ZCI_value <- calculate_ZCI(all_terms$stem, clean_terms$stem)
  
  res<-c("nPaths"=length(drp),"ZCI"=ZCI_value,"ShannonIndex"=shannon_index)
  return(res)
}))

params<-data.frame(params,div.res)

## Save results
save.image("C:/Users/danie/Desktop/WORK/DISEASOME/RData/DB_step9.RData")


p5<-ggplot(params, aes(x = percentage, y = datasets, fill = ZCI)) +
  geom_tile(color="black") +
  labs(title = "ZCI", x = "Percentage of patients",
       y = "Number of datasets",fill = "ZCI") +
  scale_fill_gradientn(
    colors = c("lightblue", "gold", "darkorange", "darkorange3", "darkred"),
    values = scales::rescale(c(0.5, 1.2)),
    limits = c(0.5, 1.2)) +
  theme(legend.text = element_text(size=8),
        axis.text.x=element_text(size=7, color = "black"),
        axis.text.y = element_text(size=7, color = "black"),
        axis.title = element_text(size=8, color = "black", face = "bold"), 
        strip.text = element_text(size=7), 
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=10,hjust=0.5,face="bold"))
plot(p5)

p5<-ggplot(params, aes(x = percentage, y = datasets, fill = ShannonIndex)) +
  geom_tile(color="black") +
  labs(title = "Shannon index", x = "Percentage of patients",
       y = "Number of datasets",fill = "Shannon index") +
  scale_fill_gradientn(
    colors = c("lightblue", "gold", "darkorange", "darkorange3", "darkred"),
    values = scales::rescale(c(3.5, 6.5)),
    limits = c(3.5, 6.5)) +
  theme(legend.text = element_text(size=8),
        axis.text.x=element_text(size=7, color = "black"),
        axis.text.y = element_text(size=7, color = "black"),
        axis.title = element_text(size=8, color = "black", face = "bold"), 
        strip.text = element_text(size=7), 
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=10,hjust=0.5,face="bold")) +
  scale_x_continuous(
    breaks = c(1,5,10,15,20,25,30),labels = paste0(c(1,5,10,15,20,25,30), " %")) +
  scale_y_continuous(breaks = seq(1, 15, by = 1),labels = seq(1, 15, by = 1))
plot(p5)



## Step 5.4: Selection of best parameters ················


## LossInfo inlfexion point
loss_matrix <- matrix(params$lossInfo, nrow = 30, ncol = 15)

# Calcular cambio relativo en la dirección de N y P
grad_N <- abs(diff(loss_matrix, lag = 1, differences = 1))  # Diferencia en N
grad_P <- abs(diff(t(loss_matrix), lag = 1, differences = 1))  # Diferencia en P


# Convertimos en un umbral relativo (ej. 10% de caída)
thresholds <- quantile(params$lossInfo, probs = c(0.9, 0.8, 0.7,0.5))  # 10%, 15%, 20%

# Asignar colores según el threshold
df <- params %>%
  mutate(critical_level = case_when(
    lossInfo <= thresholds[4] ~ "50% threshold",
    lossInfo <= thresholds[3] ~ "30% threshold",
    lossInfo <= thresholds[2] ~ "20% threshold",
    lossInfo <= thresholds[1] ~ "10% threshold",
    TRUE ~ "No threshold"
  ))

df$datasets<-as.numeric(df$datasets)

#df$nowIndex<-1:nrow(params)
df_fdr<-do.call("rbind",lapply(unique(df$percentage),function(i){
  tmp<-df[df$percentage==i,]
  tmp<-tmp[tmp$fdr<=5,]
  return(tmp[1,])
}))
#df_fdr$datasets<-df_fdr$datasets+1


# Seleccionar valores en la curva crítica (donde LossInfo baja más del umbral)
#critical_curve <- params[params$lossInfo <= threshold, ]

# Graficar la curva crítica sobre un heatmap
p6<-ggplot(df, aes(x = percentage, y = datasets, fill = ShannonIndex)) +
  geom_tile(color="black") +
  
  scale_fill_gradientn(
    colors = c("darkred","darkorange3","darkorange","gold","#86de9b"),
    values = scales::rescale(c(3.5, 6.5)),
    limits = c(3.5, 6.5)) +
  theme(legend.text = element_text(size=8),
        axis.text.x=element_text(size=7, color = "black"),
        axis.text.y = element_text(size=7, color = "black"),
        axis.title = element_text(size=8, color = "black", face = "bold"), 
        strip.text = element_text(size=7), 
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=10,hjust=0.5,face="bold"))+
  geom_point(data = df %>% filter(critical_level != "No threshold"), 
             aes(color = critical_level), shape = 4, size = 2) +
  scale_color_manual(values = c("10% threshold" = "#86cdeb",
                                "20% threshold" = "#619bc2",
                                "30% threshold" = "#296289",
                                "50% threshold" = "black")) +
  geom_line(data = df_fdr,aes(x = percentage, y = datasets), 
            color = "black", size = 0.3, alpha = 1) +
  scale_x_continuous(
    breaks = c(1,5,10,15,20,25,30),labels = paste0(c(1,5,10,15,20,25,30), " %")) +
  scale_y_continuous(breaks = seq(1, 15, by = 1),labels = seq(1, 15, by = 1)) +
  labs(title = "Parameter optimization", x = "Percentage of patients",
       y = "Number of datasets",fill = "Shannon Index") 

plot(p6)


## Number of paths

df<-params[params$fdr<5 & params$lossInfo>50,]

p7<-ggplot(df, aes(x = percentage, y = datasets, fill = nPaths)) +
  geom_tile(color="black") +
  
  scale_fill_gradientn(
    colors = c("grey","gold","darkorange","darkorange3","darkred"),
    values = scales::rescale(c(min(df$nPaths), max(df$nPaths))),
    limits = c(min(df$nPaths), max(df$nPaths))) +
  theme(legend.text = element_text(size=8),
        axis.text.x=element_text(size=7, color = "black"),
        axis.text.y = element_text(size=7, color = "black"),
        axis.title = element_text(size=8, color = "black", face = "bold"), 
        strip.text = element_text(size=7), 
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=10,hjust=0.5,face="bold"))+
  scale_x_continuous(
    breaks = c(1,5,10,15,20,25,30),labels = paste0(c(1,5,10,15,20,25,30), " %")) +
  scale_y_continuous(breaks = seq(1, 15, by = 1),labels = seq(1, 15, by = 1)) +
  labs(title = "Number of geneset selected", x = "Percentage of patients",
       y = "Number of datasets",fill = "DRPs") +
  geom_vline(xintercept = 10, color = "black", linetype = "dashed", size = 0.3) + 
  geom_hline(yintercept = 5, color = "black", linetype = "dashed", size = 0.3)

plot(p7)


ggarrange(p3,p4,p6,p7,nrow=2,ncol=2)



## Step 5.5: Filtering candidates ················

params.sel<-params[params$fdr<5 & params$lossInfo>50 & 
                     params$percentage>10 & params$datasets > 5,]

params.sel$names<-paste0("d",params.sel$datasets,"_",params.sel$percentage)

DRP<-lapply(1:nrow(params.sel),function(i){
  drp<-mScores_filterPaths(MRef = SCORES,
                           min_datasets=params.sel[i,"datasets"],
                           perc_samples=params.sel[i,"percentage"],
                           Pcutoff=0.05,plotMetrics = FALSE)
  # print(length(drp))
  return(names(drp))
})
names(DRP)<-params.sel$names


m1 = make_comb_mat(DRP,mode = "distinc")
UpSet(m1, set_order = params.sel$names,
            comb_order = order(comb_size(m1),decreasing = T))


## Lost in functional terms ...........
# 
# 
# termFreqs<-function(x){
#   
#   func_terms<-tolower(ann.db[ann.db$TraitID %in% x,"term"])
#   func_terms<-gsub("\\b([tb])\\s*cell(s?)\\b", "\\1cell", func_terms, ignore.case = TRUE)
#   
#   #func_terms<-func_terms[!func_terms %in% c("tba","tbd","undetermined")]
#   
#   func_terms <- func_terms %>% as_tibble() %>% rename(term = value) %>%
#     mutate(term = str_to_lower(term)) %>% unnest_tokens(word, term) %>%
#     anti_join(stop_words)  
#   
#   ## Words to roots
#   func_terms <- func_terms %>% mutate(stem = wordStem(word))  
#   
#   func_terms<-func_terms[!func_terms$word %in% c("process","mediated","1","2","3","4",
#                                                  "5","6","7","8","9","type","dependent",
#                                                  "regulation","signaling","binding","system",
#                                                  "activity","pathway","pathways","asembly",
#                                                  "positive","negative","response"),]
#   
#   
#   freq_table <- as.data.frame(func_terms %>% count(word, sort = TRUE))
#   freq_table[1:25,]
#   return(freq_table)
# }
# 
# all<-termFreqs(x=DRP$d6_11)
# all$n<-all$n / sum(all$n) * 100
# bgap<-termFreqs(x=setdiff(DRP$d6_11,DRP$d7_12))
# bgap$n<-bgap$n / sum(bgap$n) * 100
# 
# rownames(all)<-all$word
# rownames(bgap)<-bgap$word
# 
# m<-data.frame(all[intersect(bgap$word,all$word),"n"],
#               bgap[intersect(bgap$word,all$word),"n"])
# rownames(m)<-intersect(bgap$word,all$word)
# colnames(m)<-c("all","diff")
# 
# plot(m[,1],m[,2])


rm(DRP,m1,params.sel,drp.rp,random.db,tmp,explainedVar,i,m)

save.image("C:/Users/danie/Desktop/WORK/DISEASOME/RData/DB_step10.RData")


## ······································································ Step 6
## Get selected Genesets ---- 
#' Criteria: FDR<5%, LossInfo < 50 %, datasets > 1/3 total
#' %patietns > 10%, and lowest combo with same Shannon Index
#' SELECTION: min_datasets=7, perc_samples=12
## 7 Datasets, 12 % of patients

DRGs.DB<-mScores_filterPaths(MRef=SCORES,
                         min_datasets=7, perc_samples=12,
                         Pcutoff=0.05,plotMetrics = FALSE)


# Selected Genesets: 5212 genesets

full.DB<-custom.db

## Save R Object with only essential elements
rm(list=setdiff(ls(),c("full.DB","ann.db","DRGs.DB","DATA","SCORES")))
gc()

save.image("C:/Users/danie/Desktop/WORK/DISEASOME/RData/DB_step11.RData")



