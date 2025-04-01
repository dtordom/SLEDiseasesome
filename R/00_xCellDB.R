#' ·············································································
#' Building SLE-Diseasome 
#' R version 4.4.0 (2024-04-24 ucrt)
#' March 2025
#' ·············································································
#' Get single cell Genesets FROM xCellDB

library("readxl")
library("stringr")

## Download xCell database from https://pmc.ncbi.nlm.nih.gov/articles/PMC5688663/
data <- read_excel("C:/Users/danie/Desktop/WORK/DISEASOME/Data/xCellDB.xlsx")
data<-as.data.frame(data)
data$Celltype_Source_ID<-gsub("\\+","pos",data$Celltype_Source_ID)

data[,1:2]

uniqPaths<-data$Celltype_Source_ID
uniqPaths<-unique(stringi::stri_replace_all_regex(uniqPaths,
                                pattern = c("_1", "_2", "_3","_HPCA","_ENCODE","_FANTOM","_IRIS","_BLUEPRINT","_NOVERSHTERN"), 
                                replacement = c("","","","","","","","",""),
                                vectorize = FALSE))

db<-list()
for(i in 1:length(uniqPaths)){
  tmp<-data[str_detect(data$Celltype_Source_ID,pattern = uniqPaths[i]),]
  tmp<-tmp[,-c(1:2)]
  tmp<-unlist(tmp)
  tmp<-tmp[!is.na(tmp)]
  tmp<-unique(tmp)
  db[[uniqPaths[i]]]<-tmp
}
names(db)<-gsub(" ","_",names(db))
names(db)<-gsub("-","_",names(db))

## Add additional signatures from bibliography curated by experts
db<-append(db,list("B_cells_DN4"=c("HOPX","PDE4D","IGHE","SELL"),
                   "B_cells_DN3"= c("RHOB","VIM","CTSH"),
                   "B_cells_DN2"= c("EMP3","CIB1","PSAP","CD72","DAPP1","HCK","ZEB2","RHOB","TNFRSF1B","FCRL3","FCRL5","FGR","MPP6"),
                   "B_cells_DN1" = c("TAGLN2","IGHA2","JCHAIN","IGHA1","S100A10"),
                   "B_cells_trans"=c("VPREB3","IGHD","IIGLL5","TCL1A"),
                   "APCs" = c("IFI30","TNFAIP2","CLEC7A"),
                   "NK_cells_NKT" = c("CD3E","CD3D","CD3G","HLA-DPB1","HLA-DRB1","HLA-DQB1"),
                   
                   "NK_cells_CD56bright" = c("GPR183","ILR7","LTB","GZMK","CD62L","CCR7","CD2","KLRC1"),
                   "NK_cells_CD56dim_CD16pos" = c("FCGR3A","KIR2DL1","KIR2DL3","KIR3DL1","KIR3DL2","KIR3DL3","KIR2DL4"),
                   "NK_cells_CD56neg_CD16pos_CD7pos" = c("CCL3","CCL4","CCL5")
                   ))

## Save Object
rm(i,tmp,data,uniqPaths)
saveRDS(db,"C:/Users/danie/Desktop/WORK/DISEASOME/RData/xcellDB.rds")

