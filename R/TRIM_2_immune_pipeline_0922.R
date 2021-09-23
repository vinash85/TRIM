
###############
#immune pipeline
###############
# part1 : bulk level
# the input is tf expression and 32 known immune features
# the output is predicted immune-related tf
# the pre-calculate results are stored here: 
# liulab/xmwang/oxphos_proj/loading_data/TRIM_model/lmer.tf.dim1_opt.csv
# the script to get the result in TRIM_2_immune_pipeline.R
###############
# part2 : single-cell level
# input tf-targeted gene sets
# we defined the average expression of tf targets as tf activity
# output: the differences of tf activity between cancer cells and cd8

##################
# example to run single-cell level of immune pipeline
######
# 1. input tf target (the tf target list was extracted from cistrom db)
tf.target.list= readRDS("/liulab/xmwang/oxphos_proj/loading_data/cistrome/tf.target.list.rds")

# 2.select single-cell object
sc.file<- list.files("/tcga/tcga_ATAC_Seq/tcga_ATAC_Seq_processed_data/singlecell-data/scRNA/Data_cancer/",
                     recursive = T,full.names = T,pattern = "_res.rds"
)
sel_filename = c("HNSCC_GSE103322_cancer_res.rds","HNSCC_GSE103322_normal_res.rds",
                 "SKCM_GSE72056_cancer_res.rds","SKCM_GSE72056_normal_res.rds",
                 "NSCLC_EMTAB6149_res.rds","NSCLC_GSE117570_res.rds","NSCLC_GSE127465_res.rds" )
loc=which(basename(sc.file)%in% sel_filename )
sel_file=sc.file[loc]

# 3. calculate tf activity 
activity.list=list()

for (file in sel_file) {
  print(basename(file))
  st.obj<- readRDS(file)
  tf.activity<- CalTFactivity(st.obj = st.obj,
                            tf.target.list = tf.target.list,
                            cell_col = "assign.ident",
                            patient_col="orig.ident",
                            cell.id1 = "CD8Te",cell.id2 = "Malignant"
                            )
  activity.list[[basename(file)]]=tf.activity
}

# 4. Identify these tfs that differently activated between cancer cells and cd8
tf.diff <-IdentifyTF(activity.list,cell.id1 ="CD8")
# pre-calculated result: "/liulab/xmwang/oxphos_proj/loading_data/scRNAseq/score.diff.mul.rna.rds"


###############
# part3 : integrate the bulk-level analysis and single-cell analysis

comb.imm.list = list()
score.diff.mul=readRDS("/liulab/xmwang/oxphos_proj/loading_data/scRNAseq/score.diff.mul.rna.rds")
lmer.tf.dim1 = fread("/liulab/xmwang/oxphos_proj/loading_data/TRIM_model/lmer.tf.dim1_opt.csv",data.table = F)
merged.immune = merge(score.diff.mul,lmer.tf.dim1,by.x="TF",by.y="object")
for (i in 1:nrow(merged.immune)) {
  gene = merged.immune[i,]$TF
  coef.1 = as.numeric(merged.immune[i,2:6])
  names(coef.1)=c("Estimate","Std. Error","df","t value","Pr(>|t|)")
  coef.2 = as.numeric(merged.immune[i,7:11])
  names(coef.2)=c("Estimate","Std. Error","df","t value","Pr(>|t|)")
  tmp <- combine.lme4(coef1 = coef.1,coef2 = coef.2,n1 = 10767,n2 = 29888)
  comb.imm.list[[gene]]=tmp
}
comb.imm = do.call(rbind,comb.imm.list)
saveRDS(comb.imm, file = "./TRIM_model/comb.imm.rds")





#############
## functions in TRIM immune pipeline
#############
library(ggplot2)
library(tidyverse)
library(data.table)
library(ggpubr)
library(pheatmap)
library(ggthemes)
library(RColorBrewer)
library(cowplot)
library(lmerTest)
library(lme4)
library(Seurat)

# define the average expression of tf targets as tf activity
get_avg_expr= function(x,expr.data=st.expr_sel){
  score.tmp=expr.data[rownames(expr.data)%in% x,]
  score<- colMeans(as.matrix(score.tmp))
  return(score)
}


#integrate immune pipeline 
combine.lme4 = function(coef1, coef2, n1, n2) {
  # coef1 = numeric with names: Estimate, Std. Error, df,    t value,     Pr(>|t|)
  # coef2 = numeric with names: Estimate Std. Error df    t value     Pr(>|t|)
  # n1 = number of samples in dataset1 
  # n2 = number of samples in dataset1 
  n = n1 + n2
  var1 = n1 * (coef1["Std. Error"])^2
  var2 = n2 * (coef2["Std. Error"])^2
  Estimate = (coef1["Estimate"] + coef2["Estimate"])/2
  Var = (var1 + var2)/4
  Df = coef1["df"] + coef2["df"]
  SE = sqrt(Var/n)
  tval = Estimate/SE
  p.val = 2 * pt(abs(tval), Df, lower.tail = FALSE) 
  out = c(Estimate, SE, Df, tval, p.val)
  names(out)=c("Estimate", "SE", "Df", "tval", "p.val")
  return(out)
}



CalTFactivity<- function(st.obj,cell_col = "assign.ident",
                         patient_col="orig.ident",
                         cell.id1 = "CD8Te",cell.id2 = "Malignant",
                         tf.target.list = list()){
  
  st.expr = st.obj$RNA@assays$RNA@data
  meta.info<- st.obj$RNA@meta.data
  
  loc1 = which(grepl(meta.info[,cell_col],pattern = cell.id1))
  loc2 = which(grepl(meta.info[,cell_col],pattern = cell.id2))
  loc=c(loc1,loc2)
  st.expr_sel=st.expr[,loc]
  patient.info=meta.info[loc,patient_col]
  celltype.info = meta.info[loc,cell_col]
  common.gene<- lapply(tf.target.list, intersect, rownames(st.expr_sel)) 
  lapply(common.gene, length)->num.tmp
  tf.target.tmp=tf.target.list[which(num.tmp<30)]
  print("TF activity calculation for",names(tf.target.tmp),
        " may not accurate due to detected genes <30 in scRNAseq object!")
  
  sapply(tf.target.list, get_avg_expr,expr.data=st.expr_sel)->score.tmp
  score.tmp=as.data.frame(score.tmp)
  # integrete meta info
  score.tmp$patient = patient.info
  score.tmp$celltype= celltype.info
  return(score.tmp)
}

IdentifyTF<- function(activity.list,cell.id1 = "CD8"){
  library(lme4)
  score.rna<- do.call(rbind,activity.list)
  score.rna=as.data.frame(score.rna)
  score.rna <- score.rna %>%
    mutate(celltype=if_else(grepl(cell.id1,celltype),"0","1"),
           celltype= as.factor(celltype))
  
  score.diff.list=list()
  n=1
  ncol.num = ncol(score.rna)
  for (tf in colnames(score.rna)[-c((ncol.num-1):ncol.num)]) {
    print(paste0(n," : ",tf))
    lmer(score.rna[,tf]~celltype+(1|patient),data = score.rna)->mod.tmp
    tmp<-as.data.frame(coef(summary(mod.tmp) ))
    score.diff.list[[tf]]=tmp[2,]
    n=n+1
  }
  score.diff = do.call(rbind,score.diff.list)
  score.diff$TF = rownames(score.diff)
  return(score.diff)
}

