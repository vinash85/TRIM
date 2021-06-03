##################
#TRIM  immune module
##################

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
library(factoextra)
library(remef)
##########################
# 1. bulk level
##########################
feature.immune<- readxl::read_xlsx("./data/1-s2.0-S1074761318301213-mmc2.xlsx")
tf.exprdata = readRDS("./data/tf.exprdata.rds")


oxphos.dataset.ck<- feature.immune[,c(1:2,5:32)] %>%
  dplyr::rename(sample =`TCGA Participant Barcode`,
                cancer = `TCGA Study`) %>%
  mutate_at(colnames(.)[-c(1,2)],as.numeric) %>%
  mutate(sample = gsub(sample,pattern = "-",replacement = ".")) %>%
  column_to_rownames(var = "sample") 
#impute features
library(missMDA)
rr.impute <- imputePCA(oxphos.dataset.ck[,-1], ncp=2)
dataset.ck.impute<- as.data.frame(rr.impute$completeObs)
dataset.ck.impute$cancer =oxphos.dataset.ck$cancer

#remove the cancertype effect from the actual value of each feature
feature.rev =list()
for (feature in colnames(dataset.ck.impute)[-ncol(dataset.ck.impute)]) {
  
  dataset.ck.impute[,c(feature,"cancer")] %>%
    lmer(as.formula(paste0("`",feature,"`","~","(1|cancer)")),data = .) ->tmp.mod
  y_partial <- remef(tmp.mod, ran = "all")
  feature.rev[[feature]]=y_partial
}

feature.rev.tcga = do.call(cbind,feature.rev)
rownames(feature.rev.tcga) = rownames(dataset.ck.impute)
res.pca <- prcomp(feature.rev.tcga,scale. = T)


tf.exprdata$object = substr(tf.exprdata$object,1,12)
calc_tf_feature_lmer(res.pca = res.pca,dim = "Dim.1",
                     cancer_info = oxphos.dataset.ck$cancer,
                     tf.exprdata = tf.exprdata  )->lmer.tf.dim1




#scale the value 
# 1. get the sd of each tf
tf.sd = apply(tf.exprdata[,-1], 2, sd)
# 2. the sd of dim1
sd_pca1<- sd(res.ind$coord[,1])  
# 3. b*sd_pc/sd_tf
all(lmer.tf.dim1$object==rownames(tf.sd))
lmer.tf.dim1$sd = tf.sd
lmer.tf.dim1$Estimate_nor = (lmer.tf.dim1$Estimate*sd_pca1)/lmer.tf.dim1$sd
write.csv(lmer.tf.dim1, file = "~/OXPHOS/Figures/lmer.tf.dim1_opt.csv")

lmer.tf.dim1 = fread("~/OXPHOS/Figures/lmer.tf.dim1_opt.csv",data.table = F)
lmer.tf.dim1 = column_to_rownames(lmer.tf.dim1,var = "V1")

#####################
# 2 single cell levle 
#####################
##single cell check the degene between the tumor cell and immune cells in multiple datasets

library(Seurat)

sc.file<- list.files("/tcga/tcga_ATAC_Seq/tcga_ATAC_Seq_processed_data/singlecell-data/scRNA/Data_cancer/",
                     recursive = T,full.names = T,pattern = "_res.rds"
                      )
tf.target.list = readRDS("/liulab/xmwang/oxphos_proj/loading_data/cistrome/tf.target.list.rds")
lapply(tf.target.list, length)->num.tmp
tf.target.bcc=tf.target.list[which(num.tmp>30)]
#assign.ident
get_avg_expr= function(x,expr.data=st.expr_sel){
  score.tmp=expr.data[rownames(expr.data)%in% x,]
  score<- colMeans(as.matrix(score.tmp))
  return(score)
}
sample.score.list=list()
# note 
sel_filename = c("HNSCC_GSE103322_cancer_res.rds","HNSCC_GSE103322_normal_res.rds",
                 "SKCM_GSE72056_cancer_res.rds","SKCM_GSE72056_normal_res.rds",
                "NSCLC_EMTAB6149_res.rds","NSCLC_GSE117570_res.rds","NSCLC_GSE127465_res.rds" )
loc=which(basename(sc.file)%in% sel_filename )
sel_file=sc.file[loc]

#file=sc.file[17]
for (file in sel_file) {
  print(basename(file))
  st.obj<- readRDS(file)
  st.expr = st.obj$RNA@assays$RNA@data
  meta.info<- st.obj$RNA@meta.data
  if (!grepl(pattern = "cancer",basename(file))) {
    loc1 = which(grepl(meta.info$assign.ident,pattern = "CD8Te"))
    loc2 = which(grepl(meta.info$assign.ident,pattern = "Malignant"))
    loc=c(loc1,loc2)
    st.expr_sel=st.expr[,loc]
    patient.info=meta.info[loc,]$orig.ident
    celltype.info = meta.info[loc,]$assign.ident
  }else{
    st.expr_sel=st.expr
    patient.info=meta.info$orig.ident
    celltype.info = "Malignant"
  }
  sapply(tf.target.bcc, get_avg_expr,expr.data=st.expr_sel)->score.tmp
  score.tmp=as.data.frame(score.tmp)
  # integrete meta info
  score.tmp$patient = patient.info
  score.tmp$celltype= celltype.info
  sample.score.list[[basename(file)]]=score.tmp
  saveRDS(sample.score.list, 
          file = "/liulab/xmwang/oxphos_proj/loading_data/TRIM_model/sample.score.list.rds")
}


sample.score.list=readRDS("/liulab/xmwang/oxphos_proj/loading_data/TRIM_model/sample.score.list.rds")

# calc the degene
score.rna<- do.call(rbind,sample.score.list)
score.rna=as.data.frame(score.rna)
score.rna <- score.rna %>%
  mutate(celltype=if_else(grepl("CD8",celltype),"0","1"),
         celltype= as.factor(celltype))
library(lme4)
score.diff.list=list()
n=1
for (tf in colnames(score.rna)[-c(531:532)]) {
  print(paste0(n," : ",tf))
  lmer(score.rna[,tf]~celltype+(1|patient),data = score.rna)->mod.tmp
  tmp<-as.data.frame(coef(summary(mod.tmp) ))
  score.diff.list[[tf]]=tmp[2,]
  n=n+1
}
score.diff = do.call(rbind,score.diff.list)
score.diff$TF = rownames(score.diff)
saveRDS(score.diff,
        file = "/liulab/xmwang/oxphos_proj/loading_data/scRNAseq/score.diff.mul.rna.rds")


#############################
# 3. integrate immune pipeline 
#############################

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
comb.imm.list = list()
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



