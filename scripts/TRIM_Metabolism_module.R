
##run other signature by this pipeline

library(data.table)
library(tidyverse)
#use rp value
#h5ls("/liulab/xmwang/oxphos_proj/loading_data/cistrome/human_100kRP.hd5")

processRP<- function(RPloc="/liulab/xmwang/oxphos_proj/loading_data/cistrome/human_100kRP.hd5"){
  require(rhdf5)
  rpdata.id <- h5read(RPloc,"IDs")
  rpdata.rp <- h5read(RPloc,"RPnormalize")
  rpdata.ref <- h5read(RPloc,"RefSeq")
  #rpdata.rp2<- rpdata.rp
  rpdata.rp<- as.data.frame(rpdata.rp)
  rownames(rpdata.rp)<- rpdata.id
  colnames(rpdata.rp)<- rpdata.ref
  
  ################
  ###!!fast
  library(reshape2)
  rpdata.rp$sample_id = rownames(rpdata.rp)
  rpdata.rp.melt = melt(rpdata.rp,id= "sample_id",variable.name = "gene",value.name = "RP")
  rpdata.rp.melt$gene <-str_remove(pattern = ".*:",rpdata.rp.melt$gene) 
  data.table::setDT(rpdata.rp.melt)[, mean(RP), by = c("sample_id","gene")]->rpdata.rp.pcd
  
  return(list(rpdata.rp.pcd=rpdata.rp.pcd,
              rpdata.rp=rpdata.rp))
  
}



assemble_data<- function(cistrome.info,rpdata.rp.pcd,sg,rpdata.rp){
  #ca: atac-seq and dnase data
  cistrome.info_used <- cistrome.info[,1:4] %>%
    dplyr::filter(!(factor_type%in% c("None","not sure","other","hm","ca")))  %>%
    dplyr::mutate(DCid=as.character(DCid))
  
  cistrome.intgrated<- inner_join(rpdata.rp.pcd,cistrome.info_used[,c(1,3)],
                                  by=c("sample_id"="DCid"))   %>%
    mutate(oxphos_gene = if_else(gene%in% sg,1,0)) %>%
    dplyr::rename(RP="V1",TF="factor") 
  
  
  rpdata.rp.used = rpdata.rp[rownames(rpdata.rp)%in% cistrome.info_used$DCid,]
  rpdata.rp.used=drop_na(rpdata.rp.used)
  
  all.genes = sapply(colnames(rpdata.rp.used), function(tt) strsplit(tt, split=":")[[1]][5])
  oxphos.indicator = (all.genes %in% sg)+0
  
  out.list=list(oxphos.indicator=oxphos.indicator,
                cistrome.intgrated=cistrome.intgrated,
                rpdata.rp.used=rpdata.rp.used)
  return(out.list)
}

calc.joint.auc = function(mat, response){
  xx = apply(mat, 1, function(tt) 
    c(pROC::roc(response, as.numeric(tt), ci=T)$ci)
  )
  if(ncol(xx)== 1) return(c(xx[2], xx[2] - xx[1]))
  ci = sqrt(sum((xx[2,] - xx[1,])^2))/ncol(xx)
  auc = mean(xx[2,])
  return(c(auc, ci))
}


calc_auc_pval<- function(cistrome.intgrated, rpdata.rp.used,oxphos.indicator,
                         filename="~/OXPHOS/cistrome/auc.rds"){
  require(pROC)
  require(lmerTest)
  require(lme4)
  require(rlang)
  tf_oxphos_effect=list()
  
  n=1
  for (tf in unique(cistrome.intgrated$TF)) {
    print(paste(n,tf,sep = ":"))
    cistrome.intgrated[cistrome.intgrated$TF==tf,] ->tc.cancer
    if (length(unique(tc.cancer$sample_id))==1) {
      glm(oxphos_gene~RP,data = tc.cancer,family = binomial)->mod.tmp
    }else{
      glmer(oxphos_gene~RP+(1|sample_id),family = binomial,data = tc.cancer)->mod.tmp
    }
    
    # #get auc
    # predpr <- predict(mod.tmp,type=c("response"))
    # roc_obj <- roc( tc.cancer$oxphos_gene~ predpr)
    #get pvalue
    tmp<-as.data.frame(coef(summary(mod.tmp) ))
    #get auc 
    sample.id = tc.cancer$sample_id
    common.id= intersect(rownames(rpdata.rp.used) ,sample.id)
    if (!is_empty(common.id)) {
      
      calc.joint.auc(as.matrix(rpdata.rp.used[common.id,]),oxphos.indicator)->tmp.auc
      tmp$auc=tmp.auc[1]
      tmp$CI=tmp.auc[2]
    }else{
      next
    }
    tmp[2,]->tf_oxphos_effect[[tf]]
    saveRDS(tf_oxphos_effect,file = filename)
    n=n+1
  }
  
  return(tf_oxphos_effect)
  
}

####################
rpdata.rp.data<- processRP()
cistrome.info<- fread("/liulab/xmwang/oxphos_proj/loading_data/cistrome/DC_haveProcessed_20190506_filepath_qc.xls")
library(clusterProfiler)
c2 <- read.gmt("/liulab/xmwang/oxphos_proj/loading_data/annotation/c2.cp.kegg.v6.2.symbols.gmt")
#this pipeline only need to change genesets and output files
pro.gene = c2[grep("KEGG_LYSINE_DEGRADATION",c2$ont,ignore.case = T),]$gene

outlist<- assemble_data(cistrome.info = cistrome.info, 
                        sg = pro.gene,
                        rpdata.rp = rpdata.rp.data$rpdata.rp,
                        rpdata.rp.pcd = rpdata.rp.data$rpdata.rp.pcd )
saveRDS(outlist,file = "./outlist_pro.rds")

pro.auc<- calc_auc_pval(cistrome.intgrated = outlist$cistrome.intgrated,
                        rpdata.rp.used = outlist$rpdata.rp.used,
                        oxphos.indicator = outlist$oxphos.indicator,
                        filename = "./trim_auc_pro.rds"
                        )
print("complete Pathway AUC calculation!")



