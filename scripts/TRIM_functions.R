
# main functions in TRIM model


## metabolic pipeline
library(data.table)
library(tidyverse)
library(reshape2)
library(rhdf5)
library(pROC)
library(lmerTest)
library(lme4)
library(rlang)
library(ggrepel)
library(ggpubr)
library(ggplot2)
library(plotrix)
library(Matrix)
library(clusterProfiler)
library(ggrepel)


processRP<- function(RPloc="/liulab/xmwang/oxphos_proj/loading_data/cistrome/human_100kRP.hd5"){

  rpdata.id <- h5read(RPloc,"IDs")
  rpdata.rp <- h5read(RPloc,"RPnormalize")
  rpdata.ref <- h5read(RPloc,"RefSeq")
  rpdata.rp<- as.data.frame(rpdata.rp)
  rownames(rpdata.rp)<- rpdata.id
  colnames(rpdata.rp)<- rpdata.ref

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
  sg.indicator = (all.genes %in% sg)+0

  out.list=list(sg.indicator=sg.indicator,
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


FindSgTF<- function(cistrome.intgrated, rpdata.rp.used,sg.indicator,
                         filename="~/OXPHOS/cistrome/auc.rds"){
  library(pROC)
  library(lmerTest)
  library(lme4)
  library(rlang)
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

    #get pvalue
    tmp<-as.data.frame(coef(summary(mod.tmp) ))
    #get auc
    sample.id = tc.cancer$sample_id
    common.id= intersect(rownames(rpdata.rp.used) ,sample.id)
    if (!is_empty(common.id)) {

      calc.joint.auc(as.matrix(rpdata.rp.used[common.id,]),sg.indicator)->tmp.auc
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

### plot volcano plot
plotMetaSg<- function(data=meta_res,
                     xlab="Glycolysis regulatory potential",show.name=show.name){
  data %>%
    mutate(auc_ci = auc-CI,
           padj = -log10((pvalue)),
           padj = if_else(padj>20,20,padj),
           color_col=if_else(gene%in% show.name,"signif","not")) %>%
    ggscatter(x="auc_ci",y="padj",label = "gene",repel = T,
              xlab = xlab,
              ylab = "adjusted P value",
              label.select = show.name,color = "color_col",palette = c("grey","red"))+
    geom_hline(yintercept = 1.3,linetype="dashed",color="brown")

}



# input is the output from the metabolic pipeline
readSgdata<- function(meta_res){
  
  meta_res = do.call(rbind,meta_res)
  as.data.frame(meta_res) %>%
    rownames_to_column(var="gene") %>%
    dplyr::rename(pvalue= 'Pr(>|z|)') ->meta_res
  return(meta_res)
}

SgImmData<- function(meta_res,immune_res,thres.x = c(.15,.85),thres.y=.85){


  inner_join(immune_res,meta_res,by=c("object"="gene")) %>%
    dplyr::rename(tvalue_immune="t value",tvalue_oxphos="z value",
                  pvalue_immune = "Pr(>|t|)",pvalue_oxphos = "pvalue") %>%
    dplyr::mutate(pvalue_immune = (-log10(pvalue_immune)),
                  pvalue_auc = (-log10(pvalue_oxphos))) %>%
    dplyr::mutate(pvalue_auc = case_when(is.infinite(pvalue_auc)~max(pvalue_auc),
                                         TRUE ~pvalue_auc)) %>%
    dplyr::mutate(coef_sd =case_when(Estimate.x>0 ~ Estimate.x-1.96 *`Std. Error.x`,
                                     Estimate.x<0 ~ -(abs(Estimate.x)-1.96 *`Std. Error.x`)) )%>%
    dplyr::mutate(coef_sd = if_else(Estimate.x>0 & coef_sd <0,0,coef_sd),
                  coef_sd = if_else(Estimate.x<0 & coef_sd >0,0,coef_sd))->plot.data
  thres.y.up = quantile(plot.data$auc,thres.y,na.rm=T)
  thres.x.up = quantile(plot.data$coef_sd,thres.x[2],na.rm=T)
  thres.x.down = quantile(plot.data$coef_sd,thres.x[1],na.rm=T)
  plot.data %>%
    dplyr::mutate(color_col =if_else( pvalue_immune >1.3 &pvalue_auc>1.3 ,
                                      case_when(coef_sd>thres.x.up  & auc>thres.y.up ~ "Immune_inactive_Sg_active",
                                                coef_sd<thres.x.down  & auc>thres.y.up ~ "Immune_active_Sg_active",
                                                TRUE ~"Not_signif"),
                                      "Not_signif")) %>%
    dplyr::mutate(color_col = factor(color_col,
                                     levels = c("Immune_inactive_Sg_active",
                                                "Immune_active_Sg_active",
                                                "Not_signif")),
                  auc_ci = auc-CI) ->plot.data
  return(plot.data)
}

PlotSg<- function(merged_res,ylab="Glycolysis regulatory potential",show.name=show.name){
  

  merged_res %>%
    ggscatter(x="coef_sd",y="auc_ci",
              xlab = "Immune regulatory potential",
              ylab = ylab,
              label = "object",size = "pvalue_auc",repel = T,
              color = "color_col",alpha=0.5,
              #shape ="shape_col",
              label.select =show.name,
              legend.title = ""
    )+
    scale_color_manual(values =c("#D95F02","#7570B3", "grey","#E7298A","#1B9E77"))+
    theme(legend.position = "left")+rremove("legend")+
    geom_hline(yintercept = .5,size=.3,color="black",linetype="dashed")->p
  return(p)
}


calc_tf_feature_lmer<- function(res.pca=res.pca,dim="Dim.1",
                                cancer_info=oxphos.dataset.ck$cancer,
                                tf.exprdata=tf.exprdata){
  res.ind <- get_pca_ind(res.pca)
  res.ind.contri<- as.data.frame(res.ind$coord[,1:5])
  res.ind.contri$object <- rownames(res.ind.contri)
  res.ind.contri$cancer = cancer_info
  tf.dim.list=list()
  n=1
  colnames(tf.exprdata) = gsub(colnames(tf.exprdata),
                               pattern = "-",replacement = "_")
  
  for (tf in setdiff(colnames(tf.exprdata),"object")) {
    print(paste(n,tf,sep = ":"))
    inner_join(tf.exprdata[,c("object",tf)],
               res.ind.contri,by=c("object"="object")) %>%
      lmer(as.formula(paste0(tf,"~",dim,"+(1|cancer)")),data = .) %>%
      summary() ->tmp
    tf.dim.list[[tf]]<- tmp$coefficients[2,,drop=F] 
    n=n+1
  }
  tf.dim<- do.call(rbind,tf.dim.list)
  tf.dim<- as.data.frame(tf.dim)
  tf.dim$object = names(tf.dim.list)
  return(tf.dim)
}









