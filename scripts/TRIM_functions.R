
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
plotSgVol<- function(data=trim_auc_gly,
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
readSgdata<- function(fileloc="/liulab/xmwang/oxphos_proj/code/trim_auc_gly.rds"){
  trim_auc_gly = readRDS(fileloc)
  trim_auc_gly = do.call(rbind,trim_auc_gly)
  as.data.frame(trim_auc_gly) %>%
    rownames_to_column(var="gene") %>%
    dplyr::rename(pvalue= 'Pr(>|z|)') ->trim_auc_gly
  return(trim_auc_gly)
}

SgImmData<- function(trim_auc_gly,lmer.tf.dim1,thres.x = c(.15,.85),thres.y=.85){


  inner_join(lmer.tf.dim1,trim_auc_gly,by=c("object"="gene")) %>%
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
                                      case_when(coef_sd>thres.x.up  & auc>thres.y.up ~ "Immune_inactive_OXPHOS_active",
                                                coef_sd<thres.x.down  & auc>thres.y.up ~ "Immune_active_OXPHOS_active",
                                                TRUE ~"Not_signif"),
                                      "Not_signif")) %>%
    dplyr::mutate(color_col = factor(color_col,
                                     levels = c("Immune_inactive_OXPHOS_active",
                                                "Immune_active_OXPHOS_active",
                                                "Not_signif")),
                  auc_ci = auc-CI) ->plot.data
  return(plot.data)
}

PlotSg<- function(gly.data,ylab="Glycolysis regulatory potential",show.name=show.name){
  loc = which(gly.data$color_col!="Not_signif" )

  gly.data %>%
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







