
#' @importFrom magrittr %>%
processRP <- function(cistrome.location){
  rpdata.id <- rhdf5::h5read(cistrome.location,"IDs")
  rpdata.rp <- h5read(cistrome.location,"RPnormalize")
  rpdata.ref <- h5read(cistrome.location,"RefSeq")
  rpdata.rp<- as.data.frame(rpdata.rp)
  rownames(rpdata.rp)<- rpdata.id
  colnames(rpdata.rp)<- rpdata.ref
  
  rpdata.rp$sample_id = rownames(rpdata.rp)
  rpdata.rp.melt = reshape2::melt(rpdata.rp,id= "sample_id",variable.name = "gene",value.name = "RP")
  rpdata.rp.melt$gene <-str_remove(pattern = ".*:",rpdata.rp.melt$gene) 
  data.table::setDT(rpdata.rp.melt)[, mean(RP), by = c("sample_id","gene")]->rpdata.rp.pcd
  
  return(list(rpdata.rp.pcd=rpdata.rp.pcd,
              rpdata.rp=rpdata.rp))
  
}



assemble_data <- function(cistrome.info,rpdata.rp.pcd,pathways,rpdata.rp){
  #ca: atac-seq and dnase data
  cistrome.info_used <- cistrome.info[,1:4] %>%
    dplyr::filter(!(factor_type%in% c("None","not sure","other","hm","ca")))  %>%
    dplyr::mutate(DCid=as.character(DCid))
  
  cistrome.intgrated<- dplyr::inner_join(rpdata.rp.pcd,cistrome.info_used[,c(1,3)],
                                  by=c("sample_id"="DCid"))   %>%
    dplyr::mutate(oxphos_gene = if_else(gene%in% pathways,1,0)) %>%
    dplyr::rename(RP="V1",TF="factor") 
  
  
  rpdata.rp.used = rpdata.rp[rownames(rpdata.rp)%in% cistrome.info_used$DCid,]
  rpdata.rp.used=tidyr::drop_na(rpdata.rp.used)
  
  all.genes = sapply(colnames(rpdata.rp.used), function(tt) data.table::strsplit(tt, split=":")[[1]][5])
  pathways.indicator = (all.genes %in% pathways)+0
  
  out.list=list(pathways.indicator=pathways.indicator,
                cistrome.intgrated=cistrome.intgrated,
                rpdata.rp.used=rpdata.rp.used)
  return(out.list)
}

calc.joint.auc <- function(mat, response){
  xx = apply(mat, 1, function(tt) 
    c(pROC::roc(response, as.numeric(tt), ci=T)$ci)
  )
  if(ncol(xx)== 1) return(c(xx[2], xx[2] - xx[1]))
  ci = sqrt(sum((xx[2,] - xx[1,])^2))/ncol(xx)
  auc = mean(xx[2,])
  return(c(auc, ci))
}


calc_auc_pval<- function(cistrome.intgrated, rpdata.rp.used,pathways.indicator,
                         filename=NULL) {
  tf_oxphos_effect=list()
  
  n=1
  for (tf in unique(cistrome.intgrated$TF)) {
    print(paste(n,tf,sep = ":"))
    cistrome.intgrated[cistrome.intgrated$TF==tf,] ->tc.cancer
    if (length(unique(tc.cancer$sample_id))==1) {
      stats::glm(oxphos_gene~RP,data = tc.cancer,family = binomial)->mod.tmp
    }else{
      lmerTest::glmer(oxphos_gene~RP+(1|sample_id),family = binomial,data = tc.cancer)->mod.tmp
    }
    
    #get pvalue
    tmp<-as.data.frame(coef(summary(mod.tmp) ))
    #get auc : auc indicates the regulatory potential of a TF for this signature
    sample.id = tc.cancer$sample_id
    common.id= intersect(rownames(rpdata.rp.used) ,sample.id)
    if (!is_empty(common.id)) {
      
      calc.joint.auc(as.matrix(rpdata.rp.used[common.id,]),pathways.indicator)->tmp.auc
      tmp$auc=tmp.auc[1]
      tmp$CI=tmp.auc[2]
    }else{
      next
    }
    tmp[2,]->tf_oxphos_effect[[tf]]
    if(!is.null(filnename)) saveRDS(tf_oxphos_effect,file = filename)
    n=n+1
  }
  
  return(tf_oxphos_effect)
  
}

regulatory.pipeline <- function(pathways,
                    cistrome.location = "data/human_100kRP.hd5") {
  out_filename <- NULL 
  #cistrome.info available in R/sysdata.rda which is automatically loaded
  print("Identifying regulators of input pathways...")
  rpdata.rp.data <- processRP(cistrome.location=cistrome.location)
  outlist<- assemble_data(cistrome.info = cistrome.info, 
                          pathways = pathways,
                          rpdata.rp = rpdata.rp.data$rpdata.rp,
                          rpdata.rp.pcd = rpdata.rp.data$rpdata.rp.pcd)
  auc.res<- calc_auc_pval(cistrome.intgrated = outlist$cistrome.intgrated,
                          rpdata.rp.used = outlist$rpdata.rp.used,
                          pathways.indicator = outlist$pathways.indicator,
                          filename = out_filename)
  print("complete calculation!")
  return(auc.res)
}
