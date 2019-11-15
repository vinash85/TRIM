



#identify the oxphos-related tf using new method

#run in background

FindPathwayTF<- function(filecistrome = "/liulab/xmwang/oxphos_proj/loading_data/cistrome/NFKBpath_cis_info.rds",
                         genelist = hsa_KEGG.list$`NF-kappa B signaling pathway`,
                         file_binding = "/liulab/xmwang/oxphos_proj/loading_data/cistrome/NFKBbinding.rds",
                         file_score ="~/OXPHOS/tcga.features_oxphos.csv",
                         file_cor = "~/OXPHOS/ALL/NFKB_tf.dim.csv",
                         pathway_name ="NF-kappa B signaling pathway"){
  cistrome.intgrated<- getCistromeInfo(genelist=genelist,file = filecistrome)
  tf_path_binding <-  FindBindingTFs(cistrome.intgrated= cistrome.intgrated,file = file_binding)
  tf_path_cor <- FindCorTF(pathway_name = pathway_name,
                           tfExpr_file="/liulab/xmwang/oxphos_proj/loading_data/cistrome//tf.exprdata.rds",
                           gsva_path ='/liulab/xmwang/oxphos_proj/loading_data/TCGA/TCGA_ssgsea_full/',
                           file_cor="~/OXPHOS/ALL/NFKB_tf.dim.csv")
  #integrate 2nd and 3rd analysis
  
  return(out)
}


getCistromeInfo<- function(genelist=  genelist, 
                           file="~/OXPHOS/cistrome/cistrome.intgrated.rds"){
  require(rhdf5)
  require(data.table)
  require(tidyverse)
  #use rp value
  h5ls("/liulab/xmwang/oxphos_proj/loading_data/cistrome/human_100kRP.hd5")
  rpdata.id <- h5read("/liulab/xmwang/oxphos_proj/loading_data/cistrome/human_100kRP.hd5",
                      "IDs")
  rpdata.rp <- h5read("/liulab/xmwang/oxphos_proj/loading_data/cistrome/human_100kRP.hd5",
                      "RPnormalize")
  rpdata.ref <- h5read("/liulab/xmwang/oxphos_proj/loading_data/cistrome/human_100kRP.hd5",
                       "RefSeq")
  #rpdata.rp2<- rpdata.rp
  rpdata.rp<- as.data.frame(rpdata.rp)
  rownames(rpdata.rp)<- rpdata.id
  colnames(rpdata.rp)<- rpdata.ref
  #summarize the RP to gene level
  rpdata.rp%>% 
    rownames_to_column(var = "sample_id") %>%
    gather(-sample_id,key = "gene",value = "RP")  %>%
    mutate(gene=str_remove(pattern = ".*:",gene))  ->rpdata.rp.pcd
  data.table::setDT(rpdata.rp.pcd)[, mean(RP), by = c("sample_id","gene")]->rpdata.rp.pcd
  
  #merge the database id 
  cistrome.info<- fread("/liulab/xmwang/oxphos_proj/loading_data/cistrome/DC_haveProcessed_20190506_filepath_qc.xls")
  #ca: atac-seq and dnase data
  cistrome.info_used <- cistrome.info[,1:4] %>%
    filter(!(factor_type%in% c("None","not sure","other","hm","ca")))  %>%
    mutate(DCid=as.character(DCid))
  cistrome.intgrated<- inner_join(rpdata.rp.pcd,cistrome.info_used[,c(1,3)],
                                  by=c("sample_id"="DCid"))   %>%
    mutate(pathway_gene = if_else(gene%in% genelist,1,0)) %>%
    dplyr::rename(RP="V1",TF="factor") 
  saveRDS(cistrome.intgrated,file= file)
  return(cistrome.intgrated)
}

#tf oxphos_gene cistromeid rp
FindBindingTFs<- function(cistrome.intgrated= cistrome.intgrated,file = file){
  
  library(tidyverse)
  library(lme4)
  library(lmerTest)
  library(pROC)
  tf_path_effect=list()
  n=1
  for (tf in unique(cistrome.intgrated$TF)) {
    print(paste(n,tf,sep = ":"))
    cistrome.intgrated[cistrome.intgrated$TF==tf,] ->tc.cancer
    if (length(unique(tc.cancer$sample_id))==1) {
      glm(pathway_gene~RP,data = tc.cancer,family = binomial)->mod.tmp
    }else{
      glmer(pathway_gene~RP+(1|sample_id),family = binomial,data = tc.cancer)->mod.tmp
    }
    #get auc
    predpr <- predict(mod.tmp,type=c("response"))
    roc_obj <- roc( tc.cancer$pathway_gene~ predpr)
    #get pvalue
    tmp<-as.data.frame(coef(summary(mod.tmp) ))
    tmp$auc<-auc(roc_obj)
    tmp[2,]->tf_path_effect[[tf]]
    saveRDS(tf_path_effect,file = file)
    n=n+1
    
  }
  tf_path_effectData = do.call(rbind,tf_path_effect)
  as.data.frame(tf_path_effectData) %>%
    rownames_to_column(var="gene") %>%
    dplyr::rename(pvalue= 'Pr(>|z|)') %>%
    mutate(padj=(-log10(pvalue))) ->tf_path_effectData
  return(tf_path_effectData)
}



FindCorTF <- function(pathway_name = "NF-kappa B signaling pathway",
                      tfExpr_file="/liulab/xmwang/oxphos_proj/loading_data/scRNAseq/tf.exprdata_nor.rds",
                      gsva_path ='/liulab/xmwang/oxphos_proj/loading_data/TCGA/TCGA_ssgsea_full/',
                      file_cor="~/OXPHOS/ALL/tf.dim.csv"
){
  #1. tf expr
  tf.exprdata_nor = readRDS(file = tfExpr_file)
  #2. gsva score of the pathway
  tcga.score = ExtractScore(gsva_file=gsva_file,pathway_name = pathway_name)
  #3. correlation of tf_Expr with  pathway score in tcga
  ScoreTF<- CalTFCor(pathway_name =pathway_name,
                     file =file_cor,
                     tcga.score=tcga.score,
                     tf.exprdata_nor= tf.exprdata_nor)
  return(ScoreTF)
  
}


ExtractScore<- function(gsva_file,pathway_name = pathway_name){
  gsva_file = list.files(gsva_path,pattern = "csv",full.names = T)
  gsva.list= list()
  for (file in gsva_file) {
    cancername = gsub(pattern = "TCGA_(.*)_tpm.csv",basename(file),
                      replacement = "\\1")
    file.tmp = fread(file)
    file.tmp[file.tmp$V1 %in% pathway_name,]%>%
      column_to_rownames(var = "V1")  %>% t(.) %>%
      as.data.frame() ->gsva.list[[cancername]]
  }
  
  gsva.data= do.call(rbind,gsva.list)
  
  gsva.data %>% rownames_to_column(var = "object") %>%
    separate(col = "object",
             into = c("cancer","object"),
             sep = "[.]",extra = "merge") ->tcga.score
  
  return(tcga.score)
}



CalTFCor<- function(pathway_name =pathway_name,file ="~/OXPHOS/ALL/tf.dim.csv",
                    tcga.score = tcga.score,
                    tf.exprdata_nor= tf.exprdata_nor){
  n=1
  tf.oxphos.list=list()
  #-ncol(tf.exprdata_nor)
  for (TF in setdiff(colnames(tf.exprdata_nor),"object")) {
    print(paste(n,TF,sep = ":"))
    inner_join(tcga.score,tf.exprdata_nor[,c("object",TF)],by=c("object"="object"))%>%
      lmer(as.formula(paste0("`",TF,"`","~","`",pathway_name,"`","+(1|cancer)")),data = .) %>%
      summary() ->tmp.mod
    
    tf.oxphos.list[[TF]]<- tmp.mod$coefficients[2,,drop=F]
    n=n+1
  }
  
  tf.dim<- do.call(rbind,tf.oxphos.list)
  tf.dim<- as.data.frame(tf.dim)
  tf.dim$object = names(tf.oxphos.list)
  write.csv(tf.dim,file =file,row.names = F)
  return(tf.dim)
}

##################

# re-run the ggssea

library(GSVA)
getscore<- function(file,hsa_KEGG.list){
  cancer.tpm=readRDS(file)
  cancername = gsub(basename(file),pattern = ".rds",replacement = "")
  cancer.kegg <- gsva(as.matrix(cancer.tpm), 
                      hsa_KEGG.list, 
                      min.sz=5, max.sz=500,kcdf="Gaussian", 
                      mx.diff=TRUE, verbose=FALSE, 
                      parallel.sz=10, ssgsea.norm=TRUE)
  write.csv(cancer.kegg, 
            file = paste0("/liulab/xmwang/oxphos_proj/loading_data/TCGA/TCGA_ssgsea_full/",
                          cancername,".csv"))
  
}


  