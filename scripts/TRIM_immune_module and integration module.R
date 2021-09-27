##################
#TRIM model 
#################

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
library("factoextra")
library(remef)
###############
#Part 1: immune pipeline
###############
feature.immune<- readxl::read_xlsx("~/OXPHOS/TCGA/TCGA_paper/1-s2.0-S1074761318301213-mmc2.xlsx")
tf.exprdata = readRDS("/liulab/xmwang/oxphos_proj/loading_data/cistrome//tf.exprdata.rds")


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
str(feature.rev)
feature.rev.tcga = do.call(cbind,feature.rev)
rownames(feature.rev.tcga) = rownames(dataset.ck.impute)
write.csv(feature.rev.tcga, file = "~/OXPHOS/ALL/feature.rev.tcga.csv")
feature.rev.tcga= read.csv("~/OXPHOS/ALL/feature.rev.tcga.csv",row.names = 1)
res.pca <- prcomp(feature.rev.tcga,scale. = T)

#Check the PC1 and PC2 of cancer types in tcga
res.ind <- get_pca_ind(res.pca)
res.ind.contri<- as.data.frame(res.ind$coord[,1:10])
res.ind.contri$object <- rownames(res.ind.contri)
res.ind.contri$cancer = dataset.ck.impute$cancer
write.csv(res.ind.contri[,c(1,ncol(res.ind.contri))], row.names = F,
          file = "~/OXPHOS/Figures/Model_figure/PC1_immune.csv")

fviz_pca_var(res.pca,axes = c(1,2),
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)->p
ggsave(plot = p,filename = "~/OXPHOS/Figures/dim1_dim2.pdf")

fviz_pca_var(res.pca,axes = c(3,4),
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)->p
ggsave(plot = p,filename = "~/OXPHOS/Figures/dim3_dim4.pdf")

#check cancer type effect
my_color = c(brewer.pal(9,"Set1"),brewer.pal(8,"Dark2"),
             brewer.pal(12,"Set3"),brewer.pal(8,"Accent"))[1:length(unique(oxphos.dataset.ck$cancer))]

res.ind.contri %>%
  bind_cols(cancer = oxphos.dataset.ck$cancer) %>% 
  ggscatter(x="Dim.1",y="Dim.2",color = "cancer",xlab = "PC1",ylab = "PC2")+
  scale_color_manual(values = my_color)+xlim(c(-10,8))+
  theme(legend.position = "right")+rremove("legend")->p.cancer.no

ggsave(plot = p.cancer.no,width = 8,height = 5.5,
       filename = "~/OXPHOS/Figures/p.cancer_no_remv.pdf")

tf.exprdata$object = substr(tf.exprdata$object,1,12)

calc_tf_feature_lmer(res.pca = res.pca,dim = "Dim.1",
                     cancer_info = oxphos.dataset.ck$cancer,
                     tf.exprdata = tf.exprdata  )->lmer.tf.dim1



tf.dim.list[[tf]]<- tmp$coefficients[2,,drop=F]

as.data.frame(lmer.tf.dim1) %>%
  #rownames_to_column(var="object") %>%
  dplyr::rename(t_value="t value",pvalue="Pr(>|t|)")%>%
  mutate(Padj = (-log10(pvalue))) %>%
  ggscatter(x="Estimate",y="Padj",
            label.select =c( "ESRRA"),label = "object")->p
ggsave(plot = p,filename = "~/OXPHOS/Figures/p.cancer.lmer.pdf")
#scale the value !!
# 1. get the sd of each tf
tf.sd = apply(tf.exprdata[,-1], 2, sd)
# 2. the sd of dim1
sd(res.ind$coord[,1]) #2.459014
# 3. b*sd_pc/sd_tf
all(lmer.tf.dim1$object==rownames(tf.sd))
lmer.tf.dim1$sd = tf.sd
lmer.tf.dim1$Estimate_nor = (lmer.tf.dim1$Estimate*2.459014)/lmer.tf.dim1$sd
write.csv(lmer.tf.dim1, file = "~/OXPHOS/Figures/lmer.tf.dim1_opt.csv")

lmer.tf.dim1 = fread("~/OXPHOS/Figures/lmer.tf.dim1_opt.csv",data.table = F)
lmer.tf.dim1 = column_to_rownames(lmer.tf.dim1,var = "V1")

inner_join(tf.exprdata[,c("object","ESRRA")],
           res.ind.contri,by="object") %>%
  filter(!duplicated(object)) %>%
  column_to_rownames(var="object") %>%
  dplyr::select(-cancer)->cor.ess.dim
cor.ess.dim=cor(cor.ess.dim,method = "spearman")

ggsave(plot = corrplot(cor.ess.dim,col = rev(col2(50))),
       filename = "~/OXPHOS/Figures/lmer.tf.dim1_opt_esrra.pdf")

#plot volcano 
lmer.tf.dim1 %>%
  dplyr::rename("pvalue"="Pr(>|t|)") %>%
  mutate(log10pvalue = -log10(pvalue)) %>%
  ggscatter(x="Estimate",y = "log10pvalue",title = "Immune pipeline in Bulk RNAseq")

######################
#Part 2: integrate pipeline
#####################
trim_auc_oxphos<- readSgdata(fileloc="~/OXPHOS/cistrome/tf_oxphos_effect_cistrome_auc.rds")
auc.oxphos = readRDS("~/OXPHOS//ALL/auc.oxphos.rds")
tmp<- merge(trim_auc_oxphos,auc.oxphos,by.x=1,by.y=0)
oxphos.data=tmp[,c(1:5,7:8)]
colnames(oxphos.data)[6]="auc"
trim_auc_oxphos<-oxphos.data
rm(oxphos.data)
SgImmData(trim_auc_gly = trim_auc_oxphos,lmer.tf.dim1 = lmer.tf.dim1)->oxphos.data
write.csv(oxphos.data,"~/OXPHOS/Figures/plot.data_opt_red.csv",row.names = F)
plot.data= fread("~/OXPHOS/Figures/plot.data_opt_red.csv")
PlotSg(gly.data = oxphos.data,ylab = "OXPHOS regulatory potential")->oxphos.p
ggsave(oxphos.p,filename = "~/OXPHOS/Figures/Model_figure/tf_filter1_opt.pdf",
       height = 6,width = 10)

loc= which(plot.data$coef_sd>quantile(plot.data$coef_sd,.99)|plot.data$coef_sd<quantile(plot.data$coef_sd,.02))
x.name <- c(plot.data[loc,]$object,"IKZF1","STAT1","STAT4","IRF1","IRF4","YY1","SMAD5")


xplot<- plot.data[plot.data$object%in% x.name,] %>% 
  dplyr::mutate(color_col = if_else( coef_sd>0,
                                     "Immunosuppresive",
                                     "Immunostimulatory"))%>%
  ggdotchart( x = "object", y = "Estimate.x",
              ylab = "Immune regulatory potential",xlab = "",
              color = "color_col",                                # Color by groups
              palette = c("#7570B3","#D95F02","#1B9E77","#E7298A"), # Custom color palette
              sorting = "ascending",                       # Sort value in descending order
              add = "segments",      
              alpha=0.6,# Add segments from y = 0 to dots
              add.params = list(color = "lightgray", size = 2.5), # Change segment color and size
              group = "color_col",                                # Order by groups
              dot.size = 7,                                 # Large dot size
              label = round(.$Estimate.x,2),                        # Add mpg values as dot labels
              font.label = list(color = "white", size = 8, 
                                vjust = 0.5),               # Adjust label parameters
              ggtheme = theme_pubr()                        # ggplot2 theme
  )+
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray")+
  rremove("legend")

loc= which(plot.data$auc>quantile(plot.data$auc,.97))
y.name <- c(plot.data[loc,]$object,"MXI1","SMAD5","HOXA6")

yplot<- plot.data[plot.data$object%in% y.name,] %>% 
  dplyr::mutate(color_col=as.character(color_col),
                color_col = if_else(object%in%show.name,color_col,"Not_signif"),
                color_col=factor(color_col,
                                 levels = c("Immune_inactive_OXPHOS_active",
                                            "Immune_active_OXPHOS_active", 
                                            "Not_signif")))%>%
  ggdotchart( x = "object", y = "auc",xlab = "",
              ylab = "OXPHOS regulation potential",
              color = "color_col",                                # Color by groups
              palette = c("#D95F02","#7570B3","grey","#E7298A","#1B9E77"), # Custom color palette
              sorting = "descending",                     # Sort value in descending order
              add = "segments",      
              alpha=0.6,# Add segments from y = 0 to dots
              add.params = list(color = "lightgray", size = 2.5),
              # Change segment color and size
              rotate = T,
                                           # Order by groups
              dot.size = 7,                                 # Large dot size
              label = round(.$auc,2),                        # Add mpg values as dot labels
              font.label = list(color = "white", size = 8, 
                                vjust = 0.5),
              ggtheme = theme_pubr() 
              
  )+
  geom_hline(yintercept = .5, linetype = 2, color = "brown")+
  rremove("legend")+
  scale_y_continuous(name="OXPHOS regulatory potential", 
                     limits=c(0, .75),breaks = c(0,.5,.6,.7))


p <- p + rremove("legend")

# Arranging the plot using cowplot
library(cowplot)
library(extrafont)
#fonts()
#font_import()
plot_grid(yplot, p, NULL, xplot, ncol = 2, align = "hv", 
          rel_widths = c(1.2,2.5), rel_heights = c( 2.5,1.5),
          label_fontfamily = "Helvetica")

#p + theme(title = element_text(family = 'Helvetica'))
library(svglite)
ggsave(plot =  plot_grid(yplot, p, NULL, xplot, ncol = 2, align = "hv", 
                         rel_widths = c(1.2,2.5), rel_heights = c( 2.5,1.5),
                         label_fontfamily = "Helvetica"),
       filename = "~/OXPHOS/Figures/Model_figure/tf_filter2_opt.svg",
       width = 10,height = 8.5)

col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))


plot.data%>%
  dplyr::mutate(color_lab = if_else(object%in% show.name,"Red","B")) %>%
  ggscatter(x="coef_sd",y="auc",xlab = "Immune regulatory potential",
            alpha=0.5,
            ylab = "OXPHOS regulatory potential",add = "reg.line", conf.int = TRUE)+
  stat_cor(method = "spearman")+
  geom_point(aes(x=Estimate.x,y=Estimate.y,color=color_lab))+
  ggrepel::geom_text_repel(data=plot.data[plot.data$object%in%show.name, ], 
                  aes_string(label="object"))+
  scale_color_manual(values = c("grey","red"))+rremove("legend")->p

ggsave(plot = p, filename = "~/OXPHOS/Figures/Model_figure//cor_coef.pdf")



################
#calc in each cancer type
tf.exprdata2= tf.exprdata
tf.exprdata2$object= substr(tf.exprdata2$object,1,12)
oxphos_allcancer = read.csv("~/OXPHOS/Figures/Figure_FC/oxphos_allcancer_TCGA.csv",row.names = 1)
oxphos_allcancer$object = substr(oxphos_allcancer$object,1,12)
calc_tf_feature_cor(res.pca = res.pca,
                    oxphos_allcancer = oxphos_allcancer,
                    tf.exprdata = tf.exprdata2 )->tf.cancer.dim.cor

write.csv(tf.cancer.dim.cor,file = "~/OXPHOS/ALL/tf.cancer.dim.cor.csv")

#############
#calc the power analysis
#power analysis

# inner_join(tf.exprdata[,c("object",tf)],
#            res.ind.contri,by=c("object"="object")) %>%
#   lmer(ESRRA~Dim.1+(1|cancer),data = .) ->tmp
# library(simr)
# model_ext_class <- extend(tmp, along="Dim.1", n=1000)
# model_ext_class
# sim_treat_class <- powerSim(model_ext_class, nsim=100)
# save(sim_treat_class,tf.exprdata,res.ind.contri,
#      file = "~/OXPHOS/ALL/power_analysis_data.RData")
# p_curve_treat <- powerCurve(model_ext_class, along="Dim.1")
# 
# ggsave(plot = plot(p_curve_treat),width = 5,height = 4.5,
#        filename = "~/OXPHOS/Figures/ESRRA_power.pdf")

###################
# Part3: model validation
##################
# 3.1 check the seek results

seek.file=list.files("/liulab/xmwang/oxphos_proj/loading_data/cistrome/seek_db/",
                     pattern = ".txt",full.names = T)
enrich_result_TF<-function(x,geneset=oxphos.gene){
  require(clusterProfiler)
  require(data.table)
  print(x)
  options(stringsAsFactors = F)
  tf_gene_list=read.table(x,header=F,sep="\t",fill=T)
  oxphos.cor = tf_gene_list[tf_gene_list$V1%in% geneset, ]
  return(oxphos.cor)
}
# get.ks(crispr.data = crispr.data.lfc,geneList=new_kegg.list$mito_down,
#        change=c("greater"))->ks.enr.mito
#####
#oxphos genes:entrez id
hsa_KEGG=download_KEGG("hsa")
entrez.gene = hsa_KEGG$KEGGPATHID2EXTID
oxphos.gene =entrez.gene[entrez.gene$from=="hsa00190","to"]
fat.gene = entrez.gene[entrez.gene$from=="hsa01212","to"]
gly.gene = entrez.gene[entrez.gene$from=="hsa00010","to"]
tca.gene = entrez.gene[entrez.gene$from=="hsa00020","to"]

getSeek<- function(seek.file,geneset=oxphos.gene){
  result=apply(matrix(seek.file),1,enrich_result_TF,geneset=geneset)
  names(result)=gsub(basename(seek.file),pattern = "_res.txt",replacement = "")
  result.data = do.call(rbind,result)
  result.data %>% 
    rownames_to_column(var="TF") %>%
    mutate(TF = str_remove(pattern = "\\..*",TF)) %>%
    group_by(TF) %>%
    summarise(median_val = median(V2)) %>%
    arrange(desc(median_val))-> TF.cor.seek
  return(TF.cor.seek)
}

getSeek(seek.file = seek.file ,geneset = oxphos.gene)->seek.oxphos
getSeek(seek.file = seek.file ,geneset = fat.gene)->seek.fat
getSeek(seek.file = seek.file ,geneset = gly.gene)->seek.gly
getSeek(seek.file = seek.file ,geneset = tca.gene)->seek.tca

#chose the top oxphos tfs
#!! note consider the direction of oxphos tf regulation

#merge acu modle and cor model
path.data= read.csv("~/OXPHOS/ALL/trim_all_tfpath_cor.csv",row.names = 1)
path.data <- path.data %>%
  mutate(pathway =  case_when(pathway=="KEGG_OXIDATIVE_PHOSPHORYLATION" ~"OXPHOS",
                                pathway=="KEGG_CITRATE_CYCLE_TCA_CYCLE" ~"TCA",
                                pathway=="KEGG_FATTY_ACID_METABOLISM" ~"FA",
                                pathway=="KEGG_GLYCOLYSIS_GLUCONEOGENESIS" ~"GLY")) %>%
  mutate(final_lab = paste0(object,"(",pathway,")"))

merge(trim_auc_sum,path.data,by="final_lab") ->merged.auc.cor
top.tf<- merged.auc.cor %>% 
  group_by(pathway) %>%
  filter(auc_ci>quantile(auc_ci,0.6) &
           pvalue.x<.05 &
           Cor>quantile(Cor,.6) &
           padj<.05)
# plot.data= read.csv("~/OXPHOS/Figures/plot.data_opt.csv",row.names = 1)
# loc = which( plot.data$pvalue_auc>1.3 &  
#                plot.data$auc>quantile(plot.data$auc,0.6) &
#                plot.data$Estimate > quantile(plot.data$Estimate,.6) &
#                plot.data$`Pr...t..`<.05 
#               )
# top.oxphos = plot.data[loc,"object"] 

#plot 
plotVal<- function(seek.res=seek.oxphos,color_value=c("#D95F02"),
                   allTF =unique(merged.auc.cor$gene),
                   txt_x = .25,txt_y = .35,title_txt="OXPHOS pathway",
                   topTF=top.tf[top.tf$pathway=="OXPHOS",]$gene){
  require(ggpubr)
  require(ggplot2)
  require(tidyverse)
  as.data.frame(seek.res) %>%
    mutate(class = if_else(TF%in% allTF,
                           if_else(TF%in% topTF,
                                   "Predicted_TRs",
                                   "common_TF" ),
                           "Not")) %>% drop_na() %>%
    mutate(class = factor(class,
                          levels = c("Predicted_TRs","common_TF",
                                     "Not"))) %>%
    mutate(abs_cor = abs(median_val)) ->gs.tmp.oxphos
  gs.tmp.oxphos %>%
    filter(class!="Not") %>%
    dplyr::select(TF,median_val) %>% 
    column_to_rownames(var="TF") %>% 
    get.ks(crispr.data = .,
           geneList = gs.tmp.oxphos[gs.tmp.oxphos$class=="Predicted_TRs",]$TF,
           change = "g")->p.c
  
  plotKS(x=gs.tmp.oxphos[gs.tmp.oxphos$class=="common_TF",]$median_val,
         y=gs.tmp.oxphos[gs.tmp.oxphos$class=="Predicted_TRs",]$median_val,
         lab_y = "Predicted_TRs",lab_x = "common_TF",
         color_value = c("grey",color_value))+
         annotate("text", x = txt_x, y = txt_y, 
                  label = paste0("Predicted TRs P Value : ",
                                 format.pval(p.c$median_val$pvalue,digits = 3)))+
        labs(title = title_txt)->p1
  #barcode plot
  gs.tmp.oxphos  %>%
    #filter(median_val>(-1)) %>%
    filter(class!="Not") %>%
    ggplot() +
    theme_pubr()+
    geom_vline(aes(xintercept=median_val,color=class),size=1.5)+
    #geom_text(x=-15,y=.5,aes( label=expriment),data =tmp.data[loc,] )+
    # annotate(geom="text",  label="1",
    #          color="red")+
    #facet_grid(expriment~.,scales = "free")+
    scale_color_manual(values = c(color_value,"white"))+
    rremove("legend")+rremove("xlab")+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1))->p2
  
  p<-plot_grid(p1,p2,nrow = 2,rel_heights  = c(3.7,.8),align = "v")
  
  return(p)
  
}

plotVal(seek.res = seek.gly,topTF=top.tf[top.tf$pathway=="GLY",]$gene,
        allTF =unique(merged.auc.cor$gene),title_txt="GLY pathway",
        txt_x = .05,txt_y = .35,color_value=c("#D95F02"))->p.gly
ggsave(p.gly,
       filename = "~/OXPHOS/Figures/Model_figure/model_gly_validation.pdf",height = 6,width = 6)

plotVal(seek.res = seek.fat,topTF=top.tf[top.tf$pathway=="FA",]$gene,
        allTF =unique(merged.auc.cor$gene),title_txt="FA pathway",
        txt_x = .1,txt_y = .35,color_value=c("#4E79A7"))->p.fat
ggsave(p.fat,
       filename = "~/OXPHOS/Figures/Model_figure/model_fat_validation.pdf",height = 6,width = 6)

plotVal(seek.res = seek.tca,topTF=top.tf[top.tf$pathway=="TCA",]$gene,
        allTF =unique(merged.auc.cor$gene),title_txt="TCA pathway",
        txt_x = .25,txt_y = .35,color_value=c("#76B7B2"))->p.tca
ggsave(p.tca,
       filename = "~/OXPHOS/Figures/Model_figure/model_tca_validation.pdf",height = 6,width = 6)

plotVal(seek.res = seek.oxphos,topTF=top.tf[top.tf$pathway=="OXPHOS",]$gene,
        color_value=c("#E15759"),
        allTF =unique(merged.auc.cor$gene))->p.oxphos
ggsave(p.oxphos,
       filename = "~/OXPHOS/Figures/Model_figure/model_oxphos_validation.pdf",height = 6,width = 6)

plot_grid(p.fat,p.gly,p.oxphos,p.tca,nrow = 2)
ggsave(plot_grid(p.fat,p.gly,p.oxphos,p.tca,nrow = 2),
       filename = "~/OXPHOS/Figures/Model_figure/model_seek_validation.pdf",
       height = 10,width = 12)


my_pal<-tableau_color_pal(palette = "Tableau 10", type = c("regular"), direction = 1)
show_col(my_pal(4))
# "#4E79A7" "#F28E2B" "#E15759" "#76B7B2"



#######################
#immune tf
#chose immune genes
hsa_KEGG$KEGGPATHID2NAME[hsa_KEGG$KEGGPATHID2NAME$to%in% kegg_immune$V1,"from"]->imm.ph.id
immune.gene = hsa_KEGG$KEGGPATHID2EXTID[hsa_KEGG$KEGGPATHID2EXTID$from%in% imm.ph.id,"to"] %>%
  unique()

#extract the cor from files
result.immune=apply(matrix(seek.file),1,function(x) enrich_result_TF(x,geneset=immune.gene))
names(result.immune)=gsub(basename(seek.file),pattern = "_res.txt",replacement = "")
result.data.immune = do.call(rbind,result.immune)
result.data.immune %>% 
  rownames_to_column(var="TF") %>%
  mutate(TF = str_remove(pattern = "\\..*",TF)) %>%
  group_by(TF) %>%
  summarise(median_val = median(V2)) %>%
  arrange(desc(median_val))-> TF.cor.immune

#chose immune tf
loc = which(lmer.tf.dim1$`Pr(>|t|)`<0.01 & lmer.tf.dim1$Estimate>quantile(lmer.tf.dim1$Estimate,0.8))
top.immune= lmer.tf.dim1[loc,]$object
loc = which(lmer.tf.dim1$`Pr(>|t|)`<0.01 & lmer.tf.dim1$Estimate<quantile(lmer.tf.dim1$Estimate,0.2))
low.immune= lmer.tf.dim1[loc,]$object

#plot
as.data.frame(TF.cor.immune) %>%
  mutate(class = case_when(TF%in%top.immune ~ "predicted_pos_immune_TF",
                           TF%in%low.immune ~ "predicted_neg_immune_TF",
                           TF%in% plot.data$object ~ "common_TF" ,
                           TRUE ~ "Not")) %>% drop_na() %>%
  mutate(class = factor(class,
                        levels = c("predicted_pos_immune_TF",
                                   "predicted_neg_immune_TF",
                                   "common_TF",
                                   "Not")))->gs.tmp.immune

#p value 
gs.tmp.immune %>%
  #filter(median_val>(-1)) %>%
  filter(class!="Not") %>%
  dplyr::select(TF,median_val) %>% 
  column_to_rownames(var="TF") %>% 
  get.ks(crispr.data = .,
         #geneList = gs.tmp.immune[gs.tmp.immune$class=="predicted_neg_immune_TF",]$TF,
         geneList = gs.tmp.immune[gs.tmp.immune$class=="predicted_pos_immune_TF",]$TF,
         change = "l"
  )->p.immune.pos

 
plotKS(x=gs.tmp.immune[gs.tmp.immune$class=="common_TF","median_val"],
       y=gs.tmp.immune[gs.tmp.immune$class=="predicted_pos_immune_TF","median_val"],
       z=gs.tmp.immune[gs.tmp.immune$class=="predicted_neg_immune_TF","median_val"],
       lab_y = "predicted_pos_immune_TF",lab_z = "predicted_neg_immune_TF",
       color_value = c("grey","#8DA0CB","#FC8D62")
)-> p.imm
p.imm2<- p.imm+annotate("text", x = -.1, y = .85, 
                        label = "P Value of pos TF : 5.9e-09")+
  annotate("text", x = -.1, y = .75, 
           label = "P Value of neg TF : 3.8e-13")


gs.tmp.immune  %>%
  #filter(median_val>(-1)) %>%
  filter(class!="Not") %>%
  ggplot() +
  geom_vline(aes(xintercept=median_val,color=class),size=1.5)+
  #geom_text(x=-15,y=.5,aes( label=expriment),data =tmp.data[loc,] )+
  # annotate(geom="text",  label="1",
  #          color="red")+
  #facet_grid(expriment~.,scales = "free")+
  scale_color_manual(values = c("#FC8D62","#8DA0CB","white"))+
  theme_bw()+rremove("xlab")+rremove("legend")->p2

p2
ggsave(plot = plot_grid(p.imm2,p2,nrow = 2,rel_heights  = c(3.7,.8),align = "v"),
       filename = "~/OXPHOS/Figures/Model_figure/Seek_db_immune_comp.pdf",
       height = 6,width = 6)

#########################
# 3.2 robust analysis

feature.rev.tcga= read.csv("~/OXPHOS/ALL/feature.rev.tcga.csv",row.names = 1)
set.seed(1234)
sample(1:nrow(feature.rev.tcga),nrow(feature.rev.tcga)/2,replace = F)->loc1
setdiff(1:nrow(feature.rev.tcga),loc1)->loc2
loc = list(loc1,loc2)
lmer.tf.dim1.list=list()

for (i in 1:2) {
  res.pca <- prcomp(feature.rev.tcga[loc[[i]],],scale. = T)
  calc_tf_feature_lmer(res.pca = res.pca,dim = "Dim.1",
                       cancer_info = oxphos.dataset.ck$cancer[loc[[i]]],
                       tf.exprdata = tf.exprdata  )->lmer.tf.dim1
  lmer.tf.dim1.list[[i]]=lmer.tf.dim1
}

lmer.tf.dim1.data= do.call(rbind,lmer.tf.dim1.list)
lmer.tf.dim1.data$model = rep(c("model_1","model_2"),each=691)
# because the direction of prcomp is different between two experiment and kinda random
data.frame(coef_model1 = lmer.tf.dim1.list[[1]]$Estimate,
           coef_model2 = -lmer.tf.dim1.list[[2]]$Estimate) %>%
  ggscatter(x="coef_model1",y="coef_model2",add = "reg.line", 
            ylab = "coefficient of model 2",conf.int = TRUE,    
            xlab = "coefficient of model 1")+
  stat_cor()->p
ggsave(p, filename = "~/OXPHOS/Figures/Model_figure/robust_model.pdf",height = 4,width = 4)

############################################
# 3.3 in each cancer type check the model

for (tf in unique(cistrome.info_used$factor)) {
  print(tf)
  sample.id = cistrome.info_used[cistrome.info_used$factor==tf,]$DCid
  common.id= intersect(rownames(rpdata.rp.used) ,sample.id)
  if (!is_empty(common.id)) {
    calc.joint.auc(as.matrix(rpdata.rp.used[common.id,]),oxphos.indicator)->tmp.auc
    auc.oxphos.list[[tf]]=tmp.auc
  }else{
    next
  }
  
}

common.tf <- intersect(unique(cistrome.info_used$factor),rp.data.all$TF)
canc.auc_ci=list()
for (tf in common.tf) {
  print(tf)
  tf.allcancer = rp.data[[tf]]
  tmp.data=c()
  for (canc in unique(tf.allcancer$cancer)) {
    tf.cancer = tf.allcancer[tf.allcancer$cancer==canc,-4]
    response = (tf.cancer$gene %in% new_kegg$OXIDATIVE_PHOSPHORYLATION)+0
    if (sum(response)!=0) {
      xx<- c(pROC::roc(response, tf.cancer$RP, ci=T)$ci)
      tmp.data[canc] =  xx[1] # lowest value of 95% CI
    }else{
      next
    }
  }
  if (!is_empty(tmp.data)) {
    canc.auc_ci[[tf]] = data.frame(TF=tf,
                                   cancer = names(tmp.data),
                                   AUC = tmp.data)
  }
  
}
canc.auc_ci.data = do.call(rbind,canc.auc_ci)
write.csv(canc.auc_ci.data, file="~/OXPHOS/ALL/tf.onecancer.auc.csv")
tf.cancer.dim.cor= read.csv("~/OXPHOS/ALL/tf.cancer.dim.cor.csv",row.names = 1)
canc.auc_ci.data <- canc.auc_ci.data%>%
  dplyr::mutate(cancer=gsub(cancer,pattern = ".",fixed = T,replacement = ""),
                cancer = gsub(cancer,pattern = "_[12]",replacement = ""),
                TF_cancer = paste(TF,cancer,sep = "_"))  %>%
  group_by(TF_cancer) %>%
  summarise(AUC = mean(AUC))

  
  
#cercan.tf.featur.data = read.csv("~/OXPHOS/ALL/cercan.tf.featur.data.csv",row.names = 1)

tf.cancer.dim.cor %>%
  rownames_to_column(var="cancer") %>%
  dplyr::mutate(cancer = str_remove(cancer,pattern = "\\.(.*)") )%>%
  dplyr::mutate(TF_cancer = paste(object,cancer,sep = "_"))  %>%
  inner_join(canc.auc_ci.data,by="TF_cancer")  %>%
  dplyr::mutate(group = case_when(Cor>0.2  & AUC>0.5 ~ "Immune_inactive",
                           
                           Cor< (-0.2)  & AUC> (0.5) ~ "Immune_active",
                           TRUE ~"Not_Signif")
  )%>%
  dplyr::rename(padj.immune="Padj")%>%
  dplyr::mutate(padj.immune = if_else(is.infinite(padj.immune),260,padj.immune)) ->plot.data.canc
  # dplyr::filter(padj.immune>3 & AUC>0)
plot.data.canc %>%
  ggscatter(y="AUC",xlab = "Immune regulatory potential",
            x = "Cor",
            ylab="OXPHOS regulatory potential",
            # label = "TF_cancer",label.select = list(criteria = "object=='ESRRA'"),
            #repel = T,
            color = "group",#size = "padj.immune",
            alpha=0.6
            )+
  #stat_cor()+
  #geom_smooth(method='lm',se = FALSE,color="grey")+
  geom_hline(yintercept=c(0.5), linetype="dashed",
             color = "grey", size=1)+
  geom_vline(xintercept=c(-.2,0.2), linetype="dashed",
             color = "grey", size=1)+
  theme(legend.position = "left")+
  scale_color_manual(values =c("#7570B3","#D95F02","grey","#1B9E77","#E7298A"))+
  ggrepel::geom_label_repel(data = subset(plot.data.canc, object =="ESRRA"),
                   aes(label = TF_cancer,color=group))->p

ggsave(plot = p,
       filename = "~/OXPHOS/Figures/Model_figure/one_cancer_tf_filter.pdf",
       width = 11,height = 6)


######################
#Part4 : integrate  mulitple metabolism pathways into the TRIM 
#####################

#check singnature analysis
lmer.tf.dim1 = fread("/liulab/xmwang/oxphos_proj/loading_data/TRIM_model/lmer.tf.dim1_opt.csv",
                     data.table = F)
lmer.tf.dim1 = column_to_rownames(lmer.tf.dim1,var = "V1")


library(plotrix)
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

PlotSg<- function(gly.data,ylab="Glycolysis regulatory potential"){
  loc = which(gly.data$color_col!="Not_signif" )
  show.name<- c(gly.data[loc,]$object,"ESRRA","IKZF1","STAT4","STAT1","IRF4","IRF1","CDK7")
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



trim_auc_gly<- readSgdata(fileloc="/liulab/xmwang/oxphos_proj/code/trim_auc_gly.rds")
SgImmData(trim_auc_gly = trim_auc_gly,lmer.tf.dim1 = lmer.tf.dim1)->gly.data
gly.data[gly.data$object=="CDK7",]$color_col="Immune_inactive_OXPHOS_active"
PlotSg(gly.data = gly.data)->gly.p
gly.p

trim_auc_fat<- readSgdata(fileloc="/liulab/xmwang/oxphos_proj/code/trim_auc_fat.rds")
SgImmData(trim_auc_gly = trim_auc_fat,lmer.tf.dim1 = lmer.tf.dim1)->fat.data
fat.data[fat.data$object=="CDK7",]$color_col="Immune_inactive_OXPHOS_active"
PlotSg(gly.data = fat.data,ylab = "Fat acid regulatory potential")->fat.p
fat.p


trim_auc_tca<- readSgdata(fileloc="/liulab/xmwang/oxphos_proj/code/trim_auc_tca.rds")
SgImmData(trim_auc_gly = trim_auc_tca,lmer.tf.dim1 = lmer.tf.dim1)->tca.data
tca.data[tca.data$object=="CDK7",]$color_col="Immune_inactive_OXPHOS_active"
PlotSg(gly.data = tca.data,ylab = "TCA regulatory potential")->tca.p
tca.p

trim_auc_pro<- readSgdata(fileloc="/liulab/xmwang/oxphos_proj/code/trim_auc_pro.rds")
SgImmData(trim_auc_gly = trim_auc_pro,lmer.tf.dim1 = lmer.tf.dim1)->pro.data
PlotSg(gly.data = pro.data,ylab = "pro regulatory potential")->pro.p
pro.p
ggsave(pro.p, height = 7,width = 9,
       filename = "./OXPHOS/Figures/scRNAseq_figure/pro.sg.trim.pdf")


SgImmData(trim_auc_gly = trim_auc_oxphos,lmer.tf.dim1 = lmer.tf.dim1)->oxphos.data
oxphos.data[oxphos.data$object=="CDK7",]$color_col="Immune_inactive_OXPHOS_active"
PlotSg(gly.data = oxphos.data,ylab = "OXPHOS regulatory potential")->oxphos.p
oxphos.p

# check the overlapped genes between them

c2 <- read.gmt("/liulab/xmwang/oxphos_proj/loading_data/annotation/c2.cp.kegg.v6.2.symbols.gmt")
#KEGG_CITRATE_CYCLE_TCA_CYCLE
tca.gene = c2[grep("KEGG_CITRATE_CYCLE_TCA_CYCLE",c2$ont,ignore.case = T),]$gene

meta.sg.list<-list(tca=tca.gene,
                   fat = fat.gene,
                   gly=gly.gene,
                   oxphos=new_kegg$OXIDATIVE_PHOSPHORYLATION)
venn(meta.sg.list,zcolor = "style")

# combine auc and ci from different model

calc.joint.auc = function(mat, response){
  xx = apply(mat, 1, function(tt) 
    c(pROC::roc(response, tt, ci=T)$ci)
  )
  if(ncol(xx)== 1) return(c(xx[2], xx[2] - xx[1]))
  ci = sqrt(sum((xx[2,] - xx[1,])^2))/ncol(xx)
  auc = mean(xx[2,])
  return(c(auc, ci))
}

all(trim_auc_fat$gene==trim_auc_gly$gene)

trim_auc_oxphos.sel =trim_auc_oxphos[trim_auc_oxphos$gene%in%trim_auc_tca$gene, ]
rownames(trim_auc_oxphos.sel)=trim_auc_oxphos.sel$gene
trim_auc_oxphos.sel = trim_auc_oxphos.sel[trim_auc_tca$gene,]
all(trim_auc_tca$gene==trim_auc_oxphos.sel$gene)
auc_ci_mat<- as.matrix(data.frame(oxphos = trim_auc_oxphos.sel$auc-trim_auc_oxphos.sel$CI,
                     gly = trim_auc_gly$auc-trim_auc_gly$CI,
                     tca = trim_auc_tca$auc-trim_auc_tca$CI,
                     fat = trim_auc_fat$auc-trim_auc_fat$CI))
rownames(auc_ci_mat)=trim_auc_gly$gene

ggstripchart(expr, x = "dataset",
             y = c("GATA3", "PTEN", "XBP1"),
             combine = TRUE, 
             color = "dataset", palette = "jco",
             size = 0.1, jitter = 0.2,
             ylab = "Expression", 
             add = "median_iqr",
             add.params = list(color = "gray"))

sum.model <- data.frame(auc_ci=rowMeans(auc_ci_mat),
                        TF = trim_auc_gly$gene)

inner_join(lmer.tf.dim1,sum.model,by=c("object"="TF")) %>%
  dplyr::mutate(coef_sd =case_when(Estimate>0 ~ Estimate-1.96 *`Std. Error`,
                                   Estimate<0 ~ -(abs(Estimate)-1.96 *`Std. Error`)) )%>%
  dplyr::mutate(coef_sd = if_else(Estimate>0 & coef_sd <0,0,coef_sd),
                coef_sd = if_else(Estimate<0 & coef_sd >0,0,coef_sd))->plot.data
thres.y.up = quantile(plot.data$auc_ci,.85,na.rm=T)
thres.x.up = quantile(plot.data$coef_sd,.85,na.rm=T)
thres.x.down = quantile(plot.data$coef_sd,.15,na.rm=T)
plot.data %>%
  dplyr::rename(pvalue_immune = "Pr(>|t|)") %>%
  dplyr::mutate(color_col =if_else( pvalue_immune <.05 ,
                                    case_when(coef_sd>thres.x.up  & auc_ci>thres.y.up ~ "Immune_inactive_OXPHOS_active",
                                              coef_sd<thres.x.down  & auc_ci>thres.y.up ~ "Immune_active_OXPHOS_active",
                                              TRUE ~"Not_signif"),
                                    "Not_signif")) %>% 
  dplyr::mutate(color_col = factor(color_col,
                                   levels = c("Immune_inactive_OXPHOS_active",
                                              "Immune_active_OXPHOS_active",
                                              "Not_signif"))) ->plot.data

loc = which(plot.data$color_col!="Not_signif" )
show.name<- c(plot.data[loc,]$object,"ESRRA","IKZF1","STAT4","STAT1","IRF4","IRF1","CDK7")
#plot.data[plot.data$object=="CDK7",]$color_col = "Immune_inactive_OXPHOS_active"

tumor.tf = score.diff.mul[score.diff.mul$Estimate>.1, ]$TF
cd8.tf = score.diff.mul[score.diff.mul$Estimate< -.1, ]$TF

plot.data %>%
  mutate(class= case_when(object%in%tumor.tf ~ "Tumor TF" ,
                          object%in%cd8.tf ~ "CD8 TF" ,
                          TRUE ~ "No signif")) %>%
  ggscatter(x="coef_sd",y="auc_ci",shape = "class",
            xlab = "Immune regulatory potential",
            ylab = "Metabolism regulatory potential",
            label = "object",#size = "pvalue_auc",
            repel = T,
            color = "color_col",alpha=0.5,
            #shape ="shape_col",
            label.select =show.name,
            legend.title = ""
  )+
  scale_color_manual(values =c("#D95F02","#7570B3", "grey","#E7298A","#1B9E77"))+
  theme(legend.position = "left")+
  geom_hline(yintercept = .5,size=.3,color="black",linetype="dashed")->p
ggsave(p, filename = "~/OXPHOS/Figures/Model_figure/mul.path.pdf",height = 7,width = 11)  

show.vol=plot.data[loc,]$object
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

plotSgVol(data = trim_auc_gly,show.name = show.name)
plotSgVol(data = trim_auc_fat,show.name = show.name,
          xlab = "Fat acid metabolism regulatory potential"
          )

# plot the oxphos pipeline in one figure
# individudal
trim_auc_sum = as.data.frame(rbind(trim_auc_gly,trim_auc_fat,trim_auc_oxphos,trim_auc_tca))
trim_auc_sum$auc_ci = trim_auc_sum$auc-trim_auc_sum$CI
trim_auc_sum$label = rep(c("GLY","FA","OXPHOS","TCA"),
                         c(nrow(trim_auc_gly),nrow(trim_auc_fat),nrow(trim_auc_oxphos),nrow(trim_auc_tca)))
trim_auc_sum$final_lab = paste0(trim_auc_sum$gene,"(",trim_auc_sum$label,")")

require(ggrepel)
library(ggpubr)
require(ggplot2)
setwd("~/OXPHOS/Figures/Model_figure/")
trim_auc_sum %>%
  #filter(label=="OXPHOS") %>%
  group_by(label) %>%
  mutate(padj = p.adjust(pvalue)) %>% 
  mutate(padj= -log10(padj),
         fill=if_else(padj>2 &auc_ci>quantile(auc_ci,.8),"signif","not_signif"),
         padj = if_else(padj>50,50,padj)) ->trim.tmp
save(trim_auc_sum,sum.model, file = "~/OXPHOS/ALL/trim_auc_sum.RData")
#####
load("/liulab/xmwang/oxphos_proj/loading_data/ALL/trim_auc_sum.RData")
pathway.show = trim.tmp %>%
  filter(padj>5) %>%
  top_n(5,auc_ci) 
trim.tmp %>%
  ggplot(aes_string(x="auc_ci",y="padj",fill="label",color="label"))+
  geom_point(alpha=.2)+
  #scale_color_manual(values=c("grey","purple"))+
  labs(x="Metabolism regulatory potential",y="-log10 Padj",color="",fill="")+theme_pubr()+
  geom_vline(xintercept=.5, linetype="dashed", color = "grey")+
  geom_text_repel(data=trim.tmp[trim.tmp$final_lab%in% pathway.show$final_lab,], 
                  aes_string(label="final_lab"))+
  scale_color_tableau()+scale_fill_tableau()->p
ggsave(p, filename = "~/OXPHOS/Figures/Model_figure/integreate_mul_path.pdf",height = 6,width = 8)
ggsave(p, filename = "~/OXPHOS/Figures/Model_figure/integreate_oxphos_path.pdf",
       height = 6,width = 8)

# combined avg figure 
#combine pvalue using fisher.method
head(trim_auc_sum)
trim_auc_sum[,c("gene","pvalue","label")] %>% 
  spread(key = "label",value = "pvalue") %>%
  column_to_rownames(var = "gene") ->trim.pvalue
metaseqR::fisher.method(trim.pvalue)->trim.padj


library(metaseqR)
sum.model$rank= rank(sum.model$auc_ci) 
merge(sum.model,trim.padj,by.x="TF",by.y=0) %>%
  mutate(p.adj = -log10(p.adj)) %>%
  ggscatter(x="auc_ci",y="p.adj",label = "TF",
            label.select = list(criteria = "`p.adj`>15"))




show.name = sum.model[sum.model$rank>650,] %>%
  arrange(desc(rank)) %>%
  pull(TF)
trim.tmp[trim.tmp$gene%in%show.name,]%>%
  ggstripchart(x="gene",y="auc_ci",alpha=.3,
               ylab = "Metabolism regulatory potential",
               order = rev(show.name),add = "mean_sd",
               add.params = list(color="brown"),
               rotate=T,color = "label",
              size = 3.5)+
  #scale_shape_manual(values=c(15, 16, 17,18))+
  scale_color_tableau()+scale_fill_tableau()+labs(color="",fill="")+
  rremove("ylab")->p

ggsave(p, filename = "~/OXPHOS/Figures/Model_figure/integreate_mul_avg_path.pdf",
       height = 7,width = 9)

############################
# Part5 : integrate the three models
############################
colnames(comb.imm)=colnames(lmer.tf.dim1)[1:5]
merge(comb.imm,sum.model,by.x=0,by.y="TF") %>%
  dplyr::mutate(coef_sd =case_when(Estimate>0 ~ Estimate-1.96 *`Std. Error`,
                                   Estimate<0 ~ -(abs(Estimate)-1.96 *`Std. Error`)) )%>%
  dplyr::mutate(coef_sd = if_else(Estimate>0 & coef_sd <0,0,coef_sd),
                coef_sd = if_else(Estimate<0 & coef_sd >0,0,coef_sd))->plot.data
colnames(plot.data)[1]="object"
thres.y.up = quantile(plot.data$auc_ci,.85,na.rm=T)
thres.x.up = quantile(plot.data$coef_sd,.85,na.rm=T)
thres.x.down = quantile(plot.data$coef_sd,.15,na.rm=T)
plot.data %>%
  dplyr::rename(pvalue_immune = "Pr(>|t|)") %>%
  dplyr::mutate(color_col =if_else( pvalue_immune <.05 ,
                                    case_when(coef_sd>thres.x.up  & auc_ci>thres.y.up ~ "Immune_inactive_OXPHOS_active",
                                              coef_sd<thres.x.down  & auc_ci>thres.y.up ~ "Immune_active_OXPHOS_active",
                                              TRUE ~"Not_signif"),
                                    "Not_signif")) %>% 
  dplyr::mutate(color_col = factor(color_col,
                                   levels = c("Immune_inactive_OXPHOS_active",
                                              "Immune_active_OXPHOS_active",
                                              "Not_signif"))) ->plot.data

loc = which(plot.data$color_col!="Not_signif" )
show.name<- c(plot.data[loc,]$object,"ESRRA","IKZF1","STAT4","STAT1","IRF4","IRF1","CDK7")
#plot.data[plot.data$object=="CDK7",]$color_col = "Immune_inactive_OXPHOS_active"


plot.data %>%
  ggscatter(x="coef_sd",y="auc_ci",
            xlab = "Immune regulatory potential",
            ylab = "Metabolism regulatory potential",
            label = "object",#size = "pvalue_auc",
            repel = T,
            color = "color_col",alpha=0.5,
            #shape ="shape_col",
            label.select =show.name,
            legend.title = ""
  )+
  scale_color_manual(values =c("#D95F02","#7570B3", "grey","#E7298A","#1B9E77"))+
  theme(legend.position = "left")+
  geom_hline(yintercept = .5,size=.3,color="black",linetype="dashed")+rremove("legend")->p
ggsave(p, filename = "~/OXPHOS/Figures/Model_figure/final_integrate.pdf",height = 6.5,width = 8.5)  




