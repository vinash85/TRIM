# Transcription Regulator of Immune-Metabolism (TRIM) 
* Authors: Avinash Das Sahu, Xiaoman Wang, Keith Flaherty and X. Shiley Liu 

## Description 

Abnormal energy metabolism is a common theme for immune evasion in tumors. Targeting cancer energy metabolism can reinvigorate tumor immunity. Focusing on energy metabolism by oxidative phosphorylation (OXPHOS), TRIM identifes transcriptional regulator of immune-metabolism. In particular, TRIM analyzes 21,000 ChIP-seq experiments, and then determined their immune-modulatory potential in 11,000 tumors and 160,000 single-cells from cancer patients. 

Once a immune-metabollic regulator is identified, TRIM can evaluate immuno-modulatory effects of its targeting on cell lines, CRISPR knockout, bulk RNA-seq, immunotherapy response, and in single cell data. We are currently following up several top immune-metabolic regulators by in vitro and in vivo experiments.



## Usage

1. Metabolic pipeline : prediction potential TRs binding to user-defined signature


```Rscript

source("./script/TRIM_functions.R")
# download the two files from (link)
rpdata.rp.data<- processRP(RPloc="./data/human_100kRP.hd5")
cistrome.info<- fread("./data/DC_haveProcessed_20190506_filepath_qc.xls")

# define your signature : 
c2 <- read.gmt("/liulab/xmwang/oxphos_proj/loading_data/annotation/c2.cp.kegg.v6.2.symbols.gmt")
selected.sg = c2[grep("KEGG_OXIDATIVE_PHOSPHORYLATION",c2$ont,ignore.case = T),]$gene

outlist<- assemble_data(cistrome.info = cistrome.info, 
                        sg = selected.sg,
                        rpdata.rp = rpdata.rp.data$rpdata.rp,
                        rpdata.rp.pcd = rpdata.rp.data$rpdata.rp.pcd
                        )


meta_res<- FindSgTF(cistrome.intgrated = outlist$cistrome.intgrated,
					pdata.rp.used = outlist$rpdata.rp.used,
                    sg.indicator = outlist$sg.indicator,
                    filename = "./meta_res.rds"
                    )

# visulize the result
show.name<- meta_res[meta_res$pvalue<0.001,]$object
plotMetaSg(meta_res,xlab="Affinity to bind selected gene sets ",show.name=show.name)

```


2. Model integration 

```Rscript
# download the result of immune pipeline : link (column name )
immune_res = readRDS("./data/comb.imm.rds")

meta_res<- readSgdata(meta_res.rds)
merged_res <- SgImmData(meta_res = meta_res,immune_res = immune_res)
loc = which(merged_res$color_col!="Not_signif" )
show.name= merged_res[loc,]$object
# display the integrated result
merged.plot <- PlotSg(merged_res = merged.data,ylab="Affinity to bind selected gene sets",show.name = show.name)
merged.plot

```



## Citation



