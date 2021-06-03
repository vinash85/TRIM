# Transcription Regulator of Immune-Metabolism (TRIM) 
* Authors: Avinash Das Sahu, Xiaoman Wang, Keith Flaherty and X. Shiley Liu 

## Description 

Abnormal energy metabolism is a common theme for immune evasion in tumors. Targeting cancer energy metabolism can reinvigorate tumor immunity. Focusing on energy metabolism by oxidative phosphorylation (OXPHOS), TRIM identifes transcriptional regulator of immune-metabolism. In particular, TRIM analyzes 21,000 ChIP-seq experiments, and then determined their immune-modulatory potential in 11,000 tumors and 160,000 single-cells from cancer patients. 

Once a immune-metabollic regulator is identified, TRIM can evaluate immuno-modulatory effects of its targeting on cell lines, CRISPR knockout, bulk RNA-seq, immunotherapy response, and in single cell data. We are currently following up several top immune-metabolic regulators by in vitro and in vivo experiments.

## Installation

## Input format


## Usage

```Rscript
# metabolic pipeline : prediction potential TRs binding to user-defined signature
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
                        rpdata.rp.pcd = rpdata.rp.data$rpdata.rp.pcd )


mb.pred<- FindSgTF(cistrome.intgrated = outlist$cistrome.intgrated,
                        rpdata.rp.used = outlist$rpdata.rp.used,
                        sg.indicator = outlist$sg.indicator,
                        filename = "./trim_res.rds"
                        )

```

### Example



## Citation



