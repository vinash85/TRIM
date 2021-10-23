# Transcription Regulator of Immune-Metabolism (TRIM) 
*Authors: Avinash Das Sahu, Xiaoman Wang, X. Shirley Liu, and Keith Flaherty

## Description 

Abnormal energy metabolism is a common theme for immune evasion in tumors. Targeting cancer energy metabolism can reinvigorate tumor immunity. Focusing on energy metabolism by oxidative phosphorylation (OXPHOS), TRIM identifes transcriptional regulator of immune-metabolism. In particular, TRIM analyzes 21,000 ChIP-seq experiments, and then determined their immune-modulatory potential in 11,000 tumors and 160,000 single-cells from cancer patients. 

Once a immune-metabollic regulator is identified, TRIM can evaluate immuno-modulatory effects of its targeting on cell lines, CRISPR knockout, bulk RNA-seq, immunotherapy response, and in single cell data. We are currently following up several top immune-metabolic regulators by in vitro and in vivo experiments.


Below, we describe instruction to setup R package of BipotentR. 
The following tutorial is designed to BipotentR pacakge overview of the kinds of comparative analyses on complex cell types that are possible using the Seurat integration procedure. Here, we address a few key goals:

* Install BipotentR R package
* Download reference data
* Identify bipotent targets of CITRATE CYCLE and immune response using BipotentR 

## Install
```r
library(devtools)
install_github("vinash85/TRIM")
```

## Download reference dataset
The hdf5 file for cistrome dataset is downloaded by `download.data`.
```r
library(TRIM)
download.data()
```

## Run BipotentR
Here we describe running the package for CITRATE CYCLE  as a input pathway. The genes within TCA cycle is stored in the variable `c2.kegg`.  
The outputs will be stored in `bipotent.targets`
```r
tca.pathway =  c2.kegg[c2.kegg$term=="KEGG_CITRATE_CYCLE_TCA_CYCLE",]$gene  
## default location download.data store file
cistrome.hdf5.path =  sprintf("%s/data//human_100kRP.hd5", system.file(package = "TRIM")) 
bipotent.targets = BipotentR(tca.pathway, cistrome.hdf5.path)
```

### License

TRIM is [licensed](https://github.com/dmcable/RCTD/blob/master/LICENSE)
under the GNU General Public License v3.0.
