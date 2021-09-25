TRIM = function(pathways, cistrome.hdf5.path, single.cell.reference = NULL) {
  require(data.table)
  require(tidyverse)
	# regulatory (metabolic) module
	metabolic.result = regulatory.pipeline(pathways.list, cistrome.location=cistrome.hdf5.path)
	# immune module 
	immune.result = immune.module()
	# integration 
	result = integration(metabolic.result, immune.result)

	list(result, metabolic.result=metabolic.result, immune.module=immune.result)

}

## Create plots (used in the main figure) using TRIM results
