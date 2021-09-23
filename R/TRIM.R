TRIM = function(pathways.list, cistrome.hdf5.path, single.cell.reference = NULL) {
	# regulatory (metabolic) module
	metabolic.result = metabolic.pipeline(pathways.list, cistrome.dataset.path=cistrome.dataset.path)
	# immune module 
	immune.result = immune.module(single.cell.reference)
	# integration 
	result = integration(metabolic.result, immune.result)

	list(result, metabolic.result=metabolic.result, immune.module=immune.result)

}

## Create plots (used in the main figure) using TRIM results
