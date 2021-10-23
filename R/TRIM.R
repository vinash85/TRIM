
#' @export
BipotentR = function(pathways.list, cistrome.hdf5.path, single.cell.reference = NULL, debug = FALSE) {
  require(data.table)
  require(tidyverse)
	# regulatory (regulation) module
	regulation.result = regulatory.pipeline(pathways.list, cistrome.location=cistrome.hdf5.path, debug=debug)
	# immune module 
	immune.result = immune.module()
	# integration 
	result = integrate.modules(regulation.result, immune.result)

	list(result, regulation.result=regulation.result, immune.module=immune.result)

}

## Create plots (used in the main figure) using TRIM results
