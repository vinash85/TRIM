#' Download data for TRIM
#' @importFrom magrittr %>%
#'
#' @param download.path Folder to store downloaded data. By default the data is stored at "data/" in the package folder 
#'
#' @export
download.data = function(download.path=NULL) {
	if(is.null(download.path)){
		download.path = system.file(package="TRIM") %>% 
		paste0(., "/data")
	}
	dir.create(download.path)
	destn.file=paste0(download.path,"/human_100kRP.hd5")
	sprintf("saving file to %s", destn.file)
	download.file("https://www.dropbox.com/s/paaf5w79c73ybps/human_100kRP.hd5?dl=1",
		destn.file)
}
