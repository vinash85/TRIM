
immune.module <- function(){
  ## no input reference immune module ouputs precalculated
  immune.output = immune.module.output %>% as.data.table %>% .[,TF:=rownames(immune.module.output)]
  return(immune.output)
}

