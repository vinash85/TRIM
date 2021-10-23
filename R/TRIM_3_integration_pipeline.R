
integrate.modules <- function(regulation.result, immune.result) {
  all.result <- merge(immune.result,regulation.result,by="TF", 
    suffixes = c(".immune", ".regulation"))
  all.result = all.result %>%
  dplyr::mutate(coef_sd =case_when(Estimate.immune>0 ~ Estimate.immune-1.96 *Estimate.immune,
   Estimate.immune<0 ~ -(abs(Estimate.immune)-1.96 *SE.immune)) )%>%
  dplyr::mutate(coef_sd = if_else(Estimate.immune>0 & coef_sd <0,0,coef_sd),
    coef_sd = if_else(Estimate.immune<0 & coef_sd >0,0,coef_sd))
  thresult.y.up = quantile(all.result$auc_ci,.85,na.rm=T)
  thresult.x.up = quantile(all.result$coef_sd,.85,na.rm=T)
  thresult.x.down = quantile(all.result$coef_sd,.15,na.rm=T)
  final.result<- all.result %>%
  dplyr::mutate(color_col =if_else( p.val.immune <.05 ,
    case_when(coef_sd>thresult.x.up  & auc_ci>thresult.y.up ~ "Immune_inactive_OXPHOS_active",
      coef_sd<thresult.x.down  & auc_ci>thresult.y.up ~ "Immune_active_OXPHOS_active",
      TRUE ~"Not_signif"),
    "Not_signif")) %>% 
  dplyr::mutate(color_col = factor(color_col,
   levels = c("Immune_inactive_OXPHOS_active",
    "Immune_active_OXPHOS_active",
    "Not_signif")))
  result=list(final.result=final.result, all.TFCR.estiamte=all.result)
  return(result)
}
