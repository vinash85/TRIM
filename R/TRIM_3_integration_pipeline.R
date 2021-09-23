
###############
#integration pipeline
###############
# immune.result: comb.imm.rds in TRIM_2_immune_pipeline.R
# metabolic.result : sum.model in TRIM_1_metabolic_pipeline.R
integrate.modules = function(metabolic.result, immune.result) {
final.result <- inner_join(immune.result,metabolic.result,by=c("object"="TF")) %>%
dplyr::mutate(coef_sd =case_when(Estimate>0 ~ Estimate-1.96 *`Std. Error`,
 Estimate<0 ~ -(abs(Estimate)-1.96 *`Std. Error`)) )%>%
dplyr::mutate(coef_sd = if_else(Estimate>0 & coef_sd <0,0,coef_sd),
  coef_sd = if_else(Estimate<0 & coef_sd >0,0,coef_sd))
thresult.y.up = quantile(final.result$auc_ci,.85,na.rm=T)
thresult.x.up = quantile(final.result$coef_sd,.85,na.rm=T)
thresult.x.down = quantile(final.result$coef_sd,.15,na.rm=T)
final.result<- final.result %>%
dplyr::rename(pvalue_immune = "Pr(>|t|)") %>%
dplyr::mutate(color_col =if_else( pvalue_immune <.05 ,
  case_when(coef_sd>thresult.x.up  & auc_ci>thresult.y.up ~ "Immune_inactive_OXPHOS_active",
    coef_sd<thresult.x.down  & auc_ci>thresult.y.up ~ "Immune_active_OXPHOS_active",
    TRUE ~"Not_signif"),
  "Not_signif")) %>% 
dplyr::mutate(color_col = factor(color_col,
 levels = c("Immune_inactive_OXPHOS_active",
  "Immune_active_OXPHOS_active",
  "Not_signif")))
final.result 
}
