
HEA_pred <- function (composition = composition){
  source("descriptors_gen_dev_predict_surf.R")
  source("descriptors_gen_dev_predict_gsf.R")
  source("gbm-locfit.R")
  
  
  load("locfit_3_var_all.Rdata")
  model_gsf <- model
  load("locfit_surf_2_var_all.Rdata")
  model_surf <- model
  
  #composition=c(12,0,4,0,2,6,0,0,0,0) 
  descriptors.gsf=des_gsf(comp=composition)
  descriptors.surf=des_surf(comp=composition)
  sel=c(1:26,35:38)
  descriptors.gsf=as.matrix(descriptors.gsf[,sel])
  descriptors.surf=as.matrix(descriptors.surf[,sel])
  pre.gsf=predict.gbm.locfit(model_gsf, descriptors.gsf)
  pre.surf=predict.gbm.locfit(model_surf, descriptors.surf)
  D=pre.surf[2]/pre.gsf[2]
  result=t(as.matrix(c(pre.gsf[2],pre.surf[2],D)))
  colnames(result) <- c("GSF","Surface","D parameter")
  print(result)
  return(result)
}
