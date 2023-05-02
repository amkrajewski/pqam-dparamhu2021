init <- function (path) {
  source(paste(path,"/descriptors_gen_dev_predict_surf.R", sep = ""))
  source(paste(path,"/descriptors_gen_dev_predict_gsf.R", sep = ""))
  source(paste(path,"/gbm-locfit.R", sep = ""))
  load(paste(path,"/locfit_3_var_all.Rdata", sep = ""))
  model_gsf <<- model
  load(paste(path,"/locfit_surf_2_var_all.Rdata", sep = ""))
  model_surf <<- model
  return('init done')
}

HEA_pred <- function(composition = composition, path) {
  descriptors.gsf <- des_gsf(comp = composition, path = path)
  descriptors.surf <- des_surf(comp = composition, path = path)
  sel <- c(1:26, 35:38)
  descriptors.gsf <- as.matrix(descriptors.gsf[, sel])
  descriptors.surf <- as.matrix(descriptors.surf[, sel])
  pre.gsf <- predict.gbm.locfit(model_gsf, descriptors.gsf)
  pre.surf <- predict.gbm.locfit(model_surf, descriptors.surf)
  D <- pre.surf[2] / pre.gsf[2]
  result <- t(as.matrix(c(pre.gsf[2], pre.surf[2], D)))
  colnames(result) <- c("GSF", "Surface", "D parameter")
  return(result)
}
