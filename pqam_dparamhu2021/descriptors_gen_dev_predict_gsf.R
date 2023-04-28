des_gsf <- function(comp = comp, path) {
  comp_n <- comp / sum(comp)
  a <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
  data.1 <- data.frame(rbind(a, comp_n))

  file.3 <- paste(path,"/elemental_features.csv", sep = "")
  data.3 <- read.csv(file = file.3, header = T)
  para_ele <- as.matrix(data.3[, 2:4])

  load(file = paste(path,"/features_gsf.Rdata", sep = ""))

  feature.gsf.select <- feature.gsf.all
  n.ob <- length(data.1[, 1])
  n.feature <- length(feature.gsf.select[1, 1,])
  n.feature.ele <- length(para_ele[1,])
  n.feature.all <- n.feature + n.feature.ele

  parameter.gsf <- array(0, c(n.ob, 10, n.feature))

  comp_ele <- as.matrix(data.1)
  for (i in seq(n.feature)) {
    parameter.gsf[, , i] <- comp_ele %*% as.matrix(feature.gsf.select[, , i])
  }

  size <- 2 * n.feature.all
  predictors <- array(0, c(n.ob, size))
  for (i in seq(n.ob)) {
    nonzvec <- which(comp_ele[i,] != 0)
    avg1 <- NULL
    for (k in seq(n.feature)) {
      para.temp <- parameter.gsf[i, nonzvec, k]
      avgtemp <- para.temp %*% comp_ele[i, nonzvec]
      dvetemp <- ((1 / (1 - sum(comp_ele[i, nonzvec]^2))) * ((para.temp - avgtemp[1])^2 %*% comp_ele[i, nonzvec]))^(1 / 2)
      avg1 <- c(avg1, avgtemp, dvetemp)
    }

    avg2 <- NULL
    for (n in seq(n.feature.ele)) {
      avgtemp <- comp_ele[i, nonzvec] %*% para_ele[nonzvec, n]
      dvetemp <- ((1 / (1 - sum(comp_ele[i, nonzvec]^2))) * ((para_ele[nonzvec, n] - avgtemp[1])^2 %*% comp_ele[i, nonzvec]))^(1 / 2)
      avg2 <- c(avg2, avgtemp, dvetemp)
    }
    avg <- c(avg1, avg2)
    predictors[i,] <- avg
  }


  all <- data.frame(predictors)
  all[is.na(all)] <- 0
  return(all)
}


