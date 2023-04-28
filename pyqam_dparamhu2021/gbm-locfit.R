# Gradient Boosting Machine-Locfit: A GBM framework using local regresssion via Locfit.
# For information on Locfit see: http://ect.bell-labs.com/sl/project/locfit/index.html
#
# Copyright 2016, The Materials Project, LBNL
# Distributed under the terms of the MIT License.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the "Software"), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
# the Software, and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
# FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
# IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


# Function Name:         gbm.locfit
# Version:               V1.2
# Date:                  2015-09-16
#
# Function Description:  This function provides a cross-validation wrapper for gbm.locfit.core

gbm.locfit <- function(predictors, response, n.interactions=1, n.steps=5, learning.rate=0.15, locfit.alpha=0.7,
                       locfit.degree=1, cv.folds=10, seed.folds=NULL)
{
  # Variable checks and other preliminaries
  require(locfit)
  n.observations = length(response)
  if (is.null(n.observations) || n.observations < 2) { message("Error: response must be a vector!"); return(1) }
  predictors.dim = dim(predictors)
  if (is.null(predictors.dim)) { message("Error: predictors must be a 2-D matrix!"); return(1) }
  if (predictors.dim[1] != n.observations) { message("Error: lengths of response and predictors must match!"); return(1) }
  n.predictors = predictors.dim[2]
  if (n.predictors < 2) { message("Error: minimum of two predictors is required!"); return(1) }
  if (is.null(n.interactions) || round(n.interactions) != n.interactions || n.interactions < 1 || n.interactions > 5)
    { message("Error: n.interactions must be an integer between 1 and 5!"); return(1) }
  if (n.interactions > n.predictors)
    { message("Error: n.interactions may not exceed number of predictors!"); return(1) }
  if (is.null(n.steps) || round(n.steps) != n.steps || n.steps < 1)
    { message("Error: n.steps must be at least two!"); return(1) }
  if (is.null(learning.rate) || learning.rate <= 0.0 || learning.rate >= 1.0)
    { message("Error: learning.rate must of positive but less than 1!"); return(1) }
  if (is.null(locfit.degree) || round(locfit.degree) != locfit.degree || locfit.degree < 1)
    { message("Error: locfit.degree must be a positive integer!"); return(1) }
  if (is.null(cv.folds) || round(cv.folds) != cv.folds || cv.folds <= 1)
    { message("Error: cv.folds must be a positive integer, greater than one."); return(1) }

  # Define CV sets
  set.seed(seed.folds)
  rounds = ceiling(n.observations / cv.folds)
  sets = rep(seq(1, cv.folds), rounds)[seq(1, n.observations)]
  sets = sample(sets, n.observations, replace=FALSE)  # random permutation
  train.weights = rep(0, cv.folds);  train.mse.matrix = matrix(0, nrow=n.steps, ncol=cv.folds)
  test.weights  = rep(0, cv.folds);  test.mse.matrix  = matrix(0, nrow=n.steps, ncol=cv.folds)

  # Run CV models, dropping the test set from each run
  for ( fold in seq(1, cv.folds) ) {
    message("fold:",fold)
    train.predictors = predictors[sets != fold, ];  test.predictors = predictors[sets == fold, ]
    train.response = response[sets != fold];        test.response = response[sets == fold]
    train.weights[fold] = length(train.response);   test.weights[fold] = length(test.response)

    object = gbm.locfit.core(train.predictors, train.response, n.interactions=n.interactions, n.steps=n.steps,
              learning.rate=learning.rate, locfit.alpha=locfit.alpha, locfit.degree=locfit.degree)

    # Fill train.mse.matrix
    train.mse.matrix[ , fold] = object$mean.squared.errors

    # Fill test.mse.matrix
    test.predicted.step = rep(object$response.mean, length=length(test.predictors[, 1]) )
    for ( step in seq(1, n.steps) ) {
      cat("\r",step)
      # compute prediction error for each step, incrementally (for this fold)
      test.predicted.step = incremental.predict.gbm.locfit(object, test.predictors, test.predicted.step, step)
      test.mse.matrix[step, fold] = sum((test.response - test.predicted.step)^2) / length(test.response)
  } }

  # Normalize the train and test weights
  train.weights = train.weights / sum(train.weights)
  test.weights  = test.weights  / sum(test.weights)

  # Calculate observed.errs, prediction.errs, and prediction.stderrs
  observed.errs   = apply(train.mse.matrix, 1, weighted.mean, w=train.weights)
  prediction.errs = apply(test.mse.matrix,  1, weighted.mean, w=test.weights)
  prediction.stderrs = apply(test.mse.matrix, 1, weighted.sd, w=test.weights) / sqrt(cv.folds-1)
  ns0 = which(prediction.errs == min(prediction.errs))
  ns = min(which(prediction.errs < prediction.errs[ns0]+prediction.stderrs[ns0]))

  # Run gbm.locfit.core again with all of the data, but only to ns0 steps
  object = gbm.locfit.core(predictors, response, n.interactions=n.interactions, n.steps=ns0,
            learning.rate=learning.rate, locfit.alpha=locfit.alpha, locfit.degree=locfit.degree)

  list(locfit.list=object$locfit.list, predictor.list=object$predictor.list, learning.rate=learning.rate,
      response.mean=object$response.mean, prediction.means=object$prediction.means,
      locfit.alpha=locfit.alpha, locfit.degree=locfit.degree, cv.folds=cv.folds, observed.errs=observed.errs,
      prediction.errs=prediction.errs, prediction.stderrs=prediction.stderrs, n.steps=n.steps, ns0=ns0, ns=ns)
}


# Function Name:         gbm.locfit.core
# Version:               V1.1
# Date:                  2015-07-02
#
# Function Description:  This function provides a gbm framework for locfit.

gbm.locfit.core <- function(predictors, response, n.interactions=1, n.steps=5, learning.rate=0.15, locfit.alpha=0.7,
                            locfit.degree=1, verbose=FALSE)
{
  # Variable checks and other preliminaries
  require(locfit)
  maxk=250
  n.predictors = dim(predictors)[2]
  response.mean = mean(response)
  residual = response - response.mean  # Center response
  locfit.list = vector(mode="list", length=n.steps)
  predictor.list = vector(mode="list", length=n.steps)
  prediction.means = rep(0, length=n.steps)
  mean.squared.errors = rep(0, length=n.steps)

  for (step in seq(1, n.steps)) {
    mse.min = mean(residual^2)/learning.rate  # Set mse.min "high"
    cat("\r",step)
    if (n.interactions == 1) {
      # Run locfit on individual predictors (first level interactions only)
      for (predictor.i in seq(1, n.predictors)) {
        predictors.tmp = predictors[, predictor.i]
        locfit.i = locfit.raw(predictors.tmp, residual, alpha=locfit.alpha, deg=locfit.degree, maxk=maxk)
        predicted.i = predict(locfit.i, newdata=predictors.tmp)
        mse = mean((residual - predicted.i)^2)
        if (verbose) cat(predictor.i, ": - :", mse, "\n")
        if (mse < mse.min) {
          mse.min=mse; predicted.step = predicted.i
          predictor.list[[step]]=predictor.i; locfit.list[[step]]=locfit.i } }

    } else if (n.interactions == 2) {
      # Run locfit on pairs of predictors (first thru second level interactions)
      for (predictor.i in seq(1, n.predictors-1)) {
        for (predictor.j in seq(predictor.i+1, n.predictors)) {
          predictor.ij = c(predictor.i, predictor.j)
          predictors.tmp = predictors[, predictor.ij]
          locfit.i = locfit.raw(predictors.tmp, residual, alpha=locfit.alpha, deg=locfit.degree, maxk=maxk)
          predicted.i = predict(locfit.i, newdata=predictors.tmp)
          mse = mean((residual - predicted.i)^2)
          if (verbose) cat(predictor.i, ":", predictor.j, ":", mse, "\n")
          if (mse < mse.min) {
            mse.min=mse; predicted.step = predicted.i
            predictor.list[[step]]=predictor.ij; locfit.list[[step]]=locfit.i } } }

    } else if (n.interactions == 3) {
      # Run locfit on triples of predictors (first thru third level interactions)
      for (predictor.i in seq(1, n.predictors-2)) {
        for (predictor.j in seq(predictor.i+1, n.predictors-1)) {
          for (predictor.k in seq(predictor.j+1, n.predictors)) {
            predictor.ijk = c(predictor.i, predictor.j, predictor.k)
            predictors.tmp = predictors[, predictor.ijk]
            locfit.i = locfit.raw(predictors.tmp, residual, alpha=locfit.alpha, deg=locfit.degree, maxk=maxk)
            predicted.i = predict(locfit.i, newdata=predictors.tmp)
            mse = mean((residual - predicted.i)^2)
            if (verbose) cat(predictor.i, ":", predictor.j, ":", predictor.k, ":", mse, "\n")
            if (mse < mse.min) {
              mse.min=mse; predicted.step = predicted.i
              predictor.list[[step]]=predictor.ijk; locfit.list[[step]]=locfit.i } } } }

    } else if (n.interactions == 4) {
      # Run locfit on quads of predictors (first thru fourth level interactions)
      for (predictor.i in seq(1, n.predictors-3)) {
        for (predictor.j in seq(predictor.i+1, n.predictors-2)) {
          for (predictor.k in seq(predictor.j+1, n.predictors-1)) {
            for (predictor.l in seq(predictor.k+1, n.predictors)) {
              predictor.ijkl = c(predictor.i, predictor.j, predictor.k, predictor.l)
              predictors.tmp = predictors[, predictor.ijkl]
              locfit.i = locfit.raw(predictors.tmp, residual, alpha=locfit.alpha, deg=locfit.degree, maxk=maxk)
              predicted.i = predict(locfit.i, newdata=predictors.tmp)
              mse = mean((residual - predicted.i)^2)
              if (verbose) cat(predictor.i, ":", predictor.j, ":", predictor.k, ":", predictor.l, ":", rmse, "\n")
              if (mse < mse.min) {
                mse.min=mse; predicted.step = predicted.i
                predictor.list[[step]]=predictor.ijkl; locfit.list[[step]]=locfit.i } } } } }

    } else {
      # Run locfit on quints of predictors (first thru fifth level interactions)
      for (predictor.i in seq(1, n.predictors-4)) {
        for (predictor.j in seq(predictor.i+1, n.predictors-3)) {
          for (predictor.k in seq(predictor.j+1, n.predictors-2)) {
            for (predictor.l in seq(predictor.k+1, n.predictors-1)) {
              for (predictor.m in seq(predictor.l+1, n.predictors)) {
                predictor.ijklm = c(predictor.i, predictor.j, predictor.k, predictor.l, predictor.m)
                predictors.tmp = predictors[, predictor.ijklm]
                locfit.i = locfit.raw(predictors.tmp, residual, alpha=locfit.alpha, deg=locfit.degree, maxk=maxk)
                predicted.i = predict(locfit.i, newdata=predictors.tmp)
                mse = mean((residual - predicted.i)^2)
                if (verbose) cat(predictor.i, ":", predictor.j, ":", predictor.k, ":", predictor.l, ":", predictor.m, ":", rmse, "\n")
                if (mse < mse.min) {
                  mse.min=mse; predicted.step = predicted.i
                  predictor.list[[step]]=predictor.ijklm; locfit.list[[step]]=locfit.i } } } } } }
    }

    # Update residual
    if (verbose) cat(step, " mse.min =", mse.min, predictor.list[[step]], "\n\n")
    residual = residual - learning.rate * predicted.step
 
    prediction.means[step]    = mean(predicted.step)
    mean.squared.errors[step] = mean(residual^2)
  }
  list(locfit.list=locfit.list, predictor.list=predictor.list, learning.rate=learning.rate,
      response.mean=response.mean, prediction.means=prediction.means, mean.squared.errors=mean.squared.errors)
}


# Function Name:         plot.gbm.locfit
# Version:               V1.0
# Date:                  2015-02-28
#
# Function Description:  This function provides a plot routine for gbm.locfit.

plot.gbm.locfit <- function(object, predictors, n.steps=object$ns, i.var, n.eval=100, ...)
{
  require(locfit)
  eval.points = seq(min(predictors[, i.var]), max(predictors[, i.var]), length.out = n.eval)
  marginal.values = rep(0, n.eval)

  # Loop over evaluation points
  for (j.eval in seq(1, n.eval)) {
    tmp.predictors = predictors
    tmp.predictors[, i.var] = rep(eval.points[j.eval], dim(predictors)[1])

    marginal.values[j.eval] = mean(predict.gbm.locfit(object, tmp.predictors, n.steps))
  }
  plot(eval.points, marginal.values, ...)
}


# Function Name:         predict.gbm.locfit
# Version:               V1.0
# Date:                  2015-02-26
#
# Function Description:  This function provides a predict routine for gbm.locfit.

predict.gbm.locfit <- function(object, predictors, n.steps=object$ns)
{
  require(locfit)
  predicted = rep(object$response.mean, length=length(predictors[, 1]) )

  # Assemble predicted
  for (step in seq(1, n.steps))
    predicted = predicted + object$learning.rate *
                  predict(object$locfit.list[[step]], newdata=predictors[, object$predictor.list[[step]] ])

  predicted
}


# Function Name:         incremental.predict.gbm.locfit
# Version:               V1.0
# Date:                  2015-07-02
#
# Function Description:  This function provides a predict routine for gbm.locfit.

incremental.predict.gbm.locfit <- function(object, predictors, previous.predicted, step)
{
  require(locfit)

  incremental.predicted = previous.predicted + object$learning.rate *
                            predict(object$locfit.list[[step]], newdata=predictors[, object$predictor.list[[step]] ])

  incremental.predicted
}


# Function Name:         summary.gbm.locfit
# Version:               V1.0
# Date:                  2015-02-26
#
# Function Description:  This function provides a summary routine for gbm.locfit.

summary.gbm.locfit <- function(object, predictors, n.steps=object$ns, plot=TRUE)
{
  # Variable checks and other preliminaries
  n.observations = dim(predictors)[1]
  n.predictors = dim(predictors)[2]
  influence = vector("numeric", length=n.predictors)
  predictor.list = unique(c(object$predictor.list, recursive=TRUE))

  # Calculate predictor influence
  for ( j.pred in seq(1, n.predictors) ) {
    if (length(predictor.list[predictor.list == j.pred]) > 0) {  # predictor is used
      predictor.j.sd = sd(predictors[, j.pred])
      epsilon.j = predictor.j.sd / sqrt(n.observations)

      predictors.hi = predictors
      predictors.lo = predictors
      predictors.hi[, j.pred] = predictors.hi[, j.pred] + rep(epsilon.j, n.observations)
      predictors.lo[, j.pred] = predictors.lo[, j.pred] - rep(epsilon.j, n.observations)

      prediction.hi = predict.gbm.locfit(object, predictors.hi, n.steps)
      prediction.lo = predict.gbm.locfit(object, predictors.lo, n.steps)

      influence[j.pred] = sqrt(sum((prediction.hi - prediction.lo)^2)) / 2
    } else  # predictor is NOT used
    influence[j.pred] = 0
  }

  # Rescale influences to sum to 100 and determine order
  influence = influence / sum(influence) * 100
  order = order(influence, decreasing=TRUE)

  # Extract names
  names.predictors = colnames(predictors)

  # Report relative influence
  cat("Predictor relative influence:\n")
  for ( j.pred in seq(1, n.predictors) )
    cat(j.pred, names.predictors[order[j.pred]], influence[order[j.pred]], "\n")

  if (plot == TRUE) {
    par(las=1, mar=c(5,6,2,2))  # horizontal label text, increase y-axis margin slightly
    color = rgb(0, 1-influence[rev(order)]/max(influence), 1)
    barplot(influence[rev(order)], xlab="Relative Influence", horiz=TRUE, col=color,
        names.arg=names.predictors[rev(order)]) 
  }
  list(influence=influence, order=order)
}

weighted.sd <- function(x, w, na.rm=FALSE)
{
  if (na.rm) { w = w[i = !is.na(x)]; x = x[i] }
  sum.w = sum(w)
  sum.w2 = sum(w^2)
  mean.w = sum(x * w) / sum(w)
  sqrt((sum.w / (sum.w^2 - sum.w2)) * sum(w * (x - mean.w)^2, na.rm=na.rm))
}
