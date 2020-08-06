#mlr_helper_4_comparison_significances
#
# mlr_comparison_helper module for obtaining significance tables for comparison of model performance
#  tested under R version 3.4.0 (2017-04-21)
#
# authors: Franz Kir?ly <f.kiraly@ucl.ac.uk>




#FUNCTION pairwise.test
# obtain a matrix of *unadjusted* pairwise tests on a possibly tupled sample
#
#usage
# pairwise.test(dataframe, testfun)
#
#arguments
#  dataframe        a data frame on which the test is conducted, with N variable columns to be tested pairwise
#  testfun          the hypothesis test of signature vector x vector x params -> object
#                    where the two input vectors are paired, the output object returns a $p.value element
#                    parameters are passed to testfun by ellipsis from pairwise.test, e.g.
#                     pairwise.test(dataframe, wilcox.test, alternative = "less", paired = T)
#                    is a valid way of calling the function
#
#optional:
#  alternative      "alternative" argument is explicitly passed to testfun
#  other parameters of testfun are passed by ellipsis
#
#value
#  an (N x N) matrix of *uncorrected* p-values testing each pair of variables in data frame using testfun
#   the row variable is the first argument of testfun, the column variable is the second argument of testfun
#    the rows and columns of the matrix are named and correspond to variable names of dataframe
#
# CAVE no.1: the p-values returned are *uncorrected*, post-hoc correction is necessary
#             the right post-hoc correction depends on the statement one wants to make,
#              Bonferroni or Holm should be sufficient (though possibly too conservative for most purposes)
# CAVE no.2: while the testing is *pairwise* (by variable/column), it is not necessarily *paired* (by sample/row),
#             this depends on the behaviour of testfun, such as setting a non-default paired = T for wilcox.test
#
#Example:
# pairwise.test(mtcars, wilcox.test, alternative = "less", paired = T)
#
# author: Franz Kiraly <f.kiraly@ucl.ac.uk>
#
pairwise.test <- function(dataframe, testfun, alternative = "two.sided", ...){

  n <- ncol(dataframe)

  pwtestmat <- as.data.frame(matrix(0, ncol = n, nrow = n))
  names(pwtestmat) <- names(dataframe)
  row.names(pwtestmat) <- names(dataframe)

  if(alternative == "two.sided")
  {
    for(i in 1:n){
      for(j in i:n){

        pwtestmat[i,j] <- testfun(dataframe[,i],dataframe[,j],alternative=alternative,...)$p.value
        pwtestmat[j,i] <- pwtestmat[i,j]

      }

    }
  }
  else
  {
    for(i in 1:n){
      for(j in 1:n){

        pwtestmat[i,j] <- testfun(dataframe[,i],dataframe[,j],alternative,...)$p.value

      }

    }
  }

  attr(pwtestmat,"method") <- testfun(dataframe[,i],dataframe[,j],alternative,...)$method

  return(pwtestmat)

}


p.adjustframe <- function(x, type = "full",...){

  n <- nrow(x)

  if(type =="full"){

    xvec <- p.adjust(unlist(x),...)

    xres <- as.data.frame(matrix(xvec, nrow = n, ncol = n))

    names(xres) <- names(x)
    row.names(xres) <- row.names(x)

    x <- xres

  }


  if(type == "row"){

    for(i in 1:n){

      x[i,] <- p.adjust(x[i,],...)

    }

  }

  return(x)

}



#FUNCTION dblindexdf
# accesses elements of a data frame by double sub-index
#
#usage
# dblindexdf(df, rowinds, colinds)
#
#arguments
#  df                        a data.frame, columns (that are addressed) of the same primitive class
#  rowinds, colinds          integer vectors of equal lengths
#
#value
#  a vector of same length as rowinds, colinds
#   i-th entry is df[rowinds[i],colinds[i]]
#
#example:
#
# dblindexdf(iris,c(1,2,3),c(3,1,2))
# iris[1,3]
# iris[2,1]
# iris[3,2]
#
#author: Franz Kiraly <f.kiraly@ucl.ac.uk>
#
dblindexdf <- function(df, rowinds, colinds){

  nrows <- nrow(df)

  matinds <- rowinds + (colinds-1)*nrows

  return(as.matrix(df)[matinds])

}



#FUNCTION comparisontest
# quantitative method comparison for a benchmark results object
#
#usage
# comparisontest(bmresults)
#
#arguments
#  bmresults         an mlr bmresults object bmresults as output by the benchmark function
#
# optional:
#  testfun           a hypothesis test function for either pairwise or groupwise testing
#                      which should follow the usual signature of R hypothesis tests
#                        Default for regression and probabilistic classification:
#                            paired wilcox.test (Wilcoxon signed-rank test, a non-parametric paired test)
#                        Default for deterministic classification:
#                            (unpaired) prop.test on correctly classified fraction
#                     if pairwise: first argument should be a vector for the first sample
#                                  second argument should be a vector for the first sample
#                                  should have an argument called "alternative"
#                                  should return a list such that the element "p.value" contains the p-value
#                     if groupwise: TODO
#  type              "pairwise" or "groupwise" test; TODO, does nothing currently
#  paggrfun          how p-values from different folds are aggregated to one summary p-value
#                     default: twice the median (apply Theorem 1 in Rychlik 1994, Distributions and
#                                                 expectations of order statistics ..., with m = n/2)
#  alternative       type of alternative hypothesis, passed to testfun as argument of the same name
#                     Default: "less" = "is goodness value of row-method *smaller* than of col-method?"
#  resids            a character string that indicates which standardized residuals are used as paired test samples
#                     for regression:
#                      "abs" = absolute residuals (default)
#                      "sq" = squared residuals (abs and sq yield the same result for rank-based tests)
#                     for deterministic classification:
#                      "idty.factor" = 0/1 loss aka indicator loss, coercion to factor with levels TRUE/FALSE
#                      "idty.numeric" = 0/1 loss aka indicator loss (default)
#                      "xor" = 1 minus 0/1 loss, numeric type
#                     for probabilistic classification:
#                      "brier" = brier residuals (default)
#                      "log" = log-loss residuals (log and brier yield the same for rank-based tests, but log is numerically unstable)
#  orderby           determines the order of rows/columns in the comparison table
#                     order is increasing by the methods' goodness statistic (of the bmresults data frame)
#                      Default: if task is regr, "mae.test.mean" for aggregate MAE
#                               if task is classif, "mmce.test.mean" for aggregate mmce
#                               if task is classif.prob, "logloss.test.mean" for aggregate log-loss
#  baselines         TODO: indicates which methods are to be considered baselines
#
#value
#  a list of lists; first level distinguishes tasks
#   on the second level, for each task, a list of objects that are comparison summaries is given
#    for each task, the second level list contains, as named elements:
#      comparison: a methods times methods matrix, where rows and columns correspond to methods compared in bmresults
#                  the entry in row i and column j is the uncorrected p-value of method i compared with method j
#                   alternative: "[less/greater]" = "is goodness value of row-method [less/greater] than of col-method?"
#                  for standard setting, this is p-value of method i being better than method j
#                   as measured by over-splits aggregated p-values of wilcox.tests on the paired absolute test residuals
#
# CAVE: the p-values supplied by the method in comparison are *uncorrected*
#       in particular, multiple comparisons made need to be post-hoc adjusted
#        (though one-vs-all comparisons such as best-vs-all need not be
#           if the comparison is equivalent to best-vs-next-best)
#
#Examples:
# library(mlr)
## result descriptions below obtained from results via mlr 2.11
##  actual results may differ slightly in future mlr versions due to differences in randomisation
#
## I. deterministic classification
#
# learners <- list(makeLearner("classif.lda", id = "lda"),
#                  makeLearner("classif.randomForest", id = "RF"),
#                  makeLearner("classif.nnet", id = "NN"),
#                  makeLearner("classif.uninformed", id = "baseline"))
#
# set.seed(4242)
# bmresults <- benchmark(learners, iris.task, cv3, list(mmce, mmce.se))
#
# comparisontest(bmresults)
## at 5% s.level, all three learners have better mmce than the uninformed baseline
##  however, the three learners' mmce is not statistically distinguishable
#
#
## IIa. binary probabilistic classification
#
# learners <- list(makeLearner("classif.lda", id = "lda", predict.type = "prob"),
#                  makeLearner("classif.randomForest", id = "RF", predict.type = "prob"),
#                  makeLearner("classif.nnet", id = "NN", predict.type = "prob"),
#                  makeLearner("classif.uninformed", id = "baseline", predict.type = "prob"))
#
# set.seed(4242)
# bmresults <- benchmark(learners, pid.task, cv3, list(logloss, multiclass.brier))
#
# bmresults
# comparisontest(bmresults)
## at 5% s.level, RF and lda are not statistically distinguishable (lda < RF is not significant after m.t. correction)
##  but both RF and lda have better log/brier loss than the uninformed baseline
##  neural network's log/brier loss performance is undistinguishable from the uninformed baseline
#
#
## IIb. multiclass probabilistic classification
#
# learners <- list(makeLearner("classif.lda", id = "lda", predict.type = "prob"),
#                  makeLearner("classif.randomForest", id = "RF", predict.type = "prob"),
#                  makeLearner("classif.nnet", id = "NN", predict.type = "prob"),
#                  makeLearner("classif.uninformed", id = "baseline", predict.type = "prob"))
#
# set.seed(4242)
# bmresults <- benchmark(learners, iris.task, cv3, list(logloss, multiclass.brier))
#
# bmresults
# comparisontest(bmresults)
## at 5% level, neural network significantly outperforms lda, RF and the uninformed baseline
##  this is due to the rank-based test which is not affected by the extreme residual outliers produced by the NN
##   lda outperforms RF, all methods outperform the baseline
#
## for the unpaired Wilcoxon test (unpaired is default in wilcox.test), run
# comparisontest(bmresults, testfun = wilcox.test)
## yields the same overall picture
#
#
## III. regression
#
# learners <- list(makeLearner("regr.lm", id = "lm"),
#                  makeLearner("regr.randomForest", id = "RF"),
#                  makeLearner("regr.nnet", id = "NN"),
#                  makeLearner("regr.trainingmean", id = "baseline"))
#
# set.seed(4242)
# bmresults <- benchmark(learners, bh.task, cv3, list(rmse, rmse.se, mae, mae.se))
#
# bmresults
# comparisontest(bmresults)
## at 5% level, RF significantly outperforms lm which significantly outperforms NN and uninformed baseline
##  NN is not significantly distinguishable from the uninformed baseline
#
#
#author: Franz Kiraly <f.kiraly@ucl.ac.uk>
#   (some bugfixes by Raphael Sonabend <raphael.sonabend.15@ucl.ac.uk> and Wilbur Zhu <wilbur.zhu.16@ucl.ac.uk>)
#
comparisontest <- function(bmresults, testfun = "default", type = "pairwise",
                           alternative = "less", paggrfun = "default",
                           resids = "default", orderby = "default",
                           baselines = c("trainmean","trainmedian")){


  # set default p-value aggregation function: 2 times median
  if(paggrfun == "default") paggrfun <- function(x,...){min(2*median(x,na.rm=T,...),1)}


  # get the list of task names
  tasknames <- names(bmresults$results)

  # get the list of predictor names
  prednames <- names(bmresults$results[[1]])
  #  (should be the same for each task as all methods are tested on all tasks

  # get the kind of prediction (classif vs regr)
  classregr.type <- bmresults$learners[[1]]$type

  # get the type of prediction (deterministic or probabilistic)
  predict.type <- bmresults$learners[[1]]$predict.type

  # resultlist will be returned
  #  highest list level is tasks, a list of comparison statistics including a comparison table
  resultlist <- list()



  # case A: task is regression.
  #  (probabilistic regression currently not implemented, so predict.type is not checked)
  #  in this case, comparison by residuals

  if(classregr.type == "regr"){

    # setting default parameters to task type specific values
    if(orderby == "default") orderby <- "mae.test.mean"
    if(resids == "default") resids <- "abs"
    if(class(testfun) == "character") if(testfun == "default") testfun <- function(...){wilcox.test(paired = T,...)}

    for(taskname in tasknames){

      residfr <- as.data.frame(bmresults$results[[c(taskname,prednames[1])]]$pred)[,c("id","truth","iter")]

      iterlist <- unique(residfr$iter)

      bmresults.frame <- as.data.frame(getBMRAggrPerformances(bmresults, task.ids = taskname))
      # commented out since as.dataframe converts - to .
      #  compare parallels in case B and C
      #bmresults.frame <- as.data.frame(getBMRAggrPerformances(bmresults, task.ids = taskname))[,paste0(taskname,".",prednames)]
      #names(bmresults.frame) <- prednames
      bmresults.res <- rank(as.data.frame(t(bmresults.frame))[,orderby])
      orderedprednames <- prednames[sort(bmresults.res,index.return = T)$ix]

      # this loop creates a data frame whose columns are those common to all "pred" frames,
      #  plus paired/grouped samples of standardized residuals for each method compared
      for(predname in orderedprednames){

        newcol <- as.data.frame(bmresults$results[[c(taskname,predname)]]$pred)[,"response",drop = F]
        names(newcol) <- paste0("resids.",predname)

        residfr <- cbind(residfr,newcol)
        residfr[,paste0("resids.",predname)] <- residfr$truth - residfr[,paste0("resids.",predname)]

        residfr[,paste0("resids.",predname)] <- switch(resids,
                                                  abs = abs(residfr[,paste0("resids.",predname)]),
                                                  sq = residfr[,paste0("resids.",predname)]^2)


      }

      samplelist <- list()

      # this loop computes comparison statistics for each test sample, then aggregates by mean
      for(iternum in iterlist){

        residsubfr <- residfr[residfr$iter == iternum,paste0("resids.",orderedprednames)]
        names(residsubfr) <- orderedprednames

        testtable <- pairwise.test(residsubfr,testfun, alternative = alternative)
        #testtable <- p.adjustframe(testtable, type = "row")

        samplelist <- c(samplelist, list(testtable))

      }

      numpreds <- length(orderedprednames)
      numiters <- length(iterlist)

      arraysmplist <- array(unlist(samplelist),c(numpreds,numpreds,numiters))
      aggrptable <- apply(arraysmplist,c(1,2),paggrfun)
      rownames(aggrptable) <- orderedprednames
      colnames(aggrptable) <- orderedprednames

      teststr <- paste0(attr(testtable,"method"),"\nalternative (row compared to column): ",alternative,"\n",
                        "CAVE: table is *not* multiple testing corrected")

      resultlist <- c(resultlist,list(list(comparison = aggrptable, test = teststr)))
      #resultlist <- c(resultlist,list(list(comparison = Reduce('+',samplelist)/length(iterlist), test = "test")))


    }

  }

  # end case A: regression


  # case B: task is deterministic classification
  #  in this case, comparison by indicator losses 1/0

  if(classregr.type == "classif" && predict.type != "prob"){

    # setting default parameters to task type specific values
    if(orderby == "default") orderby <- "mmce.test.mean"
    if(resids == "default") resids <- "idty.numeric"
    if(class(testfun) == "character") if(testfun == "default") testfun <- function(x, y, alternative = alternative,...){
      result <- prop.test(c(sum(y),sum(x)),c(length(y),length(x)),p=NULL,alternative = alternative,...)
      # fix since prop.test returns NaN if all numbers happens to be equal, should be p=1 or p=1/2
      if(is.nan(result$p.value)) result$p.value <- 1/(1+(alternative !="two.sided"))
      return(result)
      }

    for(taskname in tasknames){

      residfr <- as.data.frame(bmresults$results[[c(taskname,prednames[1])]]$pred)[,c("id","truth","iter")]

      iterlist <- unique(residfr$iter)

      bmresults.frame <- as.data.frame(getBMRAggrPerformances(bmresults, task.ids = taskname))
      bmresults.res <- rank(as.data.frame(t(bmresults.frame))[,orderby])
      orderedprednames <- prednames[sort(bmresults.res,index.return = T)$ix]

      # this loop creates a data frame whose columns are those common to all "pred" frames,
      #  plus paired/grouped samples of standardized residuals for each method compared
      for(predname in orderedprednames){

        residfr[,paste0("resids.",predname)] <- residfr$truth == as.data.frame(bmresults$results[[c(taskname,predname,"pred")]])[,"response",drop = T]

        residfr[,paste0("resids.",predname)] <- switch(resids,
                                                       idty.factor = as.factor(residfr[,paste0("resids.",predname)]),
                                                       idty.numeric = as.numeric(residfr[,paste0("resids.",predname)]),
                                                       xor = 1-as.numeric(residfr[,paste0("resids.",predname)])
                                                       )


      }

      samplelist <- list()

      # this loop computes comparison statistics for each test sample, then aggregates by mean
      for(iternum in iterlist){

        residsubfr <- residfr[residfr$iter == iternum,paste0("resids.",orderedprednames)]
        names(residsubfr) <- orderedprednames

        testtable <- pairwise.test(residsubfr, testfun, alternative = alternative)
        #testtable <- p.adjustframe(testtable, type = "row")

        samplelist <- c(samplelist, list(testtable))

      }

      numpreds <- length(orderedprednames)
      numiters <- length(iterlist)

      arraysmplist <- array(unlist(samplelist),c(numpreds,numpreds,numiters))
      aggrptable <- apply(arraysmplist,c(1,2),paggrfun)
      rownames(aggrptable) <- orderedprednames
      colnames(aggrptable) <- orderedprednames

      teststr <- paste0(attr(testtable,"method"),"\nalternative (row compared to column): ",alternative,"\n",
                        "CAVE: table is *not* multiple testing corrected")

      resultlist <- c(resultlist,list(list(comparison = aggrptable, test = teststr)))
      #resultlist <- c(resultlist,list(list(comparison = Reduce('+',samplelist)/length(iterlist), test = "test")))

    }

  }

  # end case B: deterministic classification


  # case C: task is probabilistic classification
  #  in this case, comparison by logloss-residuals

  if(classregr.type == "classif" && predict.type == "prob"){

    # setting default parameters to task type specific values
    if(orderby == "default") orderby <- "logloss.test.mean"
    if(resids == "default") resids <- "brier"
    if(class(testfun) == "character") if(testfun == "default") testfun <- function(...){wilcox.test(paired = T,...)}

    for(taskname in tasknames){

      residfr <- as.data.frame(bmresults$results[[c(taskname,prednames[1])]]$pred)[,c("id","truth","iter")]

      iterlist <- unique(residfr$iter)


      bmresults.frame <- as.data.frame(getBMRAggrPerformances(bmresults, task.ids = taskname))
      names(bmresults.frame) <- prednames
      bmresults.res <- rank(as.data.frame(t(bmresults.frame))[,orderby])
      orderedprednames <- prednames[sort(bmresults.res,index.return = T)$ix]

      # this loop creates a data frame whose columns are those common to all "pred" frames,
      #  plus paired/grouped samples of standardized residuals for each method compared
      for(predname in orderedprednames){

        classnames <- bmresults$results[[c(taskname,predname)]]$pred$task.desc$class.levels
        predictions <- getPredictionProbabilities(bmresults$results[[c(taskname,predname)]]$pred,cl = classnames)
        #classnames <- colnames(predictions)
        #truth = match(as.character(as.data.frame(bmresults$results[[c(taskname,predname)]]$pred)$truth),
         #             classnames)
        truth = match(residfr$truth, classnames)

        predictedprobs <- dblindexdf(predictions,1:nrow(predictions),truth)

        residfr[,paste0("resids.",predname)] <- predictedprobs

        residfr[,paste0("resids.",predname)] <- switch(resids,
                                                       log = -log(residfr[,paste0("resids.",predname)]),
                                                       brier = (1-residfr[,paste0("resids.",predname)])^2)


      }

      samplelist <- list()

      # this loop computes comparison statistics for each test sample, then aggregates by mean
      for(iternum in iterlist){

        residsubfr <- residfr[residfr$iter == iternum,paste0("resids.",orderedprednames)]
        names(residsubfr) <- orderedprednames

        testtable <- pairwise.test(residsubfr, testfun, alternative = alternative)
        #testtable <- p.adjustframe(testtable, type = "row")

        samplelist <- c(samplelist, list(testtable))

      }

      numpreds <- length(orderedprednames)
      numiters <- length(iterlist)

      arraysmplist <- array(unlist(samplelist),c(numpreds,numpreds,numiters))
      aggrptable <- apply(arraysmplist,c(1,2),paggrfun)
      rownames(aggrptable) <- orderedprednames
      colnames(aggrptable) <- orderedprednames

      teststr <- paste0(attr(testtable,"method"),"\nalternative (row compared to column): ",alternative,"\n",
                        "CAVE: table is *not* multiple testing corrected")

      resultlist <- c(resultlist,list(list(comparison = aggrptable, test = teststr)))
      #resultlist <- c(resultlist,list(list(comparison = Reduce('+',samplelist)/length(iterlist), test = "test")))

    }

  }

  # end case C: probabilistic classification



  # add task names to the list

  names(resultlist) <- tasknames

  return(resultlist)

}

#mlr_helper_2_standarderrors
#
# mlr_comparison_helper module for standard errors based on most frequently used error metrics
#  tested under R version 3.4.0 (2017-04-21)
#
# contains mlr measures to be interpreted as standard errors for their metric counterpart
#  for example, the measure mae.se computes the standard error for the mlr-native measure mae
#
# the following mlr standard error measures are supported:
#  regression:
#   mae.se, mse.se, rmse.se, sae.se, sse.se, rsq.se
#  deterministic classification:
#   acc.se, fdr.se, fn.se, fnr.se, fp.se, fpr.se, mmce.se, npv.se, ppv.se, tn.se, tnr.se, tp.se, tpr.se
#  probabilistic classification:
#
#
# authors: Franz Kir?ly <f.kiraly@ucl.ac.uk>
#          Raphael Sonabend <raphael.sonabend.15@ucl.ac.uk>



#-------------------------------------------------------------------------
# standard errors for regression measures, implemented as custom measures
#-------------------------------------------------------------------------

mae.se = makeMeasure(
  id = "mae.se", name = "Standard Error of Mean Absolute Error",
  properties = c("regr", "req.pred", "req.truth"),
  minimize = TRUE, best = 0, worst = Inf,
  fun = function(task, model, pred, feats, extra.args) {
    measureMAE.SE(pred$data$truth,pred$data$response)
  }
)
measureMAE.SE = function(truth,response){
  N <- length(response)

  # standard error estimated as standard error of the mean of he sample of absolute residuals
  return(sd(abs(response-truth))/sqrt(N))
}


mse.se = makeMeasure(
  id = "mse.se", name = "Standard Error of Mean Squared Error",
  properties = c("regr", "req.pred", "req.truth"),
  minimize = TRUE, best = 0, worst = Inf,
  fun = function(task, model, pred, feats, extra.args) {
    measureMSE.SE(pred$data$truth,pred$data$response)
  }
)
measureMSE.SE = function(truth,response){
  N <- length(response)

  # standard error estimated as standard error of the mean of the sample of squared residuals
  return(sd((response-truth)^2)/sqrt(N))
}


rmse.se = makeMeasure(
  id = "rmse.se", name = "Standard Error of Root Mean Squared Error",
  properties = c("regr", "req.pred", "req.truth"),
  minimize = TRUE, best = 0, worst = Inf,
  fun = function(task, model, pred, feats, extra.args) {
    measureRMSE.SE(pred$data$truth,pred$data$response)
  }
)
measureRMSE.SE = function(truth,response){
  N <- length(response)
  # estimated by first order approximation of sd(result)? = Var(sqrt(MSE))
  #  first order approximation: Var(f(x)) = (df/dx)? * Var(x)
  #   substituted: Var(sqrt(MSE)) = 1/(2 sqrt(MSE))*Var(MSE)
  sqresids <- (response-truth)^2
  rmse <- sqrt(mean(sqresids))
  mse.se <- sd(sqresids)/sqrt(N)
  return(mse.se/(2*rmse))
}


sae.se = makeMeasure(
  id = "sae.se", name = "Standard Error of Sum of Absolute Errors",
  properties = c("regr", "req.pred", "req.truth"),
  minimize = TRUE, best = 0, worst = Inf,
  fun = function(task, model, pred, feats, extra.args) {
    measureSAE.SE(pred$data$truth,pred$data$response)
  }
)
measureSAE.SE = function(truth,response){

  N <- length(response)

  # standard error estimated as standard error of the sum of the sample of absolute residuals
  return((sd(abs(response-truth))*sqrt(N))*((N-1)/N))

}

sse.se = makeMeasure(
  id = "sse.se", name = "Standard Error of Sum of Squared Errors",
  properties = c("regr", "req.pred", "req.truth"),
  minimize = TRUE, best = 0, worst = Inf,
  fun = function(task, model, pred, feats, extra.args) {
    measureSSE.SE(pred$data$truth,pred$data$response)
  }
)
measureSSE.SE = function(truth,response){
  N <- length(response)

  # standard error estimated as standard error of the sum of the sample of squared residuals
  return((sd((response-truth)^2)*sqrt(N))*((N-1)/N))
}

rsq.se = makeMeasure(
  id = "rsq.se", name = "Standard Error of Coefficient of Determination",
  properties = c("regr", "req.pred", "req.truth"),
  minimize = TRUE, best = 0, worst = Inf,
  fun = function(task, model, pred, feats, extra.args) {
    measureRSQ.SE(pred$data$truth,pred$data,response)
  }
)
measureRSQ.SE = function(truth,response){
  N <- length(response)

  # fast computation of Jackknife, O(N) instead of naive O(N²)
  truthsvec <- (truth - mean(truth))^2
  residsqvec <- (response - truth)^2

  #residsqvecmean <- mean(residsqvec)
  #rsqall <- sum((truthsqvec)^2)/sum(residsqvec)

  # separately computing jackknifed values for enumerator and denominator of R-squared
  jkresids <- sum(residsqvec) - residsqvec # enumerator = RSS
  jktruths <- sum(truthsvec) - (truthsvec*N/(N-1)) # denominator = TSS; "Ns"-factor is due to mean

  # Jackknife sample = enumerator i/denominator i
  msjksmpl <- (N-1)*jkresids/jktruths
  # i-th entry of this vector contains R-squared computed from full test set,
  # minus i-th test truth/pred pair
  #  scaled by the "N-1" for the jackknife pseudo-inputs

  return(sd(msjksmpl)/sqrt(N))
}

Mape.se = makeMeasure(
  id = "Mape.se", name = "Standard Error of Mean absolute percentage error",
  properties = c("regr", "req.pred", "req.truth"),
  minimize = TRUE, best = 0, worst = Inf,
  fun = function(task, model, pred, feats, extra.args) {
    measureMape.SE(pred$data$truth,pred$data$response)
  }
)
measureMape.SE = function(truth,response){

  N <- length(response)

  p <- abs(truth-response)/truth

  return(sd(p)/sqrt(N))

}

msle.se = makeMeasure(
  id = "msle.se", name = "Standard Error of the Mean Squared Logarithmic Error",
  properties = c("regr", "req.pred", "req.truth"),
  minimize = TRUE, best = 0, worst = Inf,
  fun = function(task, model, pred, feats, extra.args) {
    measureMSLE.SE(pred$data$truth,pred$data$response)
  }
)
measureMSLE.SE = function(truth,response){

  N <- length(response)

  p <- ((log(response + 1, exp(1)) - log(truth + 1, exp(1)))^2)

  return(sd(p)/sqrt(N))
}

#-----------------------------------------------------------
# standard errors for deterministic classification measures
#-----------------------------------------------------------

acc.se = makeMeasure(
  id = "acc.se", name = "Standard Error of Accuracy",
  properties = c("classif", "classif.multi", "req.pred", "req.truth"),
  minimize = TRUE, best = 0, worst = Inf,
  fun = function(task, model, pred, feats, extra.args) {
    measureACC.SE(pred$data$truth,pred$data$response)
  }
)
measureACC.SE = function(truth,response){
  N <- length(response)

  p <- mean(response == truth)

  # standard error via binomial approximation
  #  equivalent to Jackknife error estimate
  return(sqrt(p*(1-p)/N))
}

# accuracy equals mmce up to sign and added constant
mmce.se <- acc.se
mmce.se$id <- "mmce.se"
mmce.se$name <- "Standard Error of Mean Misclassification Error"
mmce.se$fun <- function(task, model, pred, feats, extra.args){
  measureMMCE.SE(pred$data$truth,pred$data$response)
}
measureMMCE.SE <- measureACC.SE

tp.se = makeMeasure(
  id = "tp.se", name = "Standard Error of True Positives",
  properties = c("classif", "req.pred", "req.truth"),
  minimize = TRUE, best = 0, worst = Inf,
  fun = function(task, model, pred, feats, extra.args) {
    measureTP.SE(pred$data$truth,pred$data$response,pred$task.desc$positive)
  }
)
measureTP.SE <- function(truth,response,positive){
  N <- length(response)

  TP <- sum(truth == response & response == positive)
  p <- TP/N

  # standard error via binomial approximation
  return(sqrt(p*(1-p)*N))
}


fp.se = makeMeasure(
  id = "fp.se", name = "Standard Error of False Positives",
  properties = c("classif", "req.pred", "req.truth"),
  minimize = TRUE, best = 0, worst = Inf,
  fun = function(task, model, pred, feats, extra.args) {
    measureFP.SE(pred$data$truth,pred$data$response,pred$task.desc$positive)
  }
)
measureFP.SE <- function(truth,response,positive){
  N <- length(response)

  FP <- sum(truth!=response & response == positive)
  p <- FP/N

  # standard error via binomial approximation
  return(sqrt(p*(1-p)*N))
}

tn.se = makeMeasure(
  id = "tn.se", name = "Standard Error of True Negatives",
  properties = c("classif", "req.pred", "req.truth"),
  minimize = TRUE, best = 0, worst = Inf,
  fun = function(task, model, pred, feats, extra.args) {
    measureTN.SE(pred$data$truth,pred$data$response,pred$task.desc$negative)
  }
)
measureTN.SE <- function(truth,response,negative){
  N <- length(response)

  TN <- sum(truth == response & response == negative)
  p <- TN/N

  # standard error via binomial approximation
  return(sqrt(p*(1-p)*N))
}


fn.se = makeMeasure(
  id = "fn.se", name = "Standard Error of False Neatives",
  properties = c("classif", "req.pred", "req.truth"),
  minimize = TRUE, best = 0, worst = Inf,
  fun = function(task, model, pred, feats, extra.args) {
    measureFN.SE(pred$data$truth,pred$data$response,pred$task.desc$negative)

  }
)
measureFN.SE <- function(truth,response,negative){
  N <- length(response)

  FN <- sum(truth != response & response == negative)
  p <- FN/N

  # standard error via binomial approximation
  return(sqrt(p*(1-p)*N))
}


tpr.se = makeMeasure(
  id = "tpr.se", name = "Standard Error of True Positive Rate",
  properties = c("classif", "req.pred", "req.truth"),
  minimize = TRUE, best = 0, worst = Inf,
  fun = function(task, model, pred, feats, extra.args) {
    measureTPR.SE(pred$data$truth,pred$data$response,pred$task.desc$positive)
  }
)
measureTPR.SE <- function(truth,response,positive){
  N <- length(response)

  positives <- positive == truth
  positivepred <- positive == response

  p <- mean(positivepred & positives)
  q <- mean(positives)

  # standard error via binomial approximation plus first order approximation of sd(result)? = Var(p/q)
  #  first order approximation: Var(f(x,y)) = (df/dx)? * Var(x) + (df/dy)? * Var(y)
  #   substituted: Var(p/q) = 1/q? * Var(p) + p?/q^4 * Var(q); Var(p) = p*(1-p)/N; Var(q) = q*(1-p)/N by Binoimial approximation
  return(sqrt(p*(1-p)/(q^2*N) + p^2*(1-q)/(q^3*N)))

}


fpr.se = makeMeasure(
  id = "fpr.se", name = "Standard Error of False Positive Rate",
  properties = c("classif", "req.pred", "req.truth"),
  minimize = TRUE, best = 0, worst = Inf,
  fun = function(task, model, pred, feats, extra.args) {
    measureFPR.SE(pred$data$truth,pred$data$response,pred$task.desc$negative,pred$task.desc$positive)
  }
)
measureFPR.SE <- function(truth,response,negative,positive){
  N <- length(response)

  negatives <- negative == truth
  positivepred <- positive == response

  p <- mean(positivepred & negatives)
  q <- mean(negatives)

  # standard error via binomial approximation plus first order approximation of sd(result)? = Var(p/q)
  #  first order approximation: Var(f(x,y)) = (df/dx)? * Var(x) + (df/dy)? * Var(y)
  #   substituted: Var(p/q) = 1/q? * Var(p) + p?/q^4 * Var(q); Var(p) = p*(1-p)/N; Var(q) = q*(1-p)/N by Binoimial approximation
  return(sqrt(p*(1-p)/(q^2*N) + p^2*(1-q)/(q^3*N)))
}

tnr.se = makeMeasure(
  id = "tnr.se", name = "Standard Error of True Negative Rate",
  properties = c("classif", "req.pred", "req.truth"),
  minimize = TRUE, best = 0, worst = Inf,
  fun = function(task, model, pred, feats, extra.args) {
    measureTNR.SE(pred$data$truth,pred$data$response,pred$task.desc$negative)
  }
)
measureTNR.SE <- function(truth,response,negative){
  N <- length(response)

  negatives <- negative == truth
  negativepred <- negative == response

  p <- mean(negativepred & negatives)
  q <- mean(negatives)

  # standard error via binomial approximation plus first order approximation of sd(result)? = Var(p/q)
  #  first order approximation: Var(f(x,y)) = (df/dx)? * Var(x) + (df/dy)? * Var(y)
  #   substituted: Var(p/q) = 1/q? * Var(p) + p?/q^4 * Var(q); Var(p) = p*(1-p)/N; Var(q) = q*(1-p)/N by Binoimial approximation
  return(sqrt(p*(1-p)/(q^2*N) + p^2*(1-q)/(q^3*N)))
}


fnr.se = makeMeasure(
  id = "fnr.se", name = "Standard Error of False Negative Rate",
  properties = c("classif", "req.pred", "req.truth"),
  minimize = TRUE, best = 0, worst = Inf,
  fun = function(task, model, pred, feats, extra.args) {
    measureFNR.SE(pred$data$truth,pred$data$response,pred$task.desc$negative,pred$task.desc$positive)
  }
)
measureFNR.SE <- function(truth,response,negative,positive){
  N <- length(response)

  positives <- positive == truth
  negativepred <- negative == response

  p <- mean(negativepred & positives)
  q <- mean(positives)

  # standard error via binomial approximation plus first order approximation of sd(result)? = Var(p/q)
  #  first order approximation: Var(f(x,y)) = (df/dx)? * Var(x) + (df/dy)? * Var(y)
  #   substituted: Var(p/q) = 1/q? * Var(p) + p?/q^4 * Var(q); Var(p) = p*(1-p)/N; Var(q) = q*(1-p)/N by Binoimial approximation
  return(sqrt(p*(1-p)/(q^2*N) + p^2*(1-q)/(q^3*N)))
}


ppv.se = makeMeasure(
  id = "ppv.se", name = "Standard Error of Positive Predictive Value",
  properties = c("classif", "req.pred", "req.truth"),
  minimize = TRUE, best = 0, worst = Inf,
  fun = function(task, model, pred, feats, extra.args) {
    measurePPV.SE(pred$data$truth,pred$data$response,pred$task.desc$positive)
  }
)
measurePPV.SE <- function(truth,response,positive){
  N <- length(response)

  positives <- positive == truth
  positivepred <- positive == response

  p <- mean(positivepred & positives)
  q <- mean(positivepred)

  # standard error via binomial approximation plus first order approximation of sd(result)? = Var(p/q)
  #  first order approximation: Var(f(x,y)) = (df/dx)? * Var(x) + (df/dy)? * Var(y)
  #   substituted: Var(p/q) = 1/q? * Var(p) + p?/q^4 * Var(q); Var(p) = p*(1-p)/N; Var(q) = q*(1-p)/N by Binoimial approximation
  return(sqrt(p*(1-p)/(q^2*N) + p^2*(1-q)/(q^3*N)))
}


npv.se = makeMeasure(
  id = "npv.se", name = "Standard Error of Negative Predictive Value",
  properties = c("classif", "req.pred", "req.truth"),
  minimize = TRUE, best = 0, worst = Inf,
  fun = function(task, model, pred, feats, extra.args) {
    measureNPV.SE(pred$data$truth,pred$data$response,pred$task.desc$negative)
  }
)
measureNPV.SE <- function(truth,response,negative){
  N <- length(response)

  negatives <- negative == truth
  negativepred <- negative == response

  p <- mean(negativepred & negatives)
  q <- mean(negativepred)

  # standard error via binomial approximation plus first order approximation of sd(result)? = Var(p/q)
  #  first order approximation: Var(f(x,y)) = (df/dx)? * Var(x) + (df/dy)? * Var(y)
  #   substituted: Var(p/q) = 1/q? * Var(p) + p?/q^4 * Var(q); Var(p) = p*(1-p)/N; Var(q) = q*(1-p)/N by Binoimial approximation
  return(sqrt(p*(1-p)/(q^2*N) + p^2*(1-q)/(q^3*N)))
}


fdr.se = makeMeasure(
  id = "fdr.se", name = "Standard Error of False Discovery Rate",
  properties = c("classif", "req.pred", "req.truth"),
  minimize = TRUE, best = 0, worst = Inf,
  fun = function(task, model, pred, feats, extra.args) {
    measureFDR.SE(pred$data$truth,pred$data$response,pred$task.desc$positive)
  }
)
measureFDR.SE <- function(truth,response,positive){
  N <- length(response)

  negatives <- positive != truth
  positivepred <- positive == response

  p <- mean(positivepred & negatives)
  q <- mean(positivepred)

  # standard error via binomial approximation plus first order approximation of sd(result)? = Var(p/q)
  #  first order approximation: Var(f(x,y)) = (df/dx)? * Var(x) + (df/dy)? * Var(y)
  #   substituted: Var(p/q) = 1/q? * Var(p) + p?/q^4 * Var(q); Var(p) = p*(1-p)/N; Var(q) = q*(1-p)/N by Binoimial approximation
  return(sqrt(p*(1-p)/(q^2*N) + p^2*(1-q)/(q^3*N)))
}


#------------------------------------------------------------
# standard errors for probabilistic classification measures
#------------------------------------------------------------

brier.se = makeMeasure(
  id = "brier.se", name = "Standard Error of Brier Score",
  properties = c("classif", "classif.multi", "req.pred", "req.truth"),
  minimize = TRUE, best = 0, worst = Inf,
  fun = function(task, model, pred, feats, extra.args) {
    measureBrier.SE(getPredictionProbabilities(pred), pred$data$truth,
                    pred$task.desc$negative, pred$task.desc$positive)
  }
)
measureBrier.SE <- function(probabilities,truth,negative,positive){
  N <- length(truth)

  y <- as.numeric(positive == truth)
  p <- (y - probabilities)^2

  return(sd(p)/sqrt(N))
}

logloss.se = makeMeasure(
  id = "logloss.se", name = "Standard Error of Logarithmic Loss",
  properties = c("classif", "classif.multi", "req.pred", "req.truth"),
  minimize = TRUE, best = 0, worst = Inf,
  fun = function(task, model, pred, feats, extra.args) {
    measureLogloss.SE(getPredictionProbabilities(pred, cl = pred$task.desc$class.levels),
                      pred$data$truth)
  }
)
measureLogloss.SE <- function(probabilities, truth){
  eps = 1e-15
  probabilities[probabilities > 1 - eps] = 1 - eps
  probabilities[probabilities < eps] = eps
  truth = match(as.character(truth), colnames(probabilities))
  p = -log(mlr:::getRowEls(probabilities, truth))

  N <- length(probabilities)

  return(sd(p)/sqrt(N))
}


# logarithmic score equals minus log-loss
lsr.se <- logloss.se
lsr.se$id <- "lsr.se"
lsr.se$name <- "Standard Error of Logarithmic Scoring Rule"
lsr.se$fun <- function(task, model, pred, feats, extra.args){
  measureLSR.SE(getPredictionProbabilities(pred, cl = pred$task.desc$class.levels),
                pred$data$truth,pred$task.desc$positive)
}
measureLSR.SE <- measureLogloss.SE

qsr.se = makeMeasure(
  id = "qsr.se", name = "Standard Error of Quadratic Scoring Rule",
  properties = c("classif", "classif.multi", "req.prob", "req.truth"),
  minimize = TRUE, best = 0, worst = Inf,
  fun = function(task, model, pred, feats, extra.args) {
    measureQSR.SE(getPredictionProbabilities(pred, cl = pred$task.desc$class.levels),
                  pred$data$truth)
  }
)
measureQSR.SE <- function(probabilities,truth){
  N <- length(probabilities)

  if (is.null(dim(probabilities)))
    probabilities = cbind(probabilities, 1 - probabilities)
  truth = factor(truth, levels = colnames(probabilities))
  p <- rowSums((probabilities - createDummyFeatures(truth))^2)

  return(sd(p)/sqrt(N))
}

# multiclass brier equals minus qsr
multiclass.brier.se <- qsr.se
multiclass.brier.se$id <- "multiclassbrier.se"
multiclass.brier.se$name <- "Standard Error of Multiclass Brier Score"
multiclass.brier.se$fun <- function(task, model, pred, feats, extra.args){
  measureMulticlassBrier.SE(getPredictionProbabilities(pred, cl = pred$task.desc$class.levels),
                            pred$data$truth,pred$task.desc$positive)
}
measureMulticlassBrier.SE <- measureQSR.SE

