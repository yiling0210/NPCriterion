#' Marginal feature screening under the Neyman-Pearson paradigm
#'
#' @param x a design matrix
#' @param y a vector containing binary labels 0 and 1
#' @param method base classification method
#' \itemize{
#'   \item logistic: Logistic regression. \code{\link[stats]{glm}} function with family = 'binomial'
#'   \item penlog: Penalized logistic regression with LASSO penalty. \code{\link[glmnet]{glmnet}} in glmnet package
#'   \item svm: Support Vector Machines. \code{\link[e1071]{svm}} in e1071 package
#'   \item randomforest: Random Forest. \code{\link[randomForest]{randomForest}} in randomForest package
#'   \item lda: Linear Discriminant Analysis. \code{\link[MASS]{lda}} in MASS package
#'   \item slda: Sparse Linear Discriminant Analysis with LASSO penalty.
#'   \item nb: Naive Bayes. \code{\link[e1071]{naiveBayes}} in e1071 package
#'   \item nnb: Nonparametric Naive Bayes. \code{\link[naivebayes]{naive_bayes}} in naivebayes package
#'   \item ada: Ada-Boost. \code{\link[ada]{ada}} in ada package
#' }
#' @param p.adjust.methods multiple testing adjustment method. See \code{\link[stats]{p.adjust.methods}}
#' @param alpha a numeric scalar between 0 and 1 indicating the population type I error control
#' @param delta a numeric scalar between 0 and 1 indicating the violation rate. Default: 0.05
#' @param epsilon a numeric scalar between 0 and 1 indicating the significance level. Default: 0.05
#' @param N a positive integer indicating the maximum number of marginal features to be kept
#' @param l0 a numeric scalar between 0 and 1 indicating the proportion of leave-out class 0 data points. Default: 0.5
#' @param l1 a numeric scalar between 0 and 1 indicating the proportion of leave-out class 1 data points. Default: 0.5
#' @param seed random seed
#' @param ncores a positive integer that specifies the number of cores for computing. Default: number of cores - 1.
#' @param ... additional argument for base classification methods.
#' @return \code{npMarginal} returns a list with the following components:
#' \item{alpha}{user-specified population type I error control}
#' \item{delta}{user-specified violation rate}
#' \item{epsilon}{user-specified significance level}
#' \item{N}{user-specified maximum feature number}
#' \item{pval.unadj}{a vector of unadjusted p-values}
#' \item{pval.adj}{a vector of adjusted p-values}
#' \item{p.adjust.methods}{user-specified p-value adjusting method}
#' \item{feature}{features that pass marginal feature screening}
#'
#' @author Yiling Chen, \email{yiling0210@@ucla.edu}
#' @references \bold{FILL HERE}
#' @examples
#' ### example 1 ###
#' n = 1000; p = 1000
#' y = rbinom(n,size = 1, p =0.5)
#' x = matrix(NA, nrow =n,ncol =p)
#' x1 = rmvnorm(sum(y==1), mean = rep(1, p), sigma = diag(1,p))
#' x0 = rmvnorm(sum(y==0), mean = seq(from = -5, to = 4, length.out = p), sigma = diag(1,p))
#' x[y==1,] = x1
#' x[y==0,] = x0
#' table(y)
#' temp1 = npMarginal(x,y,method = "logistic",
#' N = 50,
#' p.adjust.methods = 'BH',
#' alpha =  0.05,
#' epsilon = 0.1,
#' l0 = 0.5,
#' l1 = 0.5)
#' plot(temp1$pval.unadj, ylim = c(0,1))
#' abline(h = epsilon)
#'
#' ### example 2 ###
#' n = 1000; p = 1000
#' y = rbinom(n,size = 1, p =0.5)
#' x = matrix(NA, nrow =n,ncol =p)
#' x1 = rmvnorm(sum(y==1), mean = rep(1, p), sigma = diag(1,p))
#' x0 = rmvnorm(sum(y==0), mean = seq(from = -2, to = 1, length.out = p), sigma = diag(1,p))
#' x[y==1,] = x1
#' x[y==0,] = x0
#' table(y)
#' temp2 = npMarginal(x,y,method = "logistic",
#' N = 50,
#' p.adjust.methods = 'BH',
#' alpha =  0.05,
#' epsilon = 0.1,
#' l0 = 0.5,
#' l1 = 0.5)
#' plot(temp2$pval.unadj, ylim = c(0,1))
#' abline(h = epsilon)
#'
#' @export
#' @importFrom stats glm p.adjust
#' @importFrom glmnet glmnet
#' @importFrom e1071 svm naiveBayes
#' @importFrom randomForest randomForest
#' @importFrom MASS lda
#' @importFrom naivebayes naive_bayes
#' @importFrom ada ada
#' @importFrom parallel mclapply detectCores
#' @importFrom mvtnorm rmvnorm
#'
#


npMarginal = function(x,y,
                      method = c("logistic", "penlog", "svm", "randomforest", "lda", "slda", "nb", "nnb", "ada", "tree"),
                      p.adjust.methods = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"),
                      N,
                      alpha,
                      delta = 0.05,
                      epsilon = 0.05,
                      l0 = 0.5,
                      l1 = 0.5,
                      seed = NULL,
                      ncores = detectCores() - 1,
                      ...){
  if(!is.null(seed)){
    set.seed(seed)
  }
  if(!is.matrix(x)){
    cat('x is not a matrix; converting it to a matrix')
    x = as.matrix(x, ncol=ncol(x))
  }

  if(!all(y %in% c(0,1))){
    stop('y contains value other than 0 or 1')
  }

  if(is.null(colnames(x))){
    genes = 1:ncol(x)
  }else{
    genes = colnames(x)
  }
  # set.seed(12345);method = 'logistic';


  method = match.arg(method, choices = c("logistic", "penlog", "svm", "randomforest", "lda", "slda", "nb", "nnb", "ada", "tree"))
  data = cbind.data.frame(y,x)
  n_idx = which(y == 1 )
  m_idx = which(y == 0 )
  d = ncol(x)
  leaveoutclass1 = sample(n_idx, l1 * length(n_idx))
  leaveoutclass0 = sample(m_idx, l0 * length(m_idx))
  # cat(head(leaveoutclass1))
  n2 = length(leaveoutclass1)
  m2 = length(leaveoutclass0)

  if(m2 < log(delta)/(log(1-alpha))){
    stop('the leave out class 0 sample is too small for the current alpha,epsilon and leave-out class 0 proportion.')
  }
  pval = mclapply(1:d, function(j){
    # fit = classificationScores(method,
    #                            train.x = x[-c(leaveoutclass1,leaveoutclass0),j, drop = F],
    #                            train.y = y[-c(leaveoutclass1,leaveoutclass0)],
    #                            test.x = x[c(leaveoutclass1,leaveoutclass0),j, drop = F],
    #                            test.y = y[c(leaveoutclass1,leaveoutclass0)])
    fit = classificationScores(method,
                               train.x = x[-c(leaveoutclass1,leaveoutclass0),j, drop = F],
                               train.y = y[-c(leaveoutclass1,leaveoutclass0)],
                               test.x = x[c(leaveoutclass1,leaveoutclass0),j, drop = F],
                               test.y = y[c(leaveoutclass1,leaveoutclass0)], ...)
    cutoff = estimatecutoff(score0 = fit$score0, alpha = alpha, delta = delta)

    pbinom(sum(fit$score1 <= cutoff), size = n2, prob = 1-alpha,
           lower.tail= T)
  }, mc.cores = ncores)

  pval = unlist(pval)
  adj_pval = p.adjust(pval, method =  p.adjust.methods, n = length(pval))
  candidate = genes[adj_pval <= epsilon]

  if(length(candidate) > N || (length(N) == 0 )){
    if(sum(adj_pval == min(adj_pval)) > N){
      warning('More than N features have the same smallest adjusted pvalues. Output all these features.')
      candidate = genes[which(adj_pval == min(adj_pval))]
    }else{
      candidate = genes[rank(log10(adj_pval), ties.method = 'max') <= N]
    }

  }

  return(list(
    alpha = alpha,
    delta = delta,
    epsilon = epsilon,
    N = N,
    pval.unadj = pval,
    pval.adj = adj_pval,
    p.adjust.methods = p.adjust.methods,
    feature = candidate
  ))

}


