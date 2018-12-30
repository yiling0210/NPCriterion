#' Feature selection based on Neyman-Pearson Criterion (NPC)
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
#' @param enumeration a feature set generation method, which can either be 'forward', 'backward' or 'exhaustive'. Default: 'forward'
#' @param max_feature_size an optional integer when enumeration is 'exhaustive'. When not supplied, set to be the total number of features
#' @param alpha a numeric scalar between 0 and 1 indicating the population type I error control
#' @param delta a numeric scalar between 0 and 1 indicating the violation rate. Default: 0.05
#' @param B a positive integer indicating the number of random splits. Default: 5
#' @param l0 a numeric scalar between 0 and 1 indicating the proportion of leave-out class 0 data points. Default: 0.5
#' @param l1 a numeric scalar between 0 and 1 indicating the proportion of leave-out class 1 data points. Default: 0.5
#' @param seed random seed
#' @param ncores a positive integer that specifies the number of cores for computing. Default: number of cores - 1.
#' @param ... additional argument for base classification methods.

#' @return \code{npCriterion} returns a list with the following components:
#'
#' \item{method}{the base classification method}
#' \item{alpha}{user-specified alpha value}
#' \item{delta}{user-specified delta value}
#' \item{B}{total number of random splits}
#' \item{l0}{the proportion of leave-out class 0 data points}
#' \item{l1}{the proportion of leave-out class 1 data points}
#' \item{featuresets_examined}{when 'enumeration' = 'forward', a list of size 2 whose first component is enumeration, and the second component is a vector of the features that are sequentially included;
#' when 'enumeration' = 'backward', a list of size 2 whose first component is enumeration, and the second component is a vector of the features that are sequentially excluded;
#' when 'enumeration' = 'exhaustive', a list of size 2 whose first component is enumeration, and the second component is a list of matrices whose column number ranges from 1 to max_feature_size. Rows of such a matrix represent a feature set.
#' }
#' \item{npc}{when 'enumeration' = 'forward' or 'backward', a vector of NPC values of feature sets in featuresets_examined;
#' when 'enumeration' = 'exhaustive', a list of vectors of NPC values computed on the feature sets in featuresets_examined
#' }
#' \item{npc.sd}{when 'enumeration' = 'forward' or 'backward', a vector of standard deviations of empirical type II errors of feature sets in featuresets_examined;
#' when 'enumeration' = 'exhaustive', a list of standard deviations of empirical type II errors computed on the feature sets in featuresets_examined
#' }
#' \item{npc.se}{when 'enumeration' = 'forward' or 'backward', a vector of standard errors of empirical type II errors of feature sets in featuresets_examined;
#' when 'enumeration' = 'exhaustive', a list of standard errors of empirical type II errors computed on the feature sets in featuresets_examined
#' }
#' \item{err}{when 'enumeration' = 'forward' or 'backward', a vector of CV errors of feature sets in featuresets_examined;
#' when 'enumeration' = 'exhaustive', a list of vectors of CV errors computed on the feature sets in featuresets_examined
#' }
#' \item{err.se}{when 'enumeration' = 'forward' or 'backward', a vector of standard deviations of test errors of feature sets in featuresets_examined;
#' when 'enumeration' = 'exhaustive', a list of standard deviations of test errors computed on the feature sets in featuresets_examined
#' }
#' \item{err.se}{when 'enumeration' = 'forward' or 'backward', a vector of standard errors of test errors of feature sets in featuresets_examined;
#' when 'enumeration' = 'exhaustive', a list of standard errors of test errors computed on the feature sets in featuresets_examined
#' }
#' \item{features_minNPC}{a feature set with the minimal NPC value and its corresponding NPC statistics and test errors.}
#'

#' @author Yiling Chen, \email{yiling0210@@ucla.edu}
#' @references \bold{FILL HERE}
#' @examples
#' ### Example1 #####
#' x = matrix(rnorm(20000), ncol =5)
#' y = rbinom(x%*%1:5,size = 1, p =0.5)
#' table(y)
#' temp1 = npCriterion(x,y,method = "logistic",
#' enumeration = 'forward',
#' max_feature_size = NULL,
#' alpha =  0.05,
#' delta = 0.05,
#' B = 5,
#' l0 = 0.5,
#' l1 = 0.5)
#'
#' temp2 = npCriterion(x,y,method ="svm",
#' kernel = 'radial',
#' enumeration = 'exhaustive',
#' alpha =  0.05,
#' delta = 0.05,
#' B = 5,
#' l0 = 0.5,
#' l1 = 0.5)
#' ### Example2 #####
#' y = rbinom(100000,size = 1, p =0.5)
#' x = matrix(NA, nrow =100000,ncol =2)
#' x1 = cbind(rnorm(sum(y==1),mean =1, sd =1),rnorm(sum(y==1),mean =1, sd =1))
#' x0 = cbind(rnorm(sum(y==0),mean =-1, sd =1),rnorm(sum(y==0),mean =0.5, sd =1.5))
#' pnorm(qnorm(0.95,-1,1),1,1)
#' pnorm(qnorm(0.95,0.5,1.5),1,1)
#' x[y==1,] = x1
#' x[y==0,] = x0
#' table(y)
#'
#' temp3 = npCriterion(x,y,method ="lda",
#'                     enumeration = 'exhaustive',
#'                     alpha =  0.05,
#'                     delta = 0.05,
#'                     B = 5,
#'                     l0 = 0.5,
#'                     l1 = 0.5)
#' temp3$criteria$`ell=1`
#'
#'
#' @export
#' @importFrom stats glm
#' @importFrom glmnet glmnet
#' @importFrom e1071 svm naiveBayes
#' @importFrom randomForest randomForest
#' @importFrom MASS lda
#' @importFrom naivebayes naive_bayes
#' @importFrom ada ada
#' @importFrom parallel mclapply detectCores
#'



npCriterion = function(x,y,method = c("logistic", "penlog", "svm",
                                      "randomforest", "lda", "slda", "nb", "nnb", "ada", "tree"),
                       enumeration = NULL,
                       max_feature_size = NULL,
                       alpha,
                       delta = 0.05,
                       B = 5,
                       l0 = 0.5,
                       l1 = 0.5,
                       seed = NULL,
                       ncores  = detectCores() -1,
                       ... ){
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

  n_idx = which(y == 1 )
  m_idx = which(y == 0)
  m_2 = floor(length(m_idx)*l0)
  if(m_2 < log(delta)/(log(1-alpha))){
    stop('the leave out class 0 sample is too small for the current alpha,delta and leave-out class 0 proportion.')
  }


  p = ncol(x)

  if(is.null(enumeration)){
    enumeration = 'forward'
  }

  enumeration = match.arg(enumeration, choices = c('forward', 'backward','exhaustive'))
  method = match.arg(method, choices = c("logistic", "penlog", "svm",
                                         "randomforest", "lda", "slda", "nb", "nnb", "ada", "tree"))

  obj = NULL
  obj$method = method
  obj$alpha = alpha
  obj$delta =delta
  obj$B = B
  obj$l0 = l0
  obj$l1 = l1


  if(enumeration == 'forward'){
    # forwd_it = 1
    selected_idx = NULL
    selected_npc = NULL
    unselected_idx = setdiff(1:p, selected_idx)
    while(length(unselected_idx)>0){


      stats_ls = mclapply(unselected_idx, function(j){
        stats_fixedfeat = npCriterion_fixedfeat(x[, c(selected_idx, j), drop = FALSE],
                                                y,
                                                n_idx = n_idx,
                                                m_idx = m_idx,
                                                method = method,
                                                alpha = alpha,
                                                delta = delta,
                                                B = B,
                                                l0 = l0,
                                                l1 = l1,
                                                seed = seed,
                                                ncores = 1,
                                                ...
        )
        # stats_fixedfeat = npCriterion_fixedfeat(x[, c(selected_idx, j), drop = FALSE],
        #                                         y,
        #                                         n_idx = n_idx,
        #                                         m_idx = m_idx,
        #                                         method = method,
        #                                         alpha = alpha,
        #                                         delta = delta,
        #                                         B = B,
        #                                         l0 = l0,
        #                                         l1 = l1,
        #                                         seed = seed,
        #                                         ncores = 1
        # )
        stats_fixedfeat
      },mc.cores = ncores)


      candidate_npc = sapply(stats_ls, function(x){x$NPC})
      selected_idx = c( selected_idx, unselected_idx[which.min(candidate_npc)])
      selected_npc = c(selected_npc, stats_ls[which.min(candidate_npc)])
      unselected_idx = setdiff(1:p, selected_idx)

    }

    ###### organize output from forward selection
    obj$featuresets_examined = list(enumeration = 'forward', sequential_features = selected_idx)

    criteria = matrix(unlist(selected_npc), byrow = T, ncol = length(selected_npc[[1]]))
    colnames(criteria) = names(selected_npc[[1]])
    stats_minNPC = criteria[which.min(criteria[,1]),]
    criteria = split(t(criteria), 1:ncol(criteria))
    names(criteria) = names(selected_npc[[1]])
    obj = c(obj, criteria)


    obj$features_minNPC = list(feature_idx = selected_idx[1:which.min(criteria$NPC)],
                               criteria = stats_minNPC)

    return(obj)
  }

  if(enumeration == 'backward'){
    selected_idx = 1:p
    selected_npc = NULL
    unselected_idx = setdiff(1:p, selected_idx)
    while(length(selected_idx)>1){


      stats_ls = mclapply(selected_idx, function(j){
        stats_fixedfeat = npCriterion_fixedfeat(x[, c(selected_idx, j), drop = FALSE],
                                                y,
                                                n_idx = n_idx,
                                                m_idx = m_idx,
                                                method = method,
                                                alpha = alpha,
                                                delta = delta,
                                                B = B,
                                                l0 = l0,
                                                l1 = l1,
                                                seed = seed,
                                                ncores = 1,
                                                ...
        )
        # stats_fixedfeat = npCriterion_fixedfeat(x[, setdiff(selected_idx, j), drop = FALSE],
        #                                         y,
        #                                         n_idx = n_idx,
        #                                         m_idx = m_idx,
        #                                         method = method,
        #                                         alpha = alpha,
        #                                         delta = delta,
        #                                         B = B,
        #                                         l0 = l0,
        #                                         l1 = l1,
        #                                         seed = seed,
        #                                         ncores = 1
        # )
        stats_fixedfeat
      },mc.cores = ncores)


      candidate_npc = sapply(stats_ls, function(x){x$NPC})
      unselected_idx = c(unselected_idx, selected_idx[which.min(candidate_npc)] )
      selected_idx = setdiff(selected_idx, selected_idx[which.min(candidate_npc)])

      selected_npc = c(selected_npc, stats_ls[which.min(candidate_npc)])
      # unselected_idx = setdiff(1:p, selected_idx)

    }

    obj$featuresets_examined = list(enumeration = enumeration, sequential_features = unselected_idx)

    criteria = matrix(unlist(selected_npc), byrow = T, ncol = length(selected_npc[[1]]))
    colnames(criteria) = names(selected_npc[[1]])
    stats_minNPC = criteria[which.min(criteria[,1]),]
    criteria = split(t(criteria), 1:ncol(criteria))
    names(criteria) = names(selected_npc[[1]])
    obj = c(obj, criteria)

    obj$features_minNPC= list(feature_idx = setdiff(1:p,unselected_idx[1:which.min(criteria$NPC)]),
                              stats_minNPC = stats)


    return(obj)


  }


  if(enumeration == 'exhaustive') {
    if (is.null(max_feature_size)) {
      max_feature_size = p
    }
    if (max_feature_size >= 30) {
      warning('Heavy computation due to the number of subsets needed to examine!')
    }

    stats_ls = mclapply(1:max_feature_size, function(ell) {
      subset_idx = t(combn(1:p, ell))

      stats_fixedell = lapply(1:nrow(subset_idx), function(i) {
        npCriterion_fixedfeat(x[,subset_idx[i,],drop = FALSE],
                              y,
                              n_idx = n_idx,
                              m_idx = m_idx,
                              method = method,
                              alpha = alpha,
                              delta = delta,
                              B = B,
                              l0 = l0,
                              l1 = l1,
                              seed = seed,
                              ncores = 1,
                              ...)
        # npCriterion_fixedfeat(
        #   x[, subset_idx[i, ], drop = FALSE],
        #   y,
        #   n_idx = n_idx,
        #   m_idx = m_idx,
        #   method = method,
        #   alpha = alpha,
        #   delta = delta,
        #   B = B,
        #   l0 = l0,
        #   l1 = l1,
        #   seed = seed,
        #   ncores = 1
        # )
      })

      criteria = matrix(unlist(stats_fixedell),
                        byrow = T,
                        ncol = length(stats_fixedell[[1]]))
      re = list(subset_idx = subset_idx, npc = criteria[,1],npc.sd = criteria[,2],
                npc.se = criteria[,3], err = criteria[,4], err.sd = criteria[,5], err.se = criteria[,6])

      return(re)
    }, mc.cores = ncores)

    names(stats_ls) = paste0('ell=', 1:max_feature_size)
    obj$featuresets_examined = list(enumeration = enumeration, lapply(stats_ls,function(x){x$subset_idx}))
    obj$npc = lapply(stats_ls,function(x){x$npc})
    obj$npc.sd = lapply(stats_ls,function(x){x$npc.sd})
    obj$npc.se = lapply(stats_ls,function(x){x$npc.se})
    obj$err = lapply(stats_ls,function(x){x$err})
    obj$err.sd = lapply(stats_ls,function(x){x$err.sd})
    obj$err.se = lapply(stats_ls,function(x){x$err.se})


    npc  = unlist(obj$npc)
    minNPC_ell_idx = which.min(cumsum(choose(p, 1:max_feature_size)) < which.min(unlist(npc)))
    if(minNPC_ell_idx == 1){
      minNPC_feature = stats_ls[[minNPC_ell_idx]]$subset_idx[which.min(unlist(npc))]
    }else{
      minNPC_feature = stats_ls[[minNPC_ell_idx]]$subset_idx[which.min(unlist(npc)) -  cumsum(choose(p, 1:max_feature_size))[minNPC_ell_idx - 1], ]
    }
    criteria = matrix(unlist(stats_ls[[minNPC_ell_idx]][-1]), byrow = F, nrow = nrow( stats_ls[[minNPC_ell_idx]]$subset_idx))
    minNPC_stats = criteria[which.min(criteria[,1]),]

    obj$features_minNPC = list(feature_idx = minNPC_feature, stats_minNPC = minNPC_stats)




  }


  return(obj)


}






# mysam = function(x,y){
#   sam.out = siggenes::sam(data= t(x), cl= y)
#   rank = list.siggenes(sam.out,delta = 0.01, gene.names = 1:500, order=T)
#   return(rank)
# }
