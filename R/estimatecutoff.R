estimatecutoff = function(score0, alpha, delta){
  m2  = length(score0)
  cutoff = sort(score0)[which((1-pbinom(1:m2-1, size = m2, prob = 1-alpha,
                                        lower.tail= T))<=delta)[1]]
  return(cutoff)
}
