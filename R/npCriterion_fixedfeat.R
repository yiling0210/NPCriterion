npCriterion_fixedfeat = function(x,y,
                                 n_idx,
                                 m_idx,
                                 method ,
                                 alpha,
                                 delta = 0.05,
                                 B = 5,
                                 l0 = 0.5,
                                 l1 = 0.5,
                                 seed = NULL,
                                 ncores,
                                 ...){

  re = mclapply(1:B, function(b){
    # cat('b=',b)
    leaveoutclass1 = sample(n_idx, l1 * length(n_idx))
    leaveoutclass0 = sample(m_idx, l0 * length(m_idx))
    # fit = classificationScores(method,
    #                            train.x = x[-c(leaveoutclass1,leaveoutclass0),],
    #                            train.y = y[-c(leaveoutclass1,leaveoutclass0)],
    #                            test.x = x[c(leaveoutclass1,leaveoutclass0),],
    #                            test.y = y[c(leaveoutclass1,leaveoutclass0)])
    fit = classificationScores(method,
                               train.x = x[-c(leaveoutclass1,leaveoutclass0),],
                               train.y = y[-c(leaveoutclass1,leaveoutclass0)],
                               test.x = x[c(leaveoutclass1,leaveoutclass0),],
                               test.y = y[c(leaveoutclass1,leaveoutclass0)], ...)
    ## under NP
    cutoff = estimatecutoff(score0 = fit$score0, alpha = alpha, delta = delta)
    # cat('cutoff = ', cutoff)
    type2 = estimatetype2(score1 = fit$score1, cutoff = cutoff)
    ### under classical
    pred_cl = as.numeric(c( fit$score1, fit$score0) > 0.5)
    err = mean(pred_cl != y[c(leaveoutclass1,leaveoutclass0)])


    return(list(type2= type2, err = err))
  }, mc.cores = ncores)

  NPC = unlist(re)[2*(1:B) - 1]
  cl = unlist(re)[2*(1:B)]


  list( NPC = mean(NPC), npc.se = sd(NPC), test_err = mean(cl), err.se = sd(cl))

}

