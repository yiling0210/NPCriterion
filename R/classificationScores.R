# method
# train.x = x[-c(leaveoutclass1,leaveoutclass0),];
# train.y = y[-c(leaveoutclass1,leaveoutclass0)];
# test.x = x[c(leaveoutclass1,leaveoutclass0),];
# test.y = y[c(leaveoutclass1,leaveoutclass0)]


classificationScores = function(method, train.x, train.y, test.x, test.y, ...){


  if (method == "logistic") {
    train.data = data.frame(x = train.x, y = train.y)
    fit = glm(y ~ ., data = train.data, family = "binomial")
    # colnames(fit$data)
    # # if(colnames(fit$data))
    # colnames(data.frame(x = test.x))
    test.score = predict(fit, data.frame(x = test.x), type = "response")
  }
  else if (method == "penlog") {
    fit = cv.glmnet(train.x, train.y, family = "binomial",
                    ...)
    test.score = predict(fit$glmnet.fit, newx = test.x, type = "response",
                         s = fit$lambda.min)
    test.score = as.vector(test.score)
  }
  else if (method == "svm") {
    train.y = as.factor(train.y)
    fit = svm(train.x, train.y, ...)
    test.score = attr(predict(fit, test.x, decision.values = TRUE),
                      "decision.values")[, 1]
  }
  else if (method == "randomforest") {
    if(is.null(dim(train.x))){
      train.x = as.matrix(train.x)
      test.x = as.matrix(test.x)
    }
    train.y = as.factor(train.y)
    # set.seed(12345)
    fit = randomForest(train.x, train.y, ...)
    # fit = randomForest(train.x, train.y)
    test.score = predict(fit, test.x, type = "prob")[, 2]
    # test.score[1:10]
  }
  else if (method == "lda") {
    if(is.null(dim(train.x))){
      train.x = as.matrix(train.x)
      test.x = as.matrix(test.x)
    }
    fit = lda(train.x, train.y)
    test.score = predict(fit, newdata = test.x)$posterior[, 2]
  }

  else if (method == "slda") {
    n1 = sum(train.y == 1)
    n0 = sum(train.y == 0)
    n = n1 + n0
    lda.y = train.y
    lda.y[train.y == 0] = -n/n0
    lda.y[train.y == 1] = n/n1
    fit = cv.glmnet(train.x, lda.y, ...)
    test.score = predict(fit$glmnet.fit, newx = test.x, type = "link",
                         s = fit$lambda.min)
    test.score = as.vector(test.score)
  }
  else if (method == "nb") {
    train.data = data.frame(train.x, y = train.y)
    fit <- naive_bayes(as.factor(y) ~ ., data = train.data,
                       usekernel = FALSE)
    test.score = predict(fit, data.frame(test.x), type = "prob")[,
                                                                 2]
  }
  else if (method == "nnb") {
    train.data = data.frame(train.x, y = train.y)
    fit <- naive_bayes(as.factor(y) ~ ., data = train.data,
                       usekernel = TRUE)
    test.score = predict(fit, data.frame(test.x), type = "prob")[,
                                                                 2]
  }
  else if (method == "ada") {
    train.data = data.frame(train.x, y = train.y)
    fit = ada(y ~ ., data = train.data)
    test.score = predict(fit, data.frame(test.x), type = "probs")[,
                                                                  2]
  }
  else if (method == "tree") {
    train.y = as.factor(train.y)
    train.data = data.frame(train.x, y = train.y)
    fit = tree(y ~ ., data = train.data)
    test.score = predict(fit,
                         newdata = data.frame(test.x),
                         type = "vector")
  }

  score0 = (test.score[which(test.y == 0)])
  score1 = (test.score[which(test.y == 1)])
  return(list( score0 = score0 , score1 = score1))
}
