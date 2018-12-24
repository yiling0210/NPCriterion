Yiling
# NPCriterion
Model comparison criterion under the Neyman Pearson binary classification paradigm

### Latest News

### Introduction

NPCriterion is proposed to compare models under the Neyman Pearson(NP) binary classification paradigm. Users can refer to our paper [Neyman-Pearson Criterion (NPC): A Model Selection Criterion for Asymmetric Binary Classification]() for a detailed description of the criterion and applications.

Any suggestions on the package are welcome! For technical problems, please report to Issues. For suggestions and comments on the method, please contact Yiling([yiling0210@ucla.edu](yiling0210@ucla.edu)) or Dr. Jessica Li ([jli@stat.ucla.edu](jli@stat.ucla.edu)).

### Installation
The package is not on CRAN yet. For installation please use the following codes in R
```R
install.packages("devtools")
library(devtools)

install_github("yiling0210/NPCriterion")
```


### Quick start
NPCriterion is built to integrate many common classification methods, including logistic regression, penalized logistic regression, svm, random forest, linear discriminant analysis (lds), sparse lda, naive Bayes, nonparametric Naive Bayes and adaBoost. 
Its mandatory inputs include the covariates matrix x, response y, classification method method and the type I error control alpha:
```R
npCriterion(x = x,  # covariate matrix 
            y = y,  # binary vector response
            method = "logistic", # classification method
            enumeration = "forward", # forward
            alpha = 0.05, # type I error control
            delta = 0.05, # violation rate
            B = 5, # number of random splits
            l0 = 0.5, # proportion of leave-out class 0 data
            l1 = 0.5, # proportion of leave-out class 1 data
            ncores = detectCores() - 1 # number of cores for parallel computing
            )
```

npCriterion returns a list, including features that have been examined featuresets_examined and	a feature set with the minimal NPC value features_minNPC	

For detailed usage, please refer to the package [manual]{} or vignette.
