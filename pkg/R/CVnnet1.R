`CVnnet1` <-
function(formula, data, truth, nreps=1, ri, ...)
{
    res <- numeric(length(truth))
    for(i in sort(unique(ri))) {
        for(rep in 1:nreps) {
            learn <- nnet(formula, data[ri != i, ], trace=F, ...)
            res[ri == i] <- res[ri == i] + predict(learn, data[ri == i,  ])
        }
    }
    somers2(res/nreps, truth)["C"]
}

