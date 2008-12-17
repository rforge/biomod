`CVnnet` <-
function(formula, data, truth, size=c(round(0.75*Biomod.material[["NbVar"]]), 
round(0.75*Biomod.material[["NbVar"]]), round(0.75*Biomod.material[["NbVar"]]), Biomod.material[["NbVar"]], Biomod.material[["NbVar"]], Biomod.material[["NbVar"]]), lambda=rep(
    c(0.01, 0.05, 0.1), 2), nreps=1, nifold=3, ...)
{
    choice <- numeric(length(lambda))
    ri <- sample(nifold, nrow(data), replace=T)
    for(j in seq(along=lambda)) {
        choice[j] <- CVnnet1(formula, data, truth, nreps=nreps,ri=ri, size=size[j], decay=lambda[j], ...)
    }
    cbind(size=size, decay=lambda, fit=choice)
}

