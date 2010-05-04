somers2 = function (x, y, weights = NULL, normwt = FALSE, na.rm = TRUE)
{
    if (length(y) != length(x))
        stop("y must have same length as x")
    y <- as.integer(y)
    wtpres <- length(weights)
    if (wtpres && (wtpres != length(x)))
        stop("weights must have same length as x")
    if (na.rm) {
        miss <- if (wtpres)
            is.na(x + y + weights)
        else is.na(x + y)
        nmiss <- sum(miss)
        if (nmiss > 0) {
            miss <- !miss
            x <- x[miss]
            y <- y[miss]
            if (wtpres)
                weights <- weights[miss]
        }
    }
    else nmiss <- 0
    u <- sort(unique(y))
    if (any(y %nin% 0:1))
        stop("y must be binary")
    if (wtpres) {
        if (normwt)
            weights <- length(x) * weights/sum(weights)
        n <- sum(weights)
    }
    else n <- length(x)
    if (n < 2)
        stop("must have >=2 non-missing observations")
    n1 <- if (wtpres)
        sum(weights[y == 1])
    else sum(y == 1)
    if (n1 == 0 || n1 == n)
        return(c(C = NA, Dxy = NA, n = n, Missing = nmiss))
    mean.rank <- if (wtpres)
        wtd.mean(wtd.rank(x, weights, na.rm = FALSE), weights *
            y)
    else mean(rank(x)[y == 1])
    c.index <- (mean.rank - (n1 + 1)/2)/(n - n1)
    dxy <- 2 * (c.index - 0.5)
    r <- c(c.index, dxy, n, nmiss)
    names(r) <- c("C", "Dxy", "n", "Missing")
    r
}
