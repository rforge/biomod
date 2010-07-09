Biomod.Manual<-function(manual="Biomod Manual"){

    paths <- .find.package("BIOMOD")
    paths <- paths[file_test("-d", file.path(paths, "doc"))]
    vignettes <- lapply(paths, function(dir) {
        tools::list_files_with_exts(file.path(dir, "doc"), exts=c("pdf"))
    })
    topic <- manual
    vignettes <- as.character(unlist(vignettes))
        vidx <- (tools::file_path_sans_ext(basename(vignettes)) ==
            topic)
        if (any(vidx)) {
            pdf <- sub("\\.[[:alpha:]]+$", ".pdf", vignettes)
            pidx <- file_test("-f", pdf)
            ok <- vidx & pidx
            if (any(ok)) {
                idx <- min(which(ok))
                if (sum(ok) > 1) {
                  warning(gettextf("vignette '%s' found more than once,\nusing the one found in '%s'",
                    topic, dirname(pdf[idx])), call. = FALSE,
                    domain = NA)
                }
                z <- list(file = vignettes[idx], pdf = pdf[idx])
            }
            else {
                z <- list(file = vignettes[vidx][1L], pdf = character(0L))
            }
            z$topic <- topic
            class(z) <- "vignette"
            return(z)
        }
        else warning(gettextf("vignette '%s' *not* found", topic),
            call. = FALSE, domain = NA)
}
