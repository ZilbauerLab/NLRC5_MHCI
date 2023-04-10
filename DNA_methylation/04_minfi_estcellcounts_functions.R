# modified from https://github.com/hansenlab/minfi/blob/master/R/estimateCellCounts.R
pickCompProbes <- function(beta, meta,cellTypes, numProbes = 50) {

    splitit <- function(x) {
        split(seq_along(x), x)}

    p <- as.matrix(beta)
    pd <- as.data.frame(meta)
    keep <- which(pd$CellType %in% cellTypes)
    pd <- pd[keep,]
    p <- p[,keep]

    # NOTE: Make cell type a factor
    pd$CellType <- factor(pd$CellType, levels = cellTypes)
    ffComp <- rowFtests(p, pd$CellType)
    prof <- vapply(
        X = splitit(pd$CellType),
        FUN = function(j) rowMeans2(p, cols = j),
        FUN.VALUE = numeric(nrow(p)))
    r <- rowRanges(p)
    compTable <- cbind(ffComp, prof, r, abs(r[, 1] - r[, 2]))
    names(compTable)[1] <- "Fstat"
    names(compTable)[c(-2, -1, 0) + ncol(compTable)] <-
        c("low", "high", "range")
    tIndexes <- splitit(pd$CellType)
    tstatList <- lapply(tIndexes, function(i) {
        x <- rep(0,ncol(p))
        x[i] <- 1
        return(rowttests(p, factor(x)))
    })

    probeList <- lapply(tstatList, function(x) {
        y <- x[x[, "p.value"] < 1e-8, ]
        yUp <- y[order(y[, "dm"], decreasing = TRUE), ]
        yDown <- y[order(y[, "dm"], decreasing = FALSE), ]
        c(rownames(yUp)[seq_len(numProbes)],
          rownames(yDown)[seq_len(numProbes)])})


    trainingProbes <- unique(unlist(probeList))
    p <- p[trainingProbes,]

    pMeans <- colMeans2(p)
    names(pMeans) <- pd$CellType

    form <- as.formula(
        sprintf("y ~ %s - 1", paste(levels(pd$CellType), collapse = "+")))
    phenoDF <- as.data.frame(model.matrix(~ pd$CellType - 1))
    colnames(phenoDF) <- sub("^pd\\$CellType", "", colnames(phenoDF))
    if (ncol(phenoDF) == 2) {
        # Two group solution
        X <- as.matrix(phenoDF)
        coefEsts <- t(solve(t(X) %*% X) %*% t(X) %*% t(p))
    } else {
        # > 2 groups solution
        tmp <- validationCellType(Y = p, pheno = phenoDF, modelFix = form)
        coefEsts <- tmp$coefEsts
    }

    list(
        coefEsts = coefEsts,
        compTable = compTable,
        sampleMeans = pMeans)
}



projectCellType <- function(Y, coefCellType, contrastCellType = NULL,
                            nonnegative = TRUE, lessThanOne = FALSE) {
    if (is.null(contrastCellType)) {
        Xmat <- coefCellType
    } else {
        Xmat <- tcrossprod(coefCellType, contrastCellType)
    }

    nCol <- dim(Xmat)[2]
    if (nCol == 2) {
        Dmat <- crossprod(Xmat)
        mixCoef <- t(
            apply(Y, 2, function(x) solve(Dmat, crossprod(Xmat, x))))
        colnames(mixCoef) <- colnames(Xmat)
        return(mixCoef)
    } else {
        nSubj <- dim(Y)[2]

        mixCoef <- matrix(0, nSubj, nCol)
        rownames(mixCoef) <- colnames(Y)
        colnames(mixCoef) <- colnames(Xmat)

        if (nonnegative) {
            if (lessThanOne) {
                Amat <- cbind(rep(-1, nCol), diag(nCol))
                b0vec <- c(-1, rep(0, nCol))
            } else {
                Amat <- diag(nCol)
                b0vec <- rep(0, nCol)
            }
            for (i in seq_len(nSubj)) {
                obs <- which(!is.na(Y[,i]))
                Dmat <- crossprod(Xmat[obs,])
                mixCoef[i,] <- solve.QP(
                    Dmat = Dmat,
                    dvec = crossprod(Xmat[obs,], Y[obs,i]),
                    Amat = Amat,
                    bvec = b0vec)$sol
            }
        } else {
            for (i in seq_len(nSubj)) {
                obs <- which(!is.na(Y[,i]))
                Dmat <- crossprod(Xmat[obs,])
                mixCoef[i,] <- solve(Dmat, t(Xmat[obs,]) %*% Y[obs,i])
            }
        }
        mixCoef
    }
}


projectCellTypeConstrainAlways <- function(Y, coefCellType, contrastCellType = NULL,
                            nonnegative = TRUE, lessThanOne = FALSE) {
    if (is.null(contrastCellType)) {
        Xmat <- coefCellType
    } else {
        Xmat <- tcrossprod(coefCellType, contrastCellType)
    }

    nCol <- dim(Xmat)[2]
    # if (nCol == 2) {
    #     Dmat <- crossprod(Xmat)
    #     mixCoef <- t(
    #         apply(Y, 2, function(x) solve(Dmat, crossprod(Xmat, x))))
    #     colnames(mixCoef) <- colnames(Xmat)
    #     return(mixCoef)
    # } else {
    nSubj <- dim(Y)[2]

    mixCoef <- matrix(0, nSubj, nCol)
    rownames(mixCoef) <- colnames(Y)
    colnames(mixCoef) <- colnames(Xmat)

    if (nonnegative) {
        if (lessThanOne) {
            Amat <- cbind(rep(-1, nCol), diag(nCol))
            b0vec <- c(-1, rep(0, nCol))
        } else {
            Amat <- diag(nCol)
            b0vec <- rep(0, nCol)
        }
        for (i in seq_len(nSubj)) {
            obs <- which(!is.na(Y[,i]))
            Dmat <- crossprod(Xmat[obs,])
            mixCoef[i,] <- solve.QP(
                Dmat = Dmat,
                dvec = crossprod(Xmat[obs,], Y[obs,i]),
                Amat = Amat,
                bvec = b0vec)$sol
        }
    } else {
        for (i in seq_len(nSubj)) {
            obs <- which(!is.na(Y[,i]))
            Dmat <- crossprod(Xmat[obs,])
            mixCoef[i,] <- solve(Dmat, t(Xmat[obs,]) %*% Y[obs,i])
        }
    }
    mixCoef
    # }
}

validationCellType <- function(Y, pheno, modelFix, modelBatch=NULL,
                               L.forFstat = NULL, verbose = FALSE){
    N <- dim(pheno)[1]
    pheno$y <- rep(0, N)
    xTest <- model.matrix(modelFix, pheno)
    sizeModel <- dim(xTest)[2]
    M <- dim(Y)[1]

    if (is.null(L.forFstat)) {
        # NOTE: All non-intercept coefficients
        L.forFstat <- diag(sizeModel)[-1,]
        colnames(L.forFstat) <- colnames(xTest)
        rownames(L.forFstat) <- colnames(xTest)[-1]
    }

    # Initialize various containers
    sigmaResid <- sigmaIcept <- nObserved <- nClusters <- Fstat <- rep(NA, M)
    coefEsts <- matrix(NA, M, sizeModel)
    coefVcovs <- list()

    if (verbose) cat("[validationCellType] ")
    # Loop over each CpG
    for (j in seq_len(M)) {
        # Remove missing methylation values
        ii <- !is.na(Y[j, ])
        nObserved[j] <- sum(ii)
        pheno$y <- Y[j,]

        if (j %% round(M / 10) == 0 && verbose) cat(".") # Report progress

        # Try to fit a mixed model to adjust for plate
        try({
            if (!is.null(modelBatch)) {
                fit <- try(
                    lme(modelFix, random = modelBatch, data = pheno[ii, ]))
                # NOTE: If LME can't be fit, just use OLS
                OLS <- inherits(fit, "try-error")
            } else {
                OLS <- TRUE
            }

            if (OLS) {
                fit <- lm(modelFix, data = pheno[ii, ])
                fitCoef <- fit$coef
                sigmaResid[j] <- summary(fit)$sigma
                sigmaIcept[j] <- 0
                nClusters[j] <- 0
            } else {
                fitCoef <- fit$coef$fixed
                sigmaResid[j] <- fit$sigma
                sigmaIcept[j] <- sqrt(getVarCov(fit)[1])
                nClusters[j] <- length(fit$coef$random[[1]])
            }
            coefEsts[j,] <- fitCoef
            coefVcovs[[j]] <- vcov(fit)

            useCoef <- L.forFstat %*% fitCoef
            useV <- L.forFstat %*% coefVcovs[[j]] %*% t(L.forFstat)
            Fstat[j] <- (t(useCoef) %*% solve(useV, useCoef)) / sizeModel
        })
    }
    if (verbose) cat(" done\n")

    # Name the rows so that they can be easily matched to the target data set
    rownames(coefEsts) <- rownames(Y)
    colnames(coefEsts) <- names(fitCoef)
    degFree <- nObserved - nClusters - sizeModel + 1

    # Get P values corresponding to F statistics
    Pval <- 1 - pf(Fstat, sizeModel, degFree)

    list(
        coefEsts = coefEsts,
        coefVcovs = coefVcovs,
        modelFix = modelFix,
        modelBatch = modelBatch,
        sigmaIcept = sigmaIcept,
        sigmaResid = sigmaResid,
        L.forFstat = L.forFstat,
        Pval = Pval,
        orderFstat = order(-Fstat),
        Fstat = Fstat,
        nClusters = nClusters,
        nObserved = nObserved,
        degFree = degFree)
}