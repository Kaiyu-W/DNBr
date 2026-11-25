#' Center and L2-normalize columns of a matrix
#'
#' @param x A numeric matrix.
#' @return A matrix with each column centered (mean = 0) and scaled to unit L2 norm.
#'   Columns with zero variance become \code{NaN}.
#'
center_scale <- function(x) {
    x <- t(t(x) - colMeans(x))
    t(t(x) / sqrt(colSums(x^2)))
}
#' Fast Pearson (Spearman) correlation via BLAS crossprod
#'
#' @param x A numeric matrix. Sparse matrices are converted to dense.
#' @param y An optional numeric matrix. If \code{NULL}, computes correlations
#'   among columns of \code{x}.
#' @param method Character string specifying the method: \code{"pearson"} (default)
#'   or \code{"spearman"}.
#' @param ... Currently unused.
#'
#' @return A Pearson correlation matrix. \code{NA} entries (e.g., zero-variance
#'   columns) are replaced with 1. Does not handle missing values.
#'

mycor <- function(x, y = NULL, method = c("pearson", "spearman"), ...) {
    method <- match.arg(method)
    if (inherits(x, "sparseMatrix")) x <- as.matrix(x)
    if (inherits(y, "sparseMatrix")) y <- as.matrix(y)

    if (method == "spearman") x <- apply(x, 2, rank, na.last = "keep")
    x <- center_scale(x)

    if (is.null(y)) {
        obj <- crossprod(x)
    } else {
        if (method == "spearman") y <- apply(y, 2, rank, na.last = "keep")
        y <- center_scale(y)
        obj <- crossprod(x, y)
    }

    obj[is.na(obj)] <- 1
    obj
}

#' Critical index
#'
#' @param n sample size
#' @param Sd standard deviation
#' @param pcc_in pcc inside the group
#' @param pcc_out pcc outside the group
#'
CI <- function(n, Sd, pcc_in, pcc_out) {
    ifelse(
        test = pcc_out == 0,
        yes = 0,
        no = sqrt(n) * Sd * pcc_in / pcc_out
    ) 
    # the original formula comes when n = 1
}

#' to print the process information
#'
#' @param ... parameters of function cat
#' @param quiet not execute function cat
#'
#'
mycat <- function(..., quiet = FALSE) {
    if (!quiet) {
        cat(...)

        # Ensure compatibility with Jupyter Notebook by forcing output flush
        tryCatch(
            flush.console(), 
            error = function(e) {}
        )
    }
}


#' Detect if the operating system is Windows
#'
#' @return Logical. TRUE if Windows, FALSE otherwise.
#' 
is_windows_os <- function() {
    tolower(Sys.info()["sysname"]) == "windows"
}


#' convert gene pair indices to flat vector index (for as.dist output)
#'
#' @param i integer vector representing the row indices of gene pairs
#' @param j integer vector representing the column indices of gene pairs
#' @param N the total number of genes (dimension of the original square matrix)
#'
#' @return An integer vector of positions in the flat correlation vector.
#'
idx2flat <- function(i, j, N) {
    r <- pmax(i, j)
    c <- pmin(i, j)
    as.integer((c - 1L) * N - (c * (c - 1L)) / 2L + r - c)
}
