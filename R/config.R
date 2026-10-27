#' to ignore the warning and error of correlation computation
#'
#' @param x x
#' @param y y
#' @param ... other parameters of function cor
#'
#'
mycor <- function(x, y = NULL, ...) {
    if (inherits(x, "sparseMatrix"))
        x <- as.matrix(x)
    if (inherits(y, "sparseMatrix"))
        y <- as.matrix(y)
    
    obj <- suppressWarnings(stats::cor(x = x, y = y, ...))
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