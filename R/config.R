#' to ignore the warning and error of correlation computation
#'
#' @param x x
#' @param y y
#' @param ... other parameters of function cor
#'
#' @export
#'
mycor <- function(x, y = NULL, ...) {
    if (is(x, "dgCMatrix"))
        x <- as.matrix(x)
    if (is(y, "dgCMatrix"))
        y <- as.matrix(y)
    
    obj <- suppressWarnings(stats::cor(x = x, y = y, ...))
    obj[is.na(obj)] <- 1
    obj
}

#' to print the process information
#'
#' @param ... parameters of function cat
#' @param quiet not execute function cat
#'
#' @export
#'
mycat <- function(..., quiet = F) {
    if (!quiet) cat(...)
}


