#' to ignore the warning and error of correlation computation
#'
#' @param ... parameters of function cor
#'
#' @export
#'
mycor <- function(...) {
    obj <- suppressWarnings(stats::cor(...))
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


