#' define S4 object of DNB_module
#'
#' @slot MODULE list.
#' @slot bestMODULE logical.
#'
#' @return S4 object, DNB_module
#'
mynew_module <- setClass(
    Class = "DNB_module",
    slots = list(
        MODULE = "list",
        bestMODULE = "logical"
    ),
    sealed = TRUE
)

#' define S4 object of DNB_pre, includes the numeric results of each DNB-computation and S4-DNB_module
#'
#' @slot MEAN numeric.
#' @slot SD numeric.
#' @slot CV numeric.
#' @slot MODULEs DNB_module.
#' @slot PCC_IN numeric.
#' @slot PCC_OUT numeric.
#' @slot SCORE numeric.
#' @slot rank numeric.
#' @slot resource character.
#' @slot rank_all numeric.
#'
#' @return S4 object, DNB_res
#'
mynew_res <- setClass(
    Class = "DNB_res",
    slots = list(
        MEAN = "numeric",
        SD = "numeric",
        CV = "numeric",
        MODULEs = "DNB_module",
        PCC_IN = "numeric",
        PCC_OUT = "numeric",
        SCORE = "numeric",
        rank = "numeric",
        resource = "character",
        rank_all = "numeric"
    ),
    sealed = TRUE
)

#' define S4 object of DNB_obj, includes the raw data(matrix or data.frame), the pre-results (S4-DNB_res) and the results of filtered (ntop max) modules (S4-DNB_res)
#'
#' @slot data ANY.
#' @slot used_gene ANY.
#' @slot pre_result DNB_res.
#' @slot result DNB_res.
#'
#' @return S4 object, DNB_obj
#'
mynew_dnb <- setClass(
    Class = "DNB_obj",
    slots = list(
        data = 'ANY',
        used_gene = 'ANY',
        pre_result = 'DNB_res',
        result = 'DNB_res'
    ),
    sealed = TRUE
)

#' initialize or generate the S4-DNB_res
#'
#' @param ntop the number of sample, use when initialize
#' @param ... parameters transmitting to function mynew_res() to generate S4-DNB_res after result has been computed
#'
#' @return S4 object, DNB_res
#'
mynew.DNB_res <- function(
    ntop = NULL,
    ...
) {
    if (is.null(ntop)) {
        mynew_res(...)
    } else if (is.numeric(ntop)) {
        mynew_res(
            MEAN = rep(0, ntop),
            SD = rep(0, ntop),
            CV = rep(0, ntop),
            PCC_IN = rep(0, ntop),
            PCC_OUT = rep(0, ntop),
            SCORE = rep(0, ntop),
            rank = rep(0, ntop),
            rank_all = rep(0, ntop),
            resource = rep('unknown', ntop),
            MODULEs = mynew_module(
                MODULE = as.list(rep("", ntop)),
                bestMODULE = rep(FALSE, ntop)
            )
        )
    } else stop("ERROR ntop! Should be a number!")
}

#' initialize or generate the S4-DNB_obj
#'
#' @param data raw data to compute DNB (matrix or data.frame)
#' @param ntop parameters for slot result, the number of top rank for modules
#' @param used_gene gene list of data to compute DNB, default all genes from data
#' @param ... parameters transmitting to function mynew.DNB_res() to initialize or generate the slot pre_result (S4-DNB_res)
#'
#' @return S4 object, DNB_obj
#'
mynew.DNB_obj <- function(
    data,
    ntop = NULL,
    used_gene = NULL,
    ...
) {
    mynew_dnb(
        data = data,
        used_gene = if (is.null(used_gene)) rownames(data) else used_gene,
        pre_result = mynew.DNB_res(...),
        result = mynew.DNB_res(ntop = ntop)
    )
}

#' generate the S3-DNB_output or change the attribute of a list to S3-DNB_output
#'
#' @param group group information from meta
#' @param list_obj if assigned, transfer this list into S3-DNB_output by changing its attribute
#'
#' @return S3-DNB_output
#'
mynew.DNB_output <- function(
    group,
    list_obj = NULL
) {
    if (is.null(list_obj)) {
        res <- list()
        length(res) <- length(group)
    } else {
        res <- list_obj
    }

    names(res) <- group
    class(res) <- "DNB_output"
    return(res)
}

#' show DNB_obj
#'
#' @param object DNB_obj
#'
#'
setMethod(
    f = "show",
    signature = "DNB_obj",
    definition = function(object) {
        data_dim <- dim(object@data)
        data_class <- class(object@data)[1]
        len_preres <- length(object@pre_result@rank)
        len_res <- length(object@result@rank)
        cat("A S4 object (DNB_obj), includes:\n",
            "...@data, a ", data_class, " with ", data_dim[1], " Rows and ", data_dim[2], " Cols\n",
            "...@used_gene, a list of used genes\n",
            "...@pre_result, a S4 object (DNB_res), with ", len_preres, " Modules in total\n",
            "...@result, a S4 object (DNB_res), with ", len_res, " Modules in total\n",
            sep = ""
        )
    }
)

#' show DNB_res
#'
#' @param object DNB_res
#'
#'
setMethod(
    f = "show",
    signature = "DNB_res",
    definition = function(object) {
        len <- length(object@rank)
        len2 <- length(object@MODULEs@bestMODULE)
        cat("A S4 object (DNB_res), includes: ", len, " Modules\n",
            "...@MEAN, a numeric vector of MEANs for each Module\n",
            "...@SD, a numeric vector of SDs for each Module\n",
            "...@CV, a numeric vector of CVs for each Module\n",
            "...@PCC_IN, a numeric vector of PCC_INs for each Module\n",
            "...@PCC_OUT, a numeric vector of PCC_OUTs for each Module\n",
            "...@SCORE, a numeric vector of SCORE for each Module\n",
            "...@rank, a numeric vector of rank (decreasing) for each Module in each group\n",
            "...@rank_all, a numeric vector of rank (decreasing) for each Module among groups\n",
            "...@resource, a character vector of resource (from which group) for each Module\n",
            "...@MODULEs: A S4 object (DNB_module), includes: ", len2, " Modules\n",
            "   ...@MODULE, a list of Module genes\n",
            "   ...@bestMODULE, a logical vector in order of @MODULE\n",
            sep = ""
        )
    }
)

#' show DNB_module
#'
#' @param object DNB_module
#'
#'
setMethod(
    f = "show",
    signature = "DNB_module",
    definition = function(object) {
        len <- length(object@bestMODULE)
        cat("A S4 object (DNB_module), includes: ", len, " Modules\n",
            "...@MODULE, a list of Module genes\n",
            "...@bestMODULE, a logical vector in order of @MODULE\n",
            sep = ""
        )
    }
)

#' show DNB_output
#'
#' @param x DNB_output
#' @param ... not use
#' @method print DNB_output
#' @export
#'
print.DNB_output <- function(x, ...) {
    group <- names(x)
    len <- length(group)
    cat("A S3 object (DNB_output), includes ", len, " S4 object (DNB_obj)\n",
        "The sub-objects are: ", paste(group, collapse = " "), "\n",
        "Details see DNB_output$sub-objects\n",
        sep = ""
    )
}

#' Compute data of single group for all modules from pre-result, and add the result into slot result
#'
#' @param object S4:DNB_obj
#' @param ... parameter: module_list, a list of modules, each with a gene vector, and others
#'
#' @return S4:DNB_obj
#'
setGeneric(
    "DNBcomputeSingle",
    function(object, ...)
        standardGeneric("DNBcomputeSingle")
)

#' DNBcomputeSingle for DNB_obj
#'
#' @param object S4:DNB_obj
#' @param module_list a list of module genes
#' @param force_allgene whether force to use all genes from data, default FALSE (use DNB_obj@used_gene)
#' @param size_effect whether consider the effect of sample size when compute CI of DNB, default TRUE
#'
setMethod(
    "DNBcomputeSingle", signature(object = "DNB_obj"),
    function(object, module_list, force_allgene = FALSE, size_effect = TRUE) {
        DNBcomputeSingle_DNB_obj(object = object, module_list = module_list, force_allgene = force_allgene, size_effect = size_effect)
    }
)
