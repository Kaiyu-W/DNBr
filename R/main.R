#' @title DNBcompute
#'
#' @description Compute the Dynamic Network Biomarkers(DNB) model
#'
#' @details return a S3 object includes several S4 objects
#' @details (optional) write score results into DNB_score_matrix.txt
#'
#' @param data the gene expression matrix,
#'   which can be a single-cell RNA-seq GEM with at least three group/clusters
#'   or a matrix merging bulk GEMs from at least three different sample
#' @param meta a data.frame with rownames as cell-id as well as one column of group infomation
#' @param diffgenes which genes we're interested in, or no special ones (all, default)
#' @param allgenes the whole genes that ordered in advance by expression, or the rownames of GEM (default)
#' @param meta_levels the order of meta group, default ordered by decreasing if NULL
#' @param high_method the method to select genes for the first step, by either high_cv (default) or top_gene
#' @param high_cutoff the cutoff value corresponding to the high_method,
#' 
#'  with the range between 0 - 1(all) for high_cv and 1 - #allgenes(all) for top_gene
#'
#'  or not to select highly variable genes but use all genes when -1
#' @param cutree_method the method to select numbers of tree (module) from hclust, 
#' 
#'  by either h (height, default) or k (number K)
#' @param cutree_cutoff the cutoff value corresponding to the cutree_method,
#' 
#'  with the range between 0-1 for h and a number greater than 0 for k
#' @param minModule the min number of genes of the module meeting requirements
#' @param maxModule the max number of genes of the module meeting requirements
#' @param quiet do not print output of process during calculation (against verbose), default FALSE
#' @param fastMode avoid using for loop, rathan apply-like function, default FALSE; if TRUE, quiet will be set as TRUE
#' @param writefile write results of each group into DNB_Module_information_xx.txt with tab delimiter, default FALSE
#' @param cluster_fun customized function that user design for clustering to find module, 
#' 
#'   default NULL (hierarchical by stats::hclust(d, method = "complete", members = NULL))
#'   
#'   This function should do function of clustering (e.g hclust + cutree), 
#'   
#'   with input that first arg "d" = distance matrix and output "named int vector". 
#'   
#'   If assigned, cutree_method and cutree_cutoff would be ignored
#' @param cluster_args a list of extra arguments to the cluster_fun call. 
#'   The names attribute of args gives the argument names. (same to base::do.call(args))
#' @param size_effect whether consider the effect of sample size when compute CI of DNB, default TRUE
#'
#' @return S3:DNB_output
#'
#' @examples data(data.example)
#' @examples data(meta.example)
#' @examples a <- DNBcompute(data.example, meta.example)
#' @examples a
#'
#' @author Kaiyu Wang, in ChenLab of CAS, Shanghai, China
#'
#' @export
#'
DNBcompute <- function(
    data,
    meta,
    diffgenes = NULL,
    allgenes = NULL,
    meta_levels = NULL,
    high_method = c('high_cv', 'top_gene'),
    high_cutoff = 0.6,
    cutree_method = c('h', 'k'),
    cutree_cutoff = 0.98,
    minModule = 7,
    maxModule = 60,
    quiet = FALSE,
    fastMode = FALSE,
    writefile = FALSE,
    cluster_fun = NULL,
    cluster_args = NULL,
    size_effect = TRUE
) {
    high_method <- match.arg(high_method)
    cutree_method <- match.arg(cutree_method)
    if (fastMode & !quiet)
        message("quiet will be set as TRUE if fastMode set as TRUE!")
    
    # step1. parameters config
    # if (!is.matrix(data) && !is.data.frame(data))
    if (!is(data, "Matrix") && !is.matrix(data) && !is.data.frame(data))
        stop("ERROR data input! Should be matrix or data.frame!")
    if (is.null(allgenes))
        allgenes <- rownames(data)
    if (is.null(diffgenes)) {
        diffgenes <- allgenes
    } else {
        if (!all(diffgenes %in% allgenes))
            stop("ERROR genes! Please check the allgenes or diffgenes input!")
    }
    if (!is.data.frame(meta))
        meta <- as.data.frame(meta) # if meta is a vector with names
    if (!all(colnames(data) %in% rownames(meta)))
        stop("ERROR meta! Please check the rownames of meta df (or names of meta vector) or the colnames of data input!")
    if (ncol(meta) != 1)
        stop("ERROR meta! Meta df only has one column!")
    if (!is.null(meta_levels)) {
        if (!all(meta_levels %in% meta[, 1]))
            stop("ERROR meta_levels! Should be equal to meta!")
    }

    if (high_cutoff == -1) {
        message("Use all genes when high_cutoff is -1!")
    } else {
        if (high_method == "high_cv") {
            if (high_cutoff > 1 || high_cutoff <= 0)
                stop("ERROR high_cutoff! Should be between 0 and 1 when high_method is high_cv! Support for -1 to use all genes!")
        } else {
            if (high_cutoff < 1 || high_cutoff > length(allgenes))
                stop("ERROR high_cutoff! Should be between less than length of allgenes when high_method is top_gene! Support for -1 to use all genes!")
        }
    }
    
    if (cutree_method == "h") {
        if (cutree_cutoff > 1 || cutree_cutoff <= 0)
            stop("ERROR cutree_cutoff! Should be between 0 and 1 when cutree_method is h!")
    } else {
        message("Select tree by k! Carefully choose a proper value, or error may occur later.")
        if (cutree_cutoff < 1)
            stop("ERROR cutree_cutoff! Should be between greater than 1 when cutree_method is k!")
    }

    # check customized clustering function
    if (is.null(cluster_fun)) {
        message("Modules will be computed by hierarchical clustering!")
    } else {
        if (is.function(cluster_fun)) {
            message("Modules will be computed by user customize function!")
            if (!is.null(cluster_args) && !is.list(cluster_args))
                stop("ERROR cluster_args! Should be a list or NULL!")
        } else {
            stop("ERROR cluster_fun! Should be function object!")
        }
    }

    if (is.null(meta_levels)) {
        meta <- meta[order(meta[, 1], decreasing = T), , drop = F]
        group <- unique(meta[, 1])
    } else {
        meta_tmp <- meta[1, , drop = F]
        rownames(meta_tmp) <- 'NA'
        for (i in meta_levels) {
            meta_tmp_tmp <- meta[meta[, 1] == i, , drop = F]
            meta_tmp_tmp <- meta_tmp_tmp[order(rownames(meta_tmp_tmp), decreasing = T), , drop = F]
            meta_tmp <- rbind(meta_tmp, meta_tmp_tmp)
        }
        meta <- meta_tmp[-1, , drop = F]
        group <- meta_levels
    }
    data <- data[, rownames(meta), drop = F]
    
    # step2. initiate result
    DNB_output <- mynew.DNB_output(group = group)

    # step3. process computation in order of meta
    if (fastMode) {
        iii <- 0
        pb <- utils::txtProgressBar(style = 3)
        DNB_output <- lapply(
            group,
            function(group_tmp) {
                data_tmp <- data[, rownames(meta)[meta[, 1] == group_tmp], drop = FALSE]
                # data_tmp <- data.frame(data_tmp, check.names = FALSE)
                cse <- myprocess(
                    data = data_tmp, assay_name = group_tmp, quiet = TRUE,
                    diffgenes = diffgenes, allgenes = allgenes,
                    high_method = high_method, high_cutoff = high_cutoff,
                    cutree_method = cutree_method, cutree_cutoff = cutree_cutoff,
                    minModule = minModule, maxModule = maxModule, fastMode = TRUE,
                    cluster_fun = cluster_fun, cluster_args = cluster_args, 
                    size_effect = size_effect
                )
                iii <<- iii + 1
                utils::setTxtProgressBar(pb, iii / length(group))
                return(cse)
            }
        )
        close(pb)
        DNB_output <- mynew.DNB_output(group = group, list_obj = DNB_output)
    } else {
        for (group_tmp in group) {
            cat('Now run ', group_tmp, '\n', sep = '')

            # subset data
            data_tmp <- data[, rownames(meta)[meta[, 1] == group_tmp], drop = FALSE]
            # data_tmp <- data.frame(data_tmp, check.names = FALSE)

            # process sub-data
            DNB_output[[group_tmp]] <- myprocess(
                data = data_tmp, assay_name = group_tmp, quiet = quiet,
                diffgenes = diffgenes, allgenes = allgenes,
                high_method = high_method, high_cutoff = high_cutoff,
                cutree_method = cutree_method, cutree_cutoff = cutree_cutoff,
                minModule = minModule, maxModule = maxModule, fastMode = FALSE,
                cluster_fun = cluster_fun, cluster_args = cluster_args,
                size_effect = size_effect
            )

            cat(group_tmp, ' ends up successfully!\n', sep = '')
        }
    }

    # step4. write output into files
    if (writefile)
        lapply(
            names(DNB_output),
            function(x) {
                file_names <- paste0("DNB_Module_information_", x, ".txt")
                DNB_output_mtx <- resultAllExtract(DNB_output, group = x, slot = "pre_result", mess = FALSE)
                cat("Now write output of ", x, " into ", file_names, "\n", sep = "")
                write.table(DNB_output_mtx, file = file_names, quote = F, sep = "\t")
            }
        )

    # step5. return result Object
    cat("DNB model has been computed over!\n")
    return(DNB_output)
}

#' @title DNBfilter
#'
#' @description Filter DNB model by selecting the most N modules from each group with max score,
#' @description and apply these modules to all group to recompute DNB model
#'
#' @param DNB_output S3:DNB_output, from output of DNBcompute
#' @param ntop the numbers of modules to filter
#' @param force_allgene whether force to use all genes from data, default FALSE (use the gene sets from output of DNBcompute)
#' @param size_effect whether consider the effect of sample size when compute CI of DNB, default TRUE
#' @param quiet do not message
#'
#' @return S3:DNB_output
#' @export
#'
#' @examples data(data.example)
#' @examples data(meta.example)
#' @examples a <- DNBcompute(data.example, meta.example)
#' @examples b <- DNBfilter(a, ntop = 5)
#' @examples b
#'
#' @author Kaiyu Wang, in ChenLab of CAS, Shanghai, China
#'
DNBfilter <- function(
    DNB_output,
    ntop,
    force_allgene = FALSE,
    size_effect = TRUE,
    quiet = FALSE
) {
    if (!is(DNB_output, "DNB_output"))
        stop("ERROR input! Should be S3-DNB_output!")

    group <- names(DNB_output)

    module_list_tmp <- lapply(
        DNB_output,
        function(x) {
            data <- x@pre_result
            rank_tmp <- data@rank
            index <- which(rank_tmp <= ntop)
            if (length(index) == 0) {
                NULL
            } else {
                rank_tmp <- rank_tmp[index]
                resource_tmp <- data@resource[index]
                module_tmp <- data@MODULEs@MODULE[index]
                best_tmp <- data@MODULEs@bestMODULE[index]
                names(module_tmp) <- resource_tmp
                list(
                    module = module_tmp,
                    bestMODULE = best_tmp,
                    rank = rank_tmp,
                    resource = resource_tmp
                )
            }
        }
    )

    # check filtered Modules
    tmp_Modules <- unique(unlist(lapply(module_list_tmp, function(x) if(is.null(x)) NULL else x$resource)))
    N_Module <- sum(sapply(tmp_Modules, function(x) !is.null(x)))
    if (N_Module == 0) {
        message("No Module can fit the conditions at all!")
        message("Please re-run the last step 'DNBcompute' after changing the args, such like minModule/maxModule.")
        return(DNB_output)
    } else {
        mycat("Find", N_Module, "Module(s) in total!\n", sep = " ", quiet = quiet)
    }

    # create module_list
    module_list <- list()
    length(module_list) <- 4
    names(module_list) <- c("module", "bestMODULE", "rank", "resource")
    for (i in seq(module_list_tmp)) {
        if (!is.null(module_list_tmp[[i]])) {
            module_list[['module']] <- append(module_list[['module']], module_list_tmp[[i]]$module)
            module_list[['bestMODULE']] <- append(module_list[['bestMODULE']], module_list_tmp[[i]]$bestMODULE)
            module_list[['rank']] <- append(module_list[['rank']], module_list_tmp[[i]]$rank)
            module_list[['resource']] <- append(module_list[['resource']], module_list_tmp[[i]]$resource)
        }
    }

    # de-duplicate module_list, by the unique tag: resource
    unique_index <- !duplicated(module_list[['resource']])
    for (i in seq(module_list))
        module_list[[i]] <- module_list[[i]][unique_index]

    # generate DNB_output@result
    DNB_output_new <- lapply(
        DNB_output,
        function(DNB_obj) {
            DNBcomputeSingle(
                object = DNB_obj, 
                module_list = module_list,
                force_allgene = force_allgene,
                size_effect = size_effect
            )
        }
    )
    DNB_output_new <- mynew.DNB_output(
        group = names(DNB_output_new),
        list_obj = DNB_output_new
    )

    return(DNB_output_new)
}

#' @title ScoreExtract
#'
#' @description Extract the score information from object (S3:DNB_output)
#'
#' @param object S3:DNB_output
#' @param ranking the ranking of exactly module, default 1 if NULL
#' @param group which group to select module, default random selected if NULL
#' @param resource the actual module name, ranking & group will be ignored if use
#' @param quiet do not message
#' @param ... for future use
#'
#' @return data.frame
#' @export
#'
#' @examples data(data.example)
#' @examples data(meta.example)
#' @examples a <- DNBcompute(data.example, meta.example)
#' @examples b <- DNBfilter(a, ntop = 5)
#' @examples df_score <- ScoreExtract(
#' @examples     object = b,
#' @examples     ranking = NULL,
#' @examples     group = NULL
#' @examples )
#' @examples df_score
#'
#' @author Kaiyu Wang, in ChenLab of CAS, Shanghai, China
#'
ScoreExtract <- function(
    object,
    ranking = NULL,
    group = NULL,
    resource = NULL,
    quiet = FALSE,
    ...
) {
    UseMethod("ScoreExtract")
}

#' @title resultAllExtract
#'
#' @description Get the whole result of slot pre_result or result from S4-DNB_obj and transfer that into a data.frame (matrix)
#'
#' @param object the S3-DNB_output
#' @param group which sub-group of S3-DNB_output to extract
#' @param slot which S4-DNB_obj to extract, pre_result or result
#' @param mess whether to message if slot is pre_result, default TRUE
#' @param ... for future use
#'
#' @return data.frame
#' @export
#'
#' @examples data(data.example)
#' @examples data(meta.example)
#' @examples a <- DNBcompute(data.example, meta.example)
#' @examples b <- DNBfilter(a, ntop = 5)
#' @examples df_score <- resultAllExtract(
#' @examples     object = b,
#' @examples     group = "A",
#' @examples     slot = "pre_result"
#' @examples )
#' @examples df_score
#'
#' @author Kaiyu Wang, in ChenLab of CAS, Shanghai, China
#'
resultAllExtract <- function(
    object,
    group, 
    slot,
    mess = TRUE,
    ...
) {
    UseMethod("resultAllExtract")
}

#' @title getMaxRanking
#'
#' @description Get the max Ranking of exactly group in S3:DNB_output after function DNBfilter
#'
#' @param object S3:DNB_output
#' @param group which sub-group of S3-DNB_output to extract
#' @param ... for future use
#'
#' @return numeric value of max Ranking for your group
#' @export
#'
#' @examples data(data.example)
#' @examples data(meta.example)
#' @examples a <- DNBcompute(data.example, meta.example)
#' @examples b <- DNBfilter(a, ntop = 5)
#' @examples maxRank <- getMaxRanking(b, group = "C") # get 4 instead of 5
#' @examples DNBplot(b, ranking = maxRank, group = "C", show = TRUE, save_pdf = FALSE)
#'
#' @author Kaiyu Wang, in ChenLab of CAS, Shanghai, China
#'
getMaxRanking <- function(
    object,
    group,
    ...
) {
    UseMethod("getMaxRanking")
}

#' @title getModuleGenes
#'
#' @description Get the gene list of exactly module in S3:DNB_output
#'
#' @param object S3:DNB_output
#' @param resource that is, the actual module name
#' @param ... for future use
#'
#' @return vector of genes
#' @export
#'
#' @examples data(data.example)
#' @examples data(meta.example)
#' @examples a <- DNBcompute(data.example, meta.example)
#' @examples b <- DNBfilter(a, ntop = 5)
#' @examples A_11_genelist <- getModuleGenes(b, resource = "A_11")
#' @examples print(A_11_genelist)
#'
#' @author Kaiyu Wang, in ChenLab of CAS, Shanghai, China
#'
getModuleGenes <- function(
    object,
    resource,
    ...
) {
    UseMethod("getModuleGenes")
}

#' @title DNBplot
#'
#' @description Visualize the DNB score information
#'
#' @param object S3:DNB_output or df_score (output of ScoreExtract)
#' @param ranking for DNB_output, default 1 if NULL
#' @param group for DNB_output, default random selected if NULL
#' @param resource the actual module name, ranking & group will be ignored if use 
#' @param show whether to show, default TRUE
#' @param save_pdf whether to save pdf file, default FALSE
#' @param file_prefix the file prefix if save_pdf, default NULL
#' @param meta_levels the order of meta-group in the plots, default levels(df_score$Names) if NULL
#' @param ... for future use
#'
#' @return plot or pdf
#' @export
#'
#' @examples data(data.example)
#' @examples data(meta.example)
#' @examples a <- DNBcompute(data.example, meta.example)
#' @examples b <- DNBfilter(a, ntop = 5)
#' @examples DNBplot(b, show = TRUE, save_pdf = FALSE)
#'
#' @author Kaiyu Wang, in ChenLab of CAS, Shanghai, China
#'
DNBplot <- function(
    object, 
    ranking = NULL, 
    group = NULL, 
    resource = NULL,
    show = TRUE, 
    save_pdf = FALSE, 
    file_prefix = NULL,
    meta_levels = NULL,
    ...
) {
    UseMethod("DNBplot")
}

#' @title DNBcompute_custom
#'
#' @description Compute the Dynamic Network Biomarkers(DNB) model with customized Modules
#'
#' @details input data from existing DNB_output object, or new data and meta
#' @details input customized modules from a list of gene set(s)
#' @details return a S3 object includes several S4 objects
#'
#' @param data input data, S4:DNB_output or gene expression matrix
#' @param module_list a customized list of module gene
#' @param meta a data.frame with rownames as cell-id as well as one column of group infomation, 
#'   use if data is gene expression matrix
#' @param meta_levels the order of meta group, default ordered by decreasing if NULL, 
#'   use if data is gene expression matrix
#' @param force_allgene whether force to use all genes from data, default TRUE
#' @param size_effect whether consider the effect of sample size when compute CI of DNB, default TRUE
#' @param quiet do not message
#' @param ... for future use
#' 
#' @return S3:DNB_output
#'
#' @examples data(data.example)
#' @examples data(meta.example)
#' @examples data(module_list.example)
#' @examples b <- DNBcompute_custom(data.example, module_list.example, meta.example)
#' @examples b
#' @examples DNBplot(b, group = "USER_CUSTOMIZED", ranking = 1, show = TRUE, save_pdf = FALSE)
#' @examples DNBplot(b, group = "USER_CUSTOMIZED", ranking = 1, show = TRUE, save_pdf = FALSE)
#'
#' @author Kaiyu Wang, in ChenLab of CAS, Shanghai, China
#'
#' @export
#'
DNBcompute_custom <- function(
    data,
    module_list,
    meta = NULL,
    meta_levels = NULL,
    force_allgene = TRUE,
    size_effect = TRUE,
    quiet = FALSE,
    ...
) {
    UseMethod("DNBcompute_custom")
}
