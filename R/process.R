#' Process function of each group data
#' 
#' @description the parameters come from inner DNBcompute
#'
#' @param data same to DNBcompute
#' @param assay_name group (one of meta_levels)
#' @param diffgenes same to DNBcompute
#' @param allgenes same to DNBcompute
#' @param high_method same to DNBcompute
#' @param high_cutoff same to DNBcompute
#' @param cutree_method same to DNBcompute
#' @param cutree_cutoff same to DNBcompute
#' @param minModule same to DNBcompute
#' @param maxModule same to DNBcompute
#' @param quiet same to DNBcompute
#' @param fastMode same to DNBcompute
#' @param cluster_fun same to DNBcompute
#' @param cluster_args same to DNBcompute
#' @param size_effect same to DNBcompute
#'
#'
myprocess <- function(
    data,
    assay_name,
    diffgenes,
    allgenes,
    high_method,
    high_cutoff,
    cutree_method,
    cutree_cutoff,
    minModule,
    maxModule,
    quiet,
    fastMode,
    cluster_fun,
    cluster_args,
    size_effect
) {
    mycat("\tProcess ...\n", quiet = quiet)
    data0 <- data
    data <- t(data)

    sd_all <- apply(data, 2, sd)
    mean_all <- apply(data, 2, mean)
    cv_all <- sd_all / mean_all

    # select genes
    if (high_cutoff == -1) {
        hig_gene <- allgenes
    } else {
        if (high_method == 'top_gene') {
            hig_gene <- allgenes[1:high_cutoff]
        } else if (high_method == 'high_cv') {
            hig_num <- ceiling(length(cv_all) * high_cutoff) #cv cutoff
            cv_rank <- sort(cv_all, decreasing = T, na.last = T)[1:hig_num]
            hig_gene <- names(cv_rank)
        }
    }
    if (length(hig_gene) <= 1)
        stop("ERROR number of high gene! Only find one high gene or less!")
    sle_tab <- data[, hig_gene, drop = F]

    # calculate distance (pcc)
    cor_sle <- mycor(sle_tab, method = "pearson")
    diag(cor_sle) <- 0 # useless
    dis_sle <- 1 - abs(cor_sle)
    dist <- as.dist(dis_sle)
    dist[is.na(dist)] <- max(dist, na.rm = T) * 1.01

    # clustering
    mycat("\tCluster...\n", quiet = quiet)
    if (is.null(cluster_fun)) {
        mycat("\tCluster method be stats::hclust and cutree\n", quiet = quiet)
        mycat("\tCutree method be ", cutree_method, " with cutoff value ", cutree_cutoff, "\n", sep = "", quiet = quiet)
        clu_sle <- stats::hclust(d = dist)
        if (cutree_method == 'h')
            cut_sle <- stats::cutree(clu_sle, h = cutree_cutoff) # cutoff
        else if (cutree_method == 'k')
            tryCatch(
                cut_sle <- stats::cutree(clu_sle, k = cutree_cutoff),
                error = function(x) stop("Not proper value of k!")
            )
    } else if (is.function(cluster_fun)) {
        mycat("\tCluster method be user's customized function\n", quiet = quiet)
        tryCatch(
            cut_sle <- do.call(what = cluster_fun, args = append(list(d = dist), cluster_args)),
            error = function(x) 
                stop("Unknown error! Please check the cluster_fun with name of first arg that be 'd', or cluster_args format!")
        )
    } else {
        stop("ERROR cluster_fun input!")
    }
    
    # calcute score
    size_up2_1 <- sapply(
        min(cut_sle):max(cut_sle),
        function(x)
            sum(cut_sle == x)
    ) > 1
    cut_sle_up2_1 <- (min(cut_sle):max(cut_sle))[size_up2_1]
    len_size_up2_1 <- length(cut_sle_up2_1)
    mycat("\tFind total of ", len_size_up2_1, " scores! \n", sep = "", quiet = quiet)

    # if no module fit, return null object
    if (len_size_up2_1 == 0) {
        res_obj <- mynew.DNB_obj(
            data = data0,
            used_gene = hig_gene,
            ntop = 0
        )
        return(res_obj)
    }

    score_mtx <- matrix(
        0, ncol = len_size_up2_1, nrow = 7,
        dimnames = list(
            c("mean_mean", "sd_mean", "cv_mean", "pcc_in", "pcc_out", "score", "bestMODULE"),
            cut_sle_up2_1
        )
    )

    mycat("\tCompute all modules' score ...\n", quiet = quiet)
    if (fastMode) {
        score_mtx_tmp <- sapply(
            cut_sle_up2_1,
            function(i) {
                compute_score_eachgroup(
                    cut_sle = cut_sle, group = i, diffgenes = diffgenes, hig_gene = hig_gene,
                    sle_tab = sle_tab, cv_all = cv_all, sd_all = sd_all, mean_all = mean_all,
                    minModule = minModule, maxModule = maxModule, size_effect = size_effect
                )
            }
        )
        colnames(score_mtx_tmp) <- colnames(score_mtx)
        rownames(score_mtx_tmp) <- rownames(score_mtx)
        score_mtx <- score_mtx_tmp
    } else {
        for (count in seq(cut_sle_up2_1)) {
            i <- cut_sle_up2_1[count]
            if (i %% 50 == 0)
                mycat("\t\tNow run ", i, "\n", sep = "", quiet = quiet)

            score_mtx[, count] <- compute_score_eachgroup(
                cut_sle = cut_sle, group = i, diffgenes = diffgenes, hig_gene = hig_gene,
                sle_tab = sle_tab, cv_all = cv_all, sd_all = sd_all, mean_all = mean_all,
                minModule = minModule, maxModule = maxModule, size_effect = size_effect
            )
        }
    }

    # Pre-Filter, to check all the modules
    if (!quiet) {
        mycat("\tFilter Modules by intersected diffgenes ... \n", quiet = quiet)
        fit_n <- 1
        if (len_size_up2_1 > 0) {
            for (mm in seq(cut_sle_up2_1)) {
                MODULE_names <- names(which(cut_sle == cut_sle_up2_1[mm]))
                intersectDiff <- intersect(diffgenes, MODULE_names)
                if (length(intersectDiff) > 0) {
                    mycat("\t\tFind Module_", fit_n, ": ", sep = "", quiet = quiet)
                    mycat("CI->", score_mtx[6, mm], "; ", sep = "", quiet = quiet)
                    mycat("number of gene->", length(MODULE_names), "; ", sep = "", quiet = quiet)
                    mycat("intersested genes->", paste0(intersectDiff, collapse = ","), "\n", sep = "", quiet = quiet)
                    fit_n <- fit_n + 1
                }
            }
        } else {
            mycat("\t\tFind no module!\n", quiet = quiet)
        }
    }

    # generate output
    Module_list <- lapply(
        cut_sle_up2_1,
        function(x)
            names(which(cut_sle == x))
    )

    # rank_tmp <- score_mtx[7, ]
    # rank_order1 <- -rank(score_mtx[6, ]) + length(score_mtx[6, ]) + 1
    # rank_order2 <- max(rank_order1) - rank_order1 + 1
    # rank_order3 <- rank_order2 * rank_tmp
    # rank_order <- -rank(rank_order3) + length(rank_order3) + 1
    rank_order <- rank(-1 * score_mtx[7, ] * score_mtx[6, ], ties.method = 'average')

    # Module_resource <- rep(assay_name, length(rank_order))
    Module_resource <- paste(assay_name, names(rank_order), sep = "_")
    names(Module_list) <- Module_resource

    res_obj <- mynew.DNB_obj(
        data = data0,
        used_gene = hig_gene,
        MEAN = score_mtx[1, ],
        SD = score_mtx[2, ],
        CV = score_mtx[3, ],
        PCC_IN = score_mtx[4, ],
        PCC_OUT = score_mtx[5, ],
        SCORE = score_mtx[6, ],
        rank = rank_order,
        rank_all = rep(0, length(rank_order)),
        resource = Module_resource,
        MODULEs = mynew_module(
            MODULE = Module_list,
            bestMODULE = as.logical(score_mtx[7, ])
        )
    )
    return(res_obj)
}

#' process function after modules have been defined
#'
#' @description the parameters come from inner DNBcompute
#'
#' @param cut_sle module information table (cluster, which module each gene belongs to)
#' @param group group (assay_name)
#' @param diffgenes same to DNBcompute
#' @param hig_gene same to DNBcompute
#' @param sle_tab selected table
#' @param cv_all all cv values
#' @param sd_all all sd values
#' @param mean_all all mean values
#' @param minModule same to DNBcompute
#' @param maxModule same to DNBcompute
#' @param size_effect same to DNBcompute
#'
compute_score_eachgroup <- function(
    cut_sle,
    group,
    diffgenes,
    hig_gene,
    sle_tab,
    cv_all,
    sd_all,
    mean_all,
    minModule,
    maxModule,
    size_effect
) {
    group_name <- names(which(cut_sle == group))
    group_out <- setdiff(hig_gene, group_name)
    LenMODULE <- length(group_name)
    Lendiff <- length(intersect(diffgenes, group_name))
    bestModule <- Lendiff > 0 && LenMODULE > minModule && LenMODULE < maxModule

    size <- length(group_name)

    tab_in <- sle_tab[, group_name, drop = F]
    tab_out <- sle_tab[, group_out, drop = F]

    # calculate parameters
    cv <- mean(cv_all[group_name], na.rm = T)
    Sd <- mean(sd_all[group_name])
    Mean <- mean(mean_all[group_name])
    pcc_in <- mean(abs(as.dist(mycor(tab_in))), na.rm = T)
    pcc_out <- mean(abs(mycor(tab_in, tab_out)), na.rm = T)
    score <- CI(
        n = ifelse(size_effect, size, 1), 
        Sd = Sd, pcc_in = pcc_in, pcc_out = pcc_out
    )

    res <- c(Mean, Sd, cv, pcc_in, pcc_out, score, bestModule)
    return(res)
}

#' DNBcomputeSingle for DNB_obj
#'
#' @param object S4:DNB_obj
#' @param module_list a list of module genes
#' @param force_allgene whether force to use all genes from data, default FALSE (use DNB_obj@used_gene)
#' @param size_effect whether consider the effect of sample size when compute CI of DNB, default TRUE
#'
#' @return DNB_obj
#'
DNBcomputeSingle_DNB_obj <- function(
    object, 
    module_list, 
    force_allgene = FALSE, 
    size_effect = TRUE
) {
    data <- object@data
    used_gene <- object@used_gene
    module_geneslist <- module_list$module
    bestMODULE <- module_list$bestMODULE
    Rank <- module_list$rank
    resource <- module_list$resource
    score_mtx <- singleDNB(
        data = if (force_allgene) data else data[used_gene, , drop = FALSE], 
        module_list = module_geneslist,
        size_effect = size_effect
    )
    Rank_all <- rank(-1 * score_mtx[6, ])
    
    result <- mynew.DNB_res(
        MEAN = score_mtx[1, ],
        SD = score_mtx[2, ],
        CV = score_mtx[3, ],
        PCC_IN = score_mtx[4, ],
        PCC_OUT = score_mtx[5, ],
        SCORE = score_mtx[6, ],
        rank = Rank,
        rank_all = Rank_all,
        resource = resource,
        MODULEs = mynew_module(
            MODULE = as.list(module_geneslist),
            bestMODULE = bestMODULE
        )
    )
    
    object@result <- result
    return(object)
}

#' Compute data for the modules list, and return score matrix
#'
#' @param data data
#' @param module_list modules list
#' @param cor_vec vector of correlation coefficient, when data has only one sample (col)
#' @param size_effect whether consider the effect of sample size when compute CI of DNB, default TRUE
#'
#' @return score matrix
#'
singleDNB <- function(
    data,
    module_list, # module_list is a list
    cor_vec = NULL,
    size_effect = TRUE
) {
    # prepare data
    n_sample <- ncol(data)
    all_genes <- rownames(data)
    all_module_genes <- unlist(module_list)
    if (!all(all_module_genes %in% all_genes)) {
        drop_module_genes <- setdiff(all_module_genes, all_genes)
        warning("Cannot find genes: ", paste0(drop_module_genes, collapse = ','), "! Set force_allgene to use all genes, or make sure consistency between genes from modules and data!")
        module_list <- lapply(
            module_list,
            function(x) setdiff(x, drop_module_genes)
        )
    }
    
    if (n_sample == 1) {
        if (is.null(cor_vec))
            stop("cor_vec is required when data has only one sample!")
        if (length(cor_vec) != length(all_genes))
            stop("cor_vec should be the same length of singleDNB!")
    }

    # # define result
    # score_mtx <- matrix(
    #     0, nrow = 6, ncol = length(module_list),
    #     dimnames = list(
    #         c("Mean", "Sd", "cv", "pcc_in", "pcc_out", "score"),
    #         names(module_list)
    #     )
    # )
    score_row <- c("Mean", "Sd", "cv", "pcc_in", "pcc_out", "score")
    score_col <- names(module_list)
    
    if (n_sample > 1) {
        # base computation
        data_tmp <- t(data)
        sd_all <- apply(data_tmp, 2, sd)
        mean_all <- apply(data_tmp, 2, mean)
        cv_all <- sd_all / mean_all

        # compute scores
        score_mtx <- sapply(
            module_list,
            function(module_genes) {
                dnb_in <- module_genes
                dnb_out <- setdiff(all_genes, dnb_in)
                tab_in <- data_tmp[, dnb_in, drop = F]
                tab_out <- data_tmp[, dnb_out, drop = F]
                size <- length(dnb_in)

                Mean <- mean(mean_all[dnb_in])
                Sd <- mean(sd_all[dnb_in])
                cv <- mean(cv_all[dnb_in], na.rm = T) # if mean is zero, cv is not defined, so do not consider this one.
                pcc_in <- mean(abs(as.dist(mycor(tab_in))), na.rm = T)
                pcc_out <- mean(abs(mycor(tab_in, tab_out)), na.rm = T)
                score <- CI(
                    n = ifelse(size_effect, size, 1), 
                    Sd = Sd, pcc_in = pcc_in, pcc_out = pcc_out
                )

                c(Mean, Sd, cv, pcc_in, pcc_out, score)
            }
        )
    } else {
        # only one sample
        data_tmp <- data[, 1, drop = T]
        score_mtx <- sapply(
            module_list,
            function(module_genes) {
                dnb_in <- module_genes
                dnb_out <- setdiff(all_genes, dnb_in)
                tab_in <- data_tmp[dnb_in]
                tab_out <- data_tmp[dnb_out]
                size <- length(dnb_in)

                Mean <- mean(tab_in)
                Sd <- sd(tab_in)
                cv <- if (Mean == 0) NA else Sd / Mean
                pcc_in <- mean(cor_vec[tab_in])
                pcc_out <- mean(cor_vec[tab_out])
                score <- CI(
                    n = ifelse(size_effect, size, 1), 
                    Sd = Sd, pcc_in = pcc_in, pcc_out = pcc_out
                )

                c(Mean, Sd, cv, pcc_in, pcc_out, score)
            }
        )
    }

    rownames(score_mtx) <- score_row
    colnames(score_mtx) <- score_col

    return(score_mtx)
}

#' Get the max Ranking of exactly group in S3:DNB_output after function DNBfilter
#'
#' @param object S3:DNB_output
#' @param group which sub-group of S3-DNB_output to extract
#' @param ... not use
#'
#' @return numeric value
#' @method getMaxRanking DNB_output
#' @export
#'
getMaxRanking.DNB_output <- function(
    object,
    group,
    ...
) {
    all_group <- names(object)
    if (is.null(group)) {
        stop("group must has value!")
    } else {
        if (!group %in% all_group)
            if (group != "USER_CUSTOMIZED")
                stop("ERROR group input! Should be one of meta_levels or USER_CUSTOMIZED!")
    }

    sum(grepl(paste0("^", group, "_"), object[[1]]@result@resource))
}

#' Get the gene list of exactly module in S3:DNB_output
#'
#' @param object S3:DNB_output
#' @param resource that is, the actual module name
#' @param ... not use
#'
#' @return vector of genes or character()
#' @method getModuleGenes DNB_output
#' @export
#'
getModuleGenes.DNB_output <- function(
    object,
    resource,
    ...
) {
    meta_levels <- names(object)
    for (i in meta_levels) {
        data <- object[[i]]
        genelist <- getModuleGenes(data, resource)

        if (is.null(genelist))
            next
        else
            return(genelist)
    }

    # otherwise
    cat("Cannot find Module in either @result or @pre_result\n", sep = "")
    return(character())
}

#' Get the gene list of exactly module in S4:DNB_obj
#'
#' @param object S4:DNB_obj
#' @param resource that is, the actual module name
#' @param ... not use
#'
#' @return vector of genes or NULL
#' @method getModuleGenes DNB_obj
#' @export
#'
getModuleGenes.DNB_obj <- function(
    object,
    resource,
    ...
) {
    for (xslot in c("result", "pre_result")) {
        data <- slot(object, xslot)
        index <- which(data@resource == resource)
        if (length(index) == 0) {
            next
        } else {
            cat("Find Module in @", xslot, "\n", sep = "")
            genelist <- data@MODULEs@MODULE[[index]]

            return(genelist)
        }
    }

    # otherwise
    NULL
}

#' Get the whole result of slot pre_result or result from S4-DNB_obj and transfer that into a data.frame (matrix)
#'
#' @param object the S3-DNB_output
#' @param group which sub-group of S3-DNB_output to extract
#' @param slot which S4-DNB_obj to extract, pre_result or result
#' @param mess whether to message if slot is pre_result, default TRUE
#' @param ... not use
#'
#' @return data.frame
#' @method resultAllExtract DNB_output
#' @export
#'
resultAllExtract.DNB_output <- function(
    object,
    group,
    slot,
    mess = TRUE,
    ...
) {
    slot <- match.arg(arg = slot, choices = c("pre_result", "result"), several.ok = FALSE)
    if (slot == "pre_result" & mess)
        message("Modules from pre_result may be of a large amount, please use it carefully when printing directly!")

    result <- methods::slot(object[[group]], slot)

    MEAN <- result@MEAN
    SD <- result@SD
    CV <- result@CV
    PCC_IN <- result@PCC_IN
    PCC_OUT <- result@PCC_OUT
    SCORE <- result@SCORE
    rank <- result@rank
    rank_all <- result@rank_all
    resource <- result@resource
    bestMODULE <- result@MODULEs@bestMODULE
    genes <- sapply(result@MODULEs@MODULE, function(x) paste0(x, collapse = ","))

    data.frame(
        rank = rank,
        rank_all = rank_all,
        bestMODULE = bestMODULE,
        SCORE = SCORE,
        genes = genes,
        MEAN = MEAN,
        SD = SD,
        CV = CV,
        PCC_IN = PCC_IN,
        PCC_OUT = PCC_OUT,
        resource = resource
    )
}

#' Extract the score information of @result from object (S3:DNB_output)
#'
#' @param object S3:DNB_output
#' @param ranking the ranking of exactly module, default 1 if NULL
#' @param group which group to select module, default random selected if NULL
#' @param resource the actual module name, ranking & group will be ignored if use
#' @param quiet do not message 
#' @param ... not use
#'
#' @return score data.frame
#' @method ScoreExtract DNB_output
#' @export
#'
ScoreExtract.DNB_output <- function(
    object,
    ranking = NULL,
    group = NULL,
    resource = NULL,
    quiet = FALSE,
    ...
) {
    all_group <- names(object)

    data <- object[[1]]@result
    if (length(data@rank) == 0)
        stop("No module found! Please check or run DNBfilter() !")

    if (is.null(resource)) {
        if (is.null(group)) {
            message("Randomly select group...")
            group <- sample(all_group, 1)
        } else {
            if (!group %in% all_group)
                if (group != "USER_CUSTOMIZED")
                    stop("ERROR group input! Should be one of meta_levels or USER_CUSTOMIZED!")
        }

        max_lenModule <- getMaxRanking(object, group = group)
        if (is.null(ranking)) {
            ranking <- 1
        } else {
            if (ranking < 1 | ranking > max_lenModule)
                stop("ERROR ranking input! Should be between 1 and ", max_lenModule, "!")
        }

        mycat("Use group=", group, " and ranking=", ranking, "\n", sep = "", quiet = quiet)

        index_group <- which(grepl(paste0("^", group, "_"), data@resource))
        index_rank <- which(data@rank == ranking)
        index <- intersect(index_group, index_rank)
        if (length(index) == 0) {
            rank_tmp <- object[[1]]@result@rank
            rank_tmp[data@resource != group] <- NA
            rank_tmp <- rank(-1 * rank_tmp, ties.method = "random")
            index <- which(rank_tmp == ranking)
            if (length(index) == 0)
                stop("Unknown ERROR!")
        }

        resource <- data@resource[index]
        mycat("Find resource=", resource, "\n", sep = "", quiet = quiet)
    } else {
        index <- which(data@resource == resource)
        if (length(index) == 0)
            stop("ERROR resource input! ", resource, " cannot be found in object!")
    }
    
    df_score <- data.frame(
        Names = factor(all_group, levels = all_group),
        SCORE = 0,
        PCC_IN = 0,
        PCC_OUT = 0,
        SD = 0
    )

    for (i in seq(object)) {
        index_tmp <- which(df_score$Names == names(object)[i])
        df_score$SCORE[index_tmp] <- object[[i]]@result@SCORE[index]
        df_score$PCC_IN[index_tmp] <- object[[i]]@result@PCC_IN[index]
        df_score$PCC_OUT[index_tmp] <- object[[i]]@result@PCC_OUT[index]
        df_score$SD[index_tmp] <- object[[i]]@result@SD[index]
    }

    return(df_score)
}

#' Compute the Dynamic Network Biomarkers(DNB) model with customized Modules
#'
#' @param data the gene expression matrix,
#'      which can be a single-cell RNA-seq GEM with at least three group/clusters
#'      or a matrix merging bulk GEMs from at least three different sample
#' @param module_list a customized list of module gene
#' @param meta a data.frame with rownames as cell-id as well as one column of group infomation
#' @param meta_levels the order of meta group, default ordered by decreasing if NULL
#' @param force_allgene not use
#' @param size_effect whether consider the effect of sample size when compute CI of DNB, default TRUE
#' @param quiet do not message
#' @param ... not use
#'
#' @return S3:DNB_output
#' @method DNBcompute_custom default
#' @export
#'
DNBcompute_custom.default <- function(
    data, 
    module_list, 
    meta, 
    meta_levels = NULL,
    force_allgene = TRUE,
    size_effect = TRUE,
    quiet = FALSE,
    ...
) {
    # if (!is.matrix(data) && !is.data.frame(data))
    if (!is(data, "Matrix") && !is.matrix(data) && !is.data.frame(data))
        stop("ERROR data input! Should be matrix or data.frame!")
    if (!is.data.frame(meta))
        meta <- as.data.frame(meta) # if meta is a vector with names
    if (!is.list(module_list))
        stop("ERROR module_list! Please input a list!")
    if (!all(unique(unlist(module_list)) %in% rownames(data)))
        stop("ERROR module_list! Please check the all names of module_list or the rownames of data input!")
    if (!all(colnames(data) %in% rownames(meta)))
        stop("ERROR meta! Please check the rownames of meta df (or names of meta vector) or the colnames of data input!")
    if (ncol(meta) != 1)
        stop("ERROR meta! Meta df only has one column!")
    if (!is.null(meta_levels))
        if (!all(meta_levels %in% meta[, 1]))
            stop("ERROR meta_levels! Should be equal to meta!")

    if (is.null(meta_levels)) {
        meta <- meta[order(meta[, 1], decreasing = T), , drop = F]
        meta_levels <- unique(meta[, 1])
    } else {
        meta_tmp <- meta[1, , drop = F]
        rownames(meta_tmp) <- 'NA'
        for (i in meta_levels) {
            meta_tmp_tmp <- meta[meta[, 1] == i, , drop = F]
            meta_tmp_tmp <- meta_tmp_tmp[order(rownames(meta_tmp_tmp), decreasing = T), , drop = F]
            meta_tmp <- rbind(meta_tmp, meta_tmp_tmp)
        }
        meta <- meta_tmp[-1, , drop = F]
    }
    data <- data[, rownames(meta), drop = F]

    module_len <- length(module_list)
    if (module_len == 0) 
        stop("ERROR module_list! No elements!")

    if (is.null(names(module_list))) {
        module_resource <- paste('USER_CUSTOMIZED', 1:module_len, sep = "_")
        message("The resource is named with prefix 'USER_CUSTOMIZED', so use group = 'USER_CUSTOMIZED' when DNBplot")
    } else {
        module_resource <- names(module_list)
    }

    DNB_output <- list()
    for (i in meta_levels) {
        data_tmp <- data[, rownames(meta)[meta[, 1] == i], drop = FALSE]
        data_tmp <- data.frame(data_tmp, check.names = FALSE)

        DNB_output[[i]] <- mynew_dnb(
            data = data_tmp, 
            used_gene = rownames(data_tmp),
            pre_result = mynew.DNB_res(ntop = module_len), 
            result = mynew.DNB_res(ntop = module_len)
        )
        # DNB_output[[i]]@pre_result@resource <- rep('USER_CUSTOMIZED', module_len)
        DNB_output[[i]]@pre_result@resource <- module_resource
        DNB_output[[i]]@pre_result@rank <- seq(module_list)
        DNB_output[[i]]@pre_result@MODULEs@MODULE <- module_list
        DNB_output[[i]]@pre_result@MODULEs@bestMODULE <- rep(TRUE, module_len)
    }
    DNB_output <- mynew.DNB_output(group = meta_levels, list_obj = DNB_output)

    DNB_output <- DNBfilter(DNB_output, ntop = module_len, size_effect = size_effect, quiet = quiet)

    return(DNB_output)
}

#' Compute the Dynamic Network Biomarkers(DNB) model with customized Modules
#'
#' @param data existing S3:DNB_output object
#' @param module_list a customized list of module gene
#' @param meta not use
#' @param meta_levels not use
#' @param module_list a list of module genes
#' @param force_allgene whether force to use all genes from data, default TRUE
#' @param size_effect whether consider the effect of sample size when compute CI of DNB, default TRUE
#' @param quiet do not message
#' @param ... not use
#'
#' @return S3:DNB_output
#' @method DNBcompute_custom DNB_output
#' @export
#'
DNBcompute_custom.DNB_output <- function(
    data,
    module_list,
    meta = NULL,
    meta_levels = NULL,
    force_allgene = TRUE,
    size_effect = TRUE,
    quiet = FALSE,
    ...
) {
    object <- data

    # get meta levels
    meta_levels <- names(object)

    # get data matrix and meta
    data_list <- if (force_allgene) 
        lapply(object, function(x) x@data)
    else
        lapply(object, function(x) x@data[x@used_gene, , drop = FALSE])
    data <- data.frame()
    meta_df <- data.frame()
    for (i in seq(meta_levels)) {
        data_i <- data_list[[i]]

        # get data
        data <- if (i == 1)
            rbind(data, data_i)
        else
            cbind(data, data_i)

        # get meta
        meta_df_tmp <- data.frame(a = colnames(data_i))
        colnames(meta_df_tmp) <- meta_levels[i]
        meta_df <- if (i == 1)
            meta_df_tmp
        else
            cbind(meta_df, meta_df_tmp)
    }

    # get meta
    meta_tmp <- melt(meta_df, id.var = c())
    meta <- data.frame(meta_tmp$variable, row.names = meta_tmp$value)

    # compute DNB
    object_new <- DNBcompute_custom(
        data = data, 
        meta = meta, 
        module_list = module_list, 
        meta_levels = meta_levels,
        size_effect = size_effect
    )

    # combine the new object into the old one
    for (group in meta_levels) {
        object[[group]]@result <- combine_slot(
            x = object[[group]]@result, 
            y = object_new[[group]]@result
        )
    }

    # return merged DNB_output
    return(object)
}

#' @title combine_slot
#'
#' @description Combine different slots together from same data for S4 object
#'
#' @param x slot 1
#' @param y slot 2
#'
#' @return S4:combined_slot
#'
combine_slot <- function(x, y) {
    slot_names <- slotNames(x)
    slot_names2 <- slotNames(y)
    if (!all(slot_names == slot_names2))
        stop("ERROR slot input!")

    for (i in slot_names) {
        slot1 <- slot(x, i)
        slot2 <- slot(y, i)

        slot12 <- if (isS4(slot1))
            combine_slot(slot1, slot2)
        else
            slot12 <- c(slot1, slot2)

        slot(x, i) <- slot12
    }

    return(x)
}