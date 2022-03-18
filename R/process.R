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
    cluster_args
) {
    mycat("\tProcess ...\n", quiet = quiet)
    data0 <- data
    data <- t(data)

    sd_all <- apply(data, 2, sd)
    mean_all <- apply(data, 2, mean)
    cv_all <- sd_all / mean_all

    # select genes
    if (high_method == 'top_gene') {
        hig_gene <- allgenes[1:high_cutoff]
    } else if (high_method == 'high_cv') {
        hig_num <- ceiling(length(cv_all) * high_cutoff) #cv cutoff
        cv_rank <- sort(cv_all, decreasing = T, na.last = T)[1:hig_num]
        hig_gene <- names(cv_rank)
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
        mycat("\tcutree method be ", cutree_method, " with cutoff value ", cutree_cutoff, "\n", sep = "", quiet = quiet)
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
                    minModule = minModule, maxModule = maxModule)
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
                minModule = minModule, maxModule = maxModule)
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

    res_obj <- mynew.DNB_obj(
        data = data0,
        MEAN = score_mtx[1, ],
        SD = score_mtx[2, ],
        CV = score_mtx[3, ],
        PCC_IN = score_mtx[4, ],
        PCC_OUT = score_mtx[5, ],
        SCORE = score_mtx[6, ],
        rank = rank_order,
        rank_all = rep(0, length(rank_order)),
        # resource = rep(assay_name, length(rank_order)),
        resource = paste(assay_name, names(rank_order), sep = "_"),
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
    maxModule
) {
    group_name <- names(which(cut_sle == group))
    group_out <- setdiff(hig_gene, group_name)
    LenMODULE <- length(group_name)
    Lendiff <- length(intersect(diffgenes, group_name))
    bestModule <- if (Lendiff > 0 && LenMODULE > minModule && LenMODULE < maxModule) TRUE else FALSE

    size <- length(group_name)

    tab_in <- sle_tab[, group_name, drop = F]
    tab_out <- sle_tab[, group_out, drop = F]

    # calculate parameters
    cv <- mean(cv_all[group_name], na.rm = T)
    Sd <- mean(sd_all[group_name])
    Mean <- mean(mean_all[group_name])
    pcc_in <- mean(abs(as.dist(mycor(tab_in))), na.rm = T)
    pcc_out <- mean(abs(mycor(tab_in, tab_out)), na.rm = T)

    score <- ifelse(
        pcc_out == 0,
        0,
        sqrt(size) * Sd * pcc_in / pcc_out # score for ci
    )

    res <- c(Mean, Sd, cv, pcc_in, pcc_out, score, bestModule)
    return(res)
}

#' DNBcomputeSingle for DNB_obj
#'
#' @param object S4:DNB_obj
#' @param module_list a list of module genes
#'
#' @return DNB_obj
#'
DNBcomputeSingle_DNB_obj <- function(object, module_list) {
    data <- object@data
    module_geneslist <- module_list$module
    bestMODULE <- module_list$bestMODULE
    Rank <- module_list$rank
    resource <- module_list$resource
    score_mtx <- singleDNB(data = data, module_list = module_geneslist)
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
#'
#' @return score matrix
#'
singleDNB <- function(
    data,
    module_list # module_list is a list
) {
    # prepare data
    data_tmp <- t(data)
    # data_tmp <- data.frame(data_tmp, check.names = FALSE)
    all_genes <- rownames(data)

    # base computation
    sd_all <- apply(data_tmp, 2, sd)
    mean_all <- apply(data_tmp, 2, mean)
    cv_all <- sd_all / mean_all

    # define result
    score_mtx <- matrix(
        0, nrow = 6, ncol = length(module_list),
        dimnames = list(
            c("Mean", "Sd", "cv", "pcc_in", "pcc_out", "score"),
            names(module_list)
        )
    )

    score_mtx_tmp <- sapply(
        module_list,
        function(module_genes) {
            dnb_in <- module_genes
            dnb_out <- all_genes[which(!all_genes %in% dnb_in)]
            tab_in <- data_tmp[, dnb_in, drop =F]
            tab_out <- data_tmp[, dnb_out, drop =F]
            size <- length(dnb_in)

            cv <- mean(cv_all[dnb_in], na.rm = T) # if mean is zero, cv is not defined, so do not consider this one.
            Sd <- mean(sd_all[dnb_in])
            Mean <- mean(mean_all[dnb_in])
            pcc_in <- mean(abs(as.dist(mycor(tab_in))), na.rm = T)
            pcc_out <- mean(abs(mycor(tab_in, tab_out)), na.rm = T)
            score <- ifelse(
                pcc_out == 0,
                0,
                sqrt(size) * Sd * pcc_in / pcc_out
            )

            c(Mean, Sd, cv, pcc_in, pcc_out, score)
        }
    )
    colnames(score_mtx_tmp) <- colnames(score_mtx)
    rownames(score_mtx_tmp) <- rownames(score_mtx)
    score_mtx <- score_mtx_tmp

    return(score_mtx)
}

#' Get the max Ranking of exactly group in S3:DNB_output after function DNBfilter
#'
#' @param object S3:DNB_output
#' @param group which sub-group of S3-DNB_output to extract
#' @param ... not use
#'
#' @return numeric value
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
                stop("ERROR group input! Should be one of meta_levels!")
    }

    sum(grepl(paste0("^", group, "_"), object[[1]]@result@resource))
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
#'
resultAllExtract.DNB_output <- function(
    object,
    group,
    slot,
    mess = TRUE,
    ...
) {
    match.arg(arg = slot, choices = c("pre_result", "result"), several.ok = FALSE)
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

#' Extract the score information from object (S3:DNB_output)
#'
#' @param object S3:DNB_output
#' @param ranking the ranking of exactly module, default 1 if NULL
#' @param group which group to select module, default random selected if NULL
#' @param ... not use
#'
#' @return score data.frame
#'
ScoreExtract.DNB_output <- function(
    object,
    ranking = NULL,
    group = NULL,
    ...
) {
    all_group <- names(object)
    if (is.null(group)) {
        message("Randomly select group...")
        group <- sample(all_group, 1)
    } else {
        if (!group %in% all_group)
            if (group != "USER_CUSTOMIZED")
                stop("ERROR group input! Should be one of meta_levels!")
    }

    max_lenModule <- getMaxRanking(object, group = group)
    if (is.null(ranking)) {
        ranking <- 1
    } else {
        if (ranking < 1 | ranking > max_lenModule)
            stop("ERROR ranking input! Should be between 1 and ", max_lenModule, "!")
    }

    data <- object[[1]]@result
    if (length(data@rank) == 0)
        stop("No module found! Please check or run DNBfilter() !")

    cat("Use group=", group, " and ranking=", ranking, "\n", sep = "")

    index_group <- which(data@resource == group)
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
