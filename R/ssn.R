#' define S4 object of SSNnetwork
#'
#' @slot module matrix
#' @slot network list
#' @slot scores list
#'
#' @return S4 object, SSNnetwork
#'
mynew_network <- setClass(
    Class = "SSNnetwork",
    slots = list(
        module = "matrix",
        network = 'list',
        scores = "list"
    ),
    sealed = TRUE
)

#' define S4 object of SSNscore
#'
#' @slot meta.data data.frame
#' @slot meta character
#' @slot K numeric
#' @slot freq numeric
#' @slot name character
#' @slot module list
#' @slot score data.frame
#'
#' @return S4 object, SSNscore
#'
mynew_ssnscore <- setClass(
    Class = "SSNscore",
    slots = list(
        meta.data = "data.frame",
        meta = "character",
        K = "numeric",
        freq = "numeric",
        name = "character",
        module = 'list',
        score = "data.frame"
    ),
    sealed = TRUE
)

#' define S4 object of SSN
#'
#' @slot network SSNnetwork
#' @slot score SSNscore
#'
#' @return S4 object, SSN
#'
mynew_ssn <- setClass(
    Class = "SSN",
    slots = list(
        network = "SSNnetwork",
        score = 'SSNscore'
    ),
    sealed = TRUE
)

#' sample Expression Deviation
#'
#' @param Xi gene expression from single-sample
#' @param Xref gene expression from reference-samples
#' @param scale whether to scale sED by the standard deviation of reference expression
#'
sED <- function(Xi, Xref, scale = TRUE) 
    abs(Xi - mean(Xref)) / (sd(Xref) ** scale)

#' statistical test for delta PCC
#'
#' @param delta delta PCC
#' @param r0 reference PCC
#' @param n0 reference size
#' @param two_sides test of two-sides (default) or greater?
#'
#' @return p_value
#'
deltaPCC.test <- function(delta, r0, n0, two_sides = TRUE) {
    # delta: delta r
    # r0: r of reference
    # n0: size of reference
    z <- abs(delta) / (1- r0 ^ 2) * (n0 - 1)
    p <- pnorm(z, lower.tail = FALSE) * (1 + two_sides)
    p
}

#' @title NetworkCreate
#'
#' @description Compute genes network by PCC
#'
#' @param data gene expression matrix (or data.frame), with col-sample and row-gene
#' @param method correlation method, pearson (PCC, default) or spearman
#' @param p_require whether to test and compute p_value? 
#'
#' @return a list including matrix r and/or P 
#' @export
#'
NetworkCreate <- function(
    data, 
    method = c('pearson','spearman'), 
    p_require = FALSE
) {
    match.arg(method)
    if (!is(data, "Matrix") && !is.matrix(data) && !is.data.frame(data))
        stop("ERROR data input! Should be matrix or data.frame!")

    # cor() is compute correlations between cols.
    # Here feature's network (gene's) is required.
    # Always the data - gene expression matrix, 
    # col-sample, row-gene
    # So transpose data
    data_t <- t(data)

    if (p_require) {
        # stats::cor() only output the correlation value
        # without p value
        # So use Hmisc::rcorr() to compute both two
        cor_res <- rcorr(data_t, type = method)

        feature_network <- cor_res[c('r', 'P')]
    } else {
        cor_mtx <- cor(data_t, method = method)
        feature_network <- list(r = cor_mtx)
    }

    return(feature_network) # df
}

#' @title refNetworkCreate
#'
#' @description Compute genes network from reference samples
#'
#' @param ref_gex gene expression matrix of reference
#' @param method correlation method, pearson (PCC, default) or spearman
#' @param p_thre p_value threshold for correlation 
#'
#' @return a list including matrix r, P and r_sig
#' @export
#'
refNetworkCreate <- function(ref_gex, method = 'pearson', p_thre = 1e-3) {
    refNet <- NetworkCreate(ref_gex, method, p_require = TRUE)

    cor_mtx_sig <- refNet$r
    # set not significant ones as NA 
    cor_mtx_sig[refNet$P > p_thre] <- NA

    refNet[['r_sig']] <- cor_mtx_sig

    return(refNet)
}

#' @title flattenCorrMatrix
#'
#' @description convert correlation matrix into data.frame
#' @description if the element (edge) is NA, this feature couple (nodes) will be ignored
#'
#' @param cor_mtx correlation matrix
#'
#' @return a data.frame including network nodes and edges (feature1/2, cor)
#' @export
#'
flattenCorrMatrix <- function(cor_mtx) {
    rnames <- rownames(cor_mtx)
    lt <- lower.tri(cor_mtx)
    rindex <- row(cor_mtx)[lt]
    cindex <- col(cor_mtx)[lt]
    
    df <- data.frame(
        feature1 = rnames[rindex], 
        feature2 = rnames[cindex], 
        cor = cor_mtx[lt]
    )

    # return significant elements
    df[!is.na(df$cor), ]
}

#' @title PCCsum
#'
#' @description get all PCC values of exact feature from feature-coupled network
#' @description if the element (edge) is NA, this feature couple (nodes) will be ignored
#'
#' @param feature feature name (node)
#' @param network network data.frame (after flattenCorrMatrix())
#' @param ignore_idx index where it should be ignored (use when finding second-neighbors)
#'
#' @return a list including pcc and indices
#' @export
#'
PCCsum <- function(feature, network, ignore_idx = NULL) {
    idx_in1 <- network$feature1 == feature
    idx_in2 <- network$feature2 == feature
    if (!is.null(ignore_idx)) {
        idx_in1 <- idx_in1 & !ignore_idx
        idx_in2 <- idx_in2 & !ignore_idx
    }
    idx_in <- idx_in1 | idx_in2
    pcc_in <-  sum(abs(network$cor_delta[idx_in]))

    list(
        pcc = pcc_in, 
        idx = idx_in, 
        idx1 = idx_in1, 
        idx2 = idx_in2
    )
} 

#' @title SSNcompute
#'
#' @description get all PCC values of exact feature from feature-coupled network
#' @description if the element (edge) is NA, this feature couple (nodes) will be ignored
#'
#' @param ss_gex gene expression matrix of single sample(s)
#' @param ref_gex gene expression matrix of reference
#' @param delta_p p_value threshold for delta PCC test
#' @param ref_p p_value threshold for correlation test, when building gene network from reference
#' @param nFirstOrderNeighbor the least number of first-order neighbor(s) for each gene, default 3
#' @param nSecondOrderNeighbor the least number of second-order neighbor(s) for each gene, default 1
#' @param sED_scale whether to scale sED by standard deviation
#' @param save_dir the directort for saving SSN network data.frame, default NULL (not save)
#' @param quiet do not print output of process during calculation (against verbose), default FALSE
#'
#' @return S4::SSN object, output will be stored in object@network
#' @export
#'
SSNcompute <- function(
    ss_gex, 
    ref_gex, 
    delta_p = 0.1, 
    ref_p = 0.05, 
    nFirstOrderNeighbor = 3, 
    nSecondOrderNeighbor = 1, 
    sED_scale = TRUE,
    save_dir = NULL,
    quiet = FALSE
) {
    save <- FALSE
    if (!is.null(save_dir)) {
        if (!dir.exists(save_dir))
            stop("Directory of saving files for Single Sample Network didn't exist!")
        setwd(save_dir)
        save <- TRUE
    }
    if (is.data.frame(ref_gex))
        ref_gex <- as.matrix(ref_gex)
    if (is.data.frame(ss_gex))
        ss_gex <- as.matrix(ss_gex)
    N_FON <- nFirstOrderNeighbor
    N_SON <- nSecondOrderNeighbor
    N_ref <- ncol(ref_gex)
    N_ss <- ncol(ss_gex)

    # check ss_gex
    overlapGene_withRef <- intersect(
        rownames(ss_gex), 
        rownames(ref_gex)
    )
    N_overlap <- length(overlapGene_withRef)
    if (N_overlap < 2)
        stop("No enough features from ss_gex can be found in ref_gex!")
    mycat(
        "Find total of ", N_overlap, " genes shared between References & Samples!\n", 
        sep = "", 
        quiet = quiet
    )
    ref_gex <- ref_gex[overlapGene_withRef, , drop = FALSE]
    ss_gex <- ss_gex[overlapGene_withRef, , drop = FALSE]

    # create reference network
    mycat("Now create Reference Network ...\n", quiet = quiet)
    refNet_list <- refNetworkCreate(
        ref_gex, 
        method = 'pearson', 
        p_thre = ref_p
    )
    refNet <- refNet_list[['r_sig']]
    refNet_df <- flattenCorrMatrix(refNet)

    # calculation process for single sample network
    mycat(
        "Now compute Single-Sample Networks total of ", N_ss, " ...\n", 
        sep = "", 
        quiet = quiet
    )
    if (!quiet) {
        iii <- 0
        pb <- utils::txtProgressBar(style = 3)
    }

    ssn_module <- matrix(
        0, ncol = ncol(ss_gex), nrow = nrow(ss_gex), 
        dimnames = dimnames(ss_gex)
    )
    ssn_scores <- list()
    ssn_network <- list()
    
    for (ss_name in colnames(ss_gex)) {
        # get single sample
        ss_vec <- ss_gex[, ss_name, drop = FALSE]

        # create network for reference + single-sample
        ss_ref_gex <- cbind(ref_gex, ss_vec)
        ssNet <- NetworkCreate(
            ss_ref_gex, 
            method = 'pearson', 
            p_require = FALSE
        )[['r']]

        # get delta network for single-sample
        deltaNet <- refNet_df
        deltaNet$cor_delta <- flattenCorrMatrix(
            ssNet - refNet # delta pcc
        )$cor
        # use z-test to test delta-pcc
        deltaNet$delta_pvalue <- apply(
            deltaNet[, c('cor_delta', 'cor')], 
            1, 
            function(x) deltaPCC.test(
                delta = x[1], 
                r0 = x[2], 
                n0 = N_ref
            )
        )
        # select the significant edges
        deltaNet_sig <- deltaNet[deltaNet$delta_pvalue <= delta_p, ]
        ssn_network[[ss_name]] <- deltaNet_sig

        # save the result of single sample network
        if (save) {
            file_tmp <- paste0("SSN_score_", ss_name, ".csv")
            write.csv(
                deltaNet_sig, 
                file = file_tmp, 
                col.names = T, 
                row.names = T
            )
        }

        # get the genes from significant edges
        module <- unique(c(
            deltaNet_sig$feature1, 
            deltaNet_sig$feature2
        ))
        # ssn_module[module, ss_name] <- 1

        # compute ssn score of genes from every sample based on DNB theory
        blank_ssn <- c(na = 1, pcc_in = 0, pcc_out = 0, sED_in = 0, CI = 0)
        blank_n2 <- c(pcc = 0, n = 0)
        ssn_score <- sapply(
            module, 
            function(neighbor_1) {
                # compute pcc_in (first-order neighbors)
                pcc_1n <- PCCsum(
                    feature = neighbor_1, 
                    network = deltaNet_sig
                )
                nFON <- sum(pcc_1n$idx)
                if (nFON < N_FON)
                    return(blank_ssn)
                pcc_in <- pcc_1n$pcc / nFON

                # compute pcc_out (second-order neighbors)
                neighbors_2 <- c(
                    deltaNet_sig$feature2[pcc_1n$idx1], 
                    deltaNet_sig$feature1[pcc_1n$idx2]
                )
                pcc_2n <- sapply(
                    neighbors_2, 
                    function(neighbor_2) {
                        pccout <- PCCsum(
                            feature = neighbor_2, 
                            network = deltaNet_sig, 
                            ignore_idx = pcc_1n$idx
                        )
                        nSON <- sum(pccout$idx)
                        if (nSON < N_SON)
                            blank_n2
                        else
                            c(pcc = pccout$pcc, n = nSON)
                    }
                )
                if (sum(pcc_2n['n',] == 0))
                    # if pcc_2n_n = 0, that means second-order neighbors didn't exist
                    return(blank_ssn)
                pcc_out <- sum(pcc_2n['pcc',]) / sum(pcc_2n['n',])
                
                # compute sED_in (first-order neighbors)
                sED_in <- mean(sapply(
                    c(neighbor_1, neighbors_2), 
                    function(y) sED(
                        Xi = ss_gex[y, ss_name], 
                        Xref = ref_gex[y, ],
                        scale = sED_scale
                    )
                ))

                # compute CI based on DNB theory
                CI_score <- CI(
                    n = 1, # only one core gene
                    Sd = sED_in, 
                    pcc_in = pcc_in, 
                    pcc_out = pcc_out
                )

                # return
                c(
                    na = 0, # whether to pass 
                    pcc_in = pcc_in, 
                    pcc_out = pcc_out, 
                    sED_in = sED_in, 
                    CI = CI_score
                )
            }
        )

        # remove the ones who didn't fit requirement
        ssn_score_sig <- ssn_score[-1, ssn_score['na', ] == 0, drop = FALSE]
        # generate the module matrix
        module_sig <- colnames(ssn_score_sig)
        ssn_module[module_sig, ss_name] <- 1

        # order by CI
        idx_dec <- order(ssn_score_sig['CI', ], decreasing = T)
        ssn_scores[[ss_name]] <- t(
            ssn_score_sig[, idx_dec]
        )

        if (!quiet) {
            iii <- iii + 1
            utils::setTxtProgressBar(pb, iii / N_ss)
        }
    }
    close(pb)

    ssn <- mynew_ssn(
        network = mynew_network(
            module = ssn_module, 
            network = ssn_network, 
            scores = ssn_scores
        ),
        score = mynew_ssnscore()
    )
    return(ssn)
}

#' @title selectK
#'
#' @description select top K genes for every sample sorted by CI values
#'
#' @param object S4::SSN
#' @param K the top K genes by order of CI
#' 
#' @return matrix with element 1/0 (TRUE/FALSE), colnames samples and rownames genes
#' @export
#'
selectK <- function(object, K) {
    # check input
    if (!is(object, "SSN"))
        stop("Input object must be class of SSN! Please run SSNcompute() first!")

    # get data
    ssn_list <- object@network@scores
    ssn_module0 <- object@network@module

    # define module matrix
    ssn_module <- matrix(
        0, ncol = ncol(ssn_module0), nrow = nrow(ssn_module0), 
        dimnames = dimnames(ssn_module0)
    )

    # select genes by top K values of CI
    for (ss_name in names(ssn_list)) {
        ssn_x <- ssn_list[[ss_name]]
        genes <- if (K <= nrow(ssn_x))
            rownames(ssn_x)[1:K]
        else
            rownames(ssn_x)
        ssn_module[genes, ss_name] <- 1
    }

    # return module
    return(ssn_module)
}

#' @title SSNscore_default
#'
#' @description compute the aggregated SSN score for each group 
#'
#' @param object S4::SSN
#' @param K an integer, the top K genes by order of CI
#' @param freq percentage for selection of each group's genes.
#'      That is, n_sample = freq * #repeats in one group,
#'      and n_gene = (in this group how many repeats enriched this gene).
#'      If n_gene >= n_sample, this gene is selected for this group
#' @param meta group information for all samples (repeats), 
#'      a data.frame with 1-col, or NULL (default).
#'      If NULL, object@score@meta.data should have any column and the first will be used 
#' @param method method for aggregating the gene expression to calculate the score of each sample,
#'      by average (mean, default) or Summation (sum)
#' 
#' @return S4::SSN object, result will be stored in object@score
#' @export
#'
SSNscore_default <- function(
    object, 
    K, 
    freq = 0.5, 
    meta = NULL, 
    method = c("mean", "sum")
) {
    # check input
    if (!is(object, "SSN"))
        stop("Input object must be class of SSN! Please run SSNcompute() first!")
    match.arg(method)
    if (length(K) != 1)
        stop("K should be one integer! If a vector, pls use SSNscore().")
    if (freq > 1 | freq <= 0)
        stop("freq should be between 0 and 1!")
    
    # get data
    ssn_list <- object@network@scores
    meta_data <- object@score@meta.data
    
    ssn_module0 <- object@network@module
    score_score <- object@score@score
    
    k_vec <- object@score@K
    freq_vec <- object@score@freq
    meta_vec <- object@score@meta
    name_vec <- object@score@name
    n_meta <- ncol(meta_data)

    # check meta
    if (is.null(meta)) {
        if (n_meta == 0)
            stop("meta is required! (Otherwise object@score@meta.data should exist!)")
        meta <- meta_data[, 1, drop = FALSE]
        meta_name <- colnames(meta)
        warning("Use meta.data named ", meta_name, "!")
    } else {
        if (n_meta == 0){
            object@score@meta.data <- meta
        } else {
            if (!colnames(meta) %in% colnames(meta_data))
                object@score@meta.data <- cbind(meta_data, meta)
        }
    }
    if (nrow(meta) == ncol(ssn_module0) & 
        all(rownames(meta) %in% colnames(ssn_module0))
    ) {
        meta <- meta[colnames(ssn_module0), , drop = FALSE]
    } else {
        stop("ERROR input of meta!")
    }

    meta_name <- colnames(meta)
    tmp_name <- paste(meta_name, K, freq, sep = "_")

    # create module matrix (sample-gene) by K
    module <- selectK(object, K = K)

    # create matrix of relation between meta and genes
    metaModule <- apply(module, 1, function(x)
        tapply(
            x, 
            meta[, 1], 
            function(y) mean(y) >= freq
        )
    )
    if (length(unique(meta[, 1])) == 1) 
        metaModule <- matrix(
            metaModule, 
            byrow = T, nrow = 1, 
            dimnames = list(
                unique(meta[, 1]), 
                names(metaModule)
            )
        )
    
    # compute aggregated SSN score
    SSNscore_vec <- tapply(
        ssn_list, 
        meta, 
        function(ssn_meta) {
            ssn_name <- names(ssn_meta)[1]
            meta_name <- meta[ssn_name, 1]
            genes <- metaModule[meta_name, ]
            genes <- names(genes)[genes]
            ngenes <- length(genes)

            mean(sapply(
                ssn_meta, 
                function(ssn_df) {
                    genes_overlap <- intersect(rownames(ssn_df), genes)
                    s <- sum(ssn_df[genes_overlap, 'CI'])
                    if (method == 'mean')
                        s / ngenes
                    else
                        s
                }
            ))
        }
    )
    SSNscore_df <- as.data.frame(t(data.frame(K = SSNscore_vec)))
    rownames(SSNscore_df) <- tmp_name
    
    # merge the output into raw object
    if (!tmp_name %in% name_vec) {
        object@score@module[[tmp_name]] <- module
        object@score@K <- c(k_vec, K)
        object@score@freq <- c(freq_vec, freq)
        object@score@meta <- c(meta_vec, meta_name)
        object@score@name <- c(name_vec, tmp_name)
        object@score@score <- if (length(name_vec) == 0) 
            SSNscore_df
        else
            rbind(score_score, SSNscore_df)
    }

    # return
    return(object)
}

#' @title SSNscore
#'
#' @description compute the aggregated SSN score for each group 
#'
#' @param object S4::SSN
#' @param K a vector of K, the top K genes by order of CI
#' @param freq percentage for selection of each group's genes.
#'      That is, n_sample = freq * #repeats in one group,
#'      and n_gene = (in this group how many repeats enriched this gene).
#'      If n_gene >= n_sample, this gene is selected for this group
#' @param meta group information for all samples (repeats), 
#'      a data.frame with 1-col, or NULL (default).
#'      If NULL, object@score@meta.data should have any column and the first will be used 
#' @param method method for aggregating the gene expression to calculate the score of each sample,
#'      by average (mean, default) or Summation (sum)
#' 
#' @return S4::SSN object, result will be stored in object@score
#' @export
#'
SSNscore <- function(
    object, 
    K, 
    freq = 0.5, 
    meta = NULL, 
    method = c("mean", "sum")
) {
    fun <- function(object, Kx)
        SSNscore_default(
            object = object, 
            K = Kx, 
            freq = freq, 
            meta = meta, 
            method = method
        )

    if (length(K) == 1) {
        fun(object, K)
    } else {
        for (k in K) {
            object <- fun(object, k)
        }
        object
    }
}
