#' DNBplot for S3:DNB_output
#'
#' @param object DNB_output
#' @param ranking ranking
#' @param group group
#' @param show whether to show
#' @param save_pdf whether to save pdf file
#' @param file_prefix the file prefix if save_pdf
#' @param meta_levels the order of meta-group in the plots, default levels(df_score$Names) if be null
#' @param ... not use
#'
#' @export DNBplot.DNB_output
#' @export
#'
DNBplot.DNB_output <- function(
    object,
    ranking = NULL,
    group = NULL,
    show = TRUE,
    save_pdf = FALSE,
    file_prefix = NULL,
    meta_levels = NULL
    ...
) {
    df_score <- ScoreExtract(
        object = object,
        ranking = ranking,
        group = group
    )

    plotScore(
        df_score = df_score,
        file_prefix = file_prefix,
        show = show,
        save_pdf = save_pdf,
        meta_levels = meta_levels
    )
}

#' DNBplot for data.frame
#'
#' @param object data.frame of scores
#' @param ranking not use
#' @param group not use
#' @param show whether to show
#' @param save_pdf whether to save pdf file
#' @param file_prefix the file prefix if save_pdf
#' @param meta_levels the order of meta-group in the plots, default levels(df_score$Names) if be null
#' @param ... not use
#'
#' @export DNBplot.data.frame
#' @export
#'
DNBplot.data.frame <- function(
    object,
    ranking = NULL,
    group = NULL,
    show = TRUE,
    save_pdf = FALSE,
    file_prefix = NULL,
    meta_levels = NULL
    ...
) {
    plotScore(
        df_score = object,
        file_prefix = file_prefix,
        show = show,
        save_pdf = save_pdf,
        meta_levels = meta_levels
    )
}


#' @title plotScore
#' @description main function for plot DNB, same to DNBplot
#'
#' @param df_score score data.frame, colnames includes 'SCORE','PCC_IN','PCC_OUT','SD','Names'(group)
#' @param show whether to plot in console
#' @param save_pdf whether to save into a pdf file
#' @param file_prefix if save_pdf, the file prefix name
#' @param meta_levels the order of meta-group in the plots, default levels(df_score$Names) if be null
#' @param ... parameters to ggsave()
#'
#' @return plot(if ggplot2/gridExtra installed, return ggplot object) or pdf
#' @export
#'
plotScore <- function(
    df_score,
    show = TRUE,
    save_pdf = FALSE,
    file_prefix = NULL,
    meta_levels = NULL,
    ...
) {
    if (!is.data.frame(df_score))
        stop("ERROR df_score! Should be a data.frame!")
    if (!is.null(meta_levels)) {
        if (!all(meta_levels %in% df_score$Names))
            stop("ERROR meta_levels! Should be equal to unique(df_score$Names)!")
        df_score$Names <- factor(df_score$Names, levels = meta_levels)
    }

    pp <- function(df, y) {
        Names <- "Names"
        df_score <- df[, c(Names, y)]
        colnames(df_score)[2] <- 'y'

        ggplot(df_score, aes(x = Names, y = y, group = 1)) +
        geom_point() +
        geom_line() +
        theme_classic(base_size = 15) +
        ggtitle(y) +
        theme(axis.text.x  = element_text(face ="bold", size = 12, color = "black"),
              axis.text.y  = element_text(face ="bold", size = 12, color = "black"),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              plot.title = element_text(hjust = 0.5))
    }
    pdf_file <- "DNBplot.pdf"

    p1 <- pp(df_score, 'SCORE')
    p2 <- pp(df_score, 'PCC_IN')
    p3 <- pp(df_score, 'PCC_OUT')
    p4 <- pp(df_score, 'SD')
    p <- plot_grid(p1, p2, p3, p4)

    if (save_pdf | !show) {
        pdf_file <- ifelse(
            is.null(file_prefix),
            pdf_file,
            paste(file_prefix, pdf_file, sep = "_")
        )
        cat("Plot into ", pdf_file, " ...\n", sep = "")
        ggsave(filename = pdf_file, plot = p, ...)
    }
    
    if (show) {
        return(p)
    }
}
