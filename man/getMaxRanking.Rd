% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{getMaxRanking}
\alias{getMaxRanking}
\title{getMaxRanking}
\usage{
getMaxRanking(object, group, ...)
}
\arguments{
\item{object}{S3:DNB_output}

\item{group}{which sub-group of S3-DNB_output to extract}

\item{...}{for future use}
}
\value{
numeric value of max Ranking for your group
}
\description{
Get the max Ranking of exactly group in S3:DNB_output after function DNBfilter
}
\examples{
data(data.example)
data(meta.example)
a <- DNBcompute(data.example, meta.example)
b <- DNBfilter(a, ntop = 5)
maxRank <- getMaxRanking(b, group = "C") # get 4 instead of 5
DNBplot(b, ranking = maxRank, group = "C", show = TRUE, save_pdf = FALSE)

}
\author{
Kaiyu Wang, in ChenLab of CAS, Shanghai, China
}
