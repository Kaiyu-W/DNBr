% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process.R
\name{resultAllExtract.DNB_output}
\alias{resultAllExtract.DNB_output}
\title{Get the whole result of slot pre_result or result from S4-DNB_obj and transfer that into a data.frame (matrix)}
\usage{
\method{resultAllExtract}{DNB_output}(object, group, slot, mess = TRUE, ...)
}
\arguments{
\item{object}{the S3-DNB_output}

\item{group}{which sub-group of S3-DNB_output to extract}

\item{slot}{which S4-DNB_obj to extract, pre_result or result}

\item{mess}{whether to message if slot is pre_result, default TRUE}

\item{...}{not use}
}
\value{
data.frame
}
\description{
Get the whole result of slot pre_result or result from S4-DNB_obj and transfer that into a data.frame (matrix)
}
