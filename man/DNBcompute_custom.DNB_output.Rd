% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process.R
\name{DNBcompute_custom.DNB_output}
\alias{DNBcompute_custom.DNB_output}
\title{Compute the Dynamic Network Biomarkers(DNB) model with customized Modules}
\usage{
\method{DNBcompute_custom}{DNB_output}(
  data,
  module_list,
  meta = NULL,
  meta_levels = NULL,
  quiet = FALSE,
  ...
)
}
\arguments{
\item{data}{existing S3:DNB_output object}

\item{module_list}{a customized list of module gene}

\item{meta}{not use}

\item{meta_levels}{not use}

\item{quiet}{do not message}

\item{...}{not use}
}
\value{
S3:DNB_output
}
\description{
Compute the Dynamic Network Biomarkers(DNB) model with customized Modules
}
