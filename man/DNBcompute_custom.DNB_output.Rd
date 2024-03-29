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
  force_allgene = TRUE,
  size_effect = TRUE,
  quiet = FALSE,
  ...
)
}
\arguments{
\item{data}{existing S3:DNB_output object}

\item{module_list}{a list of module genes}

\item{meta}{not use}

\item{meta_levels}{not use}

\item{force_allgene}{whether force to use all genes from data, default TRUE}

\item{size_effect}{whether consider the effect of sample size when compute CI of DNB, default TRUE}

\item{quiet}{do not message}

\item{...}{not use}
}
\value{
S3:DNB_output
}
\description{
Compute the Dynamic Network Biomarkers(DNB) model with customized Modules
}
