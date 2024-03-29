% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wTO.Complete2.R
\name{wTO.Complete}
\alias{wTO.Complete}
\title{wTO.Complete}
\usage{
wTO.Complete(
  k = 1,
  n = 100,
  Data,
  Overlap = row.names(Data),
  method = "p",
  method_resampling = "Bootstrap",
  pvalmethod = "BH",
  savecor = F,
  expected.diff = 0.2,
  lag = NULL,
  ID = NULL,
  normalize = F,
  plot = T
)
}
\arguments{
\item{k}{Number of threads to be used for computing the weight Topological Overlap. Default is set to 1.}

\item{n}{Number of resamplings, used to compute the empirical distribuitions of the links. Default is set to 100.}

\item{Data}{data.frame containing the count / expression data for the correlation.}

\item{Overlap}{Set of nodes of interest, where the Overlapping weights will be computed.}

\item{method}{Type of the correlation that should be used. "s" / "spearman" will compute the rank spearman correlation, "p" / "pearson" will compute the linear correlation. If no value is given, the default is to use "p".}

\item{method_resampling}{method of the resampling. Bootstrap, BlockBootstrap or Reshuffle. Bootstrap null hypothesis is that the wTO is random, and Reshuffle tests if the wTO is equal to zero.}

\item{pvalmethod}{method to compute the multiple test correction for the pvalue. for more information check the function \code{\link[stats]{p.adjust}}.}

\item{savecor}{T/F if need to save the correlation.}

\item{expected.diff}{Difference expected between the real wTO and resampled wTO By default, it is set to 0.2.}

\item{lag}{time dependency, lag, if you are using the BlockedBootstrap.}

\item{ID}{ID of the samples for the blocked bootstrap (for repeated measures).}

\item{normalize}{T/F Should the data be normalized?}

\item{plot}{T/F Should the diagnosis plot be plotted?}
}
\value{
a list with results.
\itemize{
\item wTO is a data.frame containig the Nodes, the wTO computed using the signed correlations, the pvalue and the adj.pvalue.
\item abs.wTO is a data.frame containig the Nodes, the wTO computed using the absolute correlations, the pvalue and the adj.pvalue.
\item Correlation is a data.frame containing the correlation between all the nodes.
\item Empirical.Quantile quantile values for the empirical distribution.
\item Quantile quantile values for the sample distribution.
}
}
\description{
Compute the wTO and also the bootstraps. Proposed at: 	arXiv:1711.04702
}
\examples{
\dontrun{
# Using spearman rank correlation and bonferroni correction for the pvalues.
wTO.Complete( k =8, n = 1000, Data = Microarray_Expression1,
 Overlap = ExampleGRF$x, method = "s", pvalmethod = "bonferroni")
 # Changing the resampling method to Reshuffle.
wTO.Complete( k =1, n = 1000, Data = Microarray_Expression1,
Overlap = ExampleGRF$x, method_resampling = "Reshuffle")
 # Changing the resampling method to BlockBootstrap, with a lag of 2.
 row.names(metagenomics_abundance) = metagenomics_abundance$OTU
 metagenomics_abundance = metagenomics_abundance[,-1]
wTO.Complete( k =1, n = 1000, Data = metagenomics_abundance, method = "s",
Overlap = row.names(metagenomics_abundance), method_resampling = "BlockBootstrap", lag = 2)
wTO.Complete( k =2, n = 1000, Data = Microarray_Expression1, method = "s",
Overlap = ExampleGRF$x, method_resampling = "BlockBootstrap", ID = rep(1:9,each = 2))
X = wTO.Complete( k =1, n = 1000, Data = Microarray_Expression1,
Overlap = ExampleGRF$x, method = "p", plot = FALSE)
 }
}
\author{
Deisy Morselli Gysi <deisy at bioinf.uni-leipzig.de>
}
