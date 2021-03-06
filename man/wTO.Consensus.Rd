% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wTO.Consensus.R
\name{wTO.Consensus}
\alias{wTO.Consensus}
\title{wTO.Consensus}
\usage{
wTO.Consensus(data)
}
\arguments{
\item{data}{list of data.frame containing the "Node.1", "Node.2" and "wTO".}
}
\description{
Consensus requires a list of data.frame containing the pair of nodes, and the wTO values for all networks that need to be joined. Reference: 	arXiv:1711.04702
}
\examples{
\dontrun{
EXAMPLE =  wTO.Complete( k =1, n = 200, Data = Microarray_Expression1,
                                      Overlap = ExampleGRF$x, method = "p")

# Constructing the consensus network
data = list(data.frame(Node.1 = EXAMPLE$wTO$Node.1,
Node.2 = EXAMPLE$wTO$Node.2,
wto_sig = EXAMPLE$wTO$wTO_sign,
pvalsig = EXAMPLE$wTO$pval_sig),
data.frame(Node.1 = EXAMPLE$wTO$Node.1,
Node.2 = EXAMPLE$wTO$Node.2,
wtoabs = EXAMPLE$wTO$wTO_abs,
pvalabs = EXAMPLE$wTO$pval_abs) )
CONS = wTO.Consensus(data)

}
}
\author{
Deisy Morselli Gysi <deisy at bioinf.uni-leipzig.de>
}
