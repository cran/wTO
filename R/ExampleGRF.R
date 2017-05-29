#' @title ExampleGRF
#' @aliases ExampleGRF
#' @description ExampleGRF data.frame containing data.frame containing names of GRFs.
#' @format data.frame 184 lines, 1 column.
#' @usage data(ExampleGRF)
#' @name ExampleGRF
#' @docType data
#' @keywords datasets

ExampleGRF = read.table("./data/ExampleGRF.txt", stringsAsFactors = F)
ExampleGRF = base::unique(ExampleGRF)
