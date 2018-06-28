#' @title metagenomics_abundance
#' @aliases metagenomics_abundance
#' @description metagenomics_abundance 
#' @format data.frame from The USC Microbial Observatory. The data is public available at <https://www.ebi.ac.uk/metagenomics/projects/ERP013549>
#' @usage data('metagenomics_abundance')
#' @name metagenomics_abundance

metagenomics_abundance  = utils::read.table("./data/metagenomics_abundance.txt", sep = '\t')

