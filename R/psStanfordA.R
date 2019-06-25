#' The Stanford pregnancy data as \code{\link{phyloseq-class}}.
#'
#' A \code{\link{phyloseq-class}} with 109 taxa from 678 biological samples.
#'
#' @format An object of \code{\link{phyloseq-class}}.
#' \describe{
#'   \item{SubjectID}{IDs of pregnant women}
#'   \item{Preterm}{A factor variable with two levels TRUE or FALSE}
#'   \item{SampleID}{Sample ID for repeated biological sample}
#'   \item{GWColl}{Gestational week at sampling}
#' }
#' @source The data use in the paper: \href{https://www.ncbi.nlm.nih.gov/pubmed/26283357}{where}
#' @references \href{https://www.ncbi.nlm.nih.gov/pubmed/26283357}{Study}
#' @examples
#' \dontrun{sample_data(psStanfordA)
#' otu_table(psStanfordA)}
"psStanfordA"
