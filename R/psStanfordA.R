#' The Stanford pregnancy data \code{'phyloseq'} object.
#'
#' A \code{'phyloseq'} class object with 109 taxa from 678 biological samples.
#'
#' @format An object of class \code{'phyloseq'}.
#' \describe{
#'   \item{SubjectID}{IDs of individuals}
#'   \item{SampleID}{IDs of repeated observations}
#'   \item{Preterm}{A factor variable with two levels TRUE or FALSE}
#'   \item{SampleID}{Sample ID for each repeated biological samples}
#'   \item{GWColl}{Gestational week at sampling}
#' }
#' @source The data use in the paper: \href{https://www.ncbi.nlm.nih.gov/pubmed/26283357}{where}
#' @references \href{https://www.ncbi.nlm.nih.gov/pubmed/26283357}{Study}
#' @examples
#' \dontrun{sample_data(psStanfordA)
#' otu_table(psStanfordA)}
"psStanfordA"
