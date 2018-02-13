#' Simulated \code{"phyloseq"} object.
#'
#' A \code{"phyloseq"} class object with 50 taxa from 300 biological samples.
#'
#' @format An object of class \code{"phyloseq"}.
#' \describe{
#'   \item{SubjectID}{IDs of individuals}
#'   \item{Preterm}{A factor variable with two levels TRUE or FALSE}
#'   \item{SampleID}{Sample ID for each repeated biological samples}
#'   \item{Time}{Time at repeated measures}
#' }
#' @source Simulated based on the pregnancy data. Paper on \href{http://www.pnas.org/content/112/35/11060.short}{where}
#' @references Pratheepa Jeganathan (Simulated dataset)
#' @examples
#' \dontrun{sample_data(pssim)
#' otu_table(pssim)}
"pssim"
