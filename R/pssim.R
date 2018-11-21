#' Simulated \code{'phyloseq'} object.
#'
#' A \code{'phyloseq'} class object with 50 taxa from 200 biological samples.
#'
#' @format An object of class \code{'phyloseq'}.
#' \describe{
#'   \item{SubjectID}{IDs of individuals}
#'   \item{Group}{A factor variable with two levels Case or Control}
#'   \item{SampleID}{Sample ID for each repeated biological sample}
#'   \item{Time}{Time at repeated biological sample}
#' }
#' @source Simulated based on the pregnancy data. Paper on \href{http://www.pnas.org/content/112/35/11060.short}{where}
#' @references Pratheepa Jeganathan (Simulated dataset)
#' @examples
#' \dontrun{sample_data(psSim)
#' otu_table(psSim)}
"psSim"
