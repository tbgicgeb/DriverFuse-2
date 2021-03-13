#' Preparing input files
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples

input_cnv = function(x) {
  library(svpluscnv)
  x = as.data.frame(x)
  x = remove.factors(x)
  cnv = validate.cnv(x)
  return(cnv)
}
