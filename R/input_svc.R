#' Preparing input files
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples

input_svc = function(x) {
  library(svpluscnv)
  x = as.data.frame(x)
  x = remove.factors(x)
  svc = validate.svc(x)
  return(svc)
}
