#' Preparing input files
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples

input = function(x) {
  library(svpluscnv)
  x = as.data.frame(x)
  x = remove.factors(x)
  svc = validate.svc(x)
  results_svc <- svc.break.annot(svc, svc.seg.size = 200000, genome.v="hg19",upstr = 100000, verbose=FALSE)
  return(results_svc)
}
