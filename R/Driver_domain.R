#' Driver domain
#' @description plotting driver fusion gene domain in cancer
#'
#' @param x Driver fusion gene
#'
#' @return plot of domain
#' @export
#'
#' @examples
Driver_domain = function(x){
  y = strsplit(x, "-")
  m = as.character(unlist(y))
  library(org.Hs.eg.db)
  library(drawProteins)
  z <- mapIds(org.Hs.eg.db, keys = m[1], keytype = "SYMBOL", column="UNIPROT")
  rel_json <- drawProteins::get_features(z)
  rel_data <- drawProteins::feature_to_dataframe(rel_json)
  p <- draw_canvas(rel_data)
  p <- draw_chains(p, rel_data, label_size = 2.5)
  p <- draw_domains(p, rel_data)
  pdf("rplot.pdf")
  plot(p)
  z1 <- mapIds(org.Hs.eg.db, keys = m[2], keytype = "SYMBOL", column="UNIPROT")
  rel_json1 <- drawProteins::get_features(z1)
  rel_data1 <- drawProteins::feature_to_dataframe(rel_json1)
  p1 <- draw_canvas(rel_data1)
  p1 <- draw_chains(p1, rel_data1, label_size = 2.5)
  p1 <- draw_domains(p1, rel_data1)
  plot(p1)
  dev.off()
}
