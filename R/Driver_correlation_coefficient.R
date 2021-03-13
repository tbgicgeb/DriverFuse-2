#' Driver obj
#' @description plotting driver fusion gene domain in cancer
#'
#' @param x Driver fusion gene
#'
#' @return plot of domain
#' @export
#'
#' @examples
Driver_correlation_coefficient = function(x){
  y = strsplit(x, "-")
  m = as.character(unlist(y))
  library(org.Hs.eg.db)
  library(drawProteins)
  gene <- m[1]
  gene.dt <- gene.track.view(symbol = gene, plot=FALSE, genome.v = "hg19")
  start <- min(gene.dt@data$txStart) - 200000
  stop <- max(gene.dt@data$txEnd) + 50000
  chr <- gene.dt@data$chrom[1]
  sampleids <- sort(results_svc@disruptSamples[[gene]])
  gene1 <- m[2]
  gene.dt1 <- gene.track.view(symbol = gene1, plot=FALSE, genome.v = "hg19")
  start1 <- min(gene.dt1@data$txStart) - 200000
  stop1 <- max(gene.dt1@data$txEnd) + 50000
  chr1 <- gene.dt1@data$chrom[1]
  cnvdat = cnv@data
  svcdat = svc@data
  genegr <- with(data.frame(chr,start,stop), GRanges(chr, IRanges(start=start, end=stop)))
  sv1gr = with(svcdat, GRanges(chrom1, IRanges(start=pos1, end=pos1)))
  sv2gr = with(svcdat, GRanges(chrom2, IRanges(start=pos2, end=pos2)))
  sv_hits1 = GenomicAlignments::findOverlaps(sv1gr,genegr)
  sv_hits2 = GenomicAlignments::findOverlaps(sv2gr,genegr)
  svtab1 <- svcdat[sort(unique(c(queryHits(sv_hits1),queryHits(sv_hits2)))),]
  genegr <- with(data.frame(chr1,start1,stop1), GRanges(chr, IRanges(start=start1, end=stop1)))
  sv1gr = with(svcdat, GRanges(chrom1, IRanges(start=pos1, end=pos1)))
  sv2gr = with(svcdat, GRanges(chrom2, IRanges(start=pos2, end=pos2)))
  sv_hits1 = GenomicAlignments::findOverlaps(sv1gr,genegr)
  sv_hits2 = GenomicAlignments::findOverlaps(sv2gr,genegr)
  svtab2 <- svcdat[sort(unique(c(queryHits(sv_hits1),queryHits(sv_hits2)))),]
  genegr <- with(data.frame(chr,start,stop), GRanges(chr, IRanges(start=start, end=stop)))
  seg1br  = with(cnvdat, GRanges(chrom, IRanges(start=start, end=start)))
  seg2br  = with(cnvdat, GRanges(chrom, IRanges(start=end, end=end)))
  seg_hits1 = GenomicAlignments::findOverlaps(seg1br,genegr)
  seg_hits2 = GenomicAlignments::findOverlaps(seg2br,genegr)
  segbrk_first <- cnvdat[sort(unique(c(queryHits(seg_hits1),queryHits(seg_hits2))))]
  genegr <- with(data.frame(chr1,start1,stop1), GRanges(chr, IRanges(start=start1, end=stop1)))
  seg1br  = with(cnvdat, GRanges(chrom, IRanges(start=start, end=start)))
  seg2br  = with(cnvdat, GRanges(chrom, IRanges(start=end, end=end)))
  seg_hits1 = GenomicAlignments::findOverlaps(seg1br,genegr)
  seg_hits2 = GenomicAlignments::findOverlaps(seg2br,genegr)
  segbrk_second <- cnvdat[sort(unique(c(queryHits(seg_hits1),queryHits(seg_hits2))))]
  cv1 <- validate.cnv(segbrk_first)
  cv2 <- validate.cnv(segbrk_second)
  sv1 <- validate.svc(svtab1)
  sv2 <- validate.svc(svtab2)
  svc_breaks  <- svc.breaks(sv1)
  cnv_breaks  <- cnv.breaks(cv1,verbose=FALSE)
  if (nrow(segbrk_first) > 0 ) {
    cnv_breaks  <- cnv.breaks(cv1,verbose=FALSE)
  } else {
    cnv_breaks <- 0
  }
  dat1 <- as.data.frame(log2(1+cbind(svc_breaks@burden,cnv_breaks@burden)))
  svc_breaks <- svc.breaks(sv2)
  if (nrow(segbrk_second) > 0 ) {
    cnv_breaks  <- cnv.breaks(cv1,verbose=FALSE)
  } else {
    cnv_breaks@burden <- 0
  }
  dat2 <- as.data.frame(log2(1+cbind(svc_breaks@burden,0)))
  z = nrow(cnvdat)/nrow(svcdat)
  dat3 = as.data.frame(c("1",z))
  dat3 = as.data.frame(t(dat3))
  corr_coff = rbind(dat1,dat2,dat3)
  head(corr_coff)
  before_correction<- c(as.numeric(corr_coff$V1))
  after_correction<- c(as.numeric(corr_coff$V2))

  #calculating correlation
  corr_coef <- abs(cor(before_correction, after_correction))


  plot1 <- ggplot2::ggplot(corr_coff, ggplot2::aes(x=before_correction, y=after_correction))+
    ggplot2::geom_point(color = "darkorange") + ggplot2::geom_smooth(method = "lm", colour = "purple") +
    ggplot2::labs(x = "Structural Variation",
                  y = "Copy Number Variation",
                  title = "Scatterplot showing Correlation between SV and CNV for fusion gene",
                  subtitle = paste("Pearson Correlation coefficient", "=", corr_coef, sep= " "))+
    ggplot2::theme(

      #adjusting axis lines and text
      axis.line.x = ggplot2::element_line(size =0.75),
      axis.line.y = ggplot2::element_line(size =0.75),
      axis.text.x = ggplot2::element_text(angle = 90, size=8.5, colour ="black"),
      axis.text.y = ggplot2::element_text(size=8.5, colour ="black"),

      #Center align the title
      plot.title = ggplot2::element_text(face = "bold"),

      # Remove panel border
      panel.border = ggplot2::element_blank(),

      # Remove panel grid lines
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),

      # Remove panel background
      panel.background = ggplot2::element_blank(),

      # Add axis line
      axis.line = ggplot2::element_line(colour = "grey")

    )

  plot(plot1)
}
