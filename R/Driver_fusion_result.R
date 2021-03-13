#' Title driver
#'
#' @param x Driver fusion gene
#'
#' @return
#' @export
#'
#' @examples
Driver_fusion_result = function(x){
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
  sampleids1 <- sort(results_svc@disruptSamples[[gene1]])
  svcdat <- svc@data
  cnvdat <- cnv@data
  results_svc <- svc.break.annot(svc, svc.seg.size = 200000, genome.v="hg19",upstr = 100000, verbose=FALSE)
  genegr <- with(data.frame(chr,start,stop), GRanges(chr, IRanges(start=start, end=stop)))
  sv1gr = with(svcdat, GRanges(chrom1, IRanges(start=pos1, end=pos1)))
  sv2gr = with(svcdat, GRanges(chrom2, IRanges(start=pos2, end=pos2)))
  sv_hits1 = GenomicAlignments::findOverlaps(sv1gr,genegr)
  sv_hits2 = GenomicAlignments::findOverlaps(sv2gr,genegr)
  svtab1 <- svcdat[sort(unique(c(queryHits(sv_hits1),queryHits(sv_hits2)))),]
  genegr1 <- with(data.frame(chr1,start1,stop1), GRanges(chr1, IRanges(start=start, end=stop)))
  sv1gr = with(svcdat, GRanges(chrom1, IRanges(start=pos1, end=pos1)))
  sv2gr = with(svcdat, GRanges(chrom2, IRanges(start=pos2, end=pos2)))
  sv_hits1 = GenomicAlignments::findOverlaps(sv1gr,genegr1)
  sv_hits2 = GenomicAlignments::findOverlaps(sv2gr,genegr1)
  svtab2 <- svcdat[sort(unique(c(queryHits(sv_hits1),queryHits(sv_hits2)))),]
  seg1br  = with(cnvdat, GRanges(chrom, IRanges(start=start, end=start)))
  seg2br  = with(cnvdat, GRanges(chrom, IRanges(start=end, end=end)))
  seg_hits1 = GenomicAlignments::findOverlaps(seg1br,genegr)
  seg_hits2 = GenomicAlignments::findOverlaps(seg2br,genegr)
  segbrk1 <- cnvdat[sort(c(queryHits(seg_hits1),queryHits(seg_hits2)))]
  seg1br  = with(cnvdat, GRanges(chrom, IRanges(start=start, end=start)))
  seg2br  = with(cnvdat, GRanges(chrom, IRanges(start=end, end=end)))
  seg_hits1 = GenomicAlignments::findOverlaps(seg1br,genegr1)
  seg_hits2 = GenomicAlignments::findOverlaps(seg2br,genegr1)
  segbrk2 <- cnvdat[sort(unique(c(queryHits(seg_hits1),queryHits(seg_hits2))))]
  if (nrow(svtab1) > 0 & nrow(svtab2) > 0 & nrow(segbrk2) & nrow(segbrk2) > 0) {
    print("Driver fusion gene")
  } else {
    print("Passenger fusion gene")
  }
}
