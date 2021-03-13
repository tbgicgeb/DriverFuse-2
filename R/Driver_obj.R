#' Driver obj
#' @description plotting driver fusion gene domain in cancer
#'
#' @param x Driver fusion gene
#'
#' @return plot of domain
#' @export
#'
#' @examples

Driver_obj = function(x){
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
  sv2gr = with(svcdat , GRanges(chrom2, IRanges(start=pos2, end=pos2)))
  sv_hits1 = GenomicAlignments::findOverlaps(sv1gr,genegr)
  sv_hits2 = GenomicAlignments::findOverlaps(sv2gr,genegr)
  svtab1 <- svcdat[sort(unique(c(queryHits(sv_hits1),queryHits(sv_hits2)))),]
  genegr <- with(data.frame(chr1,start1,stop1), GRanges(chr, IRanges(start=start1, end=stop1)))
  sv1gr = with(svcdat, GRanges(chrom1, IRanges(start=pos1, end=pos1)))
  sv2gr = with(svcdat , GRanges(chrom2, IRanges(start=pos2, end=pos2)))
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
  Driver_res1 = list(svtab1,segbrk_first,svtab2,segbrk_second)
  head(Driver_res1)
}
