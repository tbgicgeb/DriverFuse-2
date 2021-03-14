
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DriverFuse

<!-- badges: start -->

<!-- badges: end -->

The goal of DriverFuse is to become instrumental in the study of Driver
fusion genes in cancer genomics, enabling the identification of
recurrent complex rearrangements that may pinpoint disease driver
events. Here,we implement the methodological framework into an R package
incorporating multiple functionalities to integrate orthogonal data
types such as SV and CNV and characterize fusion gene in breast cancer
datasets.This toolset allows the research community to easily perform
complex analysis of high throughput profiling data and will support
extensions to further integrate analyses of other types of omics data.

input\_svc creates input structural variant file for a user sample data
that is used for mapping fusion gene to find out whether it is “Driver”
or “Passenger” fusion gene.

input\_cnv creates input copy number variation file for a user sample
data that is used for mapping fusion gene to find out whether it is
“Driver” or “Passenger” fusion gene.

input creates localization of breakpoints with respect to known genes.
The annotation identifies small segmental variants overlapping with
genes.

Driver\_result classifies input fusion gene into “Driver” or “Passenger”
fusion gene based on mapping coordinates in structural variant and copy
number variation.

Driver\_object creates object of input fusion gene based on their
mapping structural variant and copy number variation.

Driver\_correlation\_coefficient plots the correlation coefficient
between mapping structural variant and copy number variation profile for
a given input fusion gene.

Driver\_plot plots mapping structural variants and CNV profile for input
fusion transcripts. It plots the CNV profile and structural variants
that are mapping to genomic coordinates of input fusion genes.

Driver\_domain provides genomic feature annotation tools for driver
fusion gene. Exploring domain level landscape of fusion gene is
important to identify gene role in cancer progression.

## Installation

You can install the released version of DriverFuse from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("DriverFuse")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(DriverFuse)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
library(svpluscnv)
#> 
#> Warning in fun(libname, pkgname): couldn't connect to display ":0"
#> 
#> Attaching package: 'svpluscnv'
#> The following object is masked from 'package:DriverFuse':
#> 
#>     sv.model.view
csv = system.file("extdata", "sample_data.txt", package = "DriverFuse")
csv2 = system.file("extdata", "data2", package = "DriverFuse")
sample_read = function(path) {
  readr::read_csv(path)
}
data = sample_read(csv2)
#> 
#> ── Column specification ────────────────────────────────────────────────────────
#> cols(
#>   sample = col_character(),
#>   chrom = col_character(),
#>   start = col_double(),
#>   end = col_double(),
#>   probes = col_double(),
#>   segmean = col_double(),
#>   uid = col_character()
#> )
data = as.data.frame(data)
head(data)
#>   sample chrom  start    end probes   segmean          uid
#> 1  HCC_1  chr1  10000  65805   6564   1.00964 cnv_BHGaBjvp
#> 2  HCC_1  chr1  60000 110007     16 -20.00000 cnv_XmTSvjCa
#> 3  HCC_1  chr1  60000  94821     32  -7.80227 cnv_tTYjtntc
#> 4  HCC_1  chr1  65805 121611      5  -4.42384 cnv_EiGLqQyo
#> 5  HCC_1  chr1 110007 160015     15 -20.00000 cnv_LcABRqJG
#> 6  HCC_1  chr1 121611 177417     19   0.55161 cnv_EBJnrzXf
data2 = sample_read(csv)
#> 
#> ── Column specification ────────────────────────────────────────────────────────
#> cols(
#>   sample = col_character(),
#>   chrom1 = col_character(),
#>   pos1 = col_double(),
#>   strand1 = col_character(),
#>   chrom2 = col_character(),
#>   pos2 = col_double(),
#>   strand2 = col_character(),
#>   variant_type = col_character()
#> )
data2 = as.data.frame(data2)
head(data2)
#>       sample chrom1    pos1 strand1 chrom2    pos2 strand2 variant_type
#> 1 HCC_BREAST   chr1 1786636       -   chr1 1796052       -          BND
#> 2 HCC_BREAST   chr1 2584747       -   chr1 2614866       +          DUP
#> 3 HCC_BREAST   chr1 2614754       -   chr1 2615807       +          DUP
#> 4 HCC_BREAST   chr1 2583646       -   chr1 2615819       +          DUP
#> 5 HCC_BREAST   chr1 2620899       -   chr1 2622923       +          DUP
#> 6 HCC_BREAST   chr1 2583651       -   chr1 2629307       +          DUP
input_cnv = function(x) {
  x = as.data.frame(x)
  cnv = validate.cnv(x)
  return(cnv)
}
CNV7 = input_cnv(data)
input_svc = function(x) {
  x = as.data.frame(x)
  svc = validate.svc(x)
  return(svc)
}
svc2 = validate.svc(data2)
CNV7 = validate.cnv(data)
svc2 = input_svc(data2)
CNV7  = input_cnv(data)
input = function(x) {
  library(svpluscnv)
  x = as.data.frame(x)
  x = remove.factors(x)
  svc = validate.svc(x)
  results_svc <- svc.break.annot(svc, svc.seg.size = 200000, genome.v="hg19",upstr = 100000, verbose=FALSE)
  return(results_svc)
}
results_svc <- svc.break.annot(svc2, svc.seg.size = 200000, genome.v="hg19",upstr = 100000, verbose=FALSE)
#> 'select()' returned 1:1 mapping between keys and columns
library(GenomicRanges)
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: parallel
#> 
#> Attaching package: 'BiocGenerics'
#> The following objects are masked from 'package:parallel':
#> 
#>     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
#>     clusterExport, clusterMap, parApply, parCapply, parLapply,
#>     parLapplyLB, parRapply, parSapply, parSapplyLB
#> The following objects are masked from 'package:stats':
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from 'package:base':
#> 
#>     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
#>     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
#>     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
#>     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
#>     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
#>     union, unique, unsplit, which, which.max, which.min
#> Loading required package: S4Vectors
#> 
#> Attaching package: 'S4Vectors'
#> The following object is masked from 'package:base':
#> 
#>     expand.grid
#> Loading required package: IRanges
#> Loading required package: GenomeInfoDb
gene <- "PTPRD"
gene.dt <- gene.track.view(symbol = gene, plot=FALSE, genome.v = "hg19")
start <- min(gene.dt@data$txStart) - 200000
stop <- max(gene.dt@data$txEnd) + 50000
chr <- gene.dt@data$chrom[1]
genegr <- with(data.frame(chr,start,stop), GRanges(chr, IRanges(start=start, end=stop))) 
sv1gr = with(data2, GRanges(chrom1, IRanges(start=pos1, end=pos1))) 
sv2gr = with(data2 , GRanges(chrom2, IRanges(start=pos2, end=pos2))) 
sv_hits1 = GenomicAlignments::findOverlaps(sv1gr,genegr)
sv_hits2 = GenomicAlignments::findOverlaps(sv2gr,genegr)
svtab1 <- data2[sort(unique(c(queryHits(sv_hits1),queryHits(sv_hits2)))),]

Driver_result = function(x){
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
    svcdat <- data2
    cnvdat <- CNV7@data
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
    segbrk1 <- cnvdat[sort(unique(c(queryHits(seg_hits1),queryHits(seg_hits2))))]
    seg1br  = with(cnvdat, GRanges(chrom, IRanges(start=start, end=start))) 
    seg2br  = with(cnvdat, GRanges(chrom, IRanges(start=end, end=end))) 
    seg_hits1 = GenomicAlignments::findOverlaps(seg1br,genegr1)
    seg_hits2 = GenomicAlignments::findOverlaps(seg2br,genegr1)
    segbrk2 <- cnvdat[sort(unique(c(queryHits(seg_hits1),queryHits(seg_hits2))))]
    if (nrow(svtab1) > 0 & nrow(svtab2) > 0) {
        print("Driver fusion gene")
    } else {
        print("Passenger fusion gene")
    }
}
Driver_result("PTPRD-BRCA2")
#> Loading required package: AnnotationDbi
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> [1] "Passenger fusion gene"
Driver_result("MYH9-EIF3D")
#> [1] "Driver fusion gene"
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
  cnvdat = CNV7@data
  genegr <- with(data.frame(chr,start,stop), GRanges(chr, IRanges(start=start, end=stop))) 
  sv1gr = with(data2, GRanges(chrom1, IRanges(start=pos1, end=pos1))) 
  sv2gr = with(data2 , GRanges(chrom2, IRanges(start=pos2, end=pos2))) 
  sv_hits1 = GenomicAlignments::findOverlaps(sv1gr,genegr)
  sv_hits2 = GenomicAlignments::findOverlaps(sv2gr,genegr)
  svtab1 <- data2[sort(unique(c(queryHits(sv_hits1),queryHits(sv_hits2)))),]
  genegr <- with(data.frame(chr1,start1,stop1), GRanges(chr, IRanges(start=start1, end=stop1)))
  sv1gr = with(data2, GRanges(chrom1, IRanges(start=pos1, end=pos1))) 
  sv2gr = with(data2 , GRanges(chrom2, IRanges(start=pos2, end=pos2))) 
  sv_hits1 = GenomicAlignments::findOverlaps(sv1gr,genegr)
  sv_hits2 = GenomicAlignments::findOverlaps(sv2gr,genegr)
  svtab2 <- data2[sort(unique(c(queryHits(sv_hits1),queryHits(sv_hits2)))),]
  cnvdat = CNV7@data
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
Driver_obj("MYH9-EIF3D")
#> [[1]]
#>          sample chrom1     pos1 strand1 chrom2     pos2 strand2 variant_type
#> 1389 HCC_BREAST  chr22 36488264       +  chr22 36814466       -          DEL
#> 
#> [[2]]
#>    sample chrom    start      end probes   segmean            uid
#> 1: HCC_22 chr22 36473919 36523907     21  3.274240 cnv_ESYoIoCcjH
#> 2: HCC_22 chr22 36523907 36573895   3582 -5.908850 cnv_fRJiEarzuG
#> 3: HCC_22 chr22 36573895 36623884     41 -6.101530 cnv_KYGrOaHJwZ
#> 4: HCC_22 chr22 36623884 36673872      4 -9.908850 cnv_rlHBcUnJrY
#> 5: HCC_22 chr22 36673872 36723860     90 -8.323890 cnv_OMWugIPapf
#> 6: HCC_22 chr22 36723860 36773848  13753 -6.908850 cnv_ldnMrgBjcc
#> 7: HCC_22 chr22 36773848 36823836      6 -6.735850 cnv_FqBmZkxxKO
#> 8: HCC_22 chr22 36823836 36873824   1399 -0.468065 cnv_uqJNNZeToJ
#> 
#> [[3]]
#>          sample chrom1     pos1 strand1 chrom2     pos2 strand2 variant_type
#> 1389 HCC_BREAST  chr22 36488264       +  chr22 36814466       -          DEL
#> 
#> [[4]]
#>    sample chrom    start      end probes   segmean            uid
#> 1: HCC_22 chr22 36673872 36723860     90 -8.323890 cnv_OMWugIPapf
#> 2: HCC_22 chr22 36723860 36773848  13753 -6.908850 cnv_ldnMrgBjcc
#> 3: HCC_22 chr22 36773848 36823836      6 -6.735850 cnv_FqBmZkxxKO
#> 4: HCC_22 chr22 36823836 36873824   1399 -0.468065 cnv_uqJNNZeToJ

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
  cnvdat = CNV7@data
  genegr <- with(data.frame(chr,start,stop), GRanges(chr, IRanges(start=start, end=stop))) 
  sv1gr = with(data2, GRanges(chrom1, IRanges(start=pos1, end=pos1))) 
  sv2gr = with(data2 , GRanges(chrom2, IRanges(start=pos2, end=pos2))) 
  sv_hits1 = GenomicAlignments::findOverlaps(sv1gr,genegr)
  sv_hits2 = GenomicAlignments::findOverlaps(sv2gr,genegr)
  svtab1 <- data2[sort(unique(c(queryHits(sv_hits1),queryHits(sv_hits2)))),]
  genegr <- with(data.frame(chr1,start1,stop1), GRanges(chr, IRanges(start=start1, end=stop1)))
  sv1gr = with(data2, GRanges(chrom1, IRanges(start=pos1, end=pos1))) 
  sv2gr = with(data2 , GRanges(chrom2, IRanges(start=pos2, end=pos2))) 
  sv_hits1 = GenomicAlignments::findOverlaps(sv1gr,genegr)
  sv_hits2 = GenomicAlignments::findOverlaps(sv2gr,genegr)
  svtab2 <- data2[sort(unique(c(queryHits(sv_hits1),queryHits(sv_hits2)))),]
  cnvdat = CNV7@data
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
  z = nrow(CNV7@data)/nrow(data2)
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

sv.model.view <- function(svc, cnv, chr, start, stop, 
                          sampleids=NULL,
                          cnvlim=c(-2,2), 
                          addlegend='both',
                          cex.legend=1,
                          interval=NULL,
                          addtext=NULL,
                          cex.text=.8,
                          plot=TRUE,
                          summary=TRUE,
                          ...){
    

 stopifnot(!is.null(chr) && !is.null(start) && !is.null(stop))

    stopifnot(cnv@type == "cnv")
    cnvdat <- cnv@data
    
    stopifnot(svc@type == "svc")
    svcdat <- svc@data
    
    if(!is.null(sampleids)){
        missing.samples <- setdiff(sampleids,c(svcdat$sample,cnvdat$sample))
        if(length(missing.samples) == length(unique(sampleids))){
            stop("None of the samples provided were found in 'sv' and 'cnv' input data!")
        }else if(length(missing.samples) > 0){
            warning(paste("The following samples provided are not found in 'sv' and 'cnv' input data:", paste(missing.samples,collapse=" "),sep=" "))
        }
        svcdat<-svcdat[which(svcdat$sample %in% intersect(sampleids,svcdat$sample)),]
        cnvdat<-cnvdat[which(cnvdat$sample %in% intersect(sampleids,cnvdat$sample)),]
    }

    genegr <- with(data.frame(chr,start,stop), GRanges(chr, IRanges(start=start, end=stop))) 
    
    # Find samples with SV breaks within defined genomic region
    sv1gr = with(svcdat, GRanges(chrom1, IRanges(start=pos1, end=pos1))) 
    sv2gr = with(svcdat, GRanges(chrom2, IRanges(start=pos2, end=pos2))) 
    
    sv_hits1 = GenomicAlignments::findOverlaps(sv1gr,genegr)
    sv_hits2 = GenomicAlignments::findOverlaps(sv2gr,genegr)
    svtab <- svcdat[sort(unique(c(queryHits(sv_hits1),queryHits(sv_hits2)))),]
    svBreakSamples <- unique(svtab$sample)
    if(length(svBreakSamples) == 0) warning("Thre is no SV breakpoints in the defined genomic region")
        
    # obtain SVs for plotting with different colors for each svclass
    svcolormap = setNames(c("blue", "red", "orange", "black", "green","grey20"), 
                   c("DEL", "DUP", "INV", "TRA", "INS","BND"))
    svcolor <- svcolormap[svtab$svclass]
    svtab_plot <- data.table(svtab,svcolor)
    svtab_plot_seg <- svtab_plot[which(svtab_plot$svclass != "TRA")]
    svtab_plot_tra <- svtab_plot[which(svtab_plot$svclass == "TRA")]
    
    # Find samples with CNV segment breaks within defined genomic region
    seg1br  = with(cnvdat, GRanges(chrom, IRanges(start=start, end=start))) 
    seg2br  = with(cnvdat, GRanges(chrom, IRanges(start=end, end=end))) 
    seg_hits1 = GenomicAlignments::findOverlaps(seg1br,genegr)
    seg_hits2 = GenomicAlignments::findOverlaps(seg2br,genegr)
    segBreakSamples <- unique(cnvdat[sort(unique(c(queryHits(seg_hits1),queryHits(seg_hits2))))]$sample)
    if(length(segBreakSamples) == 0)    
    segbrk <- cnvdat[sort(unique(c(queryHits(seg_hits1),queryHits(seg_hits2))))]
    
    if(plot==TRUE){
        # Find overlap between all CNV segments and the defined genomic region for plotting

        seggr <- with(cnvdat, GRanges(chrom, IRanges(start=start, end=end))) 
        hits_seg = GenomicAlignments::findOverlaps(seggr,genegr)
        seg_plot <- cnvdat[queryHits(hits_seg)]
        segcolor <- map2color(seg_plot$segmean,
                  pal=colorRampPalette(c("lightblue","white","salmon"))(256),
                  limit=cnvlim)
        seg_plot <- data.table(seg_plot,segcolor)
    
        if(!is.null(sampleids)){
            sample_order <- 1:length(sampleids)
            names(sample_order) <- sampleids
        }else{
            sample_order <- 1:length(unique(c(svBreakSamples,segBreakSamples)))
            names(sample_order) <- unique(c(svBreakSamples,segBreakSamples))
        }
    
        if(!is.null(addlegend)){
            plot_ylim <- length(sample_order)*10/100+length(sample_order)
            legend_ypos <- plot_ylim - length(sample_order)*3/100 
            if(length(sample_order) < 10) plot_ylim <- length(sample_order) +1
        }else{
            plot_ylim <- length(sample_order)
        }
        
        plot(x=NULL,y=NULL,xlim=range(c(start,stop)),ylim=range(c(0,plot_ylim)),
             xaxt='n',yaxt='n',xlab='',ylab='',bty='n', ...)
    
        mtext(side=2,at=sample_order-0.5,text=names(sample_order),las=2,line = 0.5, ...)
    
        for(sid in names(sample_order)){
            ypos <- sample_order[sid]
            polygon(rbind(
                c(start-1e7,ypos+0.02),
                c(start-1e7,ypos-0.98),
                c(stop+1e7,ypos-0.98),
                c(stop+1e7,ypos+0.02)),
                col=rep(c("grey80","grey80"),length(sample_order))[ypos],border=NA)
        }
        
        for(sid in names(sample_order)){
            seg_sample_plot <- seg_plot[which(seg_plot$sample == sid),]
            ypos <- sample_order[sid]
            for(i in 1:nrow(seg_sample_plot)){
                polygon(rbind(
                    c(seg_sample_plot[i]$start,ypos),
                    c(seg_sample_plot[i]$start,ypos-1),
                    c(seg_sample_plot[i]$end,ypos-1),
                    c(seg_sample_plot[i]$end,ypos)
                ),col=seg_sample_plot[i]$segcolor,border=NA)
            }
        }
    
        for(sid in unique(svtab_plot_seg$sample)){
            svtab_plot_seg_i <- svtab_plot_seg[which(svtab_plot_seg$sample == sid)]
            ypos <- sample_order[sid]
            addrnorm <- rep(c(0,0.2,-0.2,0.1,-0.1,0.3,-0.3),nrow(svtab_plot_seg_i))
            for(i in 1:nrow(svtab_plot_seg_i)){
                polygon(rbind(
                    c(svtab_plot_seg_i[i]$pos1,ypos-0.4-addrnorm[i]),
                    c(svtab_plot_seg_i[i]$pos1,ypos-0.6-addrnorm[i]),
                    c(svtab_plot_seg_i[i]$pos2,ypos-0.6-addrnorm[i]),
                    c(svtab_plot_seg_i[i]$pos2,ypos-0.4-addrnorm[i])
                    ),col=NA,border=svtab_plot_seg_i[i]$svcolor)

                if(svtab_plot_seg_i[i]$svclass %in% addtext){
                    if(svtab_plot_seg_i[i]$pos1 < start){
                        text(start,ypos-0.5-addrnorm[i],
                             paste("<--",svtab_plot_seg_i[i]$pos1,sep=""),
                             pos=4,offset=0,cex=cex.text)
                    }
                    if(svtab_plot_seg_i[i]$pos2 > stop){
                        text(stop,ypos-0.5-addrnorm[i],
                             paste(svtab_plot_seg_i[i]$pos2,"->",sep=""),
                             pos=2,offset=0,cex=cex.text)
                    }
                }
            }
        }
    
        for(sid in unique(svtab_plot_tra$sample)){
            svtab_plot_tra_i <- svtab_plot_tra[which(svtab_plot_tra$sample == sid),]
            ypos <- sample_order[sid]
            addrnorm <- rep(c(0,0.3,-0.3,0.1,-0.1,0.2,-0.2),nrow(svtab_plot_seg_i))
            for(i in 1:nrow(svtab_plot_tra_i)){
                if(svtab_plot_tra_i[i]$chrom2 == chr){ 
                    points(svtab_plot_tra_i[i]$pos2,ypos-0.5+addrnorm[i],pch=10)
                    lines(c(svtab_plot_tra_i[i]$pos2,svtab_plot_tra_i[i]$pos2),c(ypos,ypos-1),lwd=1,lty=3)
                    if("TRA" %in% addtext){
                        text(svtab_plot_tra_i[i]$pos2,ypos-0.5+addrnorm[i],
                             paste("  ",svtab_plot_tra_i[i]$chrom1,":",svtab_plot_tra_i[i]$pos1,sep=""),
                             pos=4,offset=0,cex=cex.text)
                    }
                }            
                if(svtab_plot_tra_i[i,"chrom1"] == chr){
                    points(svtab_plot_tra_i[i]$pos1,ypos-0.5+addrnorm[i],pch=10)
                    lines(c(svtab_plot_tra_i[i]$pos1,svtab_plot_tra_i[i]$pos1),c(ypos,ypos-1),lwd=1,lty=3)
                    if("TRA" %in% addtext) {
                        text(svtab_plot_tra_i[i]$pos1,ypos-0.5+addrnorm[i],
                             paste("  ",svtab_plot_tra_i[i]$chrom2,":",svtab_plot_tra_i[i]$pos2,sep=""),
                             pos=4,offset=0,cex=cex.text)
                    }
                }
            }
        }
    
        if(is.null(interval)) interval <- round((stop - start)/5000) * 1000
        xlabs <- seq(floor(start/10000)*10000, ceiling(stop/10000)*10000,interval)
        axis(1, at = xlabs,labels=TRUE, lwd.ticks=1.5, pos=0,...)

        if(is.null(cex.legend)) cex.legend <- 1
        
        if(addlegend %in% c("sv","both")) {
            fillx <- c("white", "white", "white", "white",NA)
            borderx <- c("blue", "red","orange","grey20",NA)
            pchx <- c(NA, NA,NA, NA,10)
            names(fillx) <-names(borderx) <-names(pchx) <- c("DEL", "DUP", "INV","BND", "TRA")
            svclassin <- unique(svcdat$svclass)
            
            legend(x= start, y =legend_ypos, legend = svclassin, 
                   bty = "n", fill = fillx[svclassin], border=borderx[svclassin], 
                   pch = pchx[svclassin], horiz = TRUE, x.intersp=0.2, cex=cex.legend)
        }
        if(addlegend %in% c("cnv","both")) {
            legend(x=start + (stop-start)/2, y = legend_ypos,legend = c(paste("CNV= ",cnvlim[1],sep=""), "CNV= 0", paste("CNV= ",cnvlim[2],sep=""), "no-data"),
                   bty = "n",fill=c("lightblue","white","salmon","grey80"), border=c(NA,"black",NA,NA), 
                   horiz = TRUE, x.intersp=0.2, cex=cex.legend)
        }
    }
    if(summary){
        return(list(svbrk=svtab,segbrk=segbrk))
    }
}
Driver_plot = function(x){
  y = strsplit(x, "-")
  m = as.character(unlist(y))
  library(org.Hs.eg.db)
  library(drawProteins)
  library(svpluscnv)
  library(data.table)
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
  sv.model.view(svc2, CNV7, chr, start, stop, sampleids=sampleids, 
              addlegend = 'both', addtext=c("TRA"), cnvlim = c(-2,2), 
              cex=.7,cex.text =.8, summary = FALSE)
  sv.model.view(svc2, CNV7, chr1, start1, stop1, sampleids=sampleids1, 
              addlegend = 'both', addtext=c("TRA"), cnvlim = c(-2,2), 
              cex=.7,cex.text =.8, summary = FALSE)

}





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
  p <- draw_domains(p, rel_data, label_domains = F)
  plot(p)
  z1 <- mapIds(org.Hs.eg.db, keys = m[2], keytype = "SYMBOL", column="UNIPROT")
  rel_json1 <- drawProteins::get_features(z1)
  rel_data1 <- drawProteins::feature_to_dataframe(rel_json1)
  p1 <- draw_canvas(rel_data1)
  p1 <- draw_chains(p1, rel_data1, label_size = 2.5)
  p1 <- draw_domains(p1, rel_data1,label_domains = F)
  plot(p1)
}
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date.

You can also embed plots, for example:

    #> `geom_smooth()` using formula 'y ~ x'

<img src="man/figures/README-pressure-1.png" width="100%" />

    #> 
    #> Attaching package: 'data.table'
    #> The following object is masked from 'package:GenomicRanges':
    #> 
    #>     shift
    #> The following object is masked from 'package:IRanges':
    #> 
    #>     shift
    #> The following objects are masked from 'package:S4Vectors':
    #> 
    #>     first, second

<img src="man/figures/README-pressure-2.png" width="100%" /><img src="man/figures/README-pressure-3.png" width="100%" />

    #> 'select()' returned 1:many mapping between keys and columns
    #> [1] "Download has worked"
    #> 'select()' returned 1:1 mapping between keys and columns

<img src="man/figures/README-pressure-4.png" width="100%" />

    #> [1] "Download has worked"

<img src="man/figures/README-pressure-5.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub\!
