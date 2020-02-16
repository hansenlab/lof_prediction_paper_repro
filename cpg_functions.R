###important functions for promoter CpG density related computations
#
getOECpGRatio <- function(genomic_ranges, genome){
  oe_cpg <- vector()
  seqs <- getSeq(genome, seqnames(genomic_ranges), 
                 start(genomic_ranges), end(genomic_ranges))
  for (i in 1:length(seqs)){
    oe_cpg[i] <- length(seqs[[i]])*dinucleotideFrequency(seqs[[i]], 
                                                         as.prob = FALSE)["CG"]/(letterFrequency(seqs[[i]], "C")*letterFrequency(seqs[[i]], "G"))
  }
  oe_cpg
}


#
getDinucleotideContent <- function(dinucleotide, genomic_ranges, genome){
  dinucleotide_number <- vector()
  seqs <- getSeq(genome, seqnames(genomic_ranges), 
                 start(genomic_ranges), end(genomic_ranges))
  for (i in 1:length(seqs)){
    dinucleotide_number[i] <- dinucleotideFrequency(seqs[[i]], as.prob = FALSE)[dinucleotide]
  }
  dinucleotide_number
}


#
####the following two functions that are called to deal with overlapping promoters
keepNonOverlappingRanges <- function(genomic_ranges){ #this function selects the promoter whose gene has the lowest LOEUF
  self_overlaps <- findOverlaps(genomic_ranges, genomic_ranges, ignore.strand = TRUE)
  gr <- genomic_ranges[subjectHits(self_overlaps)]
  gr_list <- split(gr, queryHits(self_overlaps))
  new_gr <- endoapply(gr_list, function(xx) {
    xx[which.min(xx$oe_lof_upper)]
  })
  unique(unlist(new_gr))
}

keepNonOverlappingRanges2 <- function(genomic_ranges){ #this function selects the promoter that has the maximum oe cpg ratio
  self_overlaps <- findOverlaps(genomic_ranges, genomic_ranges, ignore.strand = TRUE)
  gr <- genomic_ranges[subjectHits(self_overlaps)]
  gr_list <- split(gr, queryHits(self_overlaps))
  new_gr <- endoapply(gr_list, function(xx) {
    xx[which.max(xx$oe_cpg_ratio)]
  })
  unique(unlist(new_gr))
}



###the following three functions are used to evaluate promoter overlap with structural variants
getDistributionOfOverlappingSVSizes <- function(genomic_ranges){
  width(pintersect(findOverlapPairs(genomic_ranges, structural_variants_granges)))
}


getProportionOfRangesOverlappingSVs <- function(genomic_ranges, bootstrap = c(TRUE, FALSE)){
  all_granges <- length(genomic_ranges)
  granges_with_SVs <- length(unique(queryHits(findOverlaps(genomic_ranges, structural_variants_granges))))
  granges_with_SVs / all_granges
}

getBootstrapProportionOfRangesOverlappingSVs <- function(genomic_ranges){
  all_granges <- length(genomic_ranges)
  genomic_ranges <- genomic_ranges[sample(length(genomic_ranges), replace = TRUE)]
  granges_with_SVs <- length(unique(queryHits(findOverlaps(genomic_ranges, structural_variants_granges))))
  granges_with_SVs / all_granges
}