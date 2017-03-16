##############################################
#
# General description:
#
#   The following function reads in bed files of ATAC and RNA-seq data and
#   compares coverage.

# Input:
#
#     ATAC_bedfile: path to bed file with ATAC-seq data
#     RNAseq_bedfile: path to bed file with RNA-seq data

# Output:
#   
#    : ...

# Comments:
#   
#    Requires packages csaw

##############################################

#######################################
#                                     #
#    Read in data and                 #
#   determine coverage                #
#                                     #
#######################################

# Load packages
library(GenomicRanges)
library(ShortRead)
library(csaw)
library(Rsamtools)
library(rtracklayer)
library(data.table)

# Determine paths
ATAC_bedfile   <- '/srv/gsfs0/projects/levinson/hzudohna/ATACseq/iN_deletion_3_Aligned.out.sorted.bed'
RNAseq_bedfile <- '/srv/gsfs0/projects/levinson/hzudohna/ATACseq/no_dups_iN-1133-ATAC_S7_merged.sort.bed'
RNAseq_bamfile <- '/srv/gsfs0/projects/levinson/CP/data_for_Heinrich/iN_deletion_3_Aligned.out.sorted.bam'
ATACpeaks <- import.bed('/srv/gsfs0/projects/levinson/hzudohna/ATACseq/NA_summits.bed')
ATACpeaks <- ATACpeaks[seqnames(ATACpeaks) != "chrM"]
RNAseqReadList <-lapply(1:length(ATACpeaks), function(i){
  extractReads(RNAseq_bamfile, ATACpeaks[i])
})
RNAseqReadList <-lapply(1:10, function(i){
  extractReads(RNAseq_bamfile, ATACpeaks[i])
})

unique(seqnames(ATACpeaks))

# Loop through lines and save reads to fasta files for read 1 and 2
LinesPerScan  <- 5*10^5

# Create vector of chromosome names
ChromNames <- paste("chr", c(1:22, "X", "Y"), sep = "")

#  Scan reads from ATAC-seq 
TotalReads <- 0
i <- 1
GR_ATAC <- c()
ScannedLines <- matrix(nrow = LinesPerScan)
CoverList_ATAC <- lapply(ChromNames, function(x) 0)
names(CoverList_ATAC) <- ChromNames
cat("Reading ATAC-seq data piecewise ... \n")
while (nrow(ScannedLines)  == LinesPerScan){
  ScannedLines <- scan(ATAC_bedfile, skip = (i - 1) * LinesPerScan , 
                         nlines = LinesPerScan, 
                         what = 'character', sep = '\n')
    i <- i + 1
    ScannedLines <- fread(input = paste(ScannedLines, collapse = "\n"), 
                          header = F,
      col.names = c("Chrom", "Start", "End", "Name", "Nr", "Strand"))
    GR <- GRanges(seqnames = ScannedLines$Chrom, 
            ranges = IRanges(start = ScannedLines$Start, end = ScannedLines$End),
            strand = ScannedLines$Strand)
    for (Chr in unique(seqnames(GR))){
      GRsubset <- GR[seqnames(GR) == Chr]
      CoverList_ATAC[[Chr]] <- CoverList_ATAC[[Chr]] + coverage(GR)
    }
    
    TotalReads <- TotalReads + nrow(ScannedLines)
    cat("ATAC-seq: total of", TotalReads, "reads processed\n")
  }

#  Scan reads from RNA-seq 
TotalReads <- 0
i <- 1
GR_RNAseq <- c()
ScannedLines <- matrix(nrow = LinesPerScan)
CoverList_RNAseq <- lapply(ChromNames, function(x) 0)
names(CoverList_RNAseq) <- ChromNames
cat("Reading RNA-seq data piecewise ... \n")
while (nrow(ScannedLines)  == LinesPerScan){
  ScannedLines <- scan(RNAseq_bedfile, skip = (i - 1) * LinesPerScan , 
                       nlines = LinesPerScan, 
                       what = 'character', sep = '\n')
  i <- i + 1
  ScannedLines <- fread(input = paste(ScannedLines, collapse = "\n"), header = F,
                        col.names = c("Chrom", "Start", "End", "Name", "Nr", "Strand"))
  GR <- GRanges(seqnames = ScannedLines$Chrom, 
                ranges = IRanges(start = ScannedLines$Start, end = ScannedLines$End),
                strand = ScannedLines$Strand)
  for (Chr in unique(seqnames(GR))){
    GRsubset <- GR[seqnames(GR) == Chr]
    CoverList_RNAseq[[Chr]] <- CoverList_RNAseq[[Chr]] + coverage(GR)
  }
  TotalReads <- TotalReads + nrow(ScannedLines)
  cat("RNA-seq: total of", TotalReads, "reads processed\n")
}

#######                       
# Calculate coverage                    
#######                                     

# Determine coverage
CovATAC   <- coverage(GR_ATAC)
CovRNAseq <- coverage(GR_RNAseq)

# Determine separate islands with continuous read coverage and turn islands 
# into genomic ranges
IslATAC   <- slice(CovATAC, lower = 1)
IslRNAseq <- slice(CovRNAseq, lower = 1)

#######################################
#                                     #
#    Determine 'islands' with         #
#        nonzero coverage             #
#                                     #
#######################################

   # Determine separate islands with continuous read coverage and turn islands 
   # into genomic ranges
   IslandList <- lapply(CoverList, function(x){
     Islands <- slice(x, lower = 1)
   })
   Chroms <- names(ChromLengths)
   IslandGRanges <- lapply(1:length(IslandList), function(i){
   GRanges(seqnames = Chroms[i], 
          ranges = IslandList[[i]]@listData[[1]]@ranges,
          coverTotal = viewSums(IslandList[[i]])[[1]],
          coverMax   = viewMaxs(IslandList[[i]])[[1]],
          coverMaxPos   = viewWhichMaxs(IslandList[[i]])[[1]])
   })
   IslandGRanges <- GRangesList(IslandGRanges)
   IslandGRanges <- unlist(IslandGRanges)
   cat(length(IslandGRanges), "distinct peaks\n\n")

   # Merge ranges that are less than MinGap bp apart
   IslGRanges_reduced <- reduce(IslandGRanges, min.gapwidth = MinGap,
                             with.revmap = T)
   cat(length(IslGRanges_reduced), "distinct peaks after merging peaks that")
   cat("are less than", MinGap, "apart\n\n")

   #######################################
   #                                     #
   #    Find 'islands' not overlapping   #
   #        with reference L1            #
   #                                     #
   #######################################
   
   cat("***** Finding 'islands' not overlapping with reference L1 ... *****\n")

   # Find overlaps between islands and L1HS ranges
   blnOverlapIslands_All <- overlapsAny(IslGRanges_reduced, L1GRanges)

   # Get maximum cover and position of maximum cover in reduced ranges
   maxCoverOriginal    <- IslandGRanges@elementMetadata@listData$coverMax
   maxCoverPosOriginal <- IslandGRanges@elementMetadata@listData$coverMaxPos
   maxCover <- sapply(IslGRanges_reduced@elementMetadata@listData$revmap, 
                   function(x) max(maxCoverOriginal[x]))
   maxCoverPos <- sapply(IslGRanges_reduced@elementMetadata@listData$revmap, 
                      function(x) maxCoverPosOriginal[x[which.max(maxCoverOriginal[x])]])
   
   # Get all ranges that make the maximum cover cut-off and don't overlap with
   # reference L1
   idxSuspectL1Ranges <- which(maxCover > MinMaxCover & (!blnOverlapIslands_All))
   SuspectL1Ranges    <- IslGRanges_reduced[idxSuspectL1Ranges]
   cat(length(idxSuspectL1Ranges), "peaks have maximum coverage of at least",
      MinMaxCover, "and do not overlap with reference L1\n")

   # Remove ranges of suspected L1s that are too close to reference L1
   DistToNearestL1    <- nearest(SuspectL1Ranges, L1GRanges)
   idxSuspectL1Ranges <- idxSuspectL1Ranges[DistToNearestL1 >= MinDist2L1]
   SuspectL1Ranges    <- IslGRanges_reduced[idxSuspectL1Ranges]
   maxCoverPos_SuspL1Ranges <- maxCoverPos[idxSuspectL1Ranges]
   cat(length(idxSuspectL1Ranges), "of the above peaks are at least", 
       MinDist2L1, "bp from nearest reference L1\n\n")
   
   #######################################
   #                                     #
   #    Find 'islands' overlapping       #
   #    with full-length reference L1    #
   #                                     #
   #######################################
   
   cat("***** Filtering reads overlapping with full-length reference L1 ... *****\n")

   # Find overlaps between islands and full-length L1HS ranges
   blnOverlapIslands_L1HS <- overlapsAny(IslGRanges_reduced, 
                                         L1HSFullLength_GRanges)
   
   # Get genomoc ranges of islands overlapping with full-length L1HS
   idxFullRefL1Ranges <- which(blnOverlapIslands_L1HS)
   FullRefL1Ranges    <- IslGRanges_reduced[blnOverlapIslands_L1HS]
   
   # Filter bam file to get reads in islands overlapping with full-length L1
   param <- ScanBamParam(which = FullRefL1Ranges, what = scanBamWhat())
   filterBam(file = BamFile, destination = OutBamFileFullLengthL1, param = param)
   
   #######################################################
   #                                                     #
   #    Save results                                     #
   #                                                     #
   #######################################################

   # Save results
   cat("*******  Saving results for ComparePeaksWithRefL1...   *******\n")
   save(list = c("IslGRanges_reduced", "maxCover", "maxCoverPos", 
              "idxSuspectL1Ranges", "SuspectL1Ranges", 
              "FullRefL1Ranges", "idxFullRefL1Ranges"), 
     file = OutFile)



