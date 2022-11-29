#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Process sequence
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Extract soft clip reads
#' 
#' Extract soft clip reads from alignment.
#' 
#' @param file Character, the alignment file. Supported formats include BAM and SAM.
#' @param outfq Character, the output FASTQ file.
#' @param mapq The MAPping quality, 10 by default.
#' @param length The minimum soft clip reads length, 30 by default.
#' @param tsd The max length of tsd.
#' @export
#'  

extractSoftClip <- function(file, outfq, mapq = 20, length = 30, tsd = 10) {
  if (!file.exists(file)) {
    stop("SAM/BAM file not exits.")
  }
  if (file.exists(outfq)) {
    stop("Output fastq file exists.")
  }
  processbam(bamfile = file, outfq = outfq, quantile = mapq, length = length, tsd = tsd)
}


#' Align soft clip reads
#' 
#' Align soft clip reads to the transposable elements reference.
#' 
#' @param reference The TE reference built with bowtie2
#' @param fastq The fastq file of soft clip reads, out from \bold{extractSoftClip}
#' @param bamOutput Character, name of output bam file
#' @param threads The number of cores for parallel computing.
#' 
#' @export

alignment <- function(reference, fastq, bamOutput, threads = 1) {
  samprefix <- gsub(pattern = "(\\.sam)|(\\.bam)", replacement = "", x = bamOutput)
  commands <- paste0("--local ", "--threads ", threads)
  Rbowtie2::bowtie2(bt2Index = reference, seq1 = fastq, samOutput = paste0(samprefix, ".sam"), commands)
  Rsamtools::asBam(file = paste0(samprefix, ".sam"), destination = samprefix)
  Rsamtools::sortBam(file = paste0(samprefix, ".bam"), byQname = TRUE, destination = paste0(samprefix, ".sort"))
}

#' Identify the TE insertions
#' 
#' Identify the locations of TE insertions from the newly alignment result.
#' 
#' @param file Charcater, the alignment file in bam/sam format, out from \bold{alignment}
#' @param outBed Character, the bed file for insertions.
#' @param ratio Supports the softclip reads ratio
#' @export

insertLocation <- function(file, outBed, ratio = 0.1) {
  processSam2bed(alignmentfile = file, outbedfile = outBed, ratio = ratio)
}

#' Build a Bowtie index from transposon sequences 
#' 
#' @param reference Character, the path to the fasta reference containing the
#'  transposon sequences.
#' @param index Character, the prefix of index.
#' @param overwrite Boolean, whether to overwrite existing files. Default is true.
#' default is 'TRUE'
#' @param ... Other arguments
#' 
#' @export
buildIndex <- function(reference, index, overwrite = FALSE, ...) {
  Rbowtie2::bowtie2_build(references = reference,
                          bt2Index = index,
                          overwrite = overwrite, 
                          ...)
}
