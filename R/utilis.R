#' Detection of transposon insertion sites
#' @description  Transposon insertion site search and TSD sequenc assembly
#' @param file InsertLocation function output
#' @param max.gapwidth The maximum distance between two site combined
#' @importFrom S4Vectors subjectHits
#' @importFrom S4Vectors queryHits
#' @importFrom GenomicRanges reduce
#' @importFrom GenomicRanges findOverlaps
#' @importFrom Biostrings DNAString
#' @importFrom Biostrings pairwiseAlignment
#' @importFrom Biostrings alignedSubject
#' @importFrom Biostrings subseq
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' 
#' @export

processInsertion <- function(file, max.gapwidth = 10) {
  bed <- tryCatch(
            expr = {
              file <- read.table(file = file, sep = "\t")
              bed <- GenomicRanges::GRanges(
                        seqnames = file$V1, 
                        IRanges::IRanges(start = file$V2, end = file$V3),
                        name = file$V4, 
                        sequence = file$V5, 
                        type = file$V6,
                        count = file$V7, 
                        ratio = file$V8)
             },
             error = function(...) {
                return(GenomicRanges::GRanges())
             }
          )
  if (length(bed) == 0) {
    return(bed)
  }
  insert <- GenomicRanges::reduce(bed, min.gapwidth = max.gapwidth)
  overlap <- GenomicRanges::findOverlaps(query = insert, subject = bed)
  index <- split(x = S4Vectors::subjectHits(overlap),
                 f = S4Vectors::queryHits(overlap))
  predict <- lapply(index, FUN = function(x, min.gapwidth) {
    location <- bed[x]
    insertion <- GenomicRanges::reduce(location, min.gapwidth = max.gapwidth)
    
    if (min(location$type) == 2) {
      location <- location[which.max(location$count)]
      insertion$name <- location$name
      insertion$tsd <- "unknown"
      insertion$left <- sum(location$count)
      insertion$right <- 0
      
    } else if(max(location$type) == 1) {
      location <- location[which.max(location$count)]
      insertion$name <- location$name
      insertion$tsd <- "unknown"
      insertion$left <- 0
      insertion$right <- sum(location$count)
    } else {
      right <- location[location$type == 1]
      left <- location[location$type == 2]
      
      left <- left[which.max(left$count)]
      right <- right[which.max(right$count)]
      
      dna2 <- Biostrings::DNAString(left$sequence)
      dna1 <- Biostrings::DNAString(right$sequence)
      
      pair <- Biostrings::pairwiseAlignment(pattern = dna1, subject = dna2, type = "overlap")
      p1 <- Biostrings::alignedSubject(pair)
      #insertion$name <- unique(location$name)
      insertion$name <- left$name
      
      if (Biostrings::subseq(dna1, 1, Biostrings::width(p1)) == p1 & 
          Biostrings::subseq(dna2, length(dna2)-Biostrings::width(p1)+1, length(dna2)) == p1) {
        tsd = as.character(p1)
        insertion$tsd <- tsd
      } else {
        insertion$tsd <- "unknown"
      }
      insertion$left <- left$count
      insertion$right <- right$count
      
    }
    return(insertion)
    
  }, min.gapwidth = max.gapwidth)
  
  insertion <- Reduce("c", predict)
  return(insertion)
}
