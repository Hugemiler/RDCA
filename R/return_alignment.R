return_alignment <- function(inputfile) {

  # reads alignment from inputfile, removes inserts and converts into numbers
  align_full <- read.fasta(inputfile)
  m <- length(align_full)
  ind <- sapply(align_full, length)
  n = max(ind)
  z = matrix(0, m, n)

  sequenceChars <- c("-", "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
  sequenceNumbers <- 1:21

  align_matrix <- do.call(rbind, lapply(align_full, function(x) {toupper(as.character(x))}))

  for (i in 1:m) {
    for (j in 1:n){
      try <- sequenceNumbers[which(sequenceChars == align_matrix[i,j])]
      if (length(try) == 0) { z[i,j] <- 1}
      else z[i,j] <- try
    }
  }

  q = max(z);

  return(list("m" = m,
              "n" = n,
              "q" = q,
              "align" = z))
}
