dca <- function(inputfile, outputfile){

  ##########
  # Direct Coupling Analysis
  ##########
  # INPUTS:
  #
  #   inputfile   - file containing the FASTA alignment
  #   outputfile  - file for dca results. The file is composed by N(N-1)/2
  #                 (N = length of the sequences) rows and 4 columns:
  #                 residue i (column 1), residue j (column 2),
  #                 MI(i,j) (Mutual Information between i and j), and
  #                 DI(i,j) (Direct Information between i and j).
  #                 Note: all insert columns are removed from the alignment.
  # SOME RELEVANT VARIABLES:
  # N        number of residues in each sequence (no insert)
  # M        number of sequences in the alignment
  # Meff     effective number of sequences after reweighting
  # q        equal to 21 (20 aminoacids + 1 gap)
  # align    M x N matrix containing the alignmnent
  # Pij_true N x N x q x q matrix containing the reweigthed frequency counts.
  # Pij      N x N x q x q matrix containing the reweighted frequency counts with pseudo counts.
  # C        N(q-1) x N(q-1) matrix containing the covariance matrix.
  #
  # Copyright for this implementation:
  #               2011/12 - Andrea Pagnani and Martin Weigt
  #                         andrea.pagnani@gmail.com
  #                         martin.weigt@upmc.fr
  #
  #   Permission is granted for anyone to copy, use, or modify this
  #   software and accompanying documents for any uncommercial
  #   purposes, provided this copyright notice is retained, and note is
  #   made of any changes that have been made. This software and
  #   documents are distributed without any warranty, express or
  #   implied. All use is entirely at the user's own risk.
  #
  #   Any publication resulting from applications of DCA should cite:
  #
  #      F Morcos, A Pagnani, B Lunt, A Bertolino, DS Marks, C Sander,
  #      R Zecchina, JN Onuchic, T Hwa, M Weigt (2011), Direct-coupling
  #      analysis of residue co-evolution captures native contacts across
  #      many protein families, Proc. Natl. Acad. Sci. 108:E1293-1301.

  ## Library Load
  require(seqinr)
  require(e1071)

  ## Read input and output from command line
  args <- commandArgs(trailingOnly = TRUE)
  inputfile <- args[1]
  outputfile <- args[2]

  ## set pseudocount-weight and theta
  pseudocount_weight <- 0.5; # relative weight of pseudo count
  theta <- 0.2;              # threshold for sequence id in reweighting

  ## Read FASTA input file and compute ALIGN matrix
  initial_alignment <- return_alignment(inputfile)

  ## Calculate frequencies from HUMMING DISTANCE
  frequencies <- compute_true_frequencies(initial_alignment$align,
                                          initial_alignment$m,
                                          initial_alignment$n,
                                          initial_alignment$q,
                                          theta)

  ## Add pseudocount to values
  with_pseudocount <- with_pc(frequencies$Pij_true,
                              frequencies$Pi_true,
                              pseudocount_weight,
                              initial_alignment$n,
                              initial_alignment$q);

  ## Compute covariation matrix
  covarianceMatrix <- compute_C(with_pseudocount$Pij,
                                with_pseudocount$Pi,
                                initial_alignment$n,
                                initial_alignment$q)

  ## Inversion of Covariance Matrix
  invC <- solve(covarianceMatrix)

  ## Build the output table
  result_df <- compute_Results(with_pseudocount$Pij,
                               with_pseudocount$Pi,
                               frequencies$Pij_true,
                               frequencies$Pi_true,
                               invC,
                               initial_alignment$n,
                               initial_alignment$q)

  ## Output result
  write.table(result_df,
              file = as.character(outputfile))

}
