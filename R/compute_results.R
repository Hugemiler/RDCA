compute_Results <- function(Pij,Pi,Pij_true,Pi_true,invC,n,q){
  # computes and prints the mutual and direct informations

  result_df <- data.frame()

  for (i in (1:(n-1))){
    for (j in ((i+1):n)){
      # mutual information
      calculated_mi <- calculate_mi(i,j,Pij_true,Pi_true,q)

      # direct information from mean-field
      W_mf <- returnW(invC,i,j,q);
      DI_mf_pc <- bp_link(i,j,W_mf,Pi,q);

      currentLine <- c(i, j, calculated_mi$mi, DI_mf_pc)
      result_df <- rbind(result_df, currentLine)
    }
  }
  return(result_df)

}
