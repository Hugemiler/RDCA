compute_C <- function(Pij,Pi,n,q){
  # computes correlation matrix

  cMatrix <- matrix(0, n*(q-1), n*(q-1));

  for (i in 1:n){
    for (j in 1:n){
      for (alpha in (1:(q-1))){
        for (beta in (1:(q-1))){
          cMatrix[mapkey(i,alpha,q),mapkey(j,beta,q)] <- Pij[i,j,alpha,beta] - Pi[i,alpha]*Pi[j,beta];
        }
      }
    }
  }
  return(cMatrix)
}
