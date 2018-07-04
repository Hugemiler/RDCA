compute_true_frequencies <- function(align,m,n,q,theta) {

  # computes reweighted frequency counts

  w = rep(1,m)

  if( theta > 0.0 ) {
    w = (1./(1+colSums(hamming.distance(align)/m < theta)))
  }

  Meff <- sum(w);

  Pij_true = array(data = 0, dim = c(n, n, q, q));
  Pi_true = matrix(0, n, q);

  for (j in 1:m){
    for (i in 1:n){
      Pi_true[i, align[j,i]] = Pi_true[i, align[j,i]] + w[j];
    }
  }

  Pi_true <- Pi_true/Meff;

  for (l in 1:m){
    for (i in (1:(n-1))){
      for (j in ((i+1):n)){
        Pij_true[i,j,align[l,i],align[l,j]] = Pij_true[i,j,align[l,i],align[l,j]] + w[l]
        Pij_true[j,i,align[l,j],align[l,i]] = Pij_true[i,j,align[l,i],align[l,j]]
      }
    }
  }

  Pij_true = Pij_true/Meff;

  scra = diag(q)

  for (i in 1:n){
    for (alpha in 1:q){
      for (beta in 1:q){
        Pij_true[i,i,alpha,beta] = Pi_true[i,alpha] * scra[alpha,beta];
      }
    }
  }

  return(list("Pij_true" = Pij_true,
              "Pi_true" = Pi_true,
              "Meff" = Meff))

}
