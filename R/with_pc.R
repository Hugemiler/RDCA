with_pc <- function(Pij_true, Pi_true, pseudocount_weight,n,q){

  # adds pseudocount

  Pij = (1.-pseudocount_weight)*Pij_true + pseudocount_weight/q/q*array(1, c(n,n,q,q));
  Pi = (1.-pseudocount_weight)*Pi_true + pseudocount_weight/q*array(1, c(n,q));

  scra = diag(q);

  for (i in 1:n){
    for (alpha in 1:q){
      for (beta in 1:q){
        Pij[i,i,alpha,beta] =  (1.-pseudocount_weight)*Pij_true[i,i,alpha,beta] + pseudocount_weight/q*scra[alpha,beta];
      }
    }
  }

  return(list("Pij" = Pij,
              "Pi" = Pi))
}
