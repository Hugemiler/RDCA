bp_link <- function(i,j,W,P1,q){
  # computes direct information

  mu <- compute_mu(i,j,W,P1,q);
  DI <- compute_di(i,j,W, mu$mu1,mu$mu2,P1);

  return(DI)
}
