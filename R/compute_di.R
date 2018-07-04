compute_di <- function(i,j,W,mu1,mu2,Pia){
  # computes direct information

  tiny <- 1.0e-100

  Pdir <- W * crossprod(mu1, mu2)
  Pdir <- Pdir / sum(Pdir)

  Pfac <- crossprod(Pia[i,],Pia[j,])

  tampo <- (Pdir + tiny)/as.numeric(Pfac + tiny)
  DI <- sum(diag(crossprod(Pdir, log(tampo))))

  return(DI)

}
