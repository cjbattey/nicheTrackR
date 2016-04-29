ecospat.niche.similarity.test.noplot<- function(z1, z2, rep, one.sided = F) {
  R <- length(z1$x)
  #dev.new(2, 2, pointsize = 12)
  #par(mar = c(0, 0, 0, 0))
  l <- list()
  obs.o <- ecospat.niche.overlap(z1, z2, cor = T)
  sim.o <- data.frame(matrix(nrow = rep, ncol = 2))
  names(sim.o) <- c("D", "I")
  for (k in 1:rep) {
    #plot.new()
    #text(0.5, 0.5, paste("similarity tests:", "\n", "runs to go:", 
    #                     rep - k + 1))
    if (is.null(z2$y)) {
      center <- which(z2$z.cor == 1, arr.ind = T)
      Z <- z2$Z/max(z2$Z)
      rand.center <- sample(1:R, size = 1, replace = F, 
                            prob = Z)
      xshift <- rand.center - center
      z2.sim <- z2
      z2.sim$z.cor <- rep(0, R)
      for (i in 1:R) {
        i.trans <- i + xshift
        if (i.trans > R | i.trans < 0) 
          (next)()
        z2.sim$z.cor[i.trans] <- z2$z.cor[i]
      }
      z2.sim$z.cor <- (z2$Z != 0) * 1 * z2.sim$z.cor
    }
    if (!is.null(z2$y)) {
      centroid <- which(z2$z.cor == 1, arr.ind = T)[1, ]
      Z <- z2$Z/max(z2$Z)
      rand.centroids <- which(Z > 0, arr.ind = T)
      weight <- Z[Z > 0]
      rand.centroid <- rand.centroids[sample(1:nrow(rand.centroids), 
                                             size = 1, replace = F, prob = weight), ]
      xshift <- rand.centroid[1] - centroid[1]
      yshift <- rand.centroid[2] - centroid[2]
      z2.sim <- z2
      z2.sim$z.cor <- matrix(rep(0, R * R), ncol = R, nrow = R)
      for (i in 1:R) {
        for (j in 1:R) {
          i.trans <- i + xshift
          j.trans <- j + yshift
          if (i.trans > R | i.trans < 0) 
            (next)()
          if (j.trans > R | j.trans < 0) 
            (next)()
          z2.sim$z.cor[i.trans, j.trans] <- z2$z.cor[i, j]
        }
      }
      z2.sim$z.cor <- (z2$Z != 0) * 1 * z2.sim$z.cor
    }
    o.i <- ecospat.niche.overlap(z1, z2.sim, cor = T)
    sim.o$D[k] <- o.i$D
    sim.o$I[k] <- o.i$I
  }
  #dev.off()
  l$sim <- sim.o
  l$obs <- obs.o
  if (one.sided) {
    l$p.D <- (sum(sim.o$D >= obs.o$D) + 1)/(length(sim.o$D) + 
                                              1)
    l$p.I <- (sum(sim.o$I >= obs.o$D) + 1)/(length(sim.o$I) + 
                                              1)
  }
  else {
    #trying with pnorm...
    #l$p.D <- pnorm(obs.o$I,mean(sim.o$I),sd(sim.o$I),lower.tail=F)
    #l$p.I <- pnorm(obs.o$I,mean(sim.o$I),sd(sim.o$I),lower.tail=F)
    l$p.D <- min((sum(sim.o$D <= obs.o$D) + 1), (sum(sim.o$D >= obs.o$D) + 1)) * 2/(length(sim.o$D) + 1)
    l$p.I <- min((sum(sim.o$I <= obs.o$I) + 1), (sum(sim.o$I >= obs.o$I) + 1)) * 2/(length(sim.o$I) + 1)
  }
  return(l)
}