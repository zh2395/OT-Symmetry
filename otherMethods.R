# Below are some additional codes that implement other methods mentioned in Appendix D of the supplementary file.

library(DepthProc)
library(ddalpha)

# Method: Depth-based method (DLP)
# Example C1: normal distribution
loc = (0:4)/10
cen_nor = rep(NA,5)
for (i in 1:5) {
  get_p_val = function(id) {
    set.seed(id)
    n = 200
    x <- MASS::mvrnorm(n, c(0,0), matrix(c(2,1,1,3),2,2)) + loc[i]
    ords = order(DepthProc::depthTukey(x, rbind(x,-x), threads = 1, exact = TRUE), decreasing = TRUE)
    run = function(j) {
      return(ddalpha::depth.simplicial(c(0,0), x[c(ords[j],ords[j-1],ords[j-2]),], exact = TRUE))
    }
    z_stat = (4*(sum(unlist(parallel::mclapply(3:n, run, mc.cores = 1))) + 1)-n-2)/sqrt(11*n/3)
    return(pnorm(z_stat))
  }
  cen_nor[i] = mean(unlist(parallel::mclapply(1:10000, get_p_val, mc.cores = 7)) <= 0.05)
  print(cen_nor[i])
}
# 0.0481 0.0552 0.0760 0.1209 0.2013
# Example C2: t distribution
loc = (0:4)*0.2
cen_t = rep(NA,5)
for (i in 1:5) {
  get_p_val = function(id) {
    set.seed(id)
    n = 200
    x <- MASS::mvrnorm(n, c(0,0), matrix(c(2,1,1,3),2,2))/sqrt(rchisq(n, df=1)) + loc[i]
    ords = order(DepthProc::depthTukey(x, rbind(x,-x), threads = 1, exact = TRUE), decreasing = TRUE)
    run = function(j) {
      return(ddalpha::depth.simplicial(c(0,0), x[c(ords[j],ords[j-1],ords[j-2]),], exact = TRUE))
    }
    z_stat = (4*(sum(unlist(parallel::mclapply(3:n, run, mc.cores = 1))) + 1)-n-2)/sqrt(11*n/3)
    return(pnorm(z_stat))
  }
  cen_t[i] = mean(unlist(parallel::mclapply(1:10000, get_p_val, mc.cores = 7)) <= 0.05)
  print(cen_t[i])
}
# 0.0462 0.0661 0.1225 0.2633 0.4864
# Example C3: Uniform distribution
loc = (0:4)*0.03
cen_unif = rep(NA,5)
for (i in 1:5) {
  get_p_val = function(id) {
    set.seed(id)
    n = 200
    x <- matrix(runif(n*2),n,2)*sample(c(-1,1),n,replace = TRUE) + loc[i]
    ords = order(DepthProc::depthTukey(x, rbind(x,-x), threads = 1, exact = TRUE), decreasing = TRUE)
    run = function(j) {
      return(ddalpha::depth.simplicial(c(0,0), x[c(ords[j],ords[j-1],ords[j-2]),], exact = TRUE))
    }
    z_stat = (4*(sum(unlist(parallel::mclapply(3:n, run, mc.cores = 1))) + 1)-n-2)/sqrt(11*n/3)
    return(pnorm(z_stat))
  }
  cen_unif[i] = mean(unlist(parallel::mclapply(1:10000, get_p_val, mc.cores = 7)) <= 0.05)
  print(cen_unif[i])
}
# 0.0489 0.0837 0.1473 0.2396 0.3564

# Example C4: Exponential distribution 
get_p_val = function(id) {
  set.seed(id)
  n = 100
  x <- matrix(rexp(2 * n, rate = 1) - 1, n, 2)
  ords = order(DepthProc::depthTukey(x, rbind(x,-x), threads = 1, exact = TRUE), decreasing = TRUE)
  run = function(j) {
    return(ddalpha::depth.simplicial(c(0,0), x[c(ords[j],ords[j-1],ords[j-2]),], exact = TRUE))
  }
  z_stat = (4*(sum(unlist(parallel::mclapply(3:n, run, mc.cores = 1))) + 1)-n-2)/sqrt(11*n/3)
  return(pnorm(z_stat))
}
print(mean(unlist(parallel::mclapply(1:1000, get_p_val, mc.cores = 7)) <= 0.05))
# 0.843
# Example C5: Exponential distribution (larger sample)
get_p_val = function(id) {
  set.seed(id)
  n = 200
  x <- matrix(rexp(2 * n, rate = 1) - 1, n, 2)
  ords = order(DepthProc::depthTukey(x, rbind(x,-x), threads = 1, exact = TRUE), decreasing = TRUE)
  run = function(j) {
    return(ddalpha::depth.simplicial(c(0,0), x[c(ords[j],ords[j-1],ords[j-2]),], exact = TRUE))
  }
  z_stat = (4*(sum(unlist(parallel::mclapply(3:n, run, mc.cores = 1))) + 1)-n-2)/sqrt(11*n/3)
  return(pnorm(z_stat))
}
print(mean(unlist(parallel::mclapply(1:1000, get_p_val, mc.cores = 7)) <= 0.05))
# 0.997

# Example C6: Chisq distribution
get_p_val = function(id) {
  set.seed(id)
  n = 100
  x <- matrix(rchisq(2 * n, df = 1) - 1, n, 2)
  ords = order(DepthProc::depthTukey(x, rbind(x,-x), threads = 1, exact = TRUE), decreasing = TRUE)
  run = function(j) {
    return(ddalpha::depth.simplicial(c(0,0), x[c(ords[j],ords[j-1],ords[j-2]),], exact = TRUE))
  }
  z_stat = (4*(sum(unlist(parallel::mclapply(3:n, run, mc.cores = 1))) + 1)-n-2)/sqrt(11*n/3)
  return(pnorm(z_stat))
}
print(mean(unlist(parallel::mclapply(1:1000, get_p_val, mc.cores = 7)) <= 0.05))
# 0.987

# Example C7: Pareto distribution
get_p_val = function(id) {
  set.seed(id)
  n = 200
  x <- matrix(1 / runif(2 * n) - 2, n, 2)
  ords = order(DepthProc::depthTukey(x, rbind(x,-x), threads = 1, exact = TRUE), decreasing = TRUE)
  run = function(j) {
    return(ddalpha::depth.simplicial(c(0,0), x[c(ords[j],ords[j-1],ords[j-2]),], exact = TRUE))
  }
  z_stat = (4*(sum(unlist(parallel::mclapply(3:n, run, mc.cores = 1))) + 1)-n-2)/sqrt(11*n/3)
  return(pnorm(z_stat))
}
print(mean(unlist(parallel::mclapply(1:1000, get_p_val, mc.cores = 7)) <= 0.05))
# 1



# Methods: Spatial sign (SS) and signed rank (SSR)
loc = (0:4)*0.1
rep = 10000
cores = 7
# Example S1: normal distribution
# Method: spatial signed rank
for (i in 1:5) {
  Reject_symm = function(id) {
    set.seed(id)
    n = 200
    p = 2
    x <- MASS::mvrnorm(n, c(0,0), matrix(c(2,0,0,3),2,2)) + loc[i]
    if (MNM::mv.1sample.test(x, score = "rank", method = "signchange")$p.value <= 0.05) {
      return(1)
    }
    return(0)
  }
  print(mean(unlist(parallel::mclapply(1:rep, Reject_symm, mc.cores = cores))))
}
# 0.0492 0.1882 0.6141 0.9353 0.9969
# Example S1: normal distribution
# Method: spatial sign
for (i in 1:5) {
  Reject_symm = function(id) {
    set.seed(id)
    n = 200
    p = 2
    x <- MASS::mvrnorm(n, c(0,0), matrix(c(2,0,0,3),2,2)) + loc[i]
    if (MNM::mv.1sample.test(x, score = "sign", method = "signchange")$p.value <= 0.05) {
      return(1)
    }
    return(0)
  }
  print(mean(unlist(parallel::mclapply(1:rep, Reject_symm, mc.cores = cores))))
}
# 0.052 0.1587 0.5163 0.8731 0.9879
# Example S2: t distribution
# Method: spatial signed rank
for (i in 1:5) {
  Reject_symm = function(id) {
    set.seed(id)
    n = 200
    x <- MASS::mvrnorm(n, c(0,0), matrix(c(2,0,0,3),2,2))/sqrt(rchisq(n, df=1)) + loc[i]
    if (MNM::mv.1sample.test(x, score = "rank", method = "signchange")$p.value <= 0.05) {
      return(1)
    }
    return(0)
  }
  print(mean(unlist(parallel::mclapply(1:rep, Reject_symm, mc.cores = cores))))
}
# 0.0515 0.0943 0.2398 0.4774 0.7223
# Example S2: t distribution
# Method: spatial sign
for (i in 1:5) {
  Reject_symm = function(id) {
    set.seed(id)
    n = 200
    x <- MASS::mvrnorm(n, c(0,0), matrix(c(2,0,0,3),2,2))/sqrt(rchisq(n, df=1)) + loc[i]
    if (MNM::mv.1sample.test(x, score = "sign", method = "signchange")$p.value <= 0.05) {
      return(1)
    }
    return(0)
  }
  print(mean(unlist(parallel::mclapply(1:rep, Reject_symm, mc.cores = cores))))
}
# 0.0525 0.1182 0.3447 0.6648 0.8946
# Example S3: uniform distribution
# Method: spatial signed rank
loc = (0:4)*0.025
for (i in 1:5) {
  Reject_symm = function(id) {
    set.seed(id)
    n = 200
    x <- matrix(runif(2*n,-1,1),n,2) + loc[i]
    if (MNM::mv.1sample.test(x, score = "rank", method = "signchange")$p.value <= 0.05) {
      return(1)
    }
    return(0)
  }
  print(mean(unlist(parallel::mclapply(1:rep, Reject_symm, mc.cores = cores))))
}
# 0.0577 0.1036 0.2912 0.5803 0.8328
# Example S3: uniform distribution
# Method: spatial sign
for (i in 1:5) {
  Reject_symm = function(id) {
    set.seed(id)
    n = 200
    x <- matrix(runif(2*n,-1,1),n,2) + loc[i]
    if (MNM::mv.1sample.test(x, score = "sign", method = "signchange")$p.value <= 0.05) {
      return(1)
    }
    return(0)
  }
  print(mean(unlist(parallel::mclapply(1:rep, Reject_symm, mc.cores = cores))))
}
# 0.0542 0.0828 0.1813 0.3712 0.6015
# Correlation (Example S4)
# Method: spatial signed rank
for (i in 1:1) {
  Reject_symm = function(id) {
    set.seed(id)
    n = 200
    p = 2
    x <- MASS::mvrnorm(n, c(0,0), matrix(c(2,1,1,3),2,2))
    if (MNM::mv.1sample.test(x, score = "rank", method = "signchange")$p.value <= 0.05) {
      return(1)
    }
    return(0)
  }
  print(mean(unlist(parallel::mclapply(1:rep, Reject_symm, mc.cores = cores))))
}
# 0.0493

# Correlation (Example S4)
# Method: spatial sign
for (i in 1:1) {
  Reject_symm = function(id) {
    set.seed(id)
    n = 200
    p = 2
    x <- MASS::mvrnorm(n, c(0,0), matrix(c(2,1,1,3),2,2))
    if (MNM::mv.1sample.test(x, score = "sign", method = "signchange")$p.value <= 0.05) {
      return(1)
    }
    return(0)
  }
  print(mean(unlist(parallel::mclapply(1:rep, Reject_symm, mc.cores = cores))))
}
# 0.0544

# Laplace (Example S7)
# Method: spatial signed rank
for (i in 1:1) {
  Reject_symm = function(id) {
    set.seed(id)
    n = 200
    x <- matrix(VGAM::rlaplace(2*n, 0.2, 1), n, 2)
    if (MNM::mv.1sample.test(x, score = "rank", method = "signchange")$p.value <= 0.05) {
      return(1)
    }
    return(0)
  }
  print(mean(unlist(parallel::mclapply(1:rep, Reject_symm, mc.cores = cores))))
}
# 0.8284
# Laplace (Example S7)
# Method: spatial sign
for (i in 1:1) {
  Reject_symm = function(id) {
    set.seed(id)
    n = 200
    x <- matrix(VGAM::rlaplace(2*n, 0.2, 1), n, 2)
    if (MNM::mv.1sample.test(x, score = "sign", method = "signchange")$p.value <= 0.05) {
      return(1)
    }
    return(0)
  }
  print(mean(unlist(parallel::mclapply(1:rep, Reject_symm, mc.cores = cores))))
}
# 0.877

# Exponential distribution (Example S5)
# Method: spatial signed rank
for (i in 1:1) {
  Reject_symm = function(id) {
    set.seed(id)
    n = 100
    x <- matrix(rexp(2*n)-1,n,2)
    if (MNM::mv.1sample.test(x, score = "rank", method = "signchange")$p.value <= 0.05) {
      return(1)
    }
    return(0)
  }
  print(mean(unlist(parallel::mclapply(1:rep, Reject_symm, mc.cores = cores))))
}
# 0.3372
# Exponential distribution (Example S5)
# Method: spatial sign
for (i in 1:1) {
  Reject_symm = function(id) {
    set.seed(id)
    n = 100
    x <- matrix(rexp(2*n)-1,n,2)
    if (MNM::mv.1sample.test(x, score = "sign", method = "signchange")$p.value <= 0.05) {
      return(1)
    }
    return(0)
  }
  print(mean(unlist(parallel::mclapply(1:rep, Reject_symm, mc.cores = cores))))
}
# 0.7912

# Chisq distribution (Example S6)
# Method: spatial signed rank
for (i in 1:1) {
  Reject_symm = function(id) {
    set.seed(id)
    n = 100
    x = matrix(rchisq(2*n, df=1, ncp = 0)-1,n,2)
    if (MNM::mv.1sample.test(x, score = "rank", method = "signchange")$p.value <= 0.05) {
      return(1)
    }
    return(0)
  }
  print(mean(unlist(parallel::mclapply(1:rep, Reject_symm, mc.cores = cores))))
}
# 0.5914
# Chisq distribution (Example S6)
# Method: spatial sign
for (i in 1:1) {
  Reject_symm = function(id) {
    set.seed(id)
    n = 100
    x = matrix(rchisq(2*n, df=1, ncp = 0)-1,n,2)
    if (MNM::mv.1sample.test(x, score = "sign", method = "signchange")$p.value <= 0.05) {
      return(1)
    }
    return(0)
  }
  print(mean(unlist(parallel::mclapply(1:rep, Reject_symm, mc.cores = cores))))
}
# 0.9752

# High-dim H0 (Example S8)
# Method: spatial signed rank
for (i in 1:1) {
  Reject_symm = function(id) {
    set.seed(id)
    n = 200
    p = 50
    x = matrix(rnorm(n*p),n,p)*matrix(sin(1:p),n,p,byrow = TRUE)
    if (MNM::mv.1sample.test(x, score = "rank", method = "signchange")$p.value <= 0.05) {
      return(1)
    }
    return(0)
  }
  print(mean(unlist(parallel::mclapply(1:rep, Reject_symm, mc.cores = cores))))
}
# 0.0517
# High-dim H0 (Example S8)
# Method: spatial sign
for (i in 1:1) {
  Reject_symm = function(id) {
    set.seed(id)
    n = 200
    p = 50
    x = matrix(rnorm(n*p),n,p)*matrix(sin(1:p),n,p,byrow = TRUE)
    if (MNM::mv.1sample.test(x, score = "sign", method = "signchange")$p.value <= 0.05) {
      return(1)
    }
    return(0)
  }
  print(mean(unlist(parallel::mclapply(1:rep, Reject_symm, mc.cores = cores))))
}
# 0.0514

# High-dim H1 (Example S9)
# Method: spatial signed rank
for (i in 1:1) {
  Reject_symm = function(id) {
    set.seed(id)
    n = 200
    p = 50
    x = matrix(rnorm(n*p),n,p)*matrix(sin(1:p),n,p,byrow = TRUE) + 0.003
    if (MNM::mv.1sample.test(x, score = "rank", method = "signchange")$p.value <= 0.05) {
      return(1)
    }
    return(0)
  }
  print(mean(unlist(parallel::mclapply(1:rep, Reject_symm, mc.cores = cores))))
}
# 0.6648
# High-dim H1 (Example S9)
# Method: spatial sign
for (i in 1:1) {
  Reject_symm = function(id) {
    set.seed(id)
    n = 200
    p = 50
    x = matrix(rnorm(n*p),n,p)*matrix(sin(1:p),n,p,byrow = TRUE) + 0.003
    if (MNM::mv.1sample.test(x, score = "sign", method = "signchange")$p.value <= 0.05) {
      return(1)
    }
    return(0)
  }
  print(mean(unlist(parallel::mclapply(1:rep, Reject_symm, mc.cores = cores))))
}
# 0.656

# High-dim H1 heavy-tailed (Example S10)
# Method: spatial signed rank
for (i in 1:1) {
  Reject_symm = function(id) {
    set.seed(id)
    n = 200
    p = 50
    x = matrix(rnorm(n*p),n,p)/sqrt(rchisq(n, df=1))*matrix(sin(1:p),n,p,byrow = TRUE) + 0.003
    if (MNM::mv.1sample.test(x, score = "rank", method = "signchange")$p.value <= 0.05) {
      return(1)
    }
    return(0)
  }
  print(mean(unlist(parallel::mclapply(1:rep, Reject_symm, mc.cores = cores))))
}
# 0.2912
# High-dim H1 heavy-tailed (Example S10)
# Method: spatial sign
for (i in 1:1) {
  Reject_symm = function(id) {
    set.seed(id)
    n = 200
    p = 50
    x = matrix(rnorm(n*p),n,p)/sqrt(rchisq(n, df=1))*matrix(sin(1:p),n,p,byrow = TRUE) + 0.003
    if (MNM::mv.1sample.test(x, score = "sign", method = "signchange")$p.value <= 0.05) {
      return(1)
    }
    return(0)
  }
  print(mean(unlist(parallel::mclapply(1:rep, Reject_symm, mc.cores = cores))))
}
# 0.3808

