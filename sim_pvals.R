

library(foreach)
library(doParallel)
library(tidyverse)

sim_pvals <- function(path, mm_x, ns, num_cores){
  registerDoParallel(cores=num_cores)
  # mm_x: log_min_max ratios at which to calculate p-values
  # settings we used: 
  # mm_x <- seq(0, 5, 0.01)
  
  # ns: Number of reads for which to calculate p-values
  # settings we used: 
  #ns <- seq(5, 500, 1)
  
  i <- 1
  
  print("Starting loop")
  ptm <- proc.time()
  
  df_p <- foreach(n = ns, .combine=rbind) %dopar%{
    set.seed(123)
    print(paste0("Read Count ", n))
    write(paste0("Read Count ",n), 
          file = paste0(path, "log_make_pval_table.txt"),
                        append = TRUE, sep = " ")
          
          x <- rmultinom(1e8, n, rep(0.25, 4))
          min24 <- pmin(x[2,], x[4,])
          max13 <- pmax(x[1,], x[3,])
          # add 1 to avoid 0s
          log_min_max <- log2((min24+1) / (max13+1))
          # we only care about when min_max >= 1
          log_MM_greater <- log_min_max[log_min_max >=0]
          # get the empirical CDF
          mm_ecdf <- ecdf(log_MM_greater)
          
          # The CDF (x) calculates p(X<=x).
          # For the p-value we need p(X>=x) = 1-F(x) + p(X=x).
          # We can get this by evaluating 1 - F(x - epsilon)
          
          data.frame(num_reads = n,
                     log_min_max = mm_x,
                     min_max = 2^mm_x,
                     p = 1 - mm_ecdf(mm_x - 1e-10))
  }
  saveRDS(df_p, file = paste0(path, "df_p_1e8.RDS"))
  return(df_p)
}
