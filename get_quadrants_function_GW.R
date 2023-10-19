# genomewide get reads

library(Rsamtools)
library(foreach)
library(doParallel)
library(data.table)
library(plyr)
library(tidyverse)

# my settings

# chrom_seq = c(1:22, "X")
# 
# dir <- "/aryeelab/users/corri/data/CTCF_chr_pairs/"
# prefix <- "k562_ctcf"
# results_dir <- "/aryeelab/users/corri/data/replicate_FF_results/"
# results_prefix <- "Read_count_quadrants_chr"
# 
# step = 10
# num_cores = 8

get_quadrant_reads_GW <- function(chrom_seq, dir, prefix, results_dir, results_prefix, step, 
                                  num_cores, pval_table){
  registerDoParallel(cores=num_cores)
  
  foreach ( chrom_num = chrom_seq) %dopar% {
    print(paste0("Start chr ", chrom_num))
    
    print(paste0("Read pairs: chr", chrom_num))
    pairs <- readRDS(paste0(dir, prefix, "_chr", chrom_num, ".RDS"))
    
    full_region <- data.frame(chr1 = paste0("chr",chrom_num), s1 = min(pairs$pos1), e1 = max(pairs$pos2))
    
    
    print(paste0("Index pairs: chr", chrom_num))
    
    left_anchor <- pairs %>% 
      filter(type %in% c("uu", "uU", "UR")) %>% 
      select(chr1, pos1, strand1, type) %>% 
      dplyr::rename(chr = chr1, pos = pos1, strand = strand1) 
    
    right_anchor <- pairs %>% 
      filter(type %in% c("uu", "Uu", "RU")) %>% 
      select(chr2, pos2, strand2, type) %>% 
      dplyr::rename(chr = chr2, pos = pos2, strand = strand2) 
    
    pairs <- rbind(left_anchor, right_anchor)
    
    pairs_gr <- GRanges(pairs$chr, IRanges(start=pairs$pos, end=pairs$pos), 
                        strand = pairs$strand, type = pairs$type)
    
    # get rid of things taking up memory in workspace. only need pairs_gr
    rm(left_anchor, right_anchor, pairs)
    
    # positive
    pairs_gr_POS <- pairs_gr[pairs_gr@strand == "+"]
    # negative
    pairs_gr_NEG <- pairs_gr[pairs_gr@strand == "-"]
    
    ###################################################
    # Break up regions into 150bp bins
    
    region_start_seq <- seq(from = full_region$s1[1]-75, to = full_region$e1[1]-75, by = step)
    
    nrow <- length(region_start_seq)
    
    # Q2
    window_gr_POS <- GRanges(rep(paste0("chr",chrom_num), nrow), IRanges(start=region_start_seq, 
                                                                              end= region_start_seq + 74)) 
    # 74 instead of 75 to prevent overlap with window_gr_POS_control
    
    # Q4
    window_gr_NEG <- GRanges(rep(paste0("chr",chrom_num), nrow), IRanges(start=region_start_seq + 76, 
                                                                              end= region_start_seq + 150))
    ##### control #####
    #Q1
    window_gr_POS_control <- GRanges(rep(paste0("chr",chrom_num), nrow), IRanges(start=region_start_seq + 75, 
                                                                                      end= region_start_seq + 150))
    #Q3
    window_gr_NEG_control <- GRanges(rep(paste0("chr",chrom_num), nrow), IRanges(start=region_start_seq, 
                                                                                      end= region_start_seq + 75))
    print(paste0("Find overlaps: chr", chrom_num))
    #### find overlaps ####
    ovl_anchor_POS <- findOverlaps(window_gr_POS, pairs_gr_POS)
    out_anchor_POS <- tapply(subjectHits(ovl_anchor_POS), queryHits(ovl_anchor_POS), I)
    
    ovl_anchor_NEG <- findOverlaps(window_gr_NEG, pairs_gr_NEG)
    out_anchor_NEG <- tapply(subjectHits(ovl_anchor_NEG), queryHits(ovl_anchor_NEG), I)
    
    ##### control #####
    ovl_anchor_POS_control <- findOverlaps(window_gr_POS_control, pairs_gr_POS)
    out_anchor_POS_control <- tapply(subjectHits(ovl_anchor_POS_control), queryHits(ovl_anchor_POS_control), I)
    
    ovl_anchor_NEG_control <- findOverlaps(window_gr_NEG_control, pairs_gr_NEG)
    out_anchor_NEG_control <- tapply(subjectHits(ovl_anchor_NEG_control), queryHits(ovl_anchor_NEG_control), I)
    
    
    ######################
    print(paste0("Extract Quadrant Read Counts: chr", chrom_num))
    # Find peaks
    results <- data.frame(ID = 1:nrow, 
                               chr = rep(paste0("chr",chrom_num), nrow), 
                               window_start = region_start_seq, 
                               window_mid = region_start_seq + 75, 
                               window_end = region_start_seq + 150, 
                               
                               q1 = rep(0,nrow), 
                               q2 = rep(0,nrow), 
                               q3 = rep(0, nrow),
                               q4 = rep(0, nrow))
    
    results$q2[as.numeric(names(out_anchor_POS))] <- lengths(out_anchor_POS) 
    results$q4[as.numeric(names(out_anchor_NEG))] <- lengths(out_anchor_NEG) 
    
    results$q1[as.numeric(names(out_anchor_POS_control))] <- lengths(out_anchor_POS_control) 
    results$q3[as.numeric(names(out_anchor_NEG_control))] <- lengths(out_anchor_NEG_control)
    
    results$num_reads = rowSums(results[,c("q1","q2","q3","q4")])
    
    ########################
    print(paste0("Get p-values: chr", chrom_num))
    
    
    results$min24 <- pmin(results$q2, results$q4)
    results$max13 <- pmax(results$q1, results$q3)
    results$log_min_max <- log2((results$min24 +1) / (results$max13 +1))
    results$round_log_min_max <- round_any(results$log_min_max, 0.01, f = floor)
    
    results$min_max <-  2^results$log_min_max
    results$num_reads[results$num_reads >= 500] <- 500
    
    results <- results %>% 
      mutate(label = paste(num_reads, round_log_min_max ))
    
#    df_p <- readRDS(file = "/aryeelab/users/corri/results/June02_2023_df_p_1e8.RDS")
    df_p <- readRDS(file = pval_table)
    
    df_p<- df_p %>% 
      mutate(label = paste(num_reads, log_min_max )) %>% 
      dplyr::select(p, label) %>% 
      dplyr::rename(pvalue = p)
    
    results <- left_join(results, df_p, by = "label")
    
    results$pvalue[results$min_max <= 1] <- 1
    results$pvalue[results$num_reads < 5] <- 1
    results$pvalue[is.na(results$pvalue)] <- 0
    
    
    print(paste0("Exporting results for chr", chrom_num))
    saveRDS(results, file = paste0(results_dir, results_prefix, chrom_num, ".RDS"))
  }
}
