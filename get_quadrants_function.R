library(Rsamtools)
library(foreach)
library(doParallel)
library(data.table)
library(plyr)
library(tidyverse)

get_quadrant_reads <- function(regions, step, pairs, num_cores,pval_table){
  # first subset the pairs file  
  print("Index pairs")
  
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
  print("Break up regions into 150bp bins")
  registerDoParallel(cores=num_cores)
  
  region_start_seq <- foreach(i=1:nrow(regions))%dopar%{
    data.frame(
      chr = regions$chr1[i],
      start=seq(from = regions$s1[i]-75, to = regions$e1[i]-75, by = step),
      region_id=i) 
  }
  
  region_start_seq <- bind_rows(region_start_seq)
  print("Done extracting 150bp bins")
  
  print("Make GRanges")
  
  # Q2
  window_gr_POS <- GRanges(region_start_seq$chr, IRanges(start=region_start_seq$start, 
                                                         end= region_start_seq$start + 74)) 
  # 74 instead of 75 to prevent overlap with window_gr_left_POS_control
  
  # Q4
  window_gr_NEG <- GRanges(region_start_seq$chr, IRanges(start=region_start_seq$start + 76, 
                                                         end= region_start_seq$start + 150))
  ##### control #####
  #Q1
  window_gr_POS_control <- GRanges(region_start_seq$chr, IRanges(start=region_start_seq$start + 75, 
                                                                 end= region_start_seq$start + 150))
  #Q3
  window_gr_NEG_control <- GRanges(region_start_seq$chr, IRanges(start=region_start_seq$start, 
                                                                 end= region_start_seq$start + 75))
  
  print("Find Overlaps")
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
  print("Extract Quadrant Read Counts")
  # Find peaks
  results <- data.frame(ID = 1:nrow(region_start_seq), 
                             chr = region_start_seq$chr, 
                             window_start = region_start_seq$start, 
                             window_mid = region_start_seq$start + 75, 
                             window_end = region_start_seq$start + 150,
                             region_id =region_start_seq$region_id,
                             
                             q1 = rep(0,nrow(region_start_seq)), 
                             q2 = rep(0,nrow(region_start_seq)), 
                             q3 = rep(0, nrow(region_start_seq)),
                             q4 = rep(0, nrow(region_start_seq)))
  
  
  results$q2[as.numeric(names(out_anchor_POS))] <- lengths(out_anchor_POS) 
  results$q4[as.numeric(names(out_anchor_NEG))] <- lengths(out_anchor_NEG) 
  
  results$q1[as.numeric(names(out_anchor_POS_control))] <- lengths(out_anchor_POS_control) 
  results$q3[as.numeric(names(out_anchor_NEG_control))] <- lengths(out_anchor_NEG_control)
  
  results$num_reads = rowSums(results[,c("q1","q2","q3","q4")])
  
  regions$region_id <- 1:nrow(regions)
  results <- left_join(results, regions,by = "region_id" )  
  
  results$min24 <- pmin(results$q2, results$q4)
  results$max13 <- pmax(results$q1, results$q3)
  results$log_min_max <- log2((results$min24 +1) / (results$max13 +1))
  results$round_log_min_max <- round_any(results$log_min_max, 0.01, f = floor)
  
  results$min_max <-  2^results$log_min_max
  results$num_reads[results$num_reads >= 500] <- 500
  
  results <- results %>% 
    mutate(label = paste(num_reads, round_log_min_max ))
  
  #df_p <- readRDS(file = "/aryeelab/users/corri/results/June02_2023_df_p_1e8.RDS")
  df_p <- readRDS(file = pval_table)
  
  
  df_p<- df_p %>% 
    mutate(label = paste(num_reads, log_min_max )) %>% 
    dplyr::select(p, label) %>% 
    dplyr::rename(pvalue = p)
  
  results <- left_join(results, df_p, by = "label")
  
  results$pvalue[results$min_max <= 1] <- 1
  results$pvalue[results$num_reads < 5] <- 1
  results$pvalue[is.na(results$pvalue)] <- 0
  
  print("Returning results")
  
  return(results)
  
}
