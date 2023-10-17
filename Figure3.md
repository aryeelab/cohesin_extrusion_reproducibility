Fig. 3 CTCF binding sites identified by FactorFinder with single
basepair resolution in MNase K562 CTCF HiChIP data.
================

- [Set the path](#set-the-path)
- [Figure 3A](#figure-3a)
- [Figure 3B](#figure-3b)
- [Figure 3C](#figure-3c)
- [Figure 3D](#figure-3d)
- [Figure 3E](#figure-3e)
- [Figure 3F](#figure-3f)
- [Figure 3G](#figure-3g)

``` r
library(foreach)
library(GenomicRanges)
library(readr)
library(rtracklayer)
library(scales)
library(plyranges)
library(nlme)
library(lme4)
library(viridis)
library(ggridges)
library(RColorBrewer)
library(cowplot)
library(gridExtra)
library(doParallel)
library(data.table)
library(Signac)
library(plyr)
library(dplyr)
library(ggplot2)
library(tidyverse)

num_cores = 8
registerDoParallel(cores=num_cores)
```

# Set the path

This is where we stored all those files in the protocol.

``` r
path <- "/aryeelab/users/corri/data/replicate_FF_results/"
```

# Figure 3A

**Get true positive / true negative set for genomewide pvals**

``` r
loops <-read_tsv("/aryeelab/users/corri/data/K562_CTCF_2.5kb.interactions_FitHiC_Q0.01.bed")
loops <- loops %>% 
  filter(chr1 != "chrX") %>% 
  filter(chr2 != "chrX")

LA <- loops %>% 
  dplyr::select(chr1, s1, e1) %>% 
  dplyr::rename(chr = chr1, start = s1, end = e1)

RA <- loops %>% 
  dplyr::select(chr2, s2, e2) %>% 
  dplyr::rename(chr = chr2, start = s2, end = e2)

anchors <- rbind(LA,RA) %>% 
  mutate(loc = paste0(chr, ":",start)) %>% 
  distinct(loc, .keep_all = TRUE)

anchors_gr <- makeGRangesFromDataFrame(anchors, 
                                       keep.extra.columns = TRUE,
                                       seqnames.field = "chr",
                                       start.field = "start",
                                       end.field = "end")
```

``` r
ctcf_motifs<- readRDS("/aryeelab/users/corri/data/ALL_FIMO_CTCF_hg38.RDS")
gr_motifs <- ctcf_motifs[width(ctcf_motifs)==19]
gr_motifs$motif_mid <- round( ( start(gr_motifs) + end(gr_motifs) )/2)
```

``` r
ovl <- findOverlaps(gr_motifs, anchors_gr, maxgap = 0)

out <- as.data.frame(gr_motifs[queryHits(ovl)]) %>% 
  cbind(motif_id = queryHits(ovl),
        loop_id = subjectHits(ovl),
        loop_start = start(anchors_gr)[subjectHits(ovl)],
        loop_end = end(anchors_gr)[subjectHits(ovl)]) 
```

*Filter to loop anchors with only one CTCF motif.*

``` r
uniq<- out %>% 
  dplyr::group_by(loop_id) %>% 
  dplyr::mutate(count = max(dplyr::row_number())) %>% 
  ungroup() %>% 
  filter(count == 1) %>% 
  arrange(loop_id)

uniq_gr <- makeGRangesFromDataFrame(uniq, 
                                    keep.extra.columns = FALSE,
                                    seqnames.field = "seqnames",
                                    start.field = "start",
                                    end.field = "end")
uniq_gr$motif_mid <- uniq$motif_mid
uniq_gr$loop_start <- uniq$loop_start
uniq_gr$loop_end <- uniq$loop_end
```

*The CTCF motif has to be within 30bp of a CTCF ChIP-seq peak summit.*

``` r
chip <- read_tsv("/aryeelab/users/corri/data/K562_CTCF_peaks_ENCFF736NYC.bed", col_names = FALSE)
colnames(chip) <- c("chr", "start", "end", "name", "score", "strand", "signalValue", "pval", "qval", "peak")
chip_gr <- makeGRangesFromDataFrame(chip)
chip_gr$qval <- chip$qval
chip_gr$peak_mid <- (start(chip_gr) + end(chip_gr))/2

anchors <- chip_gr %>% 
  plyranges::anchor_center() %>% 
  plyranges::mutate(width = 301)
```

``` r
ovl <- findOverlaps(uniq_gr, anchors, maxgap = 0)

out <- as.data.frame(uniq_gr[queryHits(ovl)]) %>% 
  cbind(motif_id = queryHits(ovl),
        chip_id = subjectHits(ovl),
        chip_start = start(anchors)[subjectHits(ovl)],
        chip_mid = anchors$peak_mid[subjectHits(ovl)],
        chip_end = end(anchors)[subjectHits(ovl)],
        qval_score = anchors$qval[subjectHits(ovl)]) %>% 
  mutate(dist = chip_mid - motif_mid)

out <- out %>% 
  ungroup() %>% 
  group_by(chip_id) %>% 
  arrange(abs(dist)) %>% 
  dplyr::slice(1)

out <- out %>% 
  filter(abs(dist) < 30)
```

``` r
loop_anc <- out %>%
  ungroup() %>% 
  dplyr::select(seqnames, motif_mid) %>%
  mutate(chr1 = seqnames,
         s1 = motif_mid - 100,
         e1 = motif_mid + 100)
nrow(loop_anc)
```

``` r
pairs <- readRDS("/aryeelab/users/corri/data/k562_ctcf_mapped.pairs_STRAND_TYPE_X.rds")
pairs <- pairs %>% 
  filter(chr1 != "chrX") # 386,874,029 rows
```

*took 5 minutes*

``` r
source("/aryeelab/users/corri/code/mnase-hichip/code/get_quadrants_function.R")
path <- "/aryeelab/users/corri/data/replicate_FF_results/"
pval_table <- paste0(path, "df_p_1e8.RDS")
reads <- get_quadrant_reads(regions = loop_anc, step = 1, pairs = pairs, num_cores=8,pval_table)
```

``` r
results_anc1 <- reads

results_anc1<- results_anc1 %>% 
  mutate(dist = window_mid - motif_mid)

order <- results_anc1 %>% 
  ungroup() %>% 
  dplyr::group_by(region_id) %>%
  dplyr::summarize(max_signal = max(log_min_max),
                   signal = sum(log_min_max)) %>% 
  arrange(desc(max_signal), desc(signal)) %>% 
  dplyr::mutate(order = 1:nrow(loop_anc)) # 4523

results_anc1<- left_join(results_anc1,order, by = "region_id")
```

``` r
saveRDS(results_anc1, file = paste0(path,"heatmap_sort_readcounts.RDS"))
```

``` r
results_anc1 <- readRDS(file = paste0(path,"heatmap_sort_readcounts.RDS"))

results_anc1$log_min_max_thresh <- results_anc1$log_min_max
results_anc1$log_min_max_thresh[results_anc1$log_min_max_thresh > 3] <- 3
results_anc1$log_min_max_thresh[results_anc1$log_min_max_thresh < -3] <- -3

results_anc1 %>% 
  ggplot(aes(dist, order, fill=log_min_max_thresh )) + 
  geom_tile() +
  ylab("log2 (min/max)")  +
  xlab("Distance from CTCF Motif")+
  labs(fill = "log2 (min/max)")+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "#f0f0f0",
                                        size = 0.5))+
  theme(plot.title = element_text(color = "black", family = "Times New Roman", size = 18, face = "bold"),
        axis.text.x = element_text(color = "black", family = "Times New Roman", size = 18,face = "bold"),
        axis.text.y = element_text(color = "black", family = "Times New Roman", size = 18,face = "bold"),
        axis.title.x = element_text(color = "black", family = "Times New Roman", size = 18,face = "bold"),
        axis.title.y = element_text(color = "black", family = "Times New Roman", size = 18,face = "bold"),
        axis.ticks.x=element_blank(),
        legend.text=element_text(color = "black", family = "Times New Roman", size =18,face = "bold"),
        legend.title=element_text(color = "black", family = "Times New Roman", size = 18,face = "bold"),
        strip.text.x = element_text(color = "black", family = "Times New Roman", size = 18,face = "bold"),
        strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"))+
  scale_fill_gradient2(
    low = "darkblue",
    mid = "white",
    high = "darkred")+
  geom_vline(xintercept = c(0))+
  scale_x_continuous(breaks=seq(-100,100,20))

ggsave(paste0(path,"3A_min_max_heatmap.png"), width=8, height=4)
```

<img src="Figures/3A_min_max_heatmap.png" width="100%" />

# Figure 3B

# Figure 3C

# Figure 3D

# Figure 3E

# Figure 3F

# Figure 3G
