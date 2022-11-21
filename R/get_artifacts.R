#!/usr/bin/env Rscript

get_bin_pval <- function(mmm_count, depth, mis, typ){
  ret <- binom.test(x = mmm_count,
                    n = depth,
                    p = cov_totals_mm$probs[cov_totals_mm$mismatch == mis & cov_totals_mm$type == typ],
                    alternative = "greater",
                    conf.level = 0.95)
  return(ret$p.value)
} # returns binominal probability that a given count is an artifact

filter_probable_mismatches <- function(df, p.val.thshld){
  df <- df %>% dplyr::filter(p.adj.all >= p.val.thshld)
  return(df)

} # just filters the df based on the column values

coverages_cumulative <- vector(mode = "list")
counts_all_raw <- vector(mode = "list")
counts_read_bins <- vector(mode = "list")
counts_all_updated_wide <- vector(mode = "list")

kartoteki <- list.files(full.names = F, pattern = "testdata$")

for (i in 1:length(kartoteki)){
    # go to the i-th folder
  setwd(kartoteki[i])
    # print its name
  print(getwd())
  # store sample name variable
  sample_name <- kartoteki[i]
  # and print it
  print(sample_name)

  # write mismatch_positions_RDS file
  assign(paste0(sample_name, "_mismatch_positions"), readRDS(file.path(getwd(), paste0(sample_name, "_mismatch_positions.rds")))[, -c(8,9)])

  write.table(eval(parse(text=paste0(sample_name, "_mismatch_positions"))), paste0(sample_name, "_mismatch_positions.txt"), row.names = F, col.names = F, quote = F, sep = "\t")

  # system2(command = "/Users/lecka48/Documents/R_WORK/RMLizer/R/bedscript.sh", args = sample_name, stdout = "NULL")

coverage.bed <- as_tibble(read.table(paste0(sample_name, "_coverage.bed"),
                                     header = F,
                                     as.is = T,
                                     col.names = c("chr", "start","end", "gene_id", "score", "strand", "mismatch", "is_dup", "coverage"))) %>%
    mutate(type = case_when(
      gene_id %like% "ENSMUS" ~ "Genes",
      gene_id %like% "ERCC"   ~ "ERCC",
      gene_id %like% "slam"   ~ "Spike-In",
      TRUE ~ "Other"
    ))

mismatch_positions_orig <- readRDS(file.path(getwd(), paste0(sample_name, "_mismatch_positions.rds"))) %>% dplyr::ungroup() %>% dplyr::mutate(across(start:end, as.integer)) %>% dplyr::filter(is_dup == F)

coverage.bed <-  inner_join(coverage.bed, mismatch_positions_orig)

coverage.bed <- coverage.bed %>% filter(coverage != "0", mmm_count <= coverage)

# strandedness is problematic here; probably take the original file "mismatch_positions_sorted...." (currently being deleted) and the do reduction including the coverage

coverage.bed <- coverage.bed %>% dplyr::mutate(refbase = gsub(".{1}$", "", mismatch))

# fine as initial coverage values to compute initial rates; then read filtration is changed, so that
# individual reads are compared for mismatches on matching positions, and overlaps are filtered in this way, then the
# "filtered out" count can be subtracted from the obtained value
coverage_quick_per_position <- as_tibble(read.table(paste0(sample_name,"_cov_totals.txt"),
                                                    header = F,
                                                    as.is = T,
                                                    col.names = c("chr", "refbase","gene_id", "strand", "is_dup", "coverage"))) %>%
  filter(is_dup == F) %>%
  mutate(type = case_when(
    gene_id %like% "ENSMUS" ~ "Genes",
    gene_id %like% "ERCC" ~ "ERCC",
    gene_id %like% "slam" ~ "Spike-In",
    TRUE ~ "Other"
  ))

coverage_totals_no_strand <- coverage_quick_per_position %>%
  group_by(refbase, type) %>%
  summarise(coverage = sum(coverage))

coverage_totals_no_strand_genes <- coverage_totals_no_strand %>% filter(type == "Genes")

number_of_tests_ACTG <- c(length(unique(coverage.bed$end[coverage.bed$refbase == "A"])),
                          length(unique(coverage.bed$end[coverage.bed$refbase == "C"])),
                          length(unique(coverage.bed$end[coverage.bed$refbase == "T"])),
                          length(unique(coverage.bed$end[coverage.bed$refbase == "G"]))
)

names(number_of_tests_ACTG) <- c("A","C","T","G")

# adding total coverage to the table
cov_totals_mm <- coverage.bed %>% group_by(mismatch, type, refbase) %>% summarise(mmm_count = sum(mmm_count))
cov_totals_mm <- inner_join(cov_totals_mm, coverage_totals_no_strand) %>% mutate(probs = mmm_count/coverage)
cov_totals_mm$sample_name <- sample_name
coverages_cumulative[[i]] <- cov_totals_mm
# coverage breadth needed
# %>% head(., 300000)
coverage.bed <- coverage.bed %>% rowwise() %>%
  mutate(observed_prob = mmm_count/coverage,
         p.value = get_bin_pval(mmm_count = mmm_count,
                                depth = coverage,
                                mis = mismatch,
                                typ = type)
         ) %>%
  # group_by(gene_id) %>%  # no gene specificity at this stage
  mutate(p.adj.mut = p.adjust(p.value,
                              method = "bonferroni",
                              n = number_of_tests_ACTG[refbase]),
         p.adj.all = p.adjust(p.value,
                              method = "bonferroni",
                              n = coverage_totals_no_strand$coverage[which(coverage_totals_no_strand$refbase == refbase & coverage_totals_no_strand$type == type)])
         )

coverage.bed <- coverage.bed %>% unnest(rname_lst) %>% rename(read_name = rname_lst)

write.table(coverage.bed %>% mutate(sample_name = sample_name), paste0(sample_name, "_coverage.bed.txt"), row.names = F, col.names = F, quote = F, sep = "\t")

# remove whenever coverage == 1; this may be reconsidered, though it is not possible to determine if such position is not an artifact; this can be extended to 2, or 3 (majority consensus maybe?)
#
coverage.bed.cov1out <- coverage.bed %>% dplyr::filter(mmm_count < coverage) # zeros already filtered before; this gets rid of samples with mutation with coverage of 1

coverage.bed.flowthrough <- filter_probable_mismatches(coverage.bed.cov1out, 0.1)

# this contains the reads that should be used to select good reads later
coverage.bed.filtered <- dplyr::anti_join(coverage.bed, coverage.bed.flowthrough)

######## noise regression ############

coverage.bed.filtered <- coverage.bed.filtered %>% group_by(read_name, mismatch) %>% summarise(mm.count.per.read.per.mm = n())

# coverage.bed.fltrd.unnestrn <- coverage.bed.fltrd %>% unnest(rname_lst) %>% select(read_name = rname_lst, mismatch = mismatch, mmm_count)

counts_all <- as_tibble(read.table(paste0(sample_name, "_RMLizer_counts.txt.gz"), header = T))
# counts_allx <- as_tibble(read.table("RTA0h_R1_subsampled_RMLizer_counts.txt"), header = T)

counts_all_raw[[i]] <- counts_all %>% mutate(sample_name = "sample_name")

counts_all <- (counts_all %>% filter(OVER == "a", DUP == F))[ , setdiff(colnames(counts_all), c("OVER","DUP"))] # select may be better


counts_all <- counts_all %>% pivot_longer(AA:TT, names_to = "mismatch", values_to = "count")

counts_all_mm_not_zero <- counts_all %>% filter(!mismatch %like% "AA|TT|GG|CC", count != 0)

counts_all_remainder <- dplyr::anti_join(counts_all, counts_all_mm_not_zero)

# this may be fixed in a better way to prevent negative values popping up;
counts_all_mm_not_zero_updated <- dplyr::left_join(counts_all_mm_not_zero, coverage.bed.filtered) %>%
  mutate(count = case_when(is.finite(mm.count.per.read.per.mm) & count >= mm.count.per.read.per.mm ~ count-mm.count.per.read.per.mm, TRUE ~ count)) %>%
  select(1:4)
# a <- dplyr::left_join(counts_all_mm_not_zero, coverage.bed.filtered) %>% filter(count - mm.count.per.read.per.mm < 0)
counts_all_updated <- rbind(counts_all_remainder, counts_all_mm_not_zero_updated)

counts_all_updated_wide[[i]] <- counts_all_updated %>%
  tidyr::pivot_wider(names_from = "mismatch", values_from = "count") %>%
  dplyr::mutate(sample_name = sample_name)

counts_read_bins[[i]] <- counts_all_updated %>% filter(mismatch %like% "^T.+") %>% group_by(read_name) %>% mutate(covBase = sum(count)) %>% filter(mismatch == "TC") %>% group_by(gene_id, covBase, count) %>% summarise(observed = n()) %>% dplyr::mutate(sample_name = sample_name)

setwd("../")

}




counts_read_bins <- counts_all_updated %>% filter(mismatch %like% "^T.+") %>% group_by(read_name) %>% mutate(covBase = sum(count)) %>% filter(mismatch == "TC") %>% group_by(gene_id, covBase, count) %>% summarise(observed = n()) %>% dplyr::mutate(sample_name = sample_name)

counts_read_bins <- do.call("rbind", counts_read_bins)

# maybe bin reads and not mismatches???????
# additonal filtering:
# --> at least
#
counts_all_updated_test <- counts_all_updated[sample(10000, nrow(counts_all_updated)),]





# compare expected frequency with observed frequency
counts_read_bins_no_gene <- counts_all_updated %>% #head(., 20000) %>%
    dplyr::filter(gene_id %like% "ENSMUS") %>%
    dplyr::mutate(refbase = gsub(".{1}$", "", mismatch)) %>%
  # group_by(read_name) %>%
  #   mutate(covRead = sum(count)) %>%
  group_by(read_name, refbase) %>%
    mutate(covBase = sum(count)) %>%
  group_by(covBase, count, mismatch) %>%
    summarise(observed = n())         %>%
  dplyr::mutate(sample_name = sample_name)





a <- counts_read_bins %>% group_by(gene_id) %>%
  filter(gene_id %like% "ENSMUSG") %>%
  mutate(reads_per_gene = sum(observed),
         covT_gene = sum(covBase*observed),
         sum_mis = sum(count*observed)) %>%
  mutate(TC_rate_per_gene = sum_mis/covT_gene,
         TC_prob = case_when(count == 0 ~ dbinom(2, covBase, 0.00034),
                             TRUE ~ dbinom(count, covBase, 0.00034)),
         TC_prob_one = pbinom(0, covBase, 0.00034, lower.tail = F)
         ) %>%
  mutate(total_to_remove = sum(TC_prob_one*reads_per_gene),
         obs_to_rem = sum(TC_prob*reads_per_gene))

  mutate(expected_TC = rbinom(count,covBase,non_TC_prob))



a$subtr <- a$non_TC_prob * a$observed

a <- a %>% group_by(gene_id) %>% mutate(to_subtr = sum(subtr))

    dbinom = dbinom(count, covBase, 0.00034),
         pbinom = pbinom(count, covBase, 0.00034)
         )





counts_all_updated_wide <- counts_all_updated %>% pivot_wider(names_from = "mismatch", values_from = "count")



########
########
counts_all_updated_wide_genes <- counts_all_updated_wide %>% filter(gene_id %like% "ENSMUS")


probs_per_mismatch_genes <- counts_all_updated_wide_genes %>% summarise(across(AA:GA, sum)) %>%
  dplyr::transmute(
    TA = TA / sum(c(TA, TC, TG, TT)),
    TG = TG / sum(c(TA, TC, TG, TT)),
    TC = TC / sum(c(TA, TC, TG, TT)),
    AC = AC / sum(c(AA, AC, AG, AT)),
    AG = AG / sum(c(AA, AC, AG, AT)),
    AT = AT / sum(c(AA, AC, AG, AT)),
    GA = GA / sum(c(GA, GC, GG, GT)),
    GC = GC / sum(c(GA, GC, GG, GT)),
    GT = GT / sum(c(GA, GC, GG, GT)),
    CA = CA / sum(c(CA, CC, CG, CT)),
    CG = CG / sum(c(CA, CC, CG, CT)),
    CT = CT / sum(c(CA, CC, CG, CT))
  ) %>% pivot_longer(everything())



counts_all_updated_plot <-
  counts_all_updated_wide %>% dplyr::filter(gene_id != "slam") %>% group_by(gene_id) %>% summarise(across(AA:TT, sum)) %>% dplyr::rowwise() %>%
  dplyr::transmute(
    gene_id = gene_id,
    `T->A` = TA / sum(c(TA, TC, TG, TT)) * 100,
    `T->G` = TG / sum(c(TA, TC, TG, TT)) * 100,
    `T->C` = TC / sum(c(TA, TC, TG, TT)) * 100,
    `A->C` = AC / sum(c(AA, AC, AG, AT)) * 100,
    `A->G` = AG / sum(c(AA, AC, AG, AT)) * 100,
    `A->T` = AT / sum(c(AA, AC, AG, AT)) * 100,
    `G->A` = GA / sum(c(GA, GC, GG, GT)) * 100,
    `G->C` = GC / sum(c(GA, GC, GG, GT)) * 100,
    `G->T` = GT / sum(c(GA, GC, GG, GT)) * 100,
    `C->A` = CA / sum(c(CA, CC, CG, CT)) * 100,
    `C->G` = CG / sum(c(CA, CC, CG, CT)) * 100,
    `C->T` = CT / sum(c(CA, CC, CG, CT)) * 100
  ) %>% tidyr::pivot_longer(`T->A`:`C->T`) %>% dplyr::mutate(
    highlight = case_when(gene_id %like% "ERCC" ~ "ERCC", TRUE ~ "Genes"),
    RefBase = case_when(
      name %like% "^T" ~ "T",
      name %like% "^G" ~ "G",
      name %like% "^C" ~ "C",
      name %like% "^A" ~ "A"
    )
  )

counts_all_updated_plot$highlight <- factor(counts_all_updated_plot$highlight,
                                levels = c("ERCC", "Genes"))


sample_reads <- function(x, n, t){
  z <- list()
  for (i in 1:t){
    y <- x[sample(nrow(x), n) , ]
    z[[i]] <- y %>% summarise(across(AA:TT, sum))
  }
  return(do.call("rbind", z))
}

samp_rds <- sample_reads(counts_all_updated_wide_genes, nrow(counts_all_updated_wide_genes), 1)

test <- sample_reads(counts_all_updated_wide_genes, 1000000, 5)


lm_model <- lm(TC~AG+AC+AT+GA+GC+GT+CA+CG+CT+0, data = samp_rds)

summary(lm_model)

test_data <- test %>% select(AG,AC,AT,GA,GC,GT,CA,CG,CT)

predict(lm_model, test_data)

rates_genes_plot <-
  ggplot(counts_all_updated_plot %>% dplyr::filter(highlight %like% "Genes|ERCC"),
    aes(
      x = name,
      y = value,
      fill = highlight,
      color = highlight
    )) +
  geom_quasirandom(
    dodge.width = 0.9,
    varwidth = TRUE,
    size = 0.02,
    color = "grey40"
  ) +
  geom_boxplot(outlier.shape = NA,
               lwd = 0.4,
               fatten = 0.4) +
  # facet_grid(~Replicate, scales = "fixed", space =
  #                "free") +
  ylab(expression(atop(
    "", "Mismatch frequency [%]"
  ))) +
  scale_fill_manual(
    guide = "legend",
    values = c(alpha("red", 0.5), alpha("blue", 0.5)),
    breaks = c("ERCC", "Genes"),
    labels = c("ERCC spike-in controls", "Genes")
  ) +
  scale_color_manual(guide = "none",
                     values = c("black", "black")) +
  scale_x_discrete() +
  theme(
    title = element_text(
      size = 10,
      face = "plain",
      lineheight = .8
    ),
    panel.grid.major.x = element_line(colour = "white", size = 0.4),
    panel.grid.major.y = element_line(colour = "white", size = 0.4),
    legend.position = "top",
    strip.text = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16),
    text = element_text(size = 14),
    axis.text.x = element_text(
      size = 14,
      face = "plain",
      angle = 70,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 14, face = "plain"),
    legend.text = element_text(
      angle = 0,
      hjust = 0,
      size = 12
    ),
    legend.title.align = 0.5,
    legend.key.size = unit(20, "point"),
    legend.direction = 'horizontal',
    legend.key = element_blank(),
    legend.title = element_blank()
  )    +
  coord_cartesian(ylim = c(
    0,0.1))






#### tc probability from ERCC
####

coverage.bed.by.type <- coverage.bed %>% group_by(type, refbase, mismatch) %>% summarize(mmm_count = sum(mmm_count)) %>% ungroup() %>%
  mutate(coverage = coverage_totals_no_strand$coverage[match(paste(.$refbase,
                                                                   .$type),
                                                             paste(coverage_totals_no_strand$base,
                                                                   coverage_totals_no_strand$type))]
         ) %>%
  mutate(perc_mm = mmm_count/coverage*100)

a <- coverage.bed %>% group_by(chr, end, strand) %>% summarize(mmm_count = max(coverage) - min(coverage))


a <- coverage.bed %>% group_by(chr, end, refbase, strand, mismatch, gene_id) %>% summarize(mmm_count = sum(mmm_count),
                                                             coverage = sum(coverage)) %>%
  ungroup() %>%
  mutate(coverage = coverage_totals_no_strand$coverage[match(paste(.$refbase,
                                                                   .$type),
                                                             paste(coverage_totals_no_strand$base,
                                                                   coverage_totals_no_strand$type))]
  ) %>%
  mutate(perc_mm = mmm_count/coverage*100)



# coverage must by calculated using individual reads, not pairs, because if
# insert size is > 302 and mismatch falls between the reads,
# it's erroneosly included in the coverage

################# ggplot coverage ERCC and rest ##############
ggplot(data = a %>% filter(!type=="Spike-In")) +
       aes(x = type,
           y = perc_mm,
           fill = type) +

  geom_bar(size = 1, color = "black", stat = "identity") +
  # geom_boxplot(outlier.shape = NA,
  #              lwd = 0.4,
  #              fatten = 0.4) +
  # scale_y_log10() +

  ylab(expression(atop(
    "", "% of mismatched coverage"
  ))) +

  # labs(x = expression(paste("test", italic("test")))) +
  facet_wrap(~mismatch, scales = "fixed") +
  theme(
    title = element_text(
      size = 10,
      face = "plain",
      lineheight = .8
    ),
    panel.grid.major.x = element_line(colour = "white", size = 0.4),
    panel.grid.major.y = element_line(colour = "white", size = 0.4),
    legend.position = "top",
    strip.text = element_text(size = 12),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    text = element_text(size = 14),
    axis.text.x = element_text(
      size = 12,
      face = "plain",
      angle = 70,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 14, face = "plain"),
    legend.text = element_text(
      angle = 0,
      hjust = 0,
      size = 12
    ),
    legend.title.align = 0.5,
    legend.key.size = unit(20, "point"),
    legend.direction = 'horizontal',
    legend.key = element_blank(),
    legend.title = element_blank()
  )

+
  coord_cartesian(ylim = c(0, 1))


aaa <- coverage.bed.filtered %>%
  dplyr::group_by(gene_id, mismatch, coverage, refbase, end) %>%
  dplyr::summarise(mm_count = sum(mmm_count), rname_lst = list(do.call("c", rname_lst))) %>%
  dplyr::group_by(refbase, gene_id) %>%
  dplyr::mutate(cov_per_gene_per_refbase = sum(coverage))


table(aaa$mm_count)



coverage.bed.filtered  %>% rowwise() %>%  mutate(total_reads = length(unique(rname_lst)),
                                                 reads_table = list(c(table(rname_lst))))


# add n parameter based on tota number of covered bases

ggplot(data = coverage.bed.filtered %>%
            filter(gene_id %like% "ENSMUS") %>%
            pivot_longer(p.value:p.adj, names_to = "metric", values_to = "probability"),
       aes(x = mismatch,
           y = probability,
           fill = metric)) +

  # geom_point(size = 2, color = "black") +
  # geom_boxplot(outlier.shape = NA,
  #              lwd = 0.4,
  #              fatten = 0.4) +
  # scale_y_log10() +

  geom_violin(
    trim = TRUE,
    alpha = 0.7,
    size=0.1,
    draw_quantiles = c(0.25, 0.5, 0.75),
    scale = "width"
  ) +

  geom_quasirandom(
    dodge.width = 0.9,
    varwidth = TRUE,
    size = 0.2,
    color = "black"
    # aes(color = timepoint)
  ) +

  ylab(expression(atop(
    "Probability of the observed ratio under assumption", "that global mismatch frequency is correctly estimated"
  ))) +

  # labs(x = expression(paste("test", italic("test")))) +
  # facet_wrap(~mismatch, scales = "fixed", ncol = 11) +
  theme(
    title = element_text(
      size = 10,
      face = "plain",
      lineheight = .8
    ),
    panel.grid.major.x = element_line(colour = "white", size = 0.4),
    panel.grid.major.y = element_line(colour = "white", size = 0.4),
    legend.position = "top",
    strip.text = element_text(size = 12),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    text = element_text(size = 14),
    axis.text.x = element_text(
      size = 20,
      face = "plain",
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 14, face = "plain"),
    legend.text = element_text(
      angle = 0,
      hjust = 0,
      size = 12
    ),
    legend.title.align = 0.5,
    legend.key.size = unit(20, "point"),
    legend.direction = 'horizontal',
    legend.key = element_blank(),
    legend.title = element_blank()
  ) +
  coord_cartesian(ylim = c(0, 1))


setwd("../")
}


