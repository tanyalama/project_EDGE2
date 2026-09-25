#!/usr/bin/env Rscript
# aggregate.R
# Combines all chunk RDS files -> per-species EDGE2 summary (median across trees).
# Mirrors rEDGE calculate_EDGE2_multiphylo summarise step.

suppressWarnings(suppressMessages({library(data.table)}))
WS     <- Sys.getenv("EDGE2_WS")
RUNDIR <- file.path(WS, "run", Sys.getenv("CLADE"))
source(file.path(WS, "run", "edge2_engine.R"))  # for IQR2

files <- list.files(file.path(RUNDIR, "results"), pattern="chunk_.*\\.rds$", full.names=TRUE)
cat("aggregating", length(files), "chunk files\n")
DT <- rbindlist(lapply(files, readRDS))
n_trees_done <- uniqueN(DT$tree)
cat("total rows:", nrow(DT), "| trees:", n_trees_done, "| species:", uniqueN(DT$species), "\n")

# per-species summary across trees (medians + IQRs), rEDGE-faithful
summ <- DT[, .(
  TBLmn    = mean(TBL),
  pextmed  = median(pext),
  pextiqr  = IQR2(pext),
  EDmed    = median(ED),
  EDiqr    = IQR2(ED),
  EDGEmed  = median(EDGE),
  EDGEiqr  = IQR2(EDGE),
  isEDGEsp_frac = mean(isEDGEsp),
  n_trees  = .N
), by = species]

# RLcat is fixed per species (from the input table) -> take first
rlmap <- unique(DT[, .(species, RLcat)])
summ <- merge(summ, rlmap, by="species", all.x=TRUE)

# EDGE species: over-median ED in >50% of trees AND threatened category
summ[, isEDGEsp := as.integer(isEDGEsp_frac >= 0.5)]
setorder(summ, -EDGEmed)
summ[, EDGErank := seq_len(.N)]

fwrite(summ, file.path(RUNDIR, "EDGE2_ranked_species.csv"))
ef <- list.files(file.path(RUNDIR, "results"), pattern="epd_.*\\.rds$", full.names=TRUE)
if(length(ef)) fwrite(rbindlist(lapply(ef, readRDS)), file.path(RUNDIR, "ePD_per_tree.csv"))
cat("wrote EDGE2_ranked_species.csv :", nrow(summ), "species\n")
cat("EDGE species (flagged):", sum(summ$isEDGEsp), "\n")
print(summ[1:15, .(EDGErank, species, RLcat, EDGEmed, EDmed, pextmed, isEDGEsp)])
