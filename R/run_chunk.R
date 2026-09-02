#!/usr/bin/env Rscript
# run_chunk.R  <chunk_id>
# Processes one chunk of posterior trees for the EDGE2 run.
# Reads: tree_index.txt (1 tree-file path per line, the 1000 selected trees)
#        edge_table.csv  (species, RLcat) reconciled to tips
# Env:   EDGE2_WS (scratch), CHUNK_SIZE, EXT_PROB, BASE_SEED
# Writes: results/chunk_<id>.rds  (data.table of per-tree tip results)

suppressWarnings(suppressMessages({library(data.table)}))
args <- commandArgs(TRUE)
chunk_id <- as.integer(args[1])

WS        <- Sys.getenv("EDGE2_WS")
CHUNK     <- as.integer(Sys.getenv("CHUNK_SIZE", "50"))
EXT_PROB  <- Sys.getenv("EXT_PROB", "Isaac")
BASE_SEED <- as.integer(Sys.getenv("BASE_SEED", "20240601"))
RUNDIR    <- file.path(WS, "run")
source(file.path(RUNDIR, "edge2_engine.R"))

tree_paths <- readLines(file.path(RUNDIR, "tree_index.txt"))
tab <- read.csv(file.path(RUNDIR, "edge_table.csv"), stringsAsFactors = FALSE)
stopifnot(all(c("species","RLcat") %in% names(tab)))

n_trees <- length(tree_paths)
i0 <- (chunk_id - 1L)*CHUNK + 1L
i1 <- min(chunk_id*CHUNK, n_trees)
if(i0 > n_trees){ cat("chunk", chunk_id, "empty\n"); quit(status=0) }
cat(sprintf("chunk %d : trees %d..%d of %d | ext.prob=%s\n", chunk_id, i0, i1, n_trees, EXT_PROB))

out <- vector("list", i1 - i0 + 1L)
k <- 0L
for(ti in i0:i1){
  k <- k + 1L
  tr <- ape::read.tree(tree_paths[ti])
  # drop fossil FBD backbone tips (X_ prefix) not in the table
  drop <- setdiff(tr$tip.label, tab$species)
  if(length(drop)) tr <- ape::drop.tip(tr, drop)
  tabi <- tab[tab$species %in% tr$tip.label, c("species","RLcat")]
  # per-tree seed: deterministic, reproducible
  seed_i <- BASE_SEED + ti
  res <- calculate_EDGE2(tr, tabi, ext.prob = EXT_PROB, seed = seed_i)
  dt <- as.data.table(res[, c("species","RLcat","TBL","pext","ED","EDGE","isEDGEsp")])
  dt[, tree := ti]
  out[[k]] <- dt
  if(k %% 5 == 0) cat("  done", k, "trees\n")
}
res_all <- rbindlist(out)
dir.create(file.path(RUNDIR, "results"), showWarnings = FALSE)
saveRDS(res_all, file.path(RUNDIR, "results", sprintf("chunk_%03d.rds", chunk_id)),
        compress = "xz")
cat("chunk", chunk_id, "wrote", nrow(res_all), "rows\n")
