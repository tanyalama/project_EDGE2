# edge2_engine.R
# Vendored EDGE2 engine, faithful to rEDGE (Ramos-Gutierrez & Gumbs) and the
# EDGE2 protocol (Gumbs et al. 2023, PLoS Biology 21(2):e3001991).
# Functions copied verbatim from github.com/iramosgutierrez/rEDGE except
# get_extinction_prob, which is reorganised to O(n) while preserving the
# EXACT per-row set.seed()+sample() draw sequence (identical numerical output).

suppressWarnings(suppressMessages({
  library(ape); library(phylobase); library(dplyr)
}))

# --- extinction-probability anchors per Red List category ------------------
cat_pext <- function(ext.prob = "Isaac"){
  if(!ext.prob %in% c("Isaac","IUCN50","IUCN100","IUCN500"))
    stop("ext.prob should be one of: 'Isaac'/'IUCN50'/'IUCN100'/'IUCN500'")
  if(ext.prob == "Isaac")
    cat.pext <- data.frame(rl.cat = rev(c("CR","EN","VU","NT","LC")),
                           pext   = rev(c(0.97, 0.97/2, 0.97/4, 0.97/8, 0.97/16)))
  if(ext.prob == "IUCN50")
    cat.pext <- data.frame(rl.cat = rev(c("CR","EN","VU","NT","LC")),
                           pext   = rev(c(0.97, 0.42, 0.05, 0.004, 0.00005)))
  if(ext.prob == "IUCN100")
    cat.pext <- data.frame(rl.cat = rev(c("CR","EN","VU","NT","LC")),
                           pext   = rev(c(0.999, 0.667, 0.1, 0.01, 0.0001)))
  if(ext.prob == "IUCN500")
    cat.pext <- data.frame(rl.cat = rev(c("CR","EN","VU","NT","LC")),
                           pext   = rev(c(1, 0.996, 0.39, 0.02, 0.0005)))
  cat.pext
}

# monotone logit spline through the 5 category anchors (rEDGE make_pext_curve)
make_pext_curve <- function(ext.prob){
  x_pts <- 1:5
  y_pts <- cat_pext(ext.prob)$pext
  eps <- 1e-12
  y_pts <- pmin(pmax(y_pts, eps), 1 - eps)
  logit_y <- log(y_pts/(1 - y_pts))
  phi <- splinefun(x_pts, logit_y, method = "monoH.FC")
  function(x) plogis(phi(x))
}

# GE2/pext sampling distribution per RL category (rEDGE create_pext_by_cat)
create_pext_by_cat <- function(n = 1000000, ext.prob){
  pext.dist <- data.frame(RL.num = seq(0, 6, length.out = n))
  pext.dist$pext <- make_pext_curve(ext.prob)(pext.dist$RL.num)
  d <- pext.dist
  rbind(
    data.frame(RLcat="LC", pext=d$pext[d$RL.num>=0.0 & d$RL.num<1.5]),
    data.frame(RLcat="NT", pext=d$pext[d$RL.num>=1.5 & d$RL.num<2.5]),
    data.frame(RLcat="CD", pext=d$pext[d$RL.num>=1.5 & d$RL.num<2.5]),
    data.frame(RLcat="VU", pext=d$pext[d$RL.num>=2.5 & d$RL.num<3.5]),
    data.frame(RLcat="EN", pext=d$pext[d$RL.num>=3.5 & d$RL.num<4.5]),
    data.frame(RLcat="CR", pext=d$pext[d$RL.num>=4.5 & d$RL.num<5.5]),
    data.frame(RLcat="EW", pext=d$pext[d$RL.num>=5.5 & d$RL.num<6.0])
  )
}

# O(n) reorganisation of rEDGE get_extinction_prob. cat_pext_table is
# deterministic (no RNG), so it is precomputed. Per-row draws use the same
# set.seed(seed+i) / set.seed(seed+2*i) + sample(pool,1) as the original loop,
# so numerical output is identical; only the O(n^2) which() lookup is removed.
get_extinction_prob_fast <- function(tab, ext.prob, seed){
  ct <- create_pext_by_cat(ext.prob = ext.prob)
  pool_lt999 <- ct$pext[ct$pext < 0.999]
  cats <- unique(ct$RLcat)
  pools <- lapply(cats, function(cc) ct$pext[ct$RLcat == cc]); names(pools) <- cats
  n  <- nrow(tab); rl <- as.character(tab$RLcat); pe <- numeric(n)
  for(i in seq_len(n)){
    ci <- rl[i]
    if(ci %in% c("NE","DD")){
      set.seed(seed + i); pe[i] <- sample(pool_lt999, size = 1)
    } else {
      if(ci %in% c("EX","EW")) ci <- "CR"
      if(ci == "CD")           ci <- "NT"
      set.seed(seed + 2*i);   pe[i] <- sample(pools[[ci]], size = 1)
    }
  }
  tab$pext <- pe; tab
}

into_order  <- function(tree, pext) pext[match(tree$tip.label, pext$species), ]
reorder_tree <- function(tree, ordering){
  tree@edge.length <- tree@edge.length[ordering]
  tree@edge        <- tree@edge[ordering, ]
  tree
}
IQR2 <- function(x, na.rm = TRUE)
  stats::quantile(x, 0.75, na.rm = na.rm) - stats::quantile(x, 0.25, na.rm = na.rm)

# --- core EDGE2 for one tree (rEDGE calculate_EDGE2, verbatim logic) --------
calculate_EDGE2 <- function(tree, table, ext.prob = "Isaac", seed = NULL){
  table <- table[, c("species","RLcat")]
  if(is.null(seed)) seed <- round(runif(1,1,999999999))
  table <- get_extinction_prob_fast(table, ext.prob = ext.prob, seed = seed)
  names(table) <- c("species","RLcat","pext")

  N_species <- length(tree$tip.label); N_nodes <- tree$Nnode
  N_tot <- N_species + N_nodes
  if(!identical(tree$tip.label, table$species)) table <- into_order(tree, table)
  if(!inherits(tree,"phylo")) tree <- as(tree,"phylo")

  tree_dat <- data.frame(species = as.character(tree$tip.label),
                         TBL = NA, pext = table$pext, ED = NA, EDGE = NA)
  tree_dat <- merge(table[,c("species","RLcat")], tree_dat, sort = FALSE)
  ePD <- data.frame(PD = sum(tree$edge.length), ePDloss = NA)

  tree <- as(tree,"phylo4")
  root  <- phylobase::rootNode(tree)
  nodes <- c(root, phylobase::descendants(tree, root, "all"))
  ord <- order(nodes); tree <- reorder_tree(tree, ord); nodes <- nodes[ord]
  tree_dat$TBL <- tree@edge.length[1:N_species]

  node_data <- data.frame(Node = 1:N_tot, Pext = rep(1, N_tot), Edge_Sum = NA)
  node_data[1:N_species, "Pext"] <- table[,"pext"]
  for(i in c(1:length(tree@label), N_tot:(root+1))){
    anc <- tree@edge[i,1]
    node_data[anc,"Pext"] <- node_data[anc,"Pext"] * node_data[i,"Pext"]
  }
  for(i in 1:length(nodes)) tree@edge.length[i] <- tree@edge.length[i]*node_data[i,2]
  if(is.na(tree@edge.length[root])) tree@edge.length[root] <- 0
  node_data$Edge_Sum[root] <- tree@edge.length[root]
  for(i in (root+1):N_tot){
    ans <- tree@edge[i,1]
    node_data$Edge_Sum[i] <- node_data$Edge_Sum[ans] + tree@edge.length[i]
  }
  for(i in 1:N_species){
    ans <- tree@edge[i,1]
    tree_dat$EDGE[i] <- node_data$Edge_Sum[ans] + tree@edge.length[i]
  }
  tree_dat$ED <- tree_dat$EDGE / tree_dat$pext
  EDmed <- median(tree_dat$ED)
  tree_dat$isEDGEsp <- 0
  tree_dat$isEDGEsp[tree_dat$ED >= EDmed &
                    tree_dat$RLcat %in% c("VU","EN","CR","EW","EX")] <- 1
  tree <- as(tree,"phylo"); ePD$ePDloss <- sum(tree$edge.length)
  attr(tree_dat, "ePD") <- ePD
  attr(tree_dat, "seed") <- seed
  tree_dat
}
