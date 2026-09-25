# project_EDGE2: EDGE2 rankings × VGP genome coverage across vertebrate trees

This repository holds EDGE2 (Evolutionarily Distinct and Globally Endangered, Gumbs et al. 2023) rankings for five VertLife phylogenies. Each ranking is cross-referenced against the Vertebrate Genomes Project (VGP) Ordinal List to measure how many high-priority species already have a reference genome.

**v2 (September 2026).** Mammals were re-run, and birds, squamates (with tuatara), amphibians and chondrichthyans were added. All five clades use IUCN Red List 2026-1 and the updated VGP list. The v1 mammal-only release is kept in `archive/v1_mammals/`.

## Summary

| Clade | Phylogeny (VertLife tree set; 1,000 posterior trees) | Species | EDGE spp. | EDGE with VGP genome | Expected PD loss (median) |
|---|---|---|---|---|---|
| Mammals | Upham, Esselstyn & Jetz (2019) *PLoS Biol*. `mammaltree/Completed_5911sp_topoCons_FBDasZhouEtAl.zip` (node-dated, topology-constrained, FBD-as-Zhou-et-al); every 10th of 10,000 trees | 5,911 | 610 | 38 (6.2%) | 10.1% |
| Birds | Jetz et al. (2012) *Nature*. `birdtree/Stage2/HackettStage2_0001_1000.zip` (Hackett backbone, Stage 2 full data); trees 1–1,000 | 9,993 | 531 | 14 (2.6%) | 6.7% |
| Squamates | Tonini et al. (2016) *Biol Conserv*. `squamatetree/squam_shl_new_Posterior_9755.1000-10000.trees.zip` (includes *Sphenodon punctatus*); first 1,000-tree file of 10,000 | 9,755 | 875 | 3 (0.3%) | 10.9% |
| Amphibians | Jetz & Pyron (2018) *Nat Ecol Evol*. `amphibiantree/download/amph_shl_new_Posterior_7238.1000-10000.trees.zip`; first 1,000-tree file of 10,000 | 7,238 | 1,314 | 3 (0.2%) | 15.5% |
| Chondrichthyans | Stein et al. (2018) *Nat Ecol Evol*. `sharktree/Chond.10Cal.10kTreeSet.tre` (10-calibration set); every 10th of 10,000 trees | 1,192 | 263 | 10 (3.8%) | 15.1% |

VGP coverage is several-fold higher among the most distinct and highest-ranked species than across each clade as a whole (fig 11). Mammals: 4% of all species have a genome, versus 16% of the top-25 EDGE species. Chondrichthyans: 2.9% versus 20%.

Of the 1,063 VGP species, all 637 that fall in these five clades are placed on a tree tip. The other 426 belong to lineages without a VertLife tree: ray-finned fishes, turtles, crocodilians and other chordates. See `data/vgp_species_placement.csv`.

## Interactive dashboard

**[tanyalama.github.io/project_EDGE2](https://tanyalama.github.io/project_EDGE2/)**: an R Shiny app that runs entirely in the browser (Shinylive / webR), so it needs no server. The first load takes a few seconds while R starts in the browser. It has six tabs:

| Tab | Content |
|---|---|
| Overview | VGP coverage of EDGE species and the distinctness-enrichment plot, per group |
| Explore a group | Filter by Red List category, order or family, EDGE status, genome status or name. Interactive ED2-vs-GE2 scatter, top-species plot, rank curve, order/family breakdown, and a searchable, downloadable species table |
| Sequenced EDGE species | Every EDGE species that already has a VGP genome |
| Sequencing priorities | Highest-ranked EDGE species without a genome, filterable by order or family |
| VGP placement | Every VGP name with its placement tier, tree tip and evidence |
| About | Short description of the data and definitions |

To run it locally: `shiny::runApp("app")`. After changing `data/`, rebuild the app's data with `python python/build_app_data.py`. The site is rebuilt by `.github/workflows/deploy-app.yml` on every push that touches `app/`.

## Layout

```
data/<clade>/        EDGE2_ranked_species_FULL.csv, EDGE_species_list.csv, EDGE_borderline_list.csv,
                     EDGE_DD_watchlist.csv, EDGE_species_missing_VGP_genome.csv,
                     taxonomic_summary_by_{order|family}.csv, ePD_per_tree.csv
data/mammals/comparison_vs_previous_run.csv   per-species v1 → v2 changes
data/cross_clade_summary.csv                  one row per clade
data/vgp_species_placement.csv                every VGP name: placement tier, tree tip, evidence
data/vgp_distinctness_enrichment.csv          coverage by distinctness/rank tier
data/reconciliation_all_clades.csv            tree tip → IUCN 2026-1 category and match route
figures/<clade>/     fig1–fig7 (PDF + 300-dpi PNG)
figures/cross_clade/ fig8–fig11
docs/METHODS.md      full methods, v1 → v2 changes, caveats
R/                   EDGE2 engine (vendored rEDGE), SLURM chunk/aggregate scripts
python/              table builder, figure module (Arial, Zissou1 palette), dashboard data builder
app/                 Shiny dashboard (app.R + data/), published to GitHub Pages via shinylive
archive/v1_mammals/  previous release
```

## Figures

Each clade folder contains:

| Figure | Content |
|---|---|
| fig1 | Rank curve |
| fig2 | ED2 vs GE2 |
| fig3 | Top-50 EDGE species |
| fig4 | EDGE species by order or family |
| fig5 | Top 50 by EDGE2, all categories |
| fig6 | Top 50 by ED2 |
| fig7 | Top EDGE species without a genome, and coverage by order or family |

The cross-clade folder contains:

| Figure | Content |
|---|---|
| fig8 | Threatened PD and VGP coverage |
| fig9 | Red List composition of VGP genomes |
| fig10 | Every EDGE species with a VGP genome |
| fig11 | VGP coverage vs distinctness and rank |

All figures use Arial and the Wes Anderson Zissou1 palette. `RLcat = NM` marks tips with no Red List match (a project code, not an IUCN category). Their GE2 is imputed exactly as for DD.

## Methods

### 1. Method — the EDGE2 protocol

EDGE2 (Gumbs et al. 2023, *PLoS Biology* 21(2):e3001991, doi:10.1371/journal.pbio.3001991) prioritises species by combining three quantities:

- **ED2 (Evolutionary Distinctness):** how much unique evolutionary history a species represents. It is computed on a phylogeny whose branches are weighted by the extinction probability of *other* species, so a species' ED2 rises when its close relatives are themselves at risk.
- **GE2 (Global Endangerment):** a per-species probability of extinction (`pext`) derived from its IUCN Red List category.
- **EDGE2:** the expected loss of phylogenetic diversity (PD) attributable to a species, i.e. the product of its distinctness and its extinction risk, in millions of years (Myr).

For each tree, EDGE2 is computed as follows, following the protocol and the rEDGE / `EDGE.2.calc` reference implementations:

1. Each Red List category is mapped to an extinction probability `pext` using the **Isaac et al. (2007)** model, which is the protocol and rEDGE default. CR = 0.97, halving at each step down: EN = 0.485, VU = 0.2425, NT = 0.12125, LC = 0.0606.
2. `pext` is **sampled**, not fixed. A monotone logit spline through the five category anchors defines a continuous `pext` curve, and each species draws a `pext` from the distribution for its category.
   - Data Deficient (DD) and Not Evaluated (NE) species draw from the pooled distribution across all categories (pext < 0.999).
   - Tips with no Red List match (NM, a project code rather than an IUCN category) are treated the same way.
   - This sampled distribution is the "GE2 distribution", and it carries endangerment uncertainty into the result.
3. On the tree, every internal branch is weighted by the product of the `pext` values of all species descending from it. The sum of weighted branch lengths from a tip to the root is that tip's **EDGE2** score, and **ED2 = EDGE2 / pext**.
4. A species is flagged as an **EDGE species** in a tree if two conditions hold:
   - its ED2 is at or above the median ED2 of the tree;
   - it is in a threatened category (VU, EN, CR, EW or EX).

**Uncertainty** is captured on two axes at once:

- *Phylogenetic:* the calculation is repeated across 1,000 trees from the posterior distribution of each clade.
- *Endangerment:* each species' `pext` is re-sampled for every tree.

Per-species results are summarised as the **median** and **inter-quartile range (IQR)** across the 1,000 trees. A species is a final EDGE species if it was flagged in **≥ 50 %** of trees (`isEDGEsp_frac ≥ 0.5`). The EDGE species lists, counts and VGP coverage figures in this repository report extant threatened species (VU, EN, CR); extinct (EX) and extinct-in-the-wild (EW) species stay in the trees and in the full ranked tables.

The same procedure is applied to all five clades. Per-tree total PD and expected PD loss (the sum of EDGE2 across species) are also saved (`ePD_per_tree.csv`).

### 2. Data sources & versions

| Component | Source | Version / access |
|---|---|---|
| **Phylogenies** | VertLife (data.vertlife.org) | Tree set, reference and posterior sample for each clade are given in the [Summary](#summary) table |
| **IUCN categories** | IUCN Red List of Threatened Species | v2026-1, GBIF-hosted Darwin Core archive (`hosted-datasets.gbif.org/datasets/iucn/iucn-latest.zip`); 65,076 chordate species with a global category; accessed 25 Sep 2026 |
| **Synonymy** | Mammal Diversity Database (MDD); GBIF backbone taxonomy | Current MDD master `mdd.csv` (`mammaldiversity.github.io/_data/mdd.csv`); GBIF species match API (strict, kingdom Animalia) |
| **Genome list** | Vertebrate Genomes Project Ordinal List | Sheets "VGP Phase 1+" (731 species) and "VGP Families" (332 more): 1,063 binomials, 1,022 with a GCA_/GCF_ accession |
| **EDGE2 algorithm** | rEDGE (Ramos-Gutiérrez & Gumbs) & `EDGE.2.calc` (Gumbs) | github.com/iramosgutierrez/rEDGE ; github.com/rgumbs/EDGE2 |

**Tree variant rationale:** each clade uses the VertLife completed (DNA-based plus imputed placements) posterior set, so every described species in the source taxonomy has a tip. For mammals this is the `topoCons` FBD 5,911-species set used in v1, which keeps the mammal results directly comparable with v1.

> **Note on IUCN categories.** v1 used the IUCN categories embedded in the MDD release. v2 uses the IUCN Red List v2026-1 export directly for every clade, so all five clades share one current endangerment layer. MDD is still used for mammal synonymy. The change of layer, together with new assessments since v1, accounts for most of the difference between the v1 and v2 mammal EDGE lists (`data/mammals/comparison_vs_previous_run.csv`).

### 3. Taxonomy reconciliation (tree tips → IUCN categories)

Non-species tips are dropped before the EDGE2 computation:
- Mammals: 5,987 tips = **5,911 species** plus **76 fossil FBD backbone tips** (prefixed `X_`, e.g. `X_Shuotherium`).
- Amphibians: 7,239 tips = **7,238 species** plus the *Homo sapiens* outgroup.
- Birds (9,993), squamates (9,755) and chondrichthyans (1,192): all tips are species.

Each species tip was matched to an IUCN 2026-1 assessment, taking the first route that succeeded:

1. **Direct match** of `Genus_species` to an IUCN accepted name.
2. **Synonym via MDD** (mammals only): the MDD name the tip reconciled to in v1 (MSW3 and nominal-name synonymy), then matched to IUCN.
3. **Synonym via GBIF backbone:** the GBIF accepted species for the tip name, then matched to IUCN.
4. **GBIF IUCN link:** the IUCN category GBIF attaches to that species.
5. **Unmatched → NM** ("no Red List match", flagged in `data/reconciliation_all_clades.csv`). v1 labelled these tips DD. NM tips still receive a sampled `pext` from the pooled distribution, exactly as DD and NE do, so they appear in the ranking. They are candidates for refinement with a curated synonym list.

| Clade | Direct | MDD synonym | GBIF accepted | GBIF IUCN link | NM | Reconciled |
|---|---|---|---|---|---|---|
| Mammals | 5,316 | 336 | 77 | 53 | 129 | 5,782 / 5,911 (97.8 %) |
| Birds | 7,931 | – | 1,105 | 763 | 194 | 9,799 / 9,993 (98.1 %) |
| Squamates | 8,636 | – | 612 | 268 | 239 | 9,516 / 9,755 (97.5 %) |
| Amphibians | 6,069 | – | 879 | 237 | 53 | 7,185 / 7,238 (99.3 %) |
| Chondrichthyans | 979 | – | 101 | 74 | 38 | 1,154 / 1,192 (96.8 %) |

Per-tip category counts entering the run ("conservation dependent" mapped to NT):

| Clade | LC | NT | VU | EN | CR | EW | EX | DD | NM |
|---|---|---|---|---|---|---|---|---|---|
| Mammals | 3,340 | 382 | 546 | 524 | 214 | 1 | 81 | 694 | 129 |
| Birds | 7,855 | 808 | 591 | 322 | 178 | 5 | 11 | 29 | 194 |
| Squamates | 6,126 | 508 | 531 | 687 | 309 | 2 | 21 | 1,332 | 239 |
| Amphibians | 3,562 | 397 | 714 | 1,091 | 676 | 2 | 37 | 706 | 53 |
| Chondrichthyans | 520 | 120 | 186 | 117 | 98 | 0 | 1 | 112 | 38 |

VGP names are placed on tree tips by a separate tiered procedure (exact, synonym, orthographic, split-from-tip, domestic form). Every placement and its evidence is recorded in `data/vgp_species_placement.csv`; see [`docs/METHODS.md`](docs/METHODS.md).

### 4. Compute & parameters

- **Cluster:** Smith (SLURM). Conda R 4.3.3 environment (`ape`, `phylobase`, `data.table`, `dplyr`).
- **Engine:** a vendored, self-contained R implementation (`R/edge2_engine.R`), unchanged since v1 and faithful to rEDGE. Its only change from the reference is an O(n) refactor of the per-species `pext` sampler that preserves the exact `set.seed()` + `sample()` draw sequence. In v1 it was **validated against the rEDGE source on two full 5,911-species trees: max |ΔEDGE| = 0, max |ΔED| = 0, max |Δpext| = 0, identical EDGE-species flags**, i.e. machine-precision equivalence.
- **Parallelism:** per clade, 1,000 trees are split into a 20-task SLURM array with 50 trees per task, followed by a dependent aggregation job. The five arrays ran concurrently. The slowest task per clade took 3.1 min (mammals), 5.9 min (birds), 5.0 min (squamates), 3.1 min (amphibians) and 0.8 min (chondrichthyans). The whole run took about 6 min wall time.
- **Reproducibility:** extinction model `Isaac`; base seed **20240601**; per-tree seed = `20240601 + tree_index`, which is deterministic and reproducible. The same seeds are used for every clade.

See [`docs/METHODS.md`](docs/METHODS.md) for the VGP placement rules, v1 → v2 changes and caveats.

## Reproducing

`R/run_chunk.R` is the per-tree SLURM-array driver and `R/aggregate.R` collapses per-tree results into the ranked species table; both read the clade from the `CLADE` environment variable. The steps are:

1. Download the tree sets from data.vertlife.org and the IUCN 2026-1 DwC-A (GBIF-hosted).
2. Split each tree set to one tree per file and write `run/<clade>/tree_index.txt` and `edge_table.csv`.
3. Run `sbatch --export=ALL,CLADE=<clade> R/chunk.sbatch`, followed by `R/agg.sbatch`.
4. Run `python/build_tables.py` and `python/edge2_figs.py`.

See `docs/METHODS.md` for details.

## References

Gumbs R. et al. (2023) *PLoS Biol* 21:e3001991 · Isaac N.J.B. et al. (2007) *PLoS ONE* 2:e296 · Gumbs R. et al. (2024) *Nat Commun* 15:1101 · Upham N.S. et al. (2019) *PLoS Biol* 17:e3000494 · Jetz W. et al. (2012) *Nature* 491:444 · Tonini J.F.R. et al. (2016) *Biol Conserv* 204:23 · Jetz W. & Pyron R.A. (2018) *Nat Ecol Evol* 2:850 · Stein R.W. et al. (2018) *Nat Ecol Evol* 2:288 · IUCN (2026) Red List v2026-1 · Mammal Diversity Database, mammaldiversity.org.
