# project_EDGE2: EDGE2 rankings × VGP genome coverage across vertebrate trees

This repository holds EDGE2 (Evolutionarily Distinct and Globally Endangered, Gumbs et al. 2023) rankings for five VertLife phylogenies. Each ranking is cross-referenced against the Vertebrate Genomes Project (VGP) Ordinal List to measure how many high-priority species already have a reference genome.

**v2 (September 2026).** Mammals were re-run, and birds, squamates (with tuatara), amphibians and chondrichthyans were added. All five clades use IUCN Red List 2026-1 and the updated VGP list. The v1 mammal-only release is kept in `archive/v1_mammals/`.

## Summary

| Clade | Tree | Species | EDGE spp. | EDGE with VGP genome | Expected PD loss (median) |
|---|---|---|---|---|---|
| Mammals | Upham et al. 2019 | 5,911 | 610 | 38 (6.2%) | 10.1% |
| Birds | Jetz et al. 2012 (Hackett) | 9,993 | 531 | 14 (2.6%) | 6.7% |
| Squamates | Tonini et al. 2016 | 9,755 | 875 | 3 (0.3%) | 10.9% |
| Amphibians | Jetz & Pyron 2018 | 7,238 | 1,314 | 3 (0.2%) | 15.5% |
| Chondrichthyans | Stein et al. 2018 | 1,192 | 263 | 10 (3.8%) | 15.1% |

VGP coverage is several-fold higher among the most distinct and highest-ranked species than across each clade as a whole (fig 11). Mammals: 4% of all species have a genome, versus 16% of the top-25 EDGE species. Chondrichthyans: 2.9% versus 20%.

Of the 1,063 VGP species, all 637 that fall in these five clades are placed on a tree tip. The other 426 belong to lineages without a VertLife tree: ray-finned fishes, turtles, crocodilians and other chordates. See `data/vgp_species_placement.csv`.

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
python/              table builder and figure module (Arial, Zissou1 palette)
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

### The EDGE2 protocol

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

### Data sources

- **Phylogenies (VertLife, 1,000 posterior trees each):**

  | Clade | Tree set | Trees used |
  |---|---|---|
  | Mammals | Upham, Esselstyn & Jetz (2019), completed 5,911-species node-dated set (topoCons FBD) | every 10th of 10,000 |
  | Birds | Jetz et al. (2012), Hackett backbone, Stage 2 full data | trees 1–1,000 |
  | Squamates | Tonini et al. (2016), 9,755 species, includes *Sphenodon* | first 1,000-tree file |
  | Amphibians | Jetz & Pyron (2018), 7,238 species | first 1,000-tree file |
  | Chondrichthyans | Stein et al. (2018), 1,192 species, 10-calibration set | every 10th of 10,000 |

- **Endangerment:** IUCN Red List v2026-1 global categories (GBIF-hosted Darwin Core archive), for every clade.
- **Taxonomy reconciliation:** tree tips are matched to IUCN names in this order: direct match; the Mammal Diversity Database (MDD) synonymy (mammals); the GBIF backbone accepted name; GBIF's IUCN link. Unmatched tips are coded NM.
- **Genomes:** the VGP Ordinal List (sheets "VGP Phase 1+" and "VGP Families"). A species has a genome if it has a GCA_/GCF_ assembly accession. Names are placed on tree tips by exact, synonym, orthographic or split-from-tip matching, with the evidence for each placement recorded in `data/vgp_species_placement.csv`.
- **Algorithm:** rEDGE (Ramos-Gutiérrez & Gumbs) / `EDGE.2.calc` (Gumbs). The vendored engine was validated to machine precision against the reference.

See [`docs/METHODS.md`](docs/METHODS.md) for full detail, data versions, reconciliation counts, VGP placement rules and caveats.

## Reproducing

The engine (`R/edge2_engine.R`) needs R with `ape`, `phylobase`, `data.table` and `dplyr`. `R/run_chunk.R` is the per-tree SLURM-array driver (extinction model `Isaac`, base seed 20240601), and `R/aggregate.R` collapses per-tree results into the ranked species table. Both read the clade from the `CLADE` environment variable. The steps are:

1. Download the tree sets from data.vertlife.org and the IUCN 2026-1 DwC-A (GBIF-hosted).
2. Split each tree set to one tree per file and write `run/<clade>/tree_index.txt` and `edge_table.csv`.
3. Run `sbatch --export=ALL,CLADE=<clade> R/chunk.sbatch`, followed by `R/agg.sbatch`.
4. Run `python/build_tables.py` and `python/edge2_figs.py`.

See `docs/METHODS.md` for details.

## References

Gumbs R. et al. (2023) *PLoS Biol* 21:e3001991 · Isaac N.J.B. et al. (2007) *PLoS ONE* 2:e296 · Gumbs R. et al. (2024) *Nat Commun* 15:1101 · Upham N.S. et al. (2019) *PLoS Biol* 17:e3000494 · Jetz W. et al. (2012) *Nature* 491:444 · Tonini J.F.R. et al. (2016) *Biol Conserv* 204:23 · Jetz W. & Pyron R.A. (2018) *Nat Ecol Evol* 2:850 · Stein R.W. et al. (2018) *Nat Ecol Evol* 2:288 · IUCN (2026) Red List v2026-1 · Mammal Diversity Database, mammaldiversity.org.
