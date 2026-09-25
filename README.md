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

## Reproducing

1. Download the tree sets from data.vertlife.org and the IUCN 2026-1 DwC-A (GBIF-hosted).
2. Split each tree set to one tree per file and write `run/<clade>/tree_index.txt` and `edge_table.csv`.
3. Run `sbatch --export=ALL,CLADE=<clade> R/chunk.sbatch`, followed by `R/agg.sbatch`.
4. Run `python/build_tables.py` and `python/edge2_figs.py`.

See `docs/METHODS.md` for details.

## References

Gumbs R. et al. (2023) *PLoS Biol* 21:e3001991 · Gumbs R. et al. (2024) *Nat Commun* 15:1101 · Upham N.S. et al. (2019) *PLoS Biol* 17:e3000494 · Jetz W. et al. (2012) *Nature* 491:444 · Tonini J.F.R. et al. (2016) *Biol Conserv* 204:23 · Jetz W. & Pyron R.A. (2018) *Nat Ecol Evol* 2:850 · Stein R.W. et al. (2018) *Nat Ecol Evol* 2:288 · IUCN (2026) Red List v2026-1.
