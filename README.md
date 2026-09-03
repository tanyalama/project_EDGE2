# project_edge2 — EDGE2 rankings for the Upham mammalian supertree

EDGE2 (Evolutionarily Distinct and Globally Endangered) conservation-priority
rankings for all **5,911 mammal species** in the Upham et al. (2019) supertree,
computed under the **EDGE2 protocol** (Gumbs et al. 2023) with uncertainty
propagated across **1,000 posterior trees**.

## What's here

```
data/       EDGE2 ranked lists and reconciliation tables (CSV)
figures/    Publication-quality figures (PNG, 300 dpi)
R/          The EDGE2 engine and SLURM run scripts
docs/       Full methods documentation
```

### Key outputs (`data/`)

| File | Contents |
|---|---|
| `EDGE2_ranked_species_FULL.csv` | All 5,911 species — EDGE2 rank (`EDGErank`) and ED2 rank (`EDrank`), taxonomy, Red List category, median & IQR of EDGE2 / ED2 / pext, EDGE-species flag, match provenance |
| `EDGE_species_list.csv` | 586 threatened EDGE species (VU/EN/CR, ED above median in ≥50 % of trees) — the core priority list |
| `EDGE_borderline_list.csv` | 286 near-threshold species (flagged in 25–50 % of trees) — watch list |
| `EDGE_DD_watchlist.csv` | 938 Data-Deficient / Not-Evaluated species by ED2 — assessment priority |
| `reconciliation_report.csv` | Every tree tip → matched MDD name, match type, taxonomy, IUCN category |
| `taxonomic_summary_by_order.csv` | Per-order species counts, EDGE counts, % EDGE, median EDGE2 |
| `EDGE_species_missing_VGP_genome.csv` | 552 EDGE species lacking a VGP reference genome, ranked by EDGE2 — the genome-sequencing gap (see below) |

### Top 10 EDGE species

| Rank | Species | Order | RL | EDGE2 (Myr) |
|---|---|---|---|---|
| 1 | *Burramys parvus* (mountain pygmy possum) | Diprotodontia | CR | 24.4 |
| 2 | *Daubentonia madagascariensis* (aye-aye) | Primates | EN | 20.0 |
| 3 | *Gymnobelideus leadbeateri* (Leadbeater's possum) | Diprotodontia | CR | 19.6 |
| 4 | *Myrmecobius fasciatus* (numbat) | Dasyuromorphia | EN | 14.9 |
| 5 | *Manis culionensis* (Philippine pangolin) | Pholidota | CR | 14.8 |
| 6 | *Desmana moschata* (Russian desman) | Eulipotyphla | CR | 14.3 |
| 7 | *Manis pentadactyla* (Chinese pangolin) | Pholidota | CR | 14.2 |
| 8 | *Manis javanica* (Sunda pangolin) | Pholidota | CR | 14.2 |
| 9 | *Varecia variegata* (black-and-white ruffed lemur) | Primates | CR | 12.1 |
| 10 | *Varecia rubra* (red ruffed lemur) | Primates | CR | 12.0 |

## Figures (`figures/`)

- `fig1_edge2_rank_curve.png` — EDGE2 score vs rank, with cross-tree IQR band
- `fig2_ed_vs_ge2_scatter.png` — ED2 vs GE2 (pext), coloured by Red List category
- `fig3_top50_edge_species.png` — top-50 **EDGE species** (threatened VU/EN/CR, flagged), median ± IQR
- `fig4_ordinal_summary.png` — EDGE species counts and typical EDGE2 by order
- `fig5_top50_by_edge2_allcats.png` — top-50 by **EDGE2 score**, *all* Red List categories
- `fig6_top50_by_ED.png` — top-50 by **ED2** (raw evolutionary distinctness), *all* categories
- `fig7_vgp_genome_gap.png` — EDGE species genome gap: top-30 EDGE mammals lacking a VGP genome + per-order coverage

**Three top-50 views, three questions.** `fig3` answers "which threatened species
should we act on?" (the EDGE species list). `fig5` answers "which species carry the
most EDGE2 score regardless of threat status?" — this surfaces highly distinct
non-threatened lineage relicts such as *Dromiciops gliroides* (monito del monte,
sole living microbiotherian). `fig6` answers "which species are most
evolutionarily distinct?" purely on ED2, where such relicts dominate (33 of the
top 50 are Least Concern). A species can rank high on the full EDGE2 list (e.g.
*Dromiciops*, EDGErank 39, EDrank 2) yet be absent from `fig3` because it is not
in a threatened category — this is correct EDGE2 behaviour, not an omission.

## VGP genome-sequencing gap

Cross-referencing the 586 threatened EDGE species against the Vertebrate Genomes
Project (VGP) target list (with genome presence defined as a deposited GCA/GCF
assembly accession, and names reconciled through MDD synonymy so no assembly is
missed under a taxonomic split):

- **Only 34 of 586 EDGE mammals (5.8%) have a VGP reference genome.**
- The remaining **552 have no genome — and none of them are on the VGP target
  list at all** (unplanned, not merely unsequenced).
- Coverage is uniformly low across threat categories: **CR 7/107 (6.5%),
  EN 16/239 (6.7%), VU 11/240 (4.6%)**.
- The top 3 EDGE mammals overall — and 8 of the top 10 — have no genome. The 34
  that do are largely charismatic megafauna and great apes (pangolins, orangutans,
  rhinoceroses, koala, elephants, whales, *Gorilla*, *Pan*).
- By order, the gap is dominated by **Primates (132), Rodentia (106),
  Chiroptera (78), Artiodactyla (47), and Diprotodontia (41)**.

`data/EDGE_species_missing_VGP_genome.csv` lists all 552, ranked by EDGE2, with
order/family/Red List category, EDGE2 & ED2 medians, pext, and an
`on_vgp_target_list` flag. See `figures/fig7_vgp_genome_gap.png`.

## Method (summary)

For each of 1,000 posterior trees: each species' Red List category is mapped to a
sampled extinction probability (`pext`, Isaac et al. 2007 model); internal
branches are weighted by the product of descendant `pext`; the tip-to-root sum of
weighted branches is **EDGE2**, and **ED2 = EDGE2 / pext**. Results are summarised
per species as the **median and IQR** across trees. A species is an **EDGE
species** if it is threatened (VU/EN/CR/EW/EX) and its ED2 is above the tree
median in ≥50 % of trees. See [`docs/METHODS.md`](docs/METHODS.md) for full detail,
data versions, and taxonomy reconciliation.

## Data sources

- **Phylogeny:** Upham, Esselstyn & Jetz (2019), completed 5,911-species
  node-dated credible tree set (topoCons FBD), VertLife.
- **Taxonomy + IUCN categories:** Mammal Diversity Database (current release).
- **Algorithm:** rEDGE (Ramos-Gutiérrez & Gumbs) / `EDGE.2.calc` (Gumbs);
  vendored engine validated to machine precision against the reference.

> The endangerment layer uses the IUCN categories embedded in the current MDD
> release. Substituting a fresh IUCN Red List export is a drop-in replacement of
> the category table followed by re-running `R/aggregate.R`.

## Reproducing

The engine (`R/edge2_engine.R`) needs R with `ape`, `phylobase`, `data.table`,
`dplyr`. `R/run_chunk.R` is the per-tree SLURM-array driver (extinction model
`Isaac`, base seed 20240601); `R/aggregate.R` collapses per-tree results into the
ranked species table.

## References

- Gumbs R. et al. (2023) The EDGE2 protocol. *PLoS Biology* 21(2):e3001991.
- Upham N.S., Esselstyn J.A., Jetz W. (2019) Inferring the mammal tree. *PLoS Biology* 17(12):e3000494.
- Isaac N.J.B. et al. (2007) Mammals on the EDGE. *PLoS ONE* 2(3):e296.
- Mammal Diversity Database, American Society of Mammalogists, mammaldiversity.org.
