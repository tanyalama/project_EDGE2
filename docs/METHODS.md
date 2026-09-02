# EDGE2 Rankings for the Upham et al. Mammalian Supertree — Methods

**Deliverable:** Evolutionarily Distinct and Globally Endangered (EDGE2) rankings
for every species in the Upham et al. (2019) mammal supertree, computed under the
EDGE2 protocol (Gumbs et al. 2023) with uncertainty propagated across 1,000
posterior trees.

---

## 1. Method — the EDGE2 protocol

EDGE2 (Gumbs et al. 2023, *PLoS Biology* 21(2):e3001991,
doi:10.1371/journal.pbio.3001991) prioritises species by combining:

- **ED2 (Evolutionary Distinctness):** how much unique evolutionary history a
  species represents, computed on a phylogeny with branches weighted by the
  extinction probability of *other* species — so a species' ED2 rises when its
  close relatives are themselves at risk.
- **GE2 (Global Endangerment):** a per-species probability of extinction (`pext`)
  derived from its IUCN Red List category.
- **EDGE2 = expected PD loss attributable to a species**, i.e. the product of its
  distinctness and its extinction risk, expressed in millions of years (Myr).

For each tree, EDGE2 is computed as (following the protocol and the rEDGE /
`EDGE.2.calc` reference implementations):

1. Each Red List category is mapped to an extinction probability `pext`. We use
   the **Isaac et al. (2007)** model (the protocol/rEDGE default): CR = 0.97, and
   halving per step down — EN = 0.485, VU = 0.2425, NT = 0.12125, LC = 0.0606.
2. `pext` is **sampled**, not fixed: a monotone logit spline through the five
   category anchors defines a continuous `pext` curve, and each species draws a
   `pext` from the distribution for its category. Data Deficient (DD) and Not
   Evaluated (NE) species draw from the pooled distribution across all categories
   (pext < 0.999). This is the "GE2 distribution" and is what carries endangerment
   uncertainty into the result.
3. On the tree, every internal branch is weighted by the product of the `pext`
   values of all species descending from it. Summing weighted branch lengths from
   a tip to the root gives that tip's **EDGE2** score; **ED2 = EDGE2 / pext**.
4. A species is flagged an **EDGE species** in a tree if its ED2 is at or above the
   median ED2 of the tree **and** it is in a threatened category (VU, EN, CR, EW,
   or EX).

**Uncertainty** is captured on two axes simultaneously: (a) *phylogenetic* — the
calculation is repeated across 1,000 trees from the posterior distribution; and
(b) *endangerment* — the per-species `pext` is re-sampled each tree. Per-species
results are then summarised as the **median** and **inter-quartile range (IQR)**
across the 1,000 trees. A species is a final EDGE species if it was flagged in
**≥ 50 %** of trees (`isEDGEsp_frac ≥ 0.5`).

---

## 2. Data sources & versions

| Component | Source | Version / access |
|---|---|---|
| **Phylogeny** | Upham, Esselstyn & Jetz (2019), *PLoS Biology*, "Inferring the mammal tree" | VertLife `data.vertlife.org/mammaltree`, file `Completed_5911sp_topoCons_FBDasZhouEtAl.zip` (node-dated, topology-constrained, birth-death / FBD-as-Zhou-et-al credible tree set) |
| **Posterior sample** | 1,000 trees | Every 10th tree of the 10,000-tree posterior (tree0000 … tree9990), spanning the full posterior |
| **Taxonomy + IUCN categories** | Mammal Diversity Database (MDD) | Current master `mdd.csv` from `mammaldiversity.github.io/_data/mdd.csv` (6,904 rows), field `iucnStatus` |
| **EDGE2 algorithm** | rEDGE (Ramos-Gutiérrez & Gumbs) & `EDGE.2.calc` (Gumbs) | github.com/iramosgutierrez/rEDGE ; github.com/rgumbs/EDGE2 |

**Tree variant rationale:** the `topoCons` (topology-constrained) FBD completed
5,911-species set is the node-dated credible set that includes all described
species (DNA-based + imputed placements), matching MDD taxonomy at the tips.

> **Note on IUCN categories.** This run uses the IUCN Red List categories
> **embedded in the current MDD release** as the endangerment layer. MDD tracks
> current IUCN assessments and provides the exact taxonomic namespace of the tree
> tips, so it is a complete and internally consistent baseline. If a fresh IUCN
> Red List export is supplied, swapping it in is a drop-in replacement of the
> category table (`edge_table.csv`) followed by re-running the aggregation — the
> phylogenetic computation does not change.

---

## 3. Taxonomy reconciliation (tree tips → IUCN categories)

The tree has **5,987 tips**: **5,911 real species** (`Genus_species`) plus **76
fossil FBD backbone tips** (prefixed `X_`, e.g. `X_Shuotherium`) that are dropped
before EDGE2 computation.

The 5,911 real tips were matched to MDD (and thereby to an IUCN category) by:

1. **Direct match** on `sciName` (`Genus_species`): **5,037**
2. **Synonym via MSW3** name (`MSW3_sciName`): **469**
3. **Synonym via MDD nominal-name index** (33,620 historical synonyms): **232**
4. **Unmatched: 173** — predominantly Pleistocene / recently-extinct megafauna
   (e.g. *Mammuthus*, *Smilodon*, *Ursus spelaeus*, archaic *Homo*) and a residue
   of genus-reassignment synonyms.

**Result: 5,738 / 5,911 tips (97.1 %) reconciled to an IUCN category.** The 173
unmatched tips were assigned **DD** as a placeholder (flagged in
`reconciliation_report.csv`), so they still receive a sampled `pext` and appear in
the ranking; they are candidates for refinement with a curated IUCN export.

Per-tip category counts entering the run: LC 3299, DD 826 (incl. 173 placeholders),
VU 530, EN 499, NT 365, CR 206, NE 112, EX 73, EW 1.

---

## 4. Compute & parameters

- **Cluster:** Smith (SLURM). Conda R 4.3.3 environment (`ape`, `phylobase`,
  `data.table`, `dplyr`).
- **Engine:** a vendored, self-contained R implementation
  (`edge2_engine.R`) faithful to rEDGE. Its only change from the reference is an
  O(n) refactor of the per-species `pext` sampler that preserves the exact
  `set.seed()`+`sample()` draw sequence. **Validated against the rEDGE source on
  two full 5,911-species trees: max |ΔEDGE| = 0, max |ΔED| = 0, max |Δpext| = 0,
  identical EDGE-species flags** — i.e. machine-precision equivalence.
- **Parallelism:** 1,000 trees split into a 20-task SLURM array, 50 trees per
  task (~3.3 min each; whole run ~5 min wall).
- **Reproducibility:** extinction model `Isaac`; base seed **20240601**; per-tree
  seed = `20240601 + tree_index` (deterministic and reproducible).

---

## 5. Output files

| File | Contents |
|---|---|
| `EDGE2_ranked_species_FULL.csv` | All 5,911 species: rank, taxonomy, RL category, median & IQR of EDGE2 / ED2 / pext, TBL, EDGE-species flag & fraction, match provenance |
| `EDGE_species_list.csv` | 586 threatened EDGE species (VU/EN/CR, flagged in ≥50 % of trees) — the core priority list |
| `EDGE_borderline_list.csv` | 286 threatened species near the ED-median threshold (flagged in 25–50 % of trees) — a watch list |
| `EDGE_DD_watchlist.csv` | 938 Data-Deficient / Not-Evaluated species ranked by ED2 (339 with above-median ED2) — high-priority for assessment |
| `reconciliation_report.csv` | Every tree tip: matched MDD name, match type, taxonomy, IUCN category |
| `taxonomic_summary_by_order.csv` | Per-order species counts, EDGE-species counts, % EDGE, median EDGE2 |
| `fig1_edge2_rank_curve.png` | EDGE2 score vs rank, with cross-tree IQR band |
| `fig2_ed_vs_ge2_scatter.png` | ED2 vs GE2 (pext), coloured by Red List category |
| `fig3_top50_edge_species.png` | Top-50 EDGE species, median ± IQR |
| `fig4_ordinal_summary.png` | EDGE species counts & typical EDGE2 by order |

## 6. Column dictionary (ranked list)

- `EDGErank` — rank by descending median EDGE2 (1 = highest priority)
- `EDGEmed` / `EDGEiqr` — median and IQR of EDGE2 across 1,000 trees (Myr)
- `EDmed` / `EDiqr` — median and IQR of ED2 (Myr)
- `pextmed` / `pextiqr` — median and IQR of the sampled extinction probability
- `TBLmn` — mean terminal branch length (Myr)
- `isEDGEsp` — 1 if a final EDGE species (threatened & ED above median in ≥50 % of trees)
- `isEDGEsp_frac` — fraction of trees in which the species met the EDGE-species criterion
- `match` — how the tip was reconciled (`direct`, `synonym:MSW3`, `synonym:nominalName`, `unmatched`)

---

## 7. References

- Gumbs R. et al. (2023) The EDGE2 protocol. *PLoS Biology* 21(2):e3001991.
- Upham N.S., Esselstyn J.A., Jetz W. (2019) Inferring the mammal tree. *PLoS Biology* 17(12):e3000494.
- Isaac N.J.B. et al. (2007) Mammals on the EDGE. *PLoS ONE* 2(3):e296.
- Mammal Diversity Database (American Society of Mammalogists), mammaldiversity.org.
