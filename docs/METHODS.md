# EDGE2 × VGP, v2 (September 2026): mammals re-run and four VertLife clades added

## 1. What changed from v1

| | v1 (mammals only) | v2 |
|---|---|---|
| Clades | Mammals | Mammals, birds, squamates (with *Sphenodon*), amphibians, chondrichthyans |
| Red List layer | IUCN categories embedded in MDD | IUCN Red List v2026-1 (GBIF-hosted DwC-A) for every clade |
| VGP list | earlier VGP Ordinal List | updated *VGP Ordinal List.xlsx* (sheets "VGP Phase 1+" and "VGP Families") |
| Engine | vendored rEDGE (`R/edge2_engine.R`) | unchanged: Isaac pext, base seed 20240601, 1,000 trees, 20 × 50-tree SLURM chunks; `CLADE` env var selects `run/<clade>/` |
| Extra outputs | none | per-tree total PD and expected PD loss (`ePD_per_tree.csv`); IUCN English common names |

## 2. Phylogenies (VertLife, data.vertlife.org)

| Clade | Tree set | Trees used | Tips |
|---|---|---|---|
| Mammals | Upham et al. 2019, `Completed_5911sp_topoCons_FBDasZhouEtAl` | every 10th of 10,000 (1,000) | 5,911 |
| Birds | Jetz et al. 2012, Hackett Stage 2 full-data, `HackettStage2_0001_1000` | 1,000 | 9,993 |
| Squamates | Tonini et al. 2016, `squam_shl_new_Posterior_9755` (first 1,000-tree file) | 1,000 | 9,755 (includes *Sphenodon punctatus*) |
| Amphibians | Jetz & Pyron 2018, `amph_shl_new_Posterior_7238` (first 1,000-tree file) | 1,000 | 7,238 (*Homo sapiens* outgroup removed) |
| Chondrichthyans | Stein et al. 2018, `Chond.10Cal.10kTreeSet` | every 10th of 10,000 (1,000) | 1,192 |

VertLife has no trees for ray-finned fishes, turtles or crocodilians, so those VGP lineages are not scored here.

## 3. Red List reconciliation

Each tip was assigned an IUCN 2026-1 global category using the first rule that matched:

1. Direct binomial match to an IUCN accepted name.
2. Mammals only: the MDD name from the v1 reconciliation (`reconciliation_report.csv`), then matched to IUCN.
3. The GBIF backbone accepted species for the tip name (strict match, kingdom Animalia), then matched to IUCN.
4. GBIF's IUCN Red List category link for that GBIF species.
5. Mammals only: fall back to MDD `iucnStatus`.
6. Otherwise, **NM (no Red List match)**. This is a project code, not an IUCN category: it marks tips that could not be matched to any IUCN 2026-1 assessment. v1 labelled these tips DD. The engine treats NM exactly as it treats DD and NE: GE2 is drawn at random from the full pext distribution. Relabelling therefore leaves every EDGE2 score unchanged, and genuine IUCN Data Deficient species keep the DD label.

The IUCN archive contains accepted names only and no synonyms, so rules 2 to 4 do the synonym resolution.

| Clade | direct | MDD synonym | GBIF accepted | GBIF IUCN link | NM (no match) |
|---|---|---|---|---|---|
| Mammals | 5,316 | 336 | 77 | 53 | 129 |
| Birds | 7,931 | – | 1,105 | 763 | 194 |
| Squamates | 8,636 | – | 612 | 268 | 239 |
| Amphibians | 6,069 | – | 879 | 237 | 53 |
| Chondrichthyans | 979 | – | 101 | 74 | 38 |

"Conservation dependent" was mapped to NT. EX and EW species are retained in the trees, following v1 and rEDGE. In each tree, a species is flagged if its ED2 is at or above that tree's median ED2 and it is threatened (VU, EN, CR, EW or EX). An **EDGE species** is flagged in at least 50% of trees (`isEDGEsp_frac ≥ 0.5`). Lists and coverage figures report extant threatened EDGE species (VU, EN, CR). **Borderline** species are threatened with `isEDGEsp_frac` between 0.25 and 0.5.

Order and family labels come from each clade's VertLife taxonomy file, with IUCN genus-level lookup as the fallback. Squamates are summarised by family.

## 4. VGP cross-reference

The two sheets were combined: "VGP Phase 1+" (731 species) and "VGP Families" (332 additional species), giving 1,063 binomials, of which 1,022 have a GCA_/GCF_ accession. As in v1, a species **has a genome** if its main-haplotype field contains a GCA_/GCF_ accession.

**Scope.** 426 of the 1,063 VGP names belong to lineages that have no VertLife tree, so they are out of scope rather than failed matches. These are ray-finned fishes, turtles, crocodilians, lampreys, hagfish, the coelacanth, lungfish, and non-vertebrate chordates and invertebrates. The remaining 637 names fall in the five tree-covered clades, and **all 637 are placed on a tree tip**.

**Placement tiers.** Each name is assigned the first tier it qualifies for. The per-name tier and evidence are recorded in `vgp_species_placement.csv`.

| Tier | n | Rule | Genome credited to tip |
|---|---|---|---|
| exact | 577 | VGP name identical to a tip label | yes |
| synonym_of_tip | 39 | Same taxon under another name, from GBIF backbone accepted name, MDD synonymy or IUCN accepted name. Examples: generic reassignments such as *Poecile* → *Parus* and *Mobula* → *Manta*; gender endings such as *Strigops habroptilus/-a*; VGP sheet misspellings (*Monodon monocero*, *Lophostoma evote*, *Cephalophula zebra*) and a row with swapped name columns ("Eurasian siskin" → *Spinus spinus* → *Carduelis spinus*) | yes |
| orthographic | 1 | Spelling variant of a tip: *Guaruba guaruba* / *Guaruba guarouba* | yes |
| split_from_tip | 18 | Species split after the tree was built. Evidence is an MDD "split from" or "now attributed to" note (mammals), a GBIF subspecies record under the tip in the same genus or a documented renamed genus, or a single congener in the tree. Examples: *Giraffa tippelskirchi* → *G. camelopardalis*; *Balaenoptera ricei* → *B. edeni*; *Chlamydotis macqueenii* → *C. undulata*; *Natrix helvetica* → *N. natrix* | yes, flagged |
| domestic_form_of_tip | 2 | MDD "domestic form of": *Bubalus bubalis* → *B. arnee*; *Equus asinus* → *E. africanus* | no |
| no_VertLife_tree_for_lineage | 426 | Lineage has no VertLife tree | – |

Synonym credit is given only when the VGP name has no direct tip in that tree. This prevents lumps such as *Bos taurus* → *B. indicus*. Subspecies evidence is accepted only within the same genus, because shared epithets such as *vidua* and *nanus* also occur in unrelated genera.

Five tips receive more than one VGP name:

- *Monodon monoceros*, *Lophostoma evotis* and *Cephalophus zebra*: the same taxon listed twice, once misspelled.
- *Carduelis flammea* (*Acanthis flammea* + *A. cabaret*) and *Artibeus lituratus* (+ *A. intermedius*): two taxa that the tree taxonomy lumps into one tip.

Genome counts are per tip, so these duplicates do not inflate coverage.

A deduplication bug in the first v2 draft dropped the correct row whenever a misspelled duplicate hit the same tip. That draft therefore reported 613 names placed; the fix above supersedes it.

## 5. Results summary (IUCN 2026-1)

| Clade | Species | Threatened | EDGE spp. | EDGE with VGP genome | Top-50 EDGE with genome | Expected PD loss, median % (Gy) |
|---|---|---|---|---|---|---|
| Mammals | 5,911 | 1,284 | 610 | 38 (6.2%) | 7 | 10.1% (3.0) |
| Birds | 9,993 | 1,091 | 531 | 14 (2.6%) | 4 | 6.7% (5.7) |
| Squamates | 9,755 | 1,527 | 875 | 3 (0.3%) | 1 | 10.9% (13.8) |
| Amphibians | 7,238 | 2,481 | 1,314 | 3 (0.2%) | 1 | 15.5% (21.0) |
| Chondrichthyans | 1,192 | 401 | 263 | 10 (3.8%) | 5 | 15.1% (6.1) |

**VGP coverage is concentrated among distinct lineages** (`vgp_distinctness_enrichment.csv`):

| Clade | % of all species sequenced | % of clade ED2 captured | % of top-100 ED2 species | % of top-25 EDGE species | Orders (families) with a genome |
|---|---|---|---|---|---|
| Mammals | 4.0 | 6.5 | 11 | 16 | 28/28 |
| Birds | 2.5 | 3.7 | 13 | 12 | 39/40 |
| Squamates | 0.4 | 0.8 | 4 | 4 | 19/88 families |
| Amphibians | 0.6 | 1.1 | 3 | 4 | 3/3 |
| Chondrichthyans | 2.9 | 5.2 | 10 | 20 | 12/14 |

## 6. Mammals: v2 compared with v1

- EDGE2 ranks are highly concordant (Spearman ρ = 0.967 across all 5,911 species).
- There are 610 threatened EDGE species, compared with 586 in v1: 26 gained and 2 lost.
  - Of the 26 gained, 20 were DD or NE in the MDD layer and are threatened in IUCN 2026-1, 1 moved from LC to VU, and 5 kept their category but crossed the 50%-of-trees threshold.
  - The two lost are *Myrmecobius fasciatus* (numbat), which moved from EN to NT in IUCN 2026-1, and *Habromys lepturus*, which remains CR but dropped below the threshold.
- 38 threatened EDGE mammals now have a VGP genome, compared with 34 in v1.
  - New genomes: *Nasalis larvatus*, *Pontoporia blainvillei* and *Rangifer tarandus*.
  - Newly credited through a split: *Giraffa camelopardalis*, via the *G. tippelskirchi* genome.
  - No species lost its genome.
- `data/mammals/comparison_vs_previous_run.csv` gives the category and rank change for each species.

## 7. Figures (Arial; Wes Anderson Zissou1)

The palette is `#3B9AB2 #78B7C5 #EBCC2A #E1AF00 #F21A00`, assigned to the Red List categories as follows:

| Category | Colour |
|---|---|
| LC | #3B9AB2 |
| NT | #78B7C5 |
| VU | #EBCC2A |
| EN | #E96500 |
| CR | #F21A00 |
| EX/EW | near-black |
| DD | grey (#B8B8B8) |
| NM (no Red List match) | light grey (#E6E6E6) |

"Has VGP genome" is drawn in #3B9AB2. EN uses the Zissou1Continuous ramp colour #E96500 because the discrete #E1AF00 cannot be told apart from VU at marker size.

Every figure is saved as a PDF (Type-42 embedded Arial) and a 300-dpi PNG under `figures/<clade>/`:

| Figure | Content |
|---|---|
| fig1 | EDGE2 rank curve |
| fig2 | ED2 vs GE2 |
| fig3 | Top 50 EDGE species |
| fig4 | EDGE species by order/family |
| fig5 | Top 50 by EDGE2, all categories |
| fig6 | Top 50 by ED2 |
| fig7 | VGP genome gap |

The cross-clade figures are in `figures/cross_clade/`:

| Figure | Content |
|---|---|
| fig8 | Threatened PD and VGP coverage |
| fig9 | Red List composition of VGP genomes compared with each full clade |
| fig10 | Every EDGE species that already has a VGP genome, with its within-clade EDGE2 rank |
| fig11 | VGP coverage rising from all species to top-ED2 and top-EDGE species |

NM tips are drawn in light grey (#E6E6E6), separate from IUCN DD (#B8B8B8).

## 8. Caveats

- Tree taxonomies are older than IUCN 2026-1 (BirdLife 2012, Reptile Database 2015 and similar). Recently split taxa therefore inherit their parent's category through GBIF synonymy, and NM tips (0.7% to 3.2% per clade) receive an imputed GE2.
- Squamate and amphibian posteriors were taken from the first 1,000-tree file of each 10,000-tree set rather than thinned across all 10,000.
- VGP placement is name- and taxonomy-based. Splits are credited to the broader tip; domestic forms are placed but not credited; see `vgp_species_placement.csv`.
- Reproduce: `R/` (engine, chunk/aggregate scripts, sbatch files), then `python/build_tables.py` and `python/edge2_figs.py`.
