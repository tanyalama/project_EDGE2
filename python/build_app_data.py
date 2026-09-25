"""Build the compact data files for the Shiny dashboard (app/data/) from the repo's data/ tables.
Run from the repo root:  python python/build_app_data.py"""
import os
import pandas as pd

CLADES = ["mammals", "birds", "squamates", "amphibians", "chondrichthyans"]
KEEP = ["EDGErank", "EDrank", "species", "common_name", "order", "family", "RLcat", "EDGEmed", "EDGEiqr",
        "EDmed", "EDiqr", "pextmed", "isEDGEsp_frac", "is_EDGE_species", "has_vgp_genome",
        "vgp_placement_tier", "vgp_name", "vgp_accession"]
os.makedirs("app/data", exist_ok=True)

sp = pd.concat([pd.read_csv(f"data/{cl}/EDGE2_ranked_species_FULL.csv")[KEEP].assign(clade=cl)
                for cl in CLADES], ignore_index=True)
sp = sp[["clade"] + KEEP]
sp["species"] = sp.species.str.replace("_", " ")
sp["vgp_name"] = sp.vgp_name.fillna("").str.replace("_", " ")
for c in ["EDGEmed", "EDGEiqr", "EDmed", "EDiqr"]:
    sp[c] = sp[c].round(3)
sp["pextmed"] = sp.pextmed.round(4)
sp["isEDGEsp_frac"] = sp.isEDGEsp_frac.round(3)
sp[["common_name", "vgp_placement_tier", "vgp_accession"]] = sp[["common_name", "vgp_placement_tier", "vgp_accession"]].fillna("")
sp.to_csv("app/data/species.csv", index=False)

S = pd.read_csv("data/cross_clade_summary.csv")
S[["clade", "n_species", "n_threatened", "n_EDGE", "n_borderline", "n_no_RedList_match", "n_vgp_genome",
   "n_EDGE_with_genome", "pct_EDGE_with_genome", "top50_EDGE_with_genome", "PD_Gy_median", "ePDloss_Gy_median",
   "pct_PD_at_risk", "pct_PD_at_risk_q25", "pct_PD_at_risk_q75"]].round(3).to_csv("app/data/summary.csv", index=False)

pd.read_csv("data/vgp_distinctness_enrichment.csv").to_csv("app/data/enrichment.csv", index=False)

pl = pd.read_csv("data/vgp_species_placement.csv")
for c in ["vgp_name", "vgp_name_corrected", "tree_tip"]:
    pl[c] = pl[c].fillna("").str.replace("_", " ")
pl.to_csv("app/data/vgp_placement.csv", index=False)
print(len(sp), "species;", len(pl), "VGP names")
