os.makedirs("out/data",exist_ok=True)
thr=["VU","EN","CR"]
grp_level={"mammals":"order","birds":"order","amphibians":"order","chondrichthyans":"order","squamates":"family"}
full={}; summ=[]
Rr=R.set_index(["clade","species"])
for cl in clades:
    d=res[cl].copy()
    d["EDrank"]=d.EDmed.rank(ascending=False,method="first").astype(int)
    meta=Rr.loc[cl].reindex(d.species)
    d["order"]=meta.order.values; d["family"]=meta.family.values; d["match"]=meta.match.values; d["matched_IUCN_name"]=meta.matched_IUCN_name.values
    hv=Hk.loc[cl].reindex(d.species)
    d["vgp_name"]=hv.vgp_name.values; d["vgp_accession"]=hv.vgp_accession.values; d["vgp_placement_tier"]=hv.vgp_placement_tier.values
    d["has_vgp_genome"]=hv.vgp_accession.notna().values
    onl=set(A[(A.clade==cl)&A.tree_tip.notna()].tree_tip); d["on_vgp_list"]=d.species.isin(onl).values
    nm=d.match=="unmatched_DD_placeholder"; d.loc[nm,"RLcat"]="NM"; d.loc[nm,"match"]="no_RedList_match"
    d["is_EDGE_species"]=((d.isEDGEsp==1)&d.RLcat.isin(thr)).astype(int)
    cols=["EDGErank","EDrank","species","order","family","RLcat","EDGEmed","EDGEiqr","EDmed","EDiqr","pextmed","pextiqr","TBLmn","isEDGEsp_frac","is_EDGE_species","n_trees","match","matched_IUCN_name","on_vgp_list","has_vgp_genome","vgp_placement_tier","vgp_name","vgp_accession"]
    d=d[cols].sort_values("EDGErank")
    d.insert(3,"common_name",[k2c.get(a) or k2c.get(b) for a,b in zip(d.species,d.matched_IUCN_name)])
    full[cl]=d
    p=f"out/data/{cl}"; os.makedirs(p,exist_ok=True)
    d.to_csv(f"{p}/EDGE2_ranked_species_FULL.csv",index=False)
    E=d[d.is_EDGE_species==1]; E.to_csv(f"{p}/EDGE_species_list.csv",index=False)
    d[d.RLcat.isin(thr)&(d.isEDGEsp_frac>=0.25)&(d.isEDGEsp_frac<0.5)].to_csv(f"{p}/EDGE_borderline_list.csv",index=False)
    d[d.RLcat.isin(["DD","NE","NM"])].sort_values("EDmed",ascending=False).to_csv(f"{p}/EDGE_DD_watchlist.csv",index=False)
    E[~E.has_vgp_genome].to_csv(f"{p}/EDGE_species_missing_VGP_genome.csv",index=False)
    lvl=grp_level[cl]
    t=d.groupby(lvl).agg(n_species=("species","size"),n_EDGE=("is_EDGE_species","sum"),median_EDGE2=("EDGEmed","median"),
        n_vgp_genome=("has_vgp_genome","sum")).reset_index()
    t2=E.groupby(lvl).has_vgp_genome.sum().rename("n_EDGE_with_genome"); t=t.merge(t2,on=lvl,how="left").fillna({"n_EDGE_with_genome":0})
    t["pct_EDGE"]=100*t.n_EDGE/t.n_species; t.sort_values("n_EDGE",ascending=False).to_csv(f"{p}/taxonomic_summary_by_{lvl}.csv",index=False)
    e=epd[cl]
    s=dict(clade=cl,n_species=len(d),n_threatened=d.RLcat.isin(thr).sum(),n_EDGE=len(E),n_borderline=((d.RLcat.isin(thr))&(d.isEDGEsp_frac>=0.25)&(d.isEDGEsp_frac<0.5)).sum(),
           n_no_RedList_match=(d.match=="no_RedList_match").sum(),n_on_vgp=d.on_vgp_list.sum(),n_vgp_genome=d.has_vgp_genome.sum(),
           n_EDGE_with_genome=E.has_vgp_genome.sum(),pct_EDGE_with_genome=100*E.has_vgp_genome.mean(),
           top10_EDGE_with_genome=E.head(10).has_vgp_genome.sum(), top50_EDGE_with_genome=E.head(50).has_vgp_genome.sum(),
           PD_Gy_median=e.PD.median()/1000, ePDloss_Gy_median=e.ePDloss.median()/1000, pct_PD_at_risk=100*(e.ePDloss/e.PD).median(),
           pct_PD_at_risk_q25=100*(e.ePDloss/e.PD).quantile(.25), pct_PD_at_risk_q75=100*(e.ePDloss/e.PD).quantile(.75))
    for c_ in thr: s[f"{c_}_EDGE_with_genome"]=f"{E[E.RLcat==c_].has_vgp_genome.sum()}/{(E.RLcat==c_).sum()}"
    summ.append(s)
S=pd.DataFrame(summ); S.to_csv("out/data/cross_clade_summary.csv",index=False)
