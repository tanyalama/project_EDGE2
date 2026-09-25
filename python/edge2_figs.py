# edge2_figs.py — publication figures for the EDGE2 x VGP project
# House style for this project: Arial + Wes Anderson Zissou1 palette.
# EN uses #E96500, the Zissou1Continuous ramp midpoint between #E1AF00 and #F21A00,
# because #E1AF00 is not separable from VU #EBCC2A at marker size.
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

ZISSOU = ["#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00"]
RL_ORDER = ["LC", "NT", "VU", "EN", "CR", "EW", "EX", "DD"]
RL_COL = {"LC": ZISSOU[0], "NT": ZISSOU[1], "VU": ZISSOU[2], "EN": "#E96500", "CR": ZISSOU[4],
          "EW": "#4A4A4A", "EX": "#1A1A1A", "DD": "#B8B8B8", "NE": "#B8B8B8", "NM": "#E6E6E6"}
RL_NAME = {"LC": "Least Concern", "NT": "Near Threatened", "VU": "Vulnerable", "EN": "Endangered",
           "CR": "Critically Endangered", "EW": "Extinct in the Wild", "EX": "Extinct", "DD": "Data Deficient", "NM": "No Red List match (GE2 imputed)"}
GENOME_COL = ZISSOU[0]      # has VGP genome
NOGENOME_COL = "#D9D9D9"    # no VGP genome
INK = "#2B2B2B"
CLADE_LABEL = {"mammals": "Mammals", "birds": "Birds", "squamates": "Squamates", "amphibians": "Amphibians",
               "chondrichthyans": "Chondrichthyans"}
TREE_SRC = {"mammals": "Upham et al. 2019", "birds": "Jetz et al. 2012 (Hackett backbone)",
            "squamates": "Tonini et al. 2016", "amphibians": "Jetz & Pyron 2018", "chondrichthyans": "Stein et al. 2018"}

def style():
    apply_figure_style(font="Arial", sizes=(8, 7, 6))
    mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial"], "pdf.fonttype": 42,
                         "ps.fonttype": 42, "axes.edgecolor": INK, "axes.labelcolor": INK,
                         "xtick.color": INK, "ytick.color": INK, "text.color": INK,
                         "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic",
                         "mathtext.bf": "Arial:bold"})

def sp(s):  # species label, italic
    return s.replace("_", " ")

def itm(s):  # mathtext italic species name (Arial italic via custom mathtext)
    return r"$\mathit{" + s.replace("_", r"\ ") + "}$"

def ital(ax_text):
    ax_text.set_style("italic")

def save(fig, path):
    fig.savefig(path + ".pdf", bbox_inches="tight")
    fig.savefig(path + ".png", dpi=300, bbox_inches="tight")
    plt.close(fig)

def rl_legend(ax, cats, loc="lower right", **kw):
    h = [Line2D([], [], marker="o", ls="", ms=5, mfc=RL_COL[c], mec="none", label=RL_NAME[c]) for c in cats]
    return ax.legend(handles=h, loc=loc, frameon=False, handletextpad=0.3, borderaxespad=0.2, **kw)

# ---------------- fig1: rank curve ----------------
def fig_rank_curve(d, cl, path):
    d = d.sort_values("EDGErank")
    fig, ax = plt.subplots(figsize=(3.5, 2.6))
    x = d.EDGErank.values; y = d.EDGEmed.values; q = d.EDGEiqr.values
    ax.fill_between(x, np.clip(y - q / 2, 0, None), y + q / 2, color=ZISSOU[1], alpha=.45, lw=0, label="± ½ IQR across 1,000 trees")
    ax.plot(x, y, color=ZISSOU[0], lw=1.2, label="Median EDGE2")
    E = d[d.is_EDGE_species == 1]
    ax.scatter(E.EDGErank, E.EDGEmed, s=4, color=ZISSOU[4], lw=0, zorder=3, label=f"EDGE species (n = {len(E):,})")
    ax.set_xscale("log"); ax.set_ylim(0, None)
    ax.set_xlabel("EDGE2 rank (log scale)"); ax.set_ylabel("EDGE2 score (Myr, median)")
    ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, p: f"{int(v):,}" if v >= 1 else ""))
    ax.set_title(f"{CLADE_LABEL[cl]}: EDGE2 score by rank ({len(d):,} species)", loc="left", fontsize=8)
    ax.legend(loc="upper right", frameon=False)
    ax.margins(0.03)
    save(fig, path)

# ---------------- fig2: ED2 vs GE2 ----------------
def fig_ed_ge(d, cl, path):
    fig, ax = plt.subplots(figsize=(3.5, 2.8))
    order = ["NM", "DD", "LC", "NT", "VU", "EN", "CR", "EW", "EX"]
    for c in order:
        s = d[d.RLcat == c]
        if len(s): ax.scatter(s.pextmed, s.EDmed, s=3 if c in ("LC", "DD") else 5, color=RL_COL[c], lw=0, alpha=.8, zorder=2 + order.index(c))
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("GE2: extinction probability (pext, median)"); ax.set_ylabel("ED2 (Myr, median)")
    ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, p: f"{v:g}"))
    ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, p: f"{v:g}"))
    ax.set_title(f"{CLADE_LABEL[cl]}: evolutionary distinctness vs extinction risk", loc="left", fontsize=8)
    present = [c for c in ["LC", "NT", "VU", "EN", "CR", "EX", "DD", "NM"] if (d.RLcat == c).any()]
    leg = rl_legend(ax, present, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    top = d.nlargest(3, "EDmed").sort_values("pextmed", ascending=False, kind="stable")
    lo, hi = d.EDmed.min(), d.EDmed.max(); ax.set_ylim(lo / 1.5, hi * 12)
    fig.canvas.draw()
    for i, (_, r) in enumerate(top.iterrows()):
        xa = ax.transAxes.inverted().transform(ax.transData.transform((r.pextmed, r.EDmed)))[0]
        right = xa > 0.5
        t = ax.annotate(sp(r.species), (r.pextmed, r.EDmed), xytext=(xa - 0.03 if right else xa + 0.05, 0.965 - 0.065 * i),
                        textcoords="axes fraction", fontsize=6, va="center", ha="right" if right else "left",
                        arrowprops=dict(arrowstyle="-", color="#8C8C8C", lw=0.5, shrinkA=2, shrinkB=2,
                                        relpos=(1, 0.5) if right else (0, 0.5)))
        ital(t)
    ax.margins(0.04)
    save(fig, path)

# ---------------- top-50 dot/IQR plots (fig3/5/6) ----------------
def fig_top50(d, cl, path, mode="edge"):
    if mode == "edge":
        s = d[d.is_EDGE_species == 1].sort_values("EDGEmed", ascending=False).head(50); val, iqr, xl = "EDGEmed", "EDGEiqr", "EDGE2 score (Myr; median, IQR)"
        title = f"{CLADE_LABEL[cl]}: top 50 EDGE species"
    elif mode == "all":
        s = d.sort_values("EDGEmed", ascending=False).head(50); val, iqr, xl = "EDGEmed", "EDGEiqr", "EDGE2 score (Myr; median, IQR)"
        title = f"{CLADE_LABEL[cl]}: top 50 by EDGE2, all Red List categories"
    else:
        s = d.sort_values("EDmed", ascending=False).head(50); val, iqr, xl = "EDmed", "EDiqr", "ED2 (Myr; median, IQR)"
        title = f"{CLADE_LABEL[cl]}: top 50 by evolutionary distinctness (ED2)"
    s = s.iloc[::-1]
    fig, ax = plt.subplots(figsize=(4.4, 7.2))
    yy = np.arange(len(s))
    ax.hlines(yy, s[val] - s[iqr] / 2, s[val] + s[iqr] / 2, color="#8C8C8C", lw=1.0, zorder=1)
    ax.scatter(s[val], yy, s=22, c=[RL_COL.get(c, "#B8B8B8") for c in s.RLcat], edgecolors=INK, linewidths=0.3, zorder=3)
    labels = []
    for _, r in s.iterrows():
        cn = f" ({r.common_name})" if isinstance(r.common_name, str) and len(r.common_name) < 28 else ""
        labels.append(itm(r.species) + cn.replace("$", ""))
    ax.set_yticks(yy); ax.set_yticklabels(labels)
    # genome markers in a strip right of the data
    xmax = (s[val] + s[iqr] / 2).max() * 1.06
    for i, (_, r) in enumerate(s.iterrows()):
        ax.scatter(xmax, i, marker="s", s=16, color=GENOME_COL if r.has_vgp_genome else "white",
                   edgecolors=GENOME_COL, linewidths=0.6, zorder=3, clip_on=False)
    ax.text(xmax, len(s) + 0.3, "VGP\ngenome", ha="center", va="bottom", fontsize=6)
    ax.set_xlim(0, xmax * 1.04); ax.set_ylim(-0.8, len(s) - 0.2)
    ax.set_xlabel(xl); ax.set_title(title, loc="left", fontsize=8, pad=18)
    ax.spines["left"].set_visible(False); ax.tick_params(axis="y", length=0)
    present = [c for c in ["LC", "NT", "VU", "EN", "CR", "EW", "EX", "DD", "NM"] if (s.RLcat == c).any()]
    h = [Line2D([], [], marker="o", ls="", ms=5, mfc=RL_COL[c], mec=INK, mew=.3, label=c) for c in present]
    h += [Line2D([], [], marker="s", ls="", ms=4.5, mfc=GENOME_COL, mec=GENOME_COL, label="VGP genome"),
          Line2D([], [], marker="s", ls="", ms=4.5, mfc="white", mec=GENOME_COL, label="no genome")]
    ax.legend(handles=h, loc="lower right", frameon=False, ncol=1, handletextpad=0.3, bbox_to_anchor=(0.97, 0.0))
    save(fig, path)

# ---------------- fig4: ordinal summary ----------------
def fig_ordinal(d, cl, path, level="order", topn=15):
    E = d[d.is_EDGE_species == 1]
    ct = E.groupby([level, "RLcat"]).size().unstack(fill_value=0).reindex(columns=["VU", "EN", "CR"], fill_value=0)
    ct["tot"] = ct.sum(1); ct = ct.sort_values("tot", ascending=False).head(topn).iloc[::-1]
    med = E.groupby(level).EDGEmed.median().reindex(ct.index)
    nsp = d.groupby(level).size().reindex(ct.index)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(5.6, 0.2 * len(ct) + 1.2), sharey=True, gridspec_kw={"width_ratios": [1.6, 1], "wspace": 0.06})
    left = np.zeros(len(ct)); yy = np.arange(len(ct))
    for c in ["CR", "EN", "VU"]:
        a1.barh(yy, ct[c], left=left, color=RL_COL[c], height=.72, label=RL_NAME[c], edgecolor="white", lw=.3); left += ct[c].values
    for i, (t, n) in enumerate(zip(ct.tot, nsp)):
        a1.text(t + ct.tot.max() * .01, i, f"{t}  ({100*t/n:.0f}%)", va="center", fontsize=6)
    a1.set_yticks(yy); a1.set_yticklabels(ct.index); a1.set_xlabel("EDGE species (n; % of lineage)")
    a1.set_xlim(0, ct.tot.max() * 1.25)
    a1.legend(loc="lower right", frameon=False)
    a2.hlines(yy, 0, med.values, color="#8C8C8C", lw=.8); a2.scatter(med.values, yy, s=16, color=ZISSOU[0], zorder=3)
    a2.set_xlabel("Median EDGE2 of EDGE species (Myr)"); a2.set_xlim(0, med.max() * 1.12)
    a2.tick_params(axis="y", length=0)
    lvl = "family" if level == "family" else "order"
    a1.set_title(f"{CLADE_LABEL[cl]}: EDGE species by {lvl} (top {len(ct)})", loc="left", fontsize=8)
    save(fig, path)

# ---------------- fig7: VGP genome gap ----------------
def fig_vgp_gap(d, cl, path, level="order", topn=30, ntax=15):
    E = d[d.is_EDGE_species == 1]
    miss = E[~E.has_vgp_genome].sort_values("EDGEmed", ascending=False).head(topn).iloc[::-1]
    fig = plt.figure(figsize=(7.2, 5.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.1, 1], wspace=0.55)
    a = fig.add_subplot(gs[0]); b = fig.add_subplot(gs[1])
    yy = np.arange(len(miss))
    a.hlines(yy, 0, miss.EDGEmed, color="#8C8C8C", lw=.8)
    a.scatter(miss.EDGEmed, yy, s=18, c=[RL_COL[c] for c in miss.RLcat], edgecolors=INK, linewidths=.3, zorder=3)
    a.set_yticks(yy); a.set_yticklabels([sp(s) for s in miss.species], fontstyle="italic")
    a.set_xlabel("EDGE2 score (Myr, median)"); a.set_xlim(0, miss.EDGEmed.max() * 1.08)
    a.tick_params(axis="y", length=0); a.spines["left"].set_visible(False)
    a.set_title(f"Top {len(miss)} EDGE species without a VGP genome", loc="left", fontsize=8)
    present = [c for c in ["VU", "EN", "CR"] if (miss.RLcat == c).any()]
    rl_legend(a, present, loc="lower right")
    t = E.groupby(level).agg(n=("species", "size"), g=("has_vgp_genome", "sum")).sort_values("n", ascending=False).head(ntax).iloc[::-1]
    yb = np.arange(len(t))
    b.barh(yb, t.g, color=GENOME_COL, height=.72, label="with VGP genome")
    b.barh(yb, t.n - t.g, left=t.g, color=NOGENOME_COL, height=.72, label="without genome")
    for i, (n, g) in enumerate(zip(t.n, t.g)):
        b.text(n + t.n.max() * .01, i, f"{g}/{n}", va="center", fontsize=6)
    b.set_yticks(yb); b.set_yticklabels(t.index); b.set_xlabel("EDGE species (n)"); b.set_xlim(0, t.n.max() * 1.18)
    b.tick_params(axis="y", length=0)
    b.set_title(f"VGP coverage of EDGE species by {'family' if level=='family' else 'order'}", loc="left", fontsize=8)
    b.legend(loc="lower right", frameon=False)
    ng = int(E.has_vgp_genome.sum())
    fig.suptitle(f"{CLADE_LABEL[cl]}: {ng} of {len(E):,} EDGE species ({100*ng/len(E):.1f}%) have a VGP reference genome",
                 x=0.02, ha="left", y=0.995, fontsize=8, fontweight="bold")
    save(fig, path)

# ---------------- cross-clade figures ----------------
CL_ORDER = ["mammals", "birds", "squamates", "amphibians", "chondrichthyans"]

def fig_cross_clade(full, S, epd, path):
    S = S.set_index("clade").loc[CL_ORDER]
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.5), gridspec_kw={"wspace": 0.12, "width_ratios": [1, 1.25, 1.1]}, sharey=True)
    yy = np.arange(len(CL_ORDER))[::-1]
    a = axs[0]
    a.hlines(yy, S.pct_PD_at_risk_q25, S.pct_PD_at_risk_q75, color=ZISSOU[1], lw=3)
    a.scatter(S.pct_PD_at_risk, yy, s=22, color=ZISSOU[0], zorder=3)
    for y, v, g in zip(yy, S.pct_PD_at_risk, S.ePDloss_Gy_median):
        a.text(v + 0.6, y + 0.22, f"{v:.1f}% ({g:.1f} Gy)", fontsize=6, va="bottom")
    a.set_yticks(yy); a.set_yticklabels([CLADE_LABEL[c] for c in CL_ORDER])
    a.set_xlim(0, S.pct_PD_at_risk.max() * 1.5); a.set_ylim(-0.6, len(yy) - 0.3)
    a.set_xlabel("Expected PD loss (% of total PD)")
    a.set_title("Threatened evolutionary history", loc="left", fontsize=8)
    b = axs[1]
    g = S.n_EDGE_with_genome.values; n = S.n_EDGE.values
    b.barh(yy, g, color=GENOME_COL, height=.6, label="with VGP genome")
    b.barh(yy, n - g, left=g, color=NOGENOME_COL, height=.6, label="without genome")
    for y, gi, ni in zip(yy, g, n):
        b.text(ni + n.max() * .015, y, f"{gi}/{ni:,} ({100*gi/ni:.1f}%)", va="center", fontsize=6)
    b.set_xlim(0, n.max() * 1.45); b.set_xlabel("EDGE species (n)")
    b.set_title("VGP coverage of EDGE species", loc="left", fontsize=8)
    b.legend(loc="lower right", frameon=False, bbox_to_anchor=(1.0, -0.02))
    c = axs[2]
    off = {"VU": -0.2, "EN": 0, "CR": 0.2}
    for cat in ["VU", "EN", "CR"]:
        vals = []
        for cl in CL_ORDER:
            E = full[cl][full[cl].is_EDGE_species == 1]; e = E[E.RLcat == cat]
            vals.append(100 * e.has_vgp_genome.mean() if len(e) else np.nan)
        c.scatter(vals, yy + off[cat], s=18, color=RL_COL[cat], edgecolors=INK, linewidths=.3, zorder=3, label=RL_NAME[cat])
    c.set_xlabel("EDGE species with a VGP genome (%)"); c.set_xlim(-0.4, None)
    c.set_title("Coverage by Red List category", loc="left", fontsize=8)
    c.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=3, frameon=False, columnspacing=0.8, handletextpad=0.2)
    c.margins(x=0.08)
    for ax in axs[1:]: ax.tick_params(axis="y", length=0)
    save(fig, path)

def fig_vgp_rl_bias(full, path):
    cats = ["LC", "NT", "VU", "EN", "CR", "EX", "DD", "NM"]
    fig, ax = plt.subplots(figsize=(7.2, 2.6))
    rows = []; labels = []
    for cl in CL_ORDER:
        d = full[cl]
        for lab, sub in [("all species", d), ("VGP genomes", d[d.has_vgp_genome])]:
            rc = sub.RLcat.replace({"EW": "EX", "NE": "NM"}).value_counts(normalize=True).reindex(cats, fill_value=0)
            rows.append(rc.values * 100); labels.append((cl, lab, len(sub)))
    rows = np.array(rows)
    y = []; pos = 0
    for i in range(len(CL_ORDER)):
        y += [pos, pos - 0.8]; pos -= 2.2
    y = np.array(y)
    left = np.zeros(len(rows))
    for j, cat in enumerate(cats):
        ax.barh(y, rows[:, j], left=left, height=.7, color=RL_COL[cat], edgecolor="white", lw=.3, label=RL_NAME[cat])
        left += rows[:, j]
    ax.set_yticks(y); ax.set_yticklabels([f"{CLADE_LABEL[c]} — {l} (n = {n:,})" for c, l, n in labels])
    ax.set_xlim(0, 100); ax.set_xlabel("Share of species (%)")
    ax.tick_params(axis="y", length=0)
    ax.set_title("Red List composition of VGP-sequenced species compared with each full clade", loc="left", fontsize=8)
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False)
    save(fig, path)


# ---------------- positive figures: what VGP has achieved ----------------
def fig_sequenced_edge(full, path):
    """All EDGE species that already have a VGP reference genome, by clade."""
    left = ["mammals"]; right = ["birds", "chondrichthyans", "squamates", "amphibians"]
    sel = {cl: full[cl][(full[cl].is_EDGE_species == 1) & full[cl].has_vgp_genome].sort_values("EDGEmed", ascending=False) for cl in CL_ORDER}
    nl = sum(len(sel[c]) for c in left); nr = sum(len(sel[c]) for c in right)
    fig = plt.figure(figsize=(7.2, 0.155 * max(nl, nr + 3) + 1.0))
    outer = fig.add_gridspec(1, 2, wspace=1.35)
    xmax = max(sel[c].EDGEmed.max() for c in CL_ORDER if len(sel[c])) * 1.08
    def col(spec, cls):
        gs = spec.subgridspec(len(cls), 1, height_ratios=[max(len(sel[c]), 1) + 1.2 for c in cls], hspace=0.35)
        axs = []
        for i, cl in enumerate(cls):
            ax = fig.add_subplot(gs[i]); s = sel[cl].iloc[::-1]; yy = np.arange(len(s))
            ax.hlines(yy, 0, s.EDGEmed, color="#8C8C8C", lw=.8)
            ax.scatter(s.EDGEmed, yy, s=20, c=[RL_COL[c] for c in s.RLcat], edgecolors=INK, linewidths=.3, zorder=3)
            labs = []
            for _, r in s.iterrows():
                cn = f" ({r.common_name})" if isinstance(r.common_name, str) and len(r.common_name) < 26 else ""
                labs.append(itm(r.species) + cn.replace("$", ""))
            ax.set_yticks(yy); ax.set_yticklabels(labs, fontsize=5.5)
            for y, (_, r) in zip(yy, s.iterrows()):
                ax.text(r.EDGEmed + s.EDGEmed.max() * .04, y, f"#{int(r.EDGErank)}", va="center", fontsize=5.5, color="#5A5A5A")
            xm = s.EDGEmed.max() * 1.18; ax.set_xlim(0, xm); ax.set_ylim(-0.7, len(s) - 0.3)
            ax.tick_params(axis="y", length=0); ax.spines["left"].set_visible(False)
            E = full[cl][full[cl].is_EDGE_species == 1]
            ax.set_title(f"{CLADE_LABEL[cl]}: {len(s)} of {len(E):,} EDGE species sequenced", loc="left", fontsize=7, fontweight="bold")
            axs.append(ax)
        axs[-1].set_xlabel("EDGE2 score (Myr, median; scale differs by clade)")
        return axs
    a = col(outer[0], left); b = col(outer[1], right)
    h = [Line2D([], [], marker="o", ls="", ms=5, mfc=RL_COL[c], mec=INK, mew=.3, label=RL_NAME[c]) for c in ["VU", "EN", "CR"]]
    b[-1].legend(handles=h, loc="upper center", bbox_to_anchor=(0.3, -1.1), ncol=3, frameon=False, handletextpad=0.2, columnspacing=0.8)
    fig.suptitle("EDGE species with a VGP reference genome (#n = EDGE2 rank within clade)", x=0.01, ha="left", y=0.93, fontsize=8, fontweight="bold")
    save(fig, path)

def fig_vgp_enrichment(full, M, path):
    """VGP coverage rises with evolutionary distinctness and EDGE rank."""
    M = M.set_index("clade").loc[CL_ORDER]
    fig, a = plt.subplots(figsize=(4.8, 2.9))
    yy = np.arange(len(CL_ORDER))[::-1]
    series = [("pct_species", "All species", "#B8B8B8", "o"), ("pct_ED2", "Share of clade ED2 (Myr) captured", ZISSOU[1], "o"),
              ("pct_EDGE_species", "EDGE species", "white", "s"),
              ("pct_top100_ED", "Top 100 species by ED2", ZISSOU[0], "D"), ("pct_top25_EDGE", "Top 25 EDGE species", INK, "^")]
    cols = [c for c, *_ in series]
    a.hlines(yy, M[cols].min(1), M[cols].max(1), color="#D0D0D0", lw=1, zorder=1)
    for col, lab, colr, mk in series:
        a.scatter(M[col], yy, s=28, color=colr, marker=mk, edgecolors=INK, linewidths=.4, zorder=3, label=lab)
    for y, cl in zip(yy, CL_ORDER):
        a.text(M.loc[cl, cols].max() + 0.8, y, f"{int(M.loc[cl,'n_seq'])} genomes; {int(M.loc[cl,'orders_with_genome'])}/{int(M.loc[cl,'orders'])} {'families' if cl=='squamates' else 'orders'}",
               va="center", fontsize=5.5, color="#5A5A5A")
    a.set_yticks(yy); a.set_yticklabels([CLADE_LABEL[c] for c in CL_ORDER]); a.tick_params(axis="y", length=0)
    a.set_xlabel("With a VGP reference genome (%)"); a.set_xlim(0, M[cols].max().max() * 1.45)
    a.legend(loc="upper left", bbox_to_anchor=(1.0, 1.0), frameon=False, handletextpad=0.3)
    a.set_title("VGP coverage is highest among the most distinct and highest-ranked species", loc="left", fontsize=8)
    save(fig, path)
