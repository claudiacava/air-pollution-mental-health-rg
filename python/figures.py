"""Figures 2-4 and Supplementary Figure S1.
Usage: python figures.py <results/main folder> <output folder>"""
import sys, csv, statistics as st
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

RES, OUT = sys.argv[1], sys.argv[2]
POL, MH, SES = "#1F78A8", "#D96C2C", "#7A5AA6"
INK, MUTED, GRID = "#2B2B2B", "#5F6368", "#E3E5E8"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED})
R = lambda f: list(csv.DictReader(open(f"{RES}/{f}", encoding="utf-8")))
t1, t2 = R("T1_trait_h2_intercept_ratio.csv"), R("T2_rg_pollution_mental_health.csv")
f = float
POLN = {"24003_irnt": "NO₂ (2010)", "24006_irnt": "PM2.5 (2010)", "24019_irnt": "PM10 (2007)"}
SESN = {"189_irnt": "Townsend deprivation", "6138_1": "College degree", "738": "Household income",
        "20118_5": "Urban residence"}

# ---------------- Figure 2: attenuation ratio -----------------------------------
fig, ax = plt.subplots(figsize=(7, 3.6))
rows = [(POLN[r["trait"]], f(r["ratio"]), f(r["ratio_se"]), POL) for r in t1 if r["trait"] in POLN]
rows += [(SESN[r["trait"]], f(r["ratio"]), f(r["ratio_se"]), SES) for r in t1 if r["trait"] in SESN]
rows.sort(key=lambda x: (x[3] != POL, -x[1]))
mh = [f(r["ratio"]) for r in t1 if r["domain"] in ("mental_health", "wellbeing")]
labels = [r[0] for r in rows] + ["Mental health traits (n = %d)" % len(mh)]
y = np.arange(len(labels))[::-1]
ax.axvline(0, color=MUTED, lw=0.8, ls=(0, (3, 3)))
for yi, (lab, v, se, col) in zip(y[:-1], rows):
    ax.plot([v - 1.96 * se, v + 1.96 * se], [yi, yi], color=col, lw=2, solid_capstyle="round")
    ax.scatter(v, yi, s=46, color=col, edgecolor="white", linewidth=1.5, zorder=3)
    ax.text(v + 1.96 * se + 0.02, yi, f"{v:.2f}", va="center", fontsize=8, color=INK)
rng = np.random.default_rng(1)
ax.scatter(mh, y[-1] + rng.uniform(-0.22, 0.22, len(mh)), s=14, color=MH, alpha=0.75,
           edgecolor="white", linewidth=0.6, zorder=3)
med = st.median(mh)
ax.plot([med, med], [y[-1] - 0.35, y[-1] + 0.35], color=INK, lw=2)
ax.text(med + 0.02, y[-1] + 0.42, f"median {med:.2f}", fontsize=8, color=INK)
ax.set_yticks(y); ax.set_yticklabels(labels, color=INK)
for tl, lab in zip(ax.get_yticklabels(), labels):
    tl.set_color(POL if lab in POLN.values() else SES if lab in SESN.values() else MH)
ax.set_xlabel("LDSC attenuation ratio (95% CI)")
ax.grid(axis="x", color=GRID, lw=0.6); ax.set_axisbelow(True)
for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
ax.tick_params(axis="y", length=0)
fig.tight_layout()
fig.savefig(f"{OUT}/Figure2_attenuation_ratio.png", dpi=300)
fig.savefig(f"{OUT}/Figure2_attenuation_ratio.pdf")

# ---------------- Figure 3 (single panel): heatmap of rg by subdomain ------------
SUB = [("depression", "Depression"), ("anxiety", "Anxiety"), ("trauma", "Trauma-related"),
       ("unhappiness", "Unhappinessᵃ"), ("mood", "Mood instability"),
       ("fatigue", "Fatigue"), ("distress", "General distress"), ("bipolar", "Bipolar disorder")]
def short(d, code):
    for a, b in [("Non-cancer illness code, self-reported: ", "Self-reported "),
                 ("Bipolar and major depression status: ", ""),
                 ("Mental health problems ever diagnosed by a professional: ", "Diagnosed: "),
                 ("Weight change during worst episode of depression: ", "Weight change (worst depression): "),
                 (" or was on a diet", ""), ("Both gained and lost some weight during the episode", "Both")]:
        d = d.replace(a, b)
    d = d[0].upper() + d[1:]
    return (d if len(d) <= 46 else d[:44].rstrip() + "…") + f"  [{code}]"
cols = ["24003_irnt", "24006_irnt", "24019_irnt"]
cell = {(r["trait1"], r["trait2"]): (f(r["rg"]), f(r["FDR"])) for r in t2}
traits = {}
for r in t2: traits[r["trait2"]] = (r["description2"], r["subdomain"])
order, groups = [], []
for key, name in SUB:
    tr = [t for t, (d, s) in traits.items() if s == key]
    tr.sort(key=lambda t: -np.mean([cell[(c, t)][0] for c in cols]))
    if tr: groups.append((name, len(order), len(tr))); order += tr
M = np.array([[cell[(c, t)][0] for c in cols] for t in order])
fig, ax = plt.subplots(figsize=(8.0, 0.22 * len(order) + 1.6))
im = ax.imshow(M, cmap="RdBu_r", vmin=-0.6, vmax=0.6, aspect="auto")
for i, t in enumerate(order):
    for j, c in enumerate(cols):
        rg, q = cell[(c, t)]
        ax.text(j, i, f"{rg:.2f}" + ("*" if q < 0.05 else ""), ha="center", va="center", fontsize=6.8,
                color="white" if abs(rg) > 0.42 else INK, fontweight="bold" if q < 0.05 else "normal")
ax.set_xticks(range(3)); ax.set_xticklabels([POLN[c] for c in cols], color=INK, fontsize=8.5)
ax.xaxis.tick_top()
ax.set_yticks(range(len(order))); ax.set_yticklabels([short(*traits[t][:1], t) for t in order], fontsize=6.8, color=INK)
for name, start, n in groups:
    if start: ax.axhline(start - 0.5, color="white", lw=3)
    ax.text(2.62, start + (n - 1) / 2, name, va="center", ha="left", fontsize=8, color=MH, fontweight="bold")
ax.set_xticks(np.arange(-0.5, 3, 1), minor=True); ax.grid(which="minor", axis="x", color="white", lw=2)
ax.tick_params(which="both", length=0)
for s in ax.spines.values(): s.set_visible(False)
cb = fig.colorbar(im, ax=ax, orientation="horizontal", fraction=0.025, pad=0.02, aspect=30)
cb.set_label("Genetic correlation (rg); * FDR < 0.05", fontsize=8, color=INK); cb.outline.set_visible(False)
fig.subplots_adjust(left=0.50, right=0.81, top=0.95, bottom=0.07)
fig.savefig(f"{OUT}/Figure3_rg_heatmap.png", dpi=300)
fig.savefig(f"{OUT}/Figure3_rg_heatmap.pdf")

# ---------------- Figure 4: SES genetic correlations and partial rg ---------------
t3, t5 = R("T3_rg_socioeconomic.csv"), R("T5_partial_rg_socioeconomic.csv")
PCOL = {"24003_irnt": "#2563B8", "24006_irnt": "#17A398", "24019_irnt": "#B8860B"}
fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 4.6), gridspec_kw={"width_ratios": [1, 1.5]})
ses_order = ["189_irnt", "20118_5", "738", "6138_1"]
for k, s in enumerate(ses_order):
    for j, p in enumerate(cols):
        r = [x for x in t3 if x["trait1"] == s and x["trait2"] == p][0]
        yy = len(ses_order) - 1 - k + (1 - j) * 0.22
        v, se = f(r["rg"]), f(r["se"])
        a1.plot([v - 1.96 * se, v + 1.96 * se], [yy, yy], color=PCOL[p], lw=2, solid_capstyle="round")
        a1.scatter(v, yy, s=36, color=PCOL[p], edgecolor="white", lw=1.2, zorder=3, label=POLN[p] if k == 0 else None)
a1.axvline(0, color=MUTED, lw=0.8, ls=(0, (3, 3)))
a1.set_yticks(range(len(ses_order))); a1.set_yticklabels([SESN[s] for s in ses_order][::-1], color=INK)
a1.set_xlabel("Genetic correlation with air pollution (rg, 95% CI)"); a1.set_xlim(-0.7, 1.05)
a1.set_title("A  Air pollution and socioeconomic position", loc="left", fontsize=10, color=INK, fontweight="bold")
a1.legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=3)
sigpairs = {(r["trait1"], r["trait2"]) for r in t2 if f(r["FDR"]) < 0.05 and f(r["rg"]) > 0}
conds = [("none", "Unadjusted"), ("Degree", "Degree"), ("Townsend", "Townsend"),
         ("Townsend_Degree", "Townsend +\ndegree"), ("All_SES", "All four\nSES traits")]
rng = np.random.default_rng(2)
for j, p in enumerate(cols):
    pairs = sorted(t for (pp, t) in sigpairs if pp == p)
    for k, (c, lab) in enumerate(conds):
        if c == "none": vals = [cell[(p, t)][0] for t in pairs]
        else: vals = [f(r["partial_rg"]) for r in t5 if r["adjustment"] == c and r["pollutant"] == p and r["trait"] in pairs]
        x = k + (j - 1) * 0.24
        a2.scatter(x + rng.uniform(-0.06, 0.06, len(vals)), vals, s=12, color=PCOL[p], alpha=0.7, edgecolor="white", lw=0.5, zorder=3)
        m = st.median(vals); a2.plot([x - 0.1, x + 0.1], [m, m], color=INK, lw=2, zorder=4)
a2.axhline(0, color=MUTED, lw=0.8, ls=(0, (3, 3)))
a2.set_xticks(range(len(conds))); a2.set_xticklabels([l for _, l in conds], color=INK, fontsize=8.5)
a2.set_ylabel("rg / partial rg with mental health trait")
a2.set_title("B  Pollution–mental health rg, adjusted for SES", loc="left", fontsize=10, color=INK, fontweight="bold")
for a in (a1, a2):
    a.grid(axis="x" if a is a1 else "y", color=GRID, lw=0.6); a.set_axisbelow(True)
    for s in ("top", "right"): a.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig(f"{OUT}/Figure4_SES_partial_rg.png", dpi=300)
fig.savefig(f"{OUT}/Figure4_SES_partial_rg.pdf")

# ---------------- Figure 3 (combined): A heatmap + B subdomain means -----------
t6 = R("T6_rg_by_subdomain.csv")
fig = plt.figure(figsize=(11, 0.21 * len(order) + 1.9))
gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 0.55], height_ratios=[0.55, 0.45], wspace=0.30, hspace=0.0, left=0.31, right=0.84, top=0.94, bottom=0.07)
ax = fig.add_subplot(gs[:, 0]); bx = fig.add_subplot(gs[0, 1])
im = ax.imshow(M, cmap="RdBu_r", vmin=-0.6, vmax=0.6, aspect="auto")
for i, t in enumerate(order):
    for j, c in enumerate(cols):
        rg, q = cell[(c, t)]
        ax.text(j, i, f"{rg:.2f}" + ("*" if q < 0.05 else ""), ha="center", va="center", fontsize=6.6,
                color="white" if abs(rg) > 0.42 else INK, fontweight="bold" if q < 0.05 else "normal")
ax.set_xticks(range(3)); ax.set_xticklabels([POLN[c] for c in cols], color=INK, fontsize=8); ax.xaxis.tick_top()
ax.set_yticks(range(len(order))); ax.set_yticklabels([short(traits[t][0], t) for t in order], fontsize=6.6, color=INK)
for name, start, n in groups:
    if start: ax.axhline(start - 0.5, color="white", lw=3)
    ax.text(2.58, start + (n - 1) / 2, name.replace(" disorder", "").replace("General ", "Gen. ") if n < 7 else name,
            rotation=90 if n >= 7 else 0, va="center", ha="left", fontsize=6.8 if n < 7 else 7.5, color=MH, fontweight="bold")
ax.set_xticks(np.arange(-0.5, 3, 1), minor=True); ax.grid(which="minor", axis="x", color="white", lw=2)
ax.tick_params(which="both", length=0)
for s in ax.spines.values(): s.set_visible(False)
ax.set_title("A  Trait-level rg", loc="left", fontsize=10, color=INK, fontweight="bold", pad=22)
cb = fig.colorbar(im, ax=ax, orientation="horizontal", fraction=0.02, pad=0.015, aspect=25)
cb.set_label("rg; * FDR < 0.05", fontsize=7.5, color=INK); cb.outline.set_visible(False); cb.ax.tick_params(labelsize=7)
# panel B: subdomain means
names = [g[0] for g in groups]; nsz = {g[0]: g[2] for g in groups}
for k, name in enumerate(names):
    key = [kk for kk, nm in SUB if nm == name][0]
    for j, c in enumerate(cols):
        r = [x for x in t6 if x["pollutant"] == c and x["subdomain"] == key][0]
        yy = k + (j - 1) * 0.24
        bx.plot([f(r["ic_low"]), f(r["ic_high"])], [yy, yy], color=PCOL[c], lw=2, solid_capstyle="round")
        bx.scatter(f(r["mean_rg"]), yy, s=26, color=PCOL[c], edgecolor="white", lw=1, zorder=3,
                   label=POLN[c] if k == 0 else None)
    if k: bx.axhline(k - 0.5, color=GRID, lw=0.8)
bx.set_ylim(len(names) - 0.5, -0.5)
bx.set_yticks(range(len(names))); bx.set_yticklabels([f"{n} (n={nsz[n]})" for n in names], fontsize=8, color=INK)
bx.yaxis.tick_right(); bx.tick_params(axis="y", length=0)
bx.axvline(0, color=MUTED, lw=0.8, ls=(0, (3, 3)))
bx.set_xlim(-0.35, 0.6)
bx.set_xlabel("Mean rg within subdomain (95% CI)", fontsize=8)
bx.grid(axis="x", color=GRID, lw=0.6); bx.set_axisbelow(True)
for sp in ("top", "right", "left"): bx.spines[sp].set_visible(False)
bx.legend(frameon=False, fontsize=7.5, loc="upper left", bbox_to_anchor=(0.0, -0.13), ncol=1)
bx.set_title("B  Subdomain means", loc="left", fontsize=10, color=INK, fontweight="bold", pad=22)
fig.savefig(f"{OUT}/Figure3_rg_heatmap_subdomains.png", dpi=300)
fig.savefig(f"{OUT}/Figure3_rg_heatmap_subdomains.pdf")

# ---------------- Supplementary Figure S1: LCV forest plot ------------------------
import os
if os.path.exists(f"{RES}/T7_lcv.csv"):
    t7 = R("T7_lcv.csv")
    desc = {r["trait2"]: r["description2"] for r in t2}
    mh_traits = sorted({r["trait"] for r in t7}, key=lambda t: -np.mean([f(r["gcp"]) for r in t7 if r["trait"] == t]))
    fig, axs = plt.subplots(1, 3, figsize=(11, 0.28 * len(mh_traits) + 1.4), sharey=True)
    for a, p in zip(axs, cols):
        for i, t in enumerate(mh_traits):
            r = [x for x in t7 if x["pollutant"] == p and x["trait"] == t][0]
            g, se = f(r["gcp"]), f(r["gcp_se"])
            sig_rg = f(r["rg_FDR"]) < 0.05
            a.plot([max(g - 1.96 * se, -1), min(g + 1.96 * se, 1)], [i, i], color=PCOL[p], lw=1.8, solid_capstyle="round", alpha=0.9)
            a.scatter(g, i, s=26, color=PCOL[p] if sig_rg else "white", edgecolor=PCOL[p], lw=1.3, zorder=3)
        a.axvline(0, color=MUTED, lw=0.8, ls=(0, (3, 3)))
        a.set_xlim(-1.05, 1.05); a.set_title(POLN[p], fontsize=10, color=INK, fontweight="bold")
        a.grid(axis="x", color=GRID, lw=0.6); a.set_axisbelow(True)
        for s in ("top", "right"): a.spines[s].set_visible(False)
        a.set_xlabel("GCP (95% CI)", fontsize=8.5)
    axs[0].set_yticks(range(len(mh_traits)))
    axs[0].set_yticklabels([short(desc.get(t, t), t) for t in mh_traits], fontsize=7.2, color=INK)
    axs[0].invert_yaxis(); axs[0].tick_params(axis="y", length=0)
    fig.text(0.99, 0.01, "Filled points: LDSC rg FDR < 0.05; open points: not significant.  GCP > 0: pollution → mental health; GCP < 0: mental health → pollution.",
             ha="right", fontsize=7.5, color=MUTED)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(f"{OUT}/FigureS1_LCV_GCP.png", dpi=300); fig.savefig(f"{OUT}/FigureS1_LCV_GCP.pdf")
