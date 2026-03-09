"""
Step 4: Generate figures and run report.

Outputs (all written to output_dir/figures/):
  - volcano_{species}.png            : per-species volcano plot
  - heatmap_top_conserved.png        : top 20 conserved OGs × species mean log2FC
  - conservation_score_distribution.png
  - direction_agreement.png          : pairwise species agreement on conserved OGs
  - summary_stats.png                : OG funnel chart

Also writes run_report.txt to output_dir.

Called by run_pipeline.py; can also be run standalone:
  python steps/04_report.py config.yaml output/run_dir/
"""

import sys
import os
import hashlib
import datetime
import yaml
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))


# ── figure helpers ────────────────────────────────────────────────────────────

def _save(fig, path: Path):
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {path}")


def plot_volcano(de_df: pd.DataFrame, species_name: str, display_name: str,
                 output_dir: Path, padj_threshold=0.05, log2fc_threshold=1.0):
    """Volcano plot: log2FC vs -log10(padj), colored by sig + direction."""
    df = de_df.copy()
    df["neg_log10_padj"] = -np.log10(df["padj"].clip(lower=1e-300))

    # Categories
    df["category"] = "not significant"
    sig_mask = (df["padj"] < padj_threshold) & (df["log2FC"].abs() >= log2fc_threshold)
    df.loc[sig_mask & (df["log2FC"] > 0), "category"] = "up in hard"
    df.loc[sig_mask & (df["log2FC"] < 0), "category"] = "up in easy"

    color_map = {
        "not significant": "#AAAAAA",
        "up in hard":      "#D62728",
        "up in easy":      "#1F77B4",
    }

    fig, ax = plt.subplots(figsize=(7, 5))
    for cat, color in color_map.items():
        mask = df["category"] == cat
        ax.scatter(
            df.loc[mask, "log2FC"],
            df.loc[mask, "neg_log10_padj"],
            c=color, label=f"{cat} (n={mask.sum()})",
            s=6, alpha=0.5, linewidths=0,
        )

    ax.axvline(x=log2fc_threshold,  color="black", lw=0.8, ls="--")
    ax.axvline(x=-log2fc_threshold, color="black", lw=0.8, ls="--")
    ax.axhline(y=-np.log10(padj_threshold), color="black", lw=0.8, ls="--")

    ax.set_xlabel("log2FC (Hard / Easy)", fontsize=11)
    ax.set_ylabel("-log10(padj)", fontsize=11)
    ax.set_title(f"Volcano Plot — {display_name}", fontsize=12)
    ax.legend(fontsize=8, markerscale=3)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    _save(fig, output_dir / f"volcano_{species_name}.png")


def plot_conservation_distribution(results_df: pd.DataFrame, threshold: float, output_dir: Path):
    """Histogram of conservation scores with threshold line."""
    scores = results_df["conservation_score"].dropna()
    if len(scores) == 0:
        print("  Skipping conservation distribution plot (no data).")
        return

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(scores, bins=30, color="#4878CF", edgecolor="white", linewidth=0.4)
    ax.axvline(x=threshold, color="#D62728", lw=1.5, ls="--",
               label=f"Threshold = {threshold}")
    ax.set_xlabel("Conservation Score", fontsize=11)
    ax.set_ylabel("Number of Orthogroups", fontsize=11)
    ax.set_title("Conservation Score Distribution", fontsize=12)
    ax.legend(fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    _save(fig, output_dir / "conservation_score_distribution.png")


def plot_top_conserved_heatmap(conserved_df: pd.DataFrame, active_species: list,
                                output_dir: Path, top_n: int = 20):
    """Heatmap: rows = top N conserved OGs, columns = species, values = mean log2FC."""
    cols = [f"{sp}_mean_log2fc" for sp in active_species]
    cols_present = [c for c in cols if c in conserved_df.columns]

    if not cols_present:
        print("  Skipping heatmap (no mean log2FC columns found).")
        return

    top = (conserved_df
           .sort_values("conservation_score", ascending=False)
           .head(top_n)
           .set_index("orthogroup_id")[cols_present])

    # Rename columns to species names
    top.columns = [c.replace("_mean_log2fc", "") for c in top.columns]

    fig, ax = plt.subplots(figsize=(max(4, len(active_species) * 1.5), max(6, top_n * 0.35)))
    sns.heatmap(
        top.astype(float),
        ax=ax,
        cmap="RdBu_r",
        center=0,
        linewidths=0.3,
        linecolor="#DDDDDD",
        cbar_kws={"label": "Mean log2FC (sig genes)"},
        annot=len(top) <= 30,
        fmt=".2f",
    )
    ax.set_title(f"Top {top_n} Conserved Orthogroups — Mean log2FC", fontsize=12)
    ax.set_xlabel("Species", fontsize=10)
    ax.set_ylabel("Orthogroup", fontsize=10)
    ax.tick_params(axis="y", labelsize=7)
    _save(fig, output_dir / "heatmap_top_conserved.png")


def plot_direction_agreement(conserved_df: pd.DataFrame, active_species: list, output_dir: Path):
    """Bar chart: pairwise species agreement on conserved OGs."""
    from itertools import combinations

    if len(active_species) < 2:
        print("  Skipping direction agreement (fewer than 2 species).")
        return

    pairs = list(combinations(active_species, 2))
    labels = []
    agreements = []

    for sp1, sp2 in pairs:
        col1 = f"{sp1}_direction"
        col2 = f"{sp2}_direction"
        if col1 not in conserved_df or col2 not in conserved_df:
            continue
        both = conserved_df[[col1, col2]].dropna()
        both = both[(both[col1] != "tie") & (both[col2] != "tie")]
        if len(both) == 0:
            continue
        agree = (both[col1] == both[col2]).mean() * 100
        labels.append(f"{sp1}\nvs\n{sp2}")
        agreements.append(agree)

    if not labels:
        print("  Skipping direction agreement (no valid pairs).")
        return

    fig, ax = plt.subplots(figsize=(max(4, len(labels) * 1.8), 4))
    bars = ax.bar(labels, agreements, color="#4878CF", edgecolor="white")
    ax.axhline(50, color="gray", lw=0.8, ls="--", label="Random (50%)")
    ax.set_ylim(0, 105)
    ax.set_ylabel("% Agreement on Direction", fontsize=11)
    ax.set_title("Direction Agreement Between Species\n(on Conserved Orthogroups)", fontsize=12)
    ax.legend(fontsize=8)
    for bar, val in zip(bars, agreements):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1,
                f"{val:.1f}%", ha="center", va="bottom", fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    _save(fig, output_dir / "direction_agreement.png")


def plot_summary_stats(results_df: pd.DataFrame, all_ogs: int,
                       min_species: int, output_dir: Path):
    """Funnel bar chart: total OGs -> OGs with data -> conserved."""
    ogs_with_data = (results_df["n_species_with_data"] >= min_species).sum()
    conserved     = results_df["is_conserved"].sum()

    labels = [
        f"Total OGs\n({all_ogs})",
        f"OGs with data\nin ≥{min_species} species\n({ogs_with_data})",
        f"Conserved OGs\n({conserved})",
    ]
    values = [all_ogs, ogs_with_data, conserved]
    colors = ["#AEC6CF", "#4878CF", "#D62728"]

    fig, ax = plt.subplots(figsize=(6, 4))
    bars = ax.barh(labels[::-1], values[::-1], color=colors[::-1], edgecolor="white")
    ax.set_xlabel("Number of Orthogroups", fontsize=11)
    ax.set_title("Orthogroup Filtering Summary", fontsize=12)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for bar, val in zip(bars, values[::-1]):
        ax.text(bar.get_width() + all_ogs * 0.01, bar.get_y() + bar.get_height() / 2,
                str(val), va="center", fontsize=10)
    _save(fig, output_dir / "summary_stats.png")


# ── report ────────────────────────────────────────────────────────────────────

def _file_checksum(path: Path) -> str:
    if not path.exists():
        return "MISSING"
    md5 = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            md5.update(chunk)
    return md5.hexdigest()[:8]


def write_run_report(config: dict, output_dir: Path, results_df: pd.DataFrame,
                     de_stats: dict):
    """Write a plain-text run report."""
    min_species = config["pipeline"]["min_species_for_conservation"]
    threshold   = config["pipeline"]["conservation_threshold"]

    ogs_with_data = (results_df["n_species_with_data"] >= min_species).sum() if len(results_df) else 0
    conserved     = results_df["is_conserved"].sum() if len(results_df) else 0
    total_ogs     = len(results_df)

    report_lines = [
        "=" * 60,
        "RNA-seq Cross-Species Rooting Analysis — Run Report",
        "=" * 60,
        f"Run date     : {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"Run name     : {config['pipeline']['run_name']}",
        f"Output dir   : {output_dir}",
        "",
        "--- Config Snapshot ---",
        f"  OrthoFinder results : {config['pipeline']['orthofinder_results']}",
        f"  Min species         : {min_species}",
        f"  Conservation thresh : {threshold}",
        f"  padj threshold      : {config['pipeline'].get('padj_threshold', 0.05)}",
        f"  |log2FC| threshold  : {config['pipeline'].get('log2fc_threshold', 1.0)}",
        "",
        "--- Species ---",
    ]

    for sp in config["species"]:
        report_lines.append(f"  {sp['name']}: {sp['display_name']}")
        report_lines.append(f"    expression file: {sp['expression_file']}")
        report_lines.append(f"    easy samples   : {sp['samples']['easy']}")
        report_lines.append(f"    hard samples   : {sp['samples']['hard']}")

    report_lines += [
        "",
        "--- DE Statistics (per species) ---",
    ]
    for sp_name, stats in de_stats.items():
        report_lines.append(f"  {sp_name}:")
        report_lines.append(f"    n tested   : {stats.get('n_tested', 'N/A')}")
        padj_t = config['pipeline'].get('padj_threshold', 0.05)
        log2fc_t = config['pipeline'].get('log2fc_threshold', 1.0)
        report_lines.append(f"    n sig      : {stats.get('n_sig', 'N/A')}  (padj<{padj_t}, |log2FC|>={log2fc_t})")
        report_lines.append(f"    method     : {stats.get('method', 'N/A')}")
        report_lines.append(f"    mapping pct: {stats.get('mapping_pct', 'N/A')}")

    report_lines += [
        "",
        "--- Orthogroup Conservation ---",
        f"  Total OGs evaluated       : {total_ogs}",
        f"  OGs with data (>={min_species} sp)  : {ogs_with_data}",
        f"  Conserved OGs             : {conserved}",
        "",
        "--- Output Files ---",
    ]

    output_files = [
        "orthogroup_gene_map.csv",
        "orthogroup_species_counts.csv",
        "cross_species_results.csv",
        "conserved_orthogroups.csv",
        "per_gene_de_combined.csv",
    ]
    for fname in output_files:
        p = output_dir / fname
        report_lines.append(f"  {fname:<40} md5:{_file_checksum(p)}")

    report_lines += ["", "=" * 60]

    report_path = output_dir / "run_report.txt"
    with open(report_path, "w") as f:
        f.write("\n".join(report_lines) + "\n")
    print(f"  Written: {report_path}")
    return report_path


# ── main ──────────────────────────────────────────────────────────────────────

def run_report(config: dict, output_dir: str, config_dir: str):
    """Main entry point for Step 4."""
    out_path    = Path(output_dir)
    figures_dir = out_path / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    active_species = config["species"]
    species_names  = [sp["name"] for sp in active_species]
    min_species    = config["pipeline"]["min_species_for_conservation"]
    threshold      = config["pipeline"]["conservation_threshold"]

    print("\n[Step 4] Generating figures and report...")

    # --- Load data ---
    cross_file = out_path / "cross_species_results.csv"
    cons_file  = out_path / "conserved_orthogroups.csv"

    results_df  = pd.read_csv(cross_file) if cross_file.exists() else pd.DataFrame()
    conserved_df = pd.read_csv(cons_file)  if cons_file.exists()  else pd.DataFrame()

    de_stats = {}

    # --- Read thresholds from config ---
    padj_threshold  = config["pipeline"].get("padj_threshold", 0.05)
    log2fc_threshold = config["pipeline"].get("log2fc_threshold", 1.0)

    # --- Per-species volcano plots ---
    for sp in active_species:
        de_file = Path(config_dir) / sp["de_output"]
        if not de_file.exists():
            print(f"  Skipping volcano for {sp['name']} (DE file missing)")
            continue
        de_df = pd.read_csv(de_file)
        plot_volcano(de_df, sp["name"], sp["display_name"], figures_dir,
                     padj_threshold=padj_threshold, log2fc_threshold=log2fc_threshold)

        n_tested = len(de_df)
        n_sig    = ((de_df["padj"] < padj_threshold) & (de_df["log2FC"].abs() >= log2fc_threshold)).sum()
        method   = de_df["de_method"].iloc[0] if "de_method" in de_df.columns else "unknown"
        de_stats[sp["name"]] = {
            "n_tested": n_tested,
            "n_sig": n_sig,
            "method": method,
            "mapping_pct": "see orthogroup_gene_map.csv",
        }

    # --- Conservation score distribution ---
    if not results_df.empty:
        plot_conservation_distribution(results_df, threshold, figures_dir)

    # --- Top conserved heatmap ---
    if not conserved_df.empty:
        plot_top_conserved_heatmap(conserved_df, species_names, figures_dir, top_n=20)
        plot_direction_agreement(conserved_df, species_names, figures_dir)

    # --- Summary stats ---
    if not results_df.empty:
        og_gene_map_file = out_path / "orthogroup_gene_map.csv"
        all_ogs = len(results_df)
        plot_summary_stats(results_df, all_ogs, min_species, figures_dir)

    # --- Run report ---
    write_run_report(config, out_path, results_df, de_stats)

    print("[Step 4] Complete.")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python 04_report.py config.yaml output_dir/")
        sys.exit(1)

    config_file = sys.argv[1]
    output_dir  = sys.argv[2]

    with open(config_file) as f:
        config = yaml.safe_load(f)

    config_dir = str(Path(config_file).parent)
    run_report(config, output_dir, config_dir)
