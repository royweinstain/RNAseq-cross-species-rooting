"""
Step 3: Cross-species aggregation of DE results via orthogroup membership.

For each orthogroup with expression data in >= min_species_for_conservation species:
  - Per species: flag significant genes (padj < threshold, |log2FC| >= fc_threshold)
  - Compute modal direction (up_in_hard / up_in_easy) per species
  - Majority vote across species -> conservation_score = fraction agreeing
  - is_conserved = score >= conservation_threshold AND >= min_species

Writes:
  - cross_species_results.csv   : all OGs with >=min_species data
  - conserved_orthogroups.csv   : subset where is_conserved=True
  - per_gene_de_combined.csv    : all DE results joined with OG membership

Called by run_pipeline.py; can also be run standalone:
  python steps/03_cross_species.py config.yaml output/run_dir/
"""

import sys
import yaml
import pandas as pd
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from utils import normalize_gene_id


def load_and_map_de_results(config: dict, gene_map: dict, config_dir: str) -> pd.DataFrame:
    """
    Load all species' de_results.csv files, normalize gene IDs,
    join with gene_map to add orthogroup_id.

    Returns a combined DataFrame with columns:
      gene_id, normalized_gene_id, log2FC, pvalue, padj, mean_expr,
      species, de_method, orthogroup_id, direction, is_significant
    """
    frames = []
    for sp in config["species"]:
        de_file = Path(config_dir) / sp["de_output"]
        if not de_file.exists():
            raise FileNotFoundError(
                f"DE results not found for {sp['name']}: {de_file}\n"
                f"Run Step 1 first."
            )
        df = pd.read_csv(de_file)
        df["species"] = sp["name"]
        frames.append(df)

    if not frames:
        raise ValueError("No DE result files loaded.")

    combined = pd.concat(frames, ignore_index=True)

    # Normalize gene IDs for mapping
    combined["normalized_gene_id"] = combined["gene_id"].apply(normalize_gene_id)

    # Map to orthogroups (try normalized first, then raw)
    def map_gene(row):
        og = gene_map.get(row["normalized_gene_id"])
        if og is None:
            og = gene_map.get(row["gene_id"])
        return og

    combined["orthogroup_id"] = combined.apply(map_gene, axis=1)

    # Add direction column
    combined["direction"] = np.where(combined["log2FC"] > 0, "up_in_hard", "up_in_easy")

    return combined


def compute_conservation_scores(
    combined_df: pd.DataFrame,
    config: dict,
    padj_threshold: float = 0.05,
    log2fc_threshold: float = 1.0,
) -> pd.DataFrame:
    """
    Compute per-orthogroup conservation scores across species.

    For each OG:
      1. Per species: keep genes with padj < padj_threshold AND |log2FC| >= log2fc_threshold
      2. Modal direction = whichever direction has more sig genes (ties -> skip species)
      3. conservation_score = fraction of species (with >=1 sig gene) agreeing on modal direction
      4. is_conserved = score >= conservation_threshold AND n_species_with_data >= min_species

    Returns DataFrame with one row per orthogroup.
    """
    min_species     = config["pipeline"]["min_species_for_conservation"]
    cons_threshold  = config["pipeline"]["conservation_threshold"]
    active_species  = [sp["name"] for sp in config["species"]]

    # Only work with genes that have orthogroup assignments
    df = combined_df.dropna(subset=["orthogroup_id"]).copy()

    # Flag significance
    df["is_significant"] = (
        (df["padj"] < padj_threshold) &
        (df["log2FC"].abs() >= log2fc_threshold)
    )

    rows = []

    for og_id, og_df in df.groupby("orthogroup_id"):
        row = {"orthogroup_id": og_id}
        species_directions = []  # (species, modal_direction) for species with sig genes

        for sp in active_species:
            sp_df     = og_df[og_df["species"] == sp]
            n_total   = len(sp_df)
            sig_df    = sp_df[sp_df["is_significant"]]
            n_sig     = len(sig_df)

            row[f"{sp}_n_total"] = n_total
            row[f"{sp}_n_sig"]   = n_sig

            if n_sig == 0:
                row[f"{sp}_direction"]   = None
                row[f"{sp}_consistency"] = None
                continue

            n_up_hard = (sig_df["direction"] == "up_in_hard").sum()
            n_up_easy = (sig_df["direction"] == "up_in_easy").sum()

            if n_up_hard > n_up_easy:
                modal_dir   = "up_in_hard"
                consistency = n_up_hard / n_sig
            elif n_up_easy > n_up_hard:
                modal_dir   = "up_in_easy"
                consistency = n_up_easy / n_sig
            else:
                # Exact tie — skip this species for conservation vote
                row[f"{sp}_direction"]   = "tie"
                row[f"{sp}_consistency"] = 0.5
                continue

            row[f"{sp}_direction"]   = modal_dir
            row[f"{sp}_consistency"] = consistency
            species_directions.append(modal_dir)

        n_species_with_data = len(species_directions)
        row["n_species_with_data"] = n_species_with_data

        if n_species_with_data < min_species:
            row["conservation_score"] = None
            row["modal_direction"]    = None
            row["is_conserved"]       = False
            rows.append(row)
            continue

        # Majority vote
        n_up_hard_species = sum(1 for d in species_directions if d == "up_in_hard")
        n_up_easy_species = sum(1 for d in species_directions if d == "up_in_easy")

        if n_up_hard_species > n_up_easy_species:
            modal_dir = "up_in_hard"
            n_agree   = n_up_hard_species
        elif n_up_easy_species > n_up_hard_species:
            modal_dir = "up_in_easy"
            n_agree   = n_up_easy_species
        else:
            # Exact tie across species — mark as tied, not conserved
            row["conservation_score"] = 0.5
            row["modal_direction"]    = "tie"
            row["is_conserved"]       = False
            rows.append(row)
            continue

        score = n_agree / n_species_with_data
        row["conservation_score"] = score
        row["modal_direction"]    = modal_dir
        row["is_conserved"]       = (score >= cons_threshold and
                                     n_species_with_data >= min_species)
        rows.append(row)

    result_df = pd.DataFrame(rows)

    # Add per-species mean log2FC of significant genes for heatmap use
    for sp in active_species:
        sp_sig = df[(df["species"] == sp) & df["is_significant"]]
        mean_fc = sp_sig.groupby("orthogroup_id")["log2FC"].mean().rename(f"{sp}_mean_log2fc")
        result_df = result_df.merge(mean_fc, on="orthogroup_id", how="left")

    return result_df


def run_cross_species(config: dict, gene_map: dict, output_dir: str, config_dir: str) -> pd.DataFrame:
    """Main entry point for Step 3."""
    print("\n[Step 3] Loading DE results and mapping to orthogroups...")
    combined = load_and_map_de_results(config, gene_map, config_dir)

    total_genes    = len(combined)
    mapped_genes   = combined["orthogroup_id"].notna().sum()
    mapping_rate   = mapped_genes / total_genes * 100 if total_genes > 0 else 0

    print(f"  Total genes (all species): {total_genes}")
    print(f"  Mapped to orthogroups    : {mapped_genes} ({mapping_rate:.1f}%)")

    for sp in config["species"]:
        sp_df   = combined[combined["species"] == sp["name"]]
        sp_map  = sp_df["orthogroup_id"].notna().sum()
        sp_rate = sp_map / len(sp_df) * 100 if len(sp_df) > 0 else 0
        print(f"  {sp['name']}: {sp_map}/{len(sp_df)} genes mapped ({sp_rate:.1f}%)")

    print("\n[Step 3] Computing conservation scores...")
    padj_thresh   = config["pipeline"].get("padj_threshold", 0.05)
    log2fc_thresh = config["pipeline"].get("log2fc_threshold", 1.0)
    results_df = compute_conservation_scores(combined, config, padj_thresh, log2fc_thresh)

    min_sp    = config["pipeline"]["min_species_for_conservation"]
    threshold = config["pipeline"]["conservation_threshold"]

    ogs_with_data = results_df[results_df["n_species_with_data"] >= min_sp]
    conserved     = results_df[results_df["is_conserved"] == True]

    print(f"  Total OGs evaluated           : {len(results_df)}")
    print(f"  OGs with data in >= {min_sp} species : {len(ogs_with_data)}")
    print(f"  Conserved OGs (score >= {threshold}): {len(conserved)}")

    if len(conserved) > 0:
        up_hard = (conserved["modal_direction"] == "up_in_hard").sum()
        up_easy = (conserved["modal_direction"] == "up_in_easy").sum()
        print(f"    up_in_hard: {up_hard}  |  up_in_easy: {up_easy}")

    # Write outputs
    out_path = Path(output_dir)

    # All OGs with >= min_species data
    cross_file = out_path / "cross_species_results.csv"
    ogs_with_data.sort_values("conservation_score", ascending=False, na_position="last").to_csv(
        cross_file, index=False
    )

    # Conserved only
    cons_file = out_path / "conserved_orthogroups.csv"
    conserved.sort_values("conservation_score", ascending=False).to_csv(cons_file, index=False)

    # Per-gene combined with OG membership
    per_gene_file = out_path / "per_gene_de_combined.csv"
    combined.to_csv(per_gene_file, index=False)

    print(f"\n  Written: {cross_file}")
    print(f"  Written: {cons_file}  ({len(conserved)} rows)")
    print(f"  Written: {per_gene_file}")

    return results_df


# ── standalone entry point ────────────────────────────────────────────────────

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python 03_cross_species.py config.yaml output_dir/")
        sys.exit(1)

    config_file = sys.argv[1]
    output_dir  = sys.argv[2]

    with open(config_file) as f:
        config = yaml.safe_load(f)

    # Load gene map from step 2 output
    gene_map_file = Path(output_dir) / "orthogroup_gene_map.csv"
    if not gene_map_file.exists():
        print(f"ERROR: orthogroup_gene_map.csv not found. Run Step 2 first.")
        sys.exit(1)

    gene_map_df = pd.read_csv(gene_map_file)
    gene_map = {}
    for _, row in gene_map_df.iterrows():
        gene_map[row["gene_id"]]            = row["orthogroup_id"]
        gene_map[row["normalized_gene_id"]] = row["orthogroup_id"]

    config_dir = str(Path(config_file).parent)
    run_cross_species(config, gene_map, output_dir, config_dir)
