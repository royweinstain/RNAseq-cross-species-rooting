"""
Step 2: Parse OrthoFinder output into a gene -> orthogroup map.

Supports two input formats:
  A) Custom 3-column TSV (orthofinder_custom_file in config):
       Orthogroup  arab_genes  argan_genes
       OG0000001   AT1G...     XLOC_042640_1, ...
     Used when only a subset-species OrthoFinder run is available.

  B) Standard OrthoFinder output dir (orthofinder_results in config):
       Reads Orthogroups/Orthogroups.txt  (space-separated per-OG gene lists)
       or  Orthogroups/Orthogroups.tsv   (species-as-columns matrix)

Writes to output_dir:
  - orthogroup_gene_map.csv      : gene_id, normalized_gene_id, orthogroup_id, species
  - orthogroup_species_counts.csv: orthogroup_id, <species>_count, ..., total_species

Usage: python steps/02_parse_orthogroups.py config.yaml output/run_dir/
"""

import sys
import re
import csv
import yaml
import pandas as pd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from utils import normalize_gene_id, detect_species_from_gene_id, build_prefix_map


# ── Format A: custom 3-column TSV ────────────────────────────────────────────

def parse_custom_tsv(filepath: Path, config: dict, prefix_map: dict) -> tuple[dict, dict]:
    """
    Parse the 3-column TSV produced for argan-vs-arabidopsis OrthoFinder run.

    Format:
        Orthogroup\tarab_genes\targan_genes
        OG0000001\tAT1G..., AT2G...\tXLOC_042640_1, XLOC_047655_1

    The column names used for each species are taken from config:
        species[i].orthofinder_column   (e.g. "argan_genes", "poplar_genes")
        reference_species[i].orthofinder_column  (e.g. "arab_genes")

    Returns:
        gene_map    : {normalized_gene_id -> orthogroup_id}
        og_members  : {orthogroup_id -> {species_name: [gene_ids]}}
    """
    # Build orthofinder_column -> species_name lookup
    col_to_species = {}
    for sp in config.get("species", []):
        col = sp.get("orthofinder_column")
        if col:
            col_to_species[col] = sp["name"]
    for sp in config.get("reference_species", []):
        col = sp.get("orthofinder_column")
        if col:
            col_to_species[col] = sp["name"]

    gene_map: dict[str, str] = {}
    og_members: dict[str, dict] = {}

    with open(filepath, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            # Try "Orthogroup" first, then empty-string key (Orthogroups.csv
            # from OrthoFinder has an unnamed first column)
            og_id = row.get("Orthogroup", "").strip()
            if not og_id:
                og_id = row.get("", "").strip()
            if not og_id:
                continue

            og_entry: dict[str, list] = {}

            for col, species_name in col_to_species.items():
                cell = row.get(col, "").strip()
                if not cell:
                    continue

                genes = [g.strip() for g in re.split(r",\s*", cell) if g.strip()]
                for gene in genes:
                    norm = normalize_gene_id(gene)
                    gene_map[norm] = og_id
                    gene_map[gene] = og_id   # keep raw too
                    if species_name not in og_entry:
                        og_entry[species_name] = []
                    og_entry[species_name].append(gene)

            if og_entry:
                if og_id in og_members:
                    for sp_name, genes in og_entry.items():
                        og_members[og_id].setdefault(sp_name, []).extend(genes)
                else:
                    og_members[og_id] = og_entry

    return gene_map, og_members


# ── Format B: standard OrthoFinder output ────────────────────────────────────

def find_orthogroups_file(results_dir: str) -> Path:
    """Locate Orthogroups.txt or Orthogroups.tsv inside an OrthoFinder results dir."""
    results_path = Path(results_dir)
    candidates = [
        results_path / "Orthogroups" / "Orthogroups.txt",
        results_path / "Orthogroups" / "Orthogroups.tsv",
        results_path / "Orthogroups.txt",
        results_path / "Orthogroups.tsv",
    ]
    for c in candidates:
        if c.exists():
            return c
    raise FileNotFoundError(
        f"Could not find Orthogroups.txt or Orthogroups.tsv in {results_dir}\n"
        f"Looked in: {[str(c) for c in candidates]}"
    )


def parse_orthogroups_txt(filepath: Path, prefix_map: dict) -> tuple[dict, dict]:
    """Parse OrthoFinder Orthogroups.txt (one OG per line: 'OGxxxxxxx: gene1 gene2 ...')."""
    gene_map: dict[str, str] = {}
    og_members: dict[str, dict] = {}

    with open(filepath, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or ":" not in line:
                continue
            og_id, rest = line.split(":", 1)
            og_id = og_id.strip()

            og_entry: dict[str, list] = {}
            for gene in rest.split():
                norm = normalize_gene_id(gene)
                species = (detect_species_from_gene_id(norm, prefix_map) or
                           detect_species_from_gene_id(gene, prefix_map))
                if species is None:
                    continue
                gene_map[norm] = og_id
                gene_map[gene] = og_id
                og_entry.setdefault(species, []).append(gene)

            if og_entry:
                og_members[og_id] = og_entry

    return gene_map, og_members


def parse_orthogroups_tsv(filepath: Path, config: dict, prefix_map: dict) -> tuple[dict, dict]:
    """Parse OrthoFinder Orthogroups.tsv (species as columns)."""
    col_to_species = {}
    for sp in config.get("species", []):
        col_to_species[sp["orthofinder_column"]] = sp["name"]
    for sp in config.get("reference_species", []):
        col_to_species[sp["orthofinder_column"]] = sp["name"]

    gene_map: dict[str, str] = {}
    og_members: dict[str, dict] = {}

    df = pd.read_csv(filepath, sep="\t", dtype=str).fillna("")
    og_col = df.columns[0]

    for _, row in df.iterrows():
        og_id = row[og_col].strip()
        if not og_id:
            continue

        og_entry: dict[str, list] = {}
        for col in df.columns[1:]:
            species_name = col_to_species.get(col)
            cell = row[col].strip()
            if not cell:
                continue
            genes = [g.strip() for g in re.split(r"[,\s]+", cell) if g.strip()]
            for gene in genes:
                norm = normalize_gene_id(gene)
                if species_name is None:
                    species_name = detect_species_from_gene_id(norm, prefix_map)
                if species_name is None:
                    continue
                gene_map[norm] = og_id
                gene_map[gene] = og_id
                og_entry.setdefault(species_name, []).append(gene)

        if og_entry:
            og_members[og_id] = og_entry

    return gene_map, og_members


# ── Output helpers ────────────────────────────────────────────────────────────

def build_gene_map_df(og_members: dict) -> pd.DataFrame:
    rows = []
    seen = set()
    for og_id, species_dict in og_members.items():
        for species_name, genes in species_dict.items():
            for gene in genes:
                norm = normalize_gene_id(gene)
                key = (norm, og_id)
                if key in seen:
                    continue
                seen.add(key)
                rows.append({
                    "gene_id": gene,
                    "normalized_gene_id": norm,
                    "orthogroup_id": og_id,
                    "species": species_name,
                })
    return pd.DataFrame(rows)


def build_species_counts_df(og_members: dict, active_species: list) -> pd.DataFrame:
    rows = []
    for og_id, species_dict in og_members.items():
        row = {"orthogroup_id": og_id}
        n_sp = 0
        for sp in active_species:
            count = len(species_dict.get(sp, []))
            row[f"{sp}_count"] = count
            if count > 0:
                n_sp += 1
        row["total_species_with_data"] = n_sp
        rows.append(row)
    return pd.DataFrame(rows)


# ── Main entry point ──────────────────────────────────────────────────────────

def parse_orthogroups(config: dict, output_dir: str, config_dir: str) -> tuple[dict, dict]:
    """
    Main entry point for Step 2. Detects format from config and parses.
    Returns (gene_map, og_members); writes CSV summaries to output_dir.
    """
    prefix_map     = build_prefix_map(config)
    active_species = [sp["name"] for sp in config.get("species", [])]
    custom_file    = config["pipeline"].get("orthofinder_custom_file")
    og_results     = config["pipeline"].get("orthofinder_results")

    if custom_file:
        og_file = Path(config_dir) / custom_file
        if not og_file.exists():
            raise FileNotFoundError(f"orthofinder_custom_file not found: {og_file}")
        print(f"[Step 2] Parsing custom OrthoFinder TSV: {og_file}")
        gene_map, og_members = parse_custom_tsv(og_file, config, prefix_map)
    elif og_results:
        og_dir  = Path(config_dir) / og_results
        og_file = find_orthogroups_file(str(og_dir))
        print(f"[Step 2] Parsing OrthoFinder file: {og_file}")
        if og_file.suffix == ".tsv":
            gene_map, og_members = parse_orthogroups_tsv(og_file, config, prefix_map)
        else:
            gene_map, og_members = parse_orthogroups_txt(og_file, prefix_map)
    else:
        raise ValueError(
            "Config must specify either pipeline.orthofinder_custom_file "
            "or pipeline.orthofinder_results."
        )

    print(f"  Total orthogroups parsed: {len(og_members)}")
    print(f"  Total gene->OG mappings : {len(gene_map)}")

    # Write outputs
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    gene_map_df      = build_gene_map_df(og_members)
    species_counts_df = build_species_counts_df(og_members, active_species)

    gene_map_file = out_path / "orthogroup_gene_map.csv"
    counts_file   = out_path / "orthogroup_species_counts.csv"
    gene_map_df.to_csv(gene_map_file, index=False)
    species_counts_df.to_csv(counts_file, index=False)

    print(f"  Written: {gene_map_file}  ({len(gene_map_df)} rows)")
    print(f"  Written: {counts_file}  ({len(species_counts_df)} rows)")

    # Per-species coverage
    for sp in active_species:
        n = (gene_map_df["species"] == sp).sum()
        print(f"  {sp}: {n} genes in orthogroups")

    return gene_map, og_members


# ── Standalone ────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python 02_parse_orthogroups.py config.yaml output_dir/")
        sys.exit(1)
    with open(sys.argv[1]) as f:
        config = yaml.safe_load(f)
    config_dir = str(Path(sys.argv[1]).parent)
    parse_orthogroups(config, sys.argv[2], config_dir)
