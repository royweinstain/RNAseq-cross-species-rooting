"""
Step 5: Annotate conserved orthogroups with Arabidopsis gene information.

For each conserved OG:
  - Extracts Arabidopsis gene IDs from orthogroup_gene_map.csv
  - Batch-queries MyGene.info API for gene symbols, names, GO terms
  - Writes conserved_genes_annotated.csv

Called by run_pipeline.py; can also be run standalone:
  python steps/05_annotate.py config.yaml output/run_dir/
"""

import sys
import csv
import json
import time
import urllib.request
import urllib.error
import pandas as pd
import numpy as np
from pathlib import Path


def query_mygene_batch(gene_ids: list, batch_size: int = 100) -> dict:
    """
    Query MyGene.info API for Arabidopsis gene annotations.

    Args:
        gene_ids: List of Arabidopsis gene IDs (e.g. AT1G12210.1)
        batch_size: Number of genes per API request

    Returns:
        dict mapping gene_id -> {symbol, name, go_bp, go_mf}
    """
    # Deduplicate and strip transcript suffixes (AT1G12210.1 -> AT1G12210)
    locus_map = {}  # locus_id -> set of original IDs
    for gid in gene_ids:
        locus = gid.split(".")[0] if "." in gid else gid
        locus_map.setdefault(locus, set()).add(gid)

    unique_loci = list(locus_map.keys())
    results = {}

    for i in range(0, len(unique_loci), batch_size):
        batch = unique_loci[i:i + batch_size]
        query = ",".join(batch)

        url = "https://mygene.info/v3/query"
        data = urllib.parse.urlencode({
            "q": query,
            "scopes": "symbol,entrezgene,ensembl.gene",
            "fields": "symbol,name,go.BP.term,go.MF.term",
            "species": "thale-cress",
            "size": len(batch),
        }).encode("utf-8")

        try:
            req = urllib.request.Request(url, data=data, method="POST")
            req.add_header("Content-Type", "application/x-www-form-urlencoded")
            with urllib.request.urlopen(req, timeout=30) as resp:
                response_data = json.loads(resp.read().decode("utf-8"))
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as e:
            print(f"  Warning: MyGene.info API error for batch {i//batch_size + 1}: {e}")
            continue

        # Parse response
        if isinstance(response_data, list):
            for hit in response_data:
                if not isinstance(hit, dict) or hit.get("notfound"):
                    continue
                qid = hit.get("query", "")
                symbol = hit.get("symbol", "")
                name = hit.get("name", "")

                # Extract GO terms
                go_bp = _extract_go_terms(hit.get("go", {}), "BP")
                go_mf = _extract_go_terms(hit.get("go", {}), "MF")

                # Map back to all original gene IDs for this locus
                for orig_id in locus_map.get(qid, [qid]):
                    results[orig_id] = {
                        "symbol": symbol,
                        "name": name,
                        "go_bp": go_bp,
                        "go_mf": go_mf,
                    }

        # Rate limit
        if i + batch_size < len(unique_loci):
            time.sleep(0.5)

    return results


def _extract_go_terms(go_dict: dict, category: str) -> str:
    """Extract GO term descriptions from MyGene.info go field."""
    if not go_dict:
        return ""
    terms = go_dict.get(category, [])
    if isinstance(terms, dict):
        terms = [terms]
    if not isinstance(terms, list):
        return ""
    descriptions = []
    for t in terms:
        if isinstance(t, dict) and "term" in t:
            descriptions.append(t["term"])
    # Deduplicate while preserving order
    seen = set()
    unique = []
    for d in descriptions:
        if d not in seen:
            seen.add(d)
            unique.append(d)
    return "; ".join(unique[:10])  # Limit to 10 terms


def run_annotate(config: dict, output_dir: str, config_dir: str):
    """Main entry point for Step 5."""
    out_path = Path(output_dir)

    # Load conserved orthogroups
    cons_file = out_path / "conserved_orthogroups.csv"
    if not cons_file.exists():
        raise FileNotFoundError(
            f"conserved_orthogroups.csv not found in {output_dir}\n"
            f"Run Step 3 first."
        )
    conserved_df = pd.read_csv(cons_file)
    print(f"\n[Step 5] Annotating {len(conserved_df)} conserved orthogroups...")

    # Load gene map to find arabidopsis genes per OG
    gene_map_file = out_path / "orthogroup_gene_map.csv"
    if not gene_map_file.exists():
        raise FileNotFoundError(
            f"orthogroup_gene_map.csv not found in {output_dir}\n"
            f"Run Step 2 first."
        )
    gene_map_df = pd.read_csv(gene_map_file)

    # Filter to arabidopsis genes in conserved OGs
    conserved_og_ids = set(conserved_df["orthogroup_id"])
    arab_genes = gene_map_df[
        (gene_map_df["species"] == "arabidopsis") &
        (gene_map_df["orthogroup_id"].isin(conserved_og_ids))
    ]

    # Build OG -> arab gene IDs mapping
    og_to_arab = {}
    for _, row in arab_genes.iterrows():
        og_id = row["orthogroup_id"]
        gene_id = row["gene_id"]
        og_to_arab.setdefault(og_id, []).append(gene_id)

    n_ogs_with_arab = len(og_to_arab)
    all_arab_ids = list(set(g for genes in og_to_arab.values() for g in genes))
    print(f"  {n_ogs_with_arab}/{len(conserved_og_ids)} conserved OGs have Arabidopsis orthologs")
    print(f"  {len(all_arab_ids)} unique Arabidopsis gene IDs to query")

    # Query MyGene.info
    if all_arab_ids:
        print("  Querying MyGene.info for annotations...")
        annotations = query_mygene_batch(all_arab_ids)
        print(f"  Received annotations for {len(annotations)} genes")
    else:
        annotations = {}

    # Build annotation columns for conserved_df
    arab_gene_ids_col = []
    arab_symbols_col = []
    arab_names_col = []
    go_bp_col = []
    go_mf_col = []

    for _, row in conserved_df.iterrows():
        og_id = row["orthogroup_id"]
        arab_ids = og_to_arab.get(og_id, [])

        if not arab_ids:
            arab_gene_ids_col.append("")
            arab_symbols_col.append("")
            arab_names_col.append("")
            go_bp_col.append("")
            go_mf_col.append("")
            continue

        arab_gene_ids_col.append(", ".join(arab_ids))

        symbols = []
        names = []
        bp_terms = set()
        mf_terms = set()

        for gid in arab_ids:
            ann = annotations.get(gid, {})
            if ann.get("symbol"):
                symbols.append(ann["symbol"])
            if ann.get("name"):
                names.append(ann["name"])
            if ann.get("go_bp"):
                for t in ann["go_bp"].split("; "):
                    bp_terms.add(t)
            if ann.get("go_mf"):
                for t in ann["go_mf"].split("; "):
                    mf_terms.add(t)

        # Deduplicate
        arab_symbols_col.append(", ".join(sorted(set(symbols))))
        arab_names_col.append(", ".join(sorted(set(names))))
        go_bp_col.append("; ".join(sorted(bp_terms)[:10]))
        go_mf_col.append("; ".join(sorted(mf_terms)[:10]))

    conserved_df["arab_gene_ids"] = arab_gene_ids_col
    conserved_df["arab_symbols"] = arab_symbols_col
    conserved_df["arab_names"] = arab_names_col
    conserved_df["go_bp"] = go_bp_col
    conserved_df["go_mf"] = go_mf_col

    # Write output
    annotated_file = out_path / "conserved_genes_annotated.csv"
    conserved_df.to_csv(annotated_file, index=False)
    print(f"  Written: {annotated_file}  ({len(conserved_df)} rows)")

    # Summary
    n_annotated = sum(1 for s in arab_symbols_col if s)
    print(f"  {n_annotated}/{len(conserved_df)} OGs have Arabidopsis annotations")
    print("[Step 5] Complete.")


# ── standalone entry point ────────────────────────────────────────────────────

if __name__ == "__main__":
    import yaml

    if len(sys.argv) < 3:
        print("Usage: python 05_annotate.py config.yaml output_dir/")
        sys.exit(1)

    config_file = sys.argv[1]
    output_dir = sys.argv[2]

    with open(config_file) as f:
        config = yaml.safe_load(f)

    config_dir = str(Path(config_file).parent)
    run_annotate(config, output_dir, config_dir)
