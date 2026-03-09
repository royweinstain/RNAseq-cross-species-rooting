#!/usr/bin/env python3
"""
RNA-seq Cross-Species Rooting Analysis Pipeline — Orchestrator

Usage:
  python run_pipeline.py config.yaml
  python run_pipeline.py config.yaml --dry-run
  python run_pipeline.py config.yaml --steps 2,3,4
  python run_pipeline.py config.yaml --species argan
  python run_pipeline.py config.yaml --steps 1 --species argan

Steps:
  1  Per-species DE analysis (R: DESeq2 or limma-voom)
  2  Parse OrthoFinder orthogroups
  3  Cross-species conservation scoring
  4  Figures and run report
  5  Arabidopsis annotation of conserved OGs (MyGene.info)

All output is written to output/<YYYY-MM-DD_HHMMSS>/ (or --output-dir).
"""

import sys
import argparse
import datetime
import subprocess
import importlib.util
import yaml
from pathlib import Path

PIPELINE_DIR = Path(__file__).parent
STEPS_DIR    = PIPELINE_DIR / "steps"


def _load_step(filename: str):
    """Load a step module from steps/ by filename (handles numeric prefixes)."""
    module_path = STEPS_DIR / filename
    spec   = importlib.util.spec_from_file_location(filename.replace(".py", ""), module_path)
    module = importlib.util.module_from_spec(spec)
    # Make steps/ importable inside the module (for `from utils import ...`)
    if str(STEPS_DIR) not in sys.path:
        sys.path.insert(0, str(STEPS_DIR))
    spec.loader.exec_module(module)
    return module


# ── validation ────────────────────────────────────────────────────────────────

REQUIRED_PIPELINE_KEYS = [
    "run_name", "orthofinder_results", "output_dir",
    "min_species_for_conservation", "conservation_threshold",
]
REQUIRED_SPECIES_KEYS = [
    "name", "display_name", "expression_file",
    "gene_id_column", "gene_id_prefixes", "orthofinder_column",
    "samples", "de_output",
]


def validate_config(config: dict, config_dir: str) -> list:
    errors = []

    if "pipeline" not in config:
        errors.append("Missing top-level 'pipeline' key in config.")
        return errors

    for key in REQUIRED_PIPELINE_KEYS:
        if key not in config["pipeline"]:
            errors.append(f"pipeline.{key} is missing.")

    if "species" not in config or not config["species"]:
        errors.append("No species defined under 'species'.")
        return errors

    for sp in config["species"]:
        name = sp.get("name", "<unnamed>")
        for key in REQUIRED_SPECIES_KEYS:
            if key not in sp:
                errors.append(f"species '{name}': missing key '{key}'.")

        if "samples" in sp:
            if not sp["samples"].get("easy"):
                errors.append(f"species '{name}': samples.easy is empty.")
            if not sp["samples"].get("hard"):
                errors.append(f"species '{name}': samples.hard is empty.")

        if "expression_file" in sp:
            expr_path = Path(config_dir) / sp["expression_file"]
            if not expr_path.exists():
                errors.append(f"species '{name}': expression_file not found: {expr_path}")

    # Check OrthoFinder source — either custom file or results dir must be set
    custom_file = config["pipeline"].get("orthofinder_custom_file")
    og_results  = config["pipeline"].get("orthofinder_results")
    if custom_file:
        custom_path = Path(config_dir) / custom_file
        if not custom_path.exists():
            errors.append(f"pipeline.orthofinder_custom_file not found: {custom_path}")
    elif og_results:
        og_dir = Path(config_dir) / og_results
        if not og_dir.exists():
            errors.append(f"pipeline.orthofinder_results directory not found: {og_dir}")
    else:
        errors.append(
            "pipeline must specify either orthofinder_custom_file or orthofinder_results."
        )

    return errors


# ── step runners ──────────────────────────────────────────────────────────────

def _find_rscript() -> str:
    """Return path to Rscript, checking PATH then common Windows install locations."""
    import shutil
    found = shutil.which("Rscript")
    if found:
        return found
    candidates = [
        r"C:\Program Files\R\R-4.5.2\bin\Rscript.exe",
        r"C:\Program Files\R\R-4.5.1\bin\Rscript.exe",
        r"C:\Program Files\R\R-4.4.3\bin\Rscript.exe",
        r"C:\Program Files\R\R-4.4.2\bin\Rscript.exe",
    ]
    # Also try globbing for any R version under Program Files
    import glob
    candidates += glob.glob(r"C:\Program Files\R\R-*\bin\Rscript.exe")
    for c in candidates:
        if Path(c).exists():
            return c
    raise FileNotFoundError(
        "Rscript not found on PATH or in C:\\Program Files\\R\\.\n"
        "Install R from https://cran.r-project.org/bin/windows/base/"
    )


def run_step1(config: dict, config_file: str, species_filter):
    print("\n" + "=" * 60)
    print("STEP 1: Per-species differential expression (R)")
    print("=" * 60)

    rscript  = _find_rscript()
    r_script = STEPS_DIR / "01_de_analysis.R"
    cmd = [rscript, str(r_script), config_file]
    if species_filter:
        cmd += ["--species", species_filter]

    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"R script failed with exit code {result.returncode}")

    # Validate outputs
    config_dir = str(Path(config_file).parent)
    species_to_check = config["species"]
    if species_filter:
        species_to_check = [s for s in species_to_check if s["name"] == species_filter]

    missing = [
        str(Path(config_dir) / sp["de_output"])
        for sp in species_to_check
        if not (Path(config_dir) / sp["de_output"]).exists()
    ]
    if missing:
        raise RuntimeError(
            "Step 1 completed but DE output files are missing:\n" +
            "\n".join(f"  {p}" for p in missing)
        )
    print("Step 1 complete.")


def run_step2(config: dict, output_dir: str, config_dir: str):
    print("\n" + "=" * 60)
    print("STEP 2: Parse OrthoFinder orthogroups")
    print("=" * 60)
    mod = _load_step("02_parse_orthogroups.py")
    return mod.parse_orthogroups(config, output_dir, config_dir)


def run_step3(config: dict, gene_map: dict, output_dir: str, config_dir: str):
    print("\n" + "=" * 60)
    print("STEP 3: Cross-species conservation scoring")
    print("=" * 60)
    mod = _load_step("03_cross_species.py")
    return mod.run_cross_species(config, gene_map, output_dir, config_dir)


def run_step4(config: dict, output_dir: str, config_dir: str):
    print("\n" + "=" * 60)
    print("STEP 4: Figures and report")
    print("=" * 60)
    mod = _load_step("04_report.py")
    mod.run_report(config, output_dir, config_dir)


def run_step5(config: dict, output_dir: str, config_dir: str):
    print("\n" + "=" * 60)
    print("STEP 5: Arabidopsis annotation")
    print("=" * 60)
    mod = _load_step("05_annotate.py")
    mod.run_annotate(config, output_dir, config_dir)


def load_gene_map_from_csv(output_dir: str) -> dict:
    """Reload gene map from CSV (when skipping Step 2)."""
    import pandas as pd
    p = Path(output_dir) / "orthogroup_gene_map.csv"
    if not p.exists():
        raise FileNotFoundError(
            f"orthogroup_gene_map.csv not found in {output_dir}\n"
            f"Run Step 2 first, or include step 2 in --steps."
        )
    df = pd.read_csv(p)
    gene_map = {}
    for _, row in df.iterrows():
        gene_map[str(row["gene_id"])]             = row["orthogroup_id"]
        gene_map[str(row["normalized_gene_id"])]  = row["orthogroup_id"]
    print(f"  Loaded {len(gene_map)} gene->OG mappings from {p}")
    return gene_map


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="RNA-seq cross-species rooting analysis pipeline"
    )
    parser.add_argument("config", help="Path to config.yaml")
    parser.add_argument("--dry-run", action="store_true",
                        help="Validate config and paths, then exit")
    parser.add_argument("--steps", default="1,2,3,4,5",
                        help="Comma-separated steps to run (default: 1,2,3,4,5)")
    parser.add_argument("--species", default=None,
                        help="Run Step 1 for one species only (by name)")
    parser.add_argument("--output-dir", default=None,
                        help="Override output directory")
    args = parser.parse_args()

    config_file = str(Path(args.config).resolve())
    config_dir  = str(Path(config_file).parent)

    if not Path(config_file).exists():
        print(f"ERROR: Config file not found: {config_file}")
        sys.exit(1)

    with open(config_file) as f:
        config = yaml.safe_load(f)

    errors = validate_config(config, config_dir)

    if args.dry_run:
        if errors:
            print("Config validation FAILED:")
            for e in errors:
                print(f"  - {e}")
            sys.exit(1)
        print("Config is valid.")
        print(f"  Run name : {config['pipeline']['run_name']}")
        print(f"  Species  : {[s['name'] for s in config['species']]}")
        print(f"  Steps    : {args.steps}")
        sys.exit(0)

    if errors:
        print("Config validation errors:")
        for e in errors:
            print(f"  - {e}")
        print("\nAborting. Run with --dry-run for details, or fix errors above.")
        sys.exit(1)

    try:
        steps = [int(s.strip()) for s in args.steps.split(",")]
    except ValueError:
        print(f"ERROR: --steps must be comma-separated integers, got: {args.steps}")
        sys.exit(1)

    # Build output directory
    if args.output_dir:
        output_dir = args.output_dir
    else:
        timestamp  = datetime.datetime.now().strftime("%Y-%m-%d_%H%M%S")
        base_out   = Path(config_dir) / config["pipeline"]["output_dir"]
        output_dir = str(base_out / timestamp)

    Path(output_dir).mkdir(parents=True, exist_ok=True)

    print(f"\nOutput directory : {output_dir}")
    print(f"Run name         : {config['pipeline']['run_name']}")
    print(f"Steps            : {steps}")
    if args.species:
        print(f"Species filter   : {args.species}")

    gene_map = {}

    try:
        if 1 in steps:
            run_step1(config, config_file, args.species)

        if 2 in steps:
            gene_map, _ = run_step2(config, output_dir, config_dir)
        elif any(s in steps for s in [3, 4, 5]):
            gene_map = load_gene_map_from_csv(output_dir)

        if 3 in steps:
            run_step3(config, gene_map, output_dir, config_dir)

        if 4 in steps:
            run_step4(config, output_dir, config_dir)

        if 5 in steps:
            run_step5(config, output_dir, config_dir)

    except Exception as e:
        print(f"\nPIPELINE ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

    print(f"\nPipeline complete. Results in: {output_dir}")


if __name__ == "__main__":
    main()
