"""
Shared utilities for the RNA-seq cross-species rooting analysis pipeline.
"""

import re


def normalize_gene_id(gene_id: str) -> str:
    """
    Strip version/isoform suffixes from gene IDs so they match across sources.

    Handles:
      AT1G01010.1                        -> AT1G01010       (Arabidopsis .version suffix)
      Potri.001G000100.1                 -> Potri.001G000100 (Poplar .version suffix)
      Potri.001G000400.0                 -> Potri.001G000400 (Poplar .0 suffix)
      Potri.001G000100.1|PACid_27040687  -> Potri.001G000100 (Poplar OrthoFinder format)
      XLOC_042640_1                      -> XLOC_042640      (Argan _N suffix from OrthoFinder)
      g19068.t1                          -> g19068           (Argan .tN transcript suffix)

    Safe for XLOC_000001 because the trailing segment is 6 digits (>3), not matched.
    """
    gid = str(gene_id).strip()
    # Strip |PACid_... suffix (OrthoFinder poplar format: Potri.001G000100.1|PACid_27040687)
    gid = gid.split('|')[0]
    # Strip .tN transcript suffix (e.g. g19068.t1 -> g19068, common in MAKER/Augustus annotations)
    gid = re.sub(r'\.t\d+$', '', gid)
    # Strip dot-version suffix first (e.g. .1, .2, .10)
    gid = re.sub(r'\.\d{1,3}$', '', gid)
    # Strip underscore-isoform suffix only if 1-3 digits (avoids touching XLOC_000001)
    gid = re.sub(r'_\d{1,3}$', '', gid)
    return gid


def detect_species_from_gene_id(gene_id: str, species_prefixes_dict: dict) -> str | None:
    """
    Map a gene ID to its species name by prefix matching.

    Args:
        gene_id: A single gene identifier string.
        species_prefixes_dict: Dict mapping species_name -> list of prefixes.
                               E.g. {"argan": ["XLOC_"], "poplar": ["Potri."]}

    Returns:
        Species name string if matched, else None.
    """
    normalized = normalize_gene_id(gene_id)
    for species_name, prefixes in species_prefixes_dict.items():
        for prefix in prefixes:
            if normalized.startswith(prefix):
                return species_name
    return None


def build_prefix_map(config: dict) -> dict:
    """
    Build a species_name -> [prefixes] dict from config for use with
    detect_species_from_gene_id().

    Includes both expression species and reference species.
    """
    prefix_map = {}
    for sp in config.get("species", []):
        prefix_map[sp["name"]] = sp.get("gene_id_prefixes", [])
    for sp in config.get("reference_species", []):
        prefix_map[sp["name"]] = sp.get("gene_id_prefixes", [])
    return prefix_map
