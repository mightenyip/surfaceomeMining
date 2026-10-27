"""In silico pharma target stress-test and Pugh matrix for the 15 validated/likely ischemia surface targets."""

from __future__ import annotations

from pathlib import Path
from surfaceome_mining import get_source

CANDIDATES_15 = [
    "BCAM",
    "SLC2A4",
    "DDR1",
    "CD46",
    "KCNQ4",
    "IFNGR2",
    "PIK3IP1",
    "CACNA1S",
    "MOSMO",
    "PDIA3",
    "CTSD",
    "FREM2",
    "TMEM131L",
    "ADAM32",
    "CATSPER2",
]


def main() -> None:
    source = get_source("deliverome", cache_dir=Path(".scratch"))
    header = (
        f"{'Gene':8} | {'Tier':9} | {'TM':2} | {'ECD(aa)':7} | {'pLDDT':5} | "
        f"{'Disord%':7} | {'SBind':5} | {'PDBs':4} | {'Shed/Sec':9} | {'Intern.':8} | "
        f"{'Paralog%':8} | {'MouseECD%':9}"
    )
    print(header)
    print("-" * len(header))

    for sym in CANDIDATES_15:
        bundle = source.gene(sym, include=("record",))
        if not bundle or not bundle.record:
            continue
        rec = bundle.record
        f = rec.get("filters") or {}
        det = rec.get("deterministic_features") or {}
        topo = det.get("canonical_topology") or {}
        struct = det.get("structure") or {}
        sbind = det.get("surface_bind") or {}

        tier = rec.get("deep_dive_tier") or "—"
        tm = topo.get("tm_helix_count") or 0
        ecd_len = topo.get("ecd_length_residues") or 0
        plddt = struct.get("ecd_mean_plddt") or 0.0
        disord = (struct.get("ecd_disordered_fraction") or 0.0) * 100.0
        n_sbind = sbind.get("n_sites") or 0
        n_pdbs = len(sbind.get("pdbs") or [])
        shed = "Shed" if f.get("has_shed_form") else ("Secr" if f.get("has_secreted_form") else "Clean")
        if f.get("has_shed_form") and f.get("has_secreted_form"):
            shed = "Shed+Sec"
        intern = (bundle.catalog.internalization if bundle.catalog else None) or "none"
        paralog = f.get("max_paralog_ecd_pct_identity")
        paralog_str = f"{paralog:7.1f}%" if paralog is not None else "   0.0% "
        mouse_ecd = f.get("mouse_ortholog_ecd_pct_identity")
        mouse_str = f"{mouse_ecd:8.1f}%" if mouse_ecd is not None else "    None "

        print(
            f"{sym:8} | {tier:9} | {tm:2d} | {ecd_len:7d} | {plddt:5.1f} | "
            f"{disord:6.1f}% | {n_sbind:5d} | {n_pdbs:4d} | {shed:9} | {intern:8} | "
            f"{paralog_str:8} | {mouse_str:9}"
        )


if __name__ == "__main__":
    main()
