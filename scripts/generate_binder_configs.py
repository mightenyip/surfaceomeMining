"""Generate ready-to-run binder design configurations (BindCraft, RFdiffusion, Chai-1) for BCAM and SLC2A4."""

from __future__ import annotations

import json
from pathlib import Path

CACHE_DIR = Path(".scratch/deliverome/genes")
OUT_DIR = Path("binder_design")


def load_gene_record(symbol: str) -> dict:
    return json.loads((CACHE_DIR / f"{symbol}.json").read_text())


def write_fasta(path: Path, entries: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []
    for header, seq in entries:
        lines.append(f">{header}")
        for i in range(0, len(seq), 80):
            lines.append(seq[i : i + 80])
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    bcam = load_gene_record("BCAM")["deterministic_features"]
    slc2a4 = load_gene_record("SLC2A4")["deterministic_features"]

    bcam_h = bcam["canonical_topology"]["sequence"]
    bcam_m = bcam["orthologs"]["mouse"][0]["sequence"]
    bcam_c = bcam["orthologs"]["cynomolgus"][0]["sequence"]

    slc_h = slc2a4["canonical_topology"]["sequence"]
    slc_m = slc2a4["orthologs"]["mouse"][0]["sequence"]
    slc_c = slc2a4["orthologs"]["cynomolgus"][0]["sequence"]

    # 1. Write FASTA files (full-length + domain crops for docking)
    write_fasta(
        OUT_DIR / "sequences" / "bcam_constructs.fasta",
        [
            ("BCAM_HUMAN_P50895_full_1_628", bcam_h),
            ("BCAM_HUMAN_P50895_ECD_31_550", bcam_h[30:550]),
            ("BCAM_HUMAN_P50895_Ig1_Ig2_31_255_Site1_target", bcam_h[30:255]),
            ("BCAM_HUMAN_P50895_Ig3_Ig5_256_550_AntiShed_stalk", bcam_h[255:550]),
            ("Bcam_MOUSE_Q9R069_full_1_622", bcam_m),
            ("BCAM_CYNO_full_1_659", bcam_c),
        ],
    )

    write_fasta(
        OUT_DIR / "sequences" / "slc2a4_constructs.fasta",
        [
            ("SLC2A4_HUMAN_P14672_full_1_509", slc_h),
            ("Slc2a4_MOUSE_P14142_full_1_509", slc_m),
            ("SLC2A4_CYNO_full_1_509", slc_c),
        ],
    )

    # 2. Write BindCraft target settings JSONs (https://github.com/martinpacesa/BindCraft)
    bindcraft_dir = OUT_DIR / "bindcraft"
    bindcraft_dir.mkdir(parents=True, exist_ok=True)

    bindcraft_configs = {
        "bcam_ig2_site1.json": {
            "design_path": "./results/bcam_ig2_site1",
            "binder_name": "BCAM_Ig2_Site1",
            "starting_pdb": "./structures/AF-P50895-F1-model_v6.pdb",
            "chains": "A",
            "target_hotspot_residues": "A168,A179,A181,A183,A185,A187,A190,A193",
            "lengths": [65, 120],
            "number_of_final_designs": 50,
            "_notes": (
                "Targets Deliverome SURFACE-Bind Site 1 (anchor V181 in Ig2 domain). "
                "93.3% identical between Human and Mouse across P181-P195 (PKITWYRNGQRLEVP)."
            ),
        },
        "bcam_ig5_antished_stalk.json": {
            "design_path": "./results/bcam_ig5_stalk",
            "binder_name": "BCAM_Ig5_Stalk",
            "starting_pdb": "./structures/AF-P50895-F1-model_v6.pdb",
            "chains": "A",
            "target_hotspot_residues": "A341,A345,A348,A507,A509,A511,A513,A535,A537",
            "lengths": [65, 110],
            "number_of_final_designs": 50,
            "_notes": (
                "Targets membrane-proximal Ig3 linker (D341-L350, pLDDT=90.3) and Ig5 beta-sheet stalk "
                "(K507-V539) to selectively anchor to intact, uncleaved BCAM on healthy border-zone cells."
            ),
        },
        "slc2a4_exofacial_vestibule.json": {
            "design_path": "./results/slc2a4_vestibule",
            "binder_name": "SLC2A4_GLUT4_Vestibule",
            "starting_pdb": "./structures/AF-P14672-F1-model_v6.pdb",
            "chains": "A",
            "target_hotspot_residues": "A66,A68,A77,A193,A195,A198,A311,A314,A317",
            "lengths": [55, 95],
            "number_of_final_designs": 50,
            "_notes": (
                "Targets the glycan-clear exofacial vestibule of GLUT4 (EL1 G65-W81 + EL3 G191-L201 + EL4 Y309-P322), "
                "staying >12 A clear of the N57 glycosylation site. 90%+ conserved in Mouse Slc2a4."
            ),
        },
    }

    for fname, payload in bindcraft_configs.items():
        (bindcraft_dir / fname).write_text(json.dumps(payload, indent=2) + "\n")

    # 3. Write RFdiffusion runner script (https://github.com/RosettaCommons/RFdiffusion)
    rfd_dir = OUT_DIR / "rfdiffusion"
    rfd_dir.mkdir(parents=True, exist_ok=True)
    rfd_script = """#!/usr/bin/env bash
# Ready-to-run RFdiffusion commands for BCAM and SLC2A4 (GLUT4)
# Prerequisite: download AlphaFold models into ./structures/
set -euo pipefail

mkdir -p structures outputs/bcam_ig2 outputs/bcam_stalk outputs/slc2a4_vestibule

# Download canonical AlphaFold v6 PDB structures if not present
if [ ! -f structures/AF-P50895-F1-model_v6.pdb ]; then
  curl -sL https://alphafold.ebi.ac.uk/files/AF-P50895-F1-model_v6.pdb -o structures/AF-P50895-F1-model_v6.pdb
fi
if [ ! -f structures/AF-P14672-F1-model_v6.pdb ]; then
  curl -sL https://alphafold.ebi.ac.uk/files/AF-P14672-F1-model_v6.pdb -o structures/AF-P14672-F1-model_v6.pdb
fi

# 1. BCAM Option A: Conserved Ig2 Pocket (SURFACE-Bind Site 1, V181)
# Cropped to Ig1-Ig2 (residues A31-255) for fast diffusion
run_inference.py \\
  inference.input_pdb=structures/AF-P50895-F1-model_v6.pdb \\
  inference.output_prefix=outputs/bcam_ig2/bcam_ig2_binder \\
  'contigmap.contigs=[A31-255/0 70-110]' \\
  'ppi.hotspot_res=[A168,A181,A183,A185,A187,A193]' \\
  inference.num_designs=100 \\
  denoiser.noise_scale_ca=0.5 \\
  denoiser.noise_scale_frame=0.5

# 2. BCAM Option B: Anti-Shedding Membrane-Proximal Stalk (Ig3-Ig5, residues A259-550)
run_inference.py \\
  inference.input_pdb=structures/AF-P50895-F1-model_v6.pdb \\
  inference.output_prefix=outputs/bcam_stalk/bcam_stalk_binder \\
  'contigmap.contigs=[A259-550/0 70-110]' \\
  'ppi.hotspot_res=[A341,A345,A507,A509,A511,A535,A537]' \\
  inference.num_designs=100 \\
  denoiser.noise_scale_ca=0.5 \\
  denoiser.noise_scale_frame=0.5

# 3. SLC2A4 (GLUT4): Exofacial Vestibule (EL1 G65-W81 + EL3 G191-L201 + EL4 Y309-P322)
run_inference.py \\
  inference.input_pdb=structures/AF-P14672-F1-model_v6.pdb \\
  inference.output_prefix=outputs/slc2a4_vestibule/slc2a4_binder \\
  'contigmap.contigs=[A1-509/0 55-90]' \\
  'ppi.hotspot_res=[A66,A77,A193,A195,A311,A314,A317]' \\
  inference.num_designs=100 \\
  denoiser.noise_scale_ca=0.5 \\
  denoiser.noise_scale_frame=0.5
"""
    rfd_path = rfd_dir / "run_rfdiffusion.sh"
    rfd_path.write_text(rfd_script)
    rfd_path.chmod(0o755)

    # 4. Write Chai-1 cross-validation & paralog/ortholog counter-screen script
    chai_dir = OUT_DIR / "chai1"
    chai_dir.mkdir(parents=True, exist_ok=True)
    chai_script = '''"""Orthogonal in silico validation and counter-screening with Chai-1 (Apache-2.0 open source).

Install:
    pip install chai_lab

Usage:
    python3 binder_design/chai1/validate_with_chai1.py --binder-seq "MKT..." --target BCAM
"""

from __future__ import annotations

import argparse
from pathlib import Path

TARGETS = {
    "BCAM": {
        "human_on_target": ("BCAM_HUMAN_ECD", 31, 550, "bcam_constructs.fasta", "BCAM_HUMAN_P50895_ECD_31_550"),
        "mouse_ortholog": ("Bcam_MOUSE_full", 1, 622, "bcam_constructs.fasta", "Bcam_MOUSE_Q9R069_full_1_622"),
    },
    "SLC2A4": {
        "human_on_target": ("SLC2A4_HUMAN", 1, 509, "slc2a4_constructs.fasta", "SLC2A4_HUMAN_P14672_full_1_509"),
        "mouse_ortholog": ("Slc2a4_MOUSE", 1, 509, "slc2a4_constructs.fasta", "Slc2a4_MOUSE_P14142_full_1_509"),
    },
}


def read_fasta_entry(fasta_path: Path, entry_id: str) -> str:
    seq_lines: list[str] = []
    recording = False
    for line in fasta_path.read_text().splitlines():
        if line.startswith(">"):
            recording = line[1:].strip() == entry_id
            continue
        if recording:
            seq_lines.append(line.strip())
    if not seq_lines:
        raise KeyError(f"Entry {entry_id} not found in {fasta_path}")
    return "".join(seq_lines)


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate a candidate binder in Chai-1 against Human + Mouse targets.")
    parser.add_argument("--binder-seq", required=True, help="Amino acid sequence of the designed binder")
    parser.add_argument("--target", choices=["BCAM", "SLC2A4"], required=True)
    parser.add_argument("--out-dir", type=Path, default=Path("binder_design/chai1/runs"))
    args = parser.parse_args()

    seq_dir = Path(__file__).resolve().parent.parent / "sequences"
    spec = TARGETS[args.target]

    for label, (name, _, _, fasta_file, entry_id) in spec.items():
        target_seq = read_fasta_entry(seq_dir / fasta_file, entry_id)
        pair_fasta = args.out_dir / f"{args.target}_{label}.fasta"
        pair_fasta.parent.mkdir(parents=True, exist_ok=True)
        pair_fasta.write_text(
            f">protein|name=target_{name}\\n{target_seq}\\n"
            f">protein|name=designed_binder\\n{args.binder_seq.strip()}\\n"
        )
        print(f"Wrote Chai-1 input FASTA: {pair_fasta}")
        print(f"  Run command: chai fold {pair_fasta} {args.out_dir / label}")


if __name__ == "__main__":
    main()
'''
    (chai_dir / "validate_with_chai1.py").write_text(chai_script)
    print("Generated binder_design/ package successfully!")


if __name__ == "__main__":
    main()
