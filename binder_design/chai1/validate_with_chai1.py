"""Orthogonal in silico validation and counter-screening with Chai-1 (Apache-2.0 open source).

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
            f">protein|name=target_{name}\n{target_seq}\n"
            f">protein|name=designed_binder\n{args.binder_seq.strip()}\n"
        )
        print(f"Wrote Chai-1 input FASTA: {pair_fasta}")
        print(f"  Run command: chai fold {pair_fasta} {args.out_dir / label}")


if __name__ == "__main__":
    main()
