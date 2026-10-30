#!/usr/bin/env bash
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
run_inference.py \
  inference.input_pdb=structures/AF-P50895-F1-model_v6.pdb \
  inference.output_prefix=outputs/bcam_ig2/bcam_ig2_binder \
  'contigmap.contigs=[A31-255/0 70-110]' \
  'ppi.hotspot_res=[A168,A181,A183,A185,A187,A193]' \
  inference.num_designs=100 \
  denoiser.noise_scale_ca=0.5 \
  denoiser.noise_scale_frame=0.5

# 2. BCAM Option B: Anti-Shedding Membrane-Proximal Stalk (Ig3-Ig5, residues A259-550)
run_inference.py \
  inference.input_pdb=structures/AF-P50895-F1-model_v6.pdb \
  inference.output_prefix=outputs/bcam_stalk/bcam_stalk_binder \
  'contigmap.contigs=[A259-550/0 70-110]' \
  'ppi.hotspot_res=[A341,A345,A507,A509,A511,A535,A537]' \
  inference.num_designs=100 \
  denoiser.noise_scale_ca=0.5 \
  denoiser.noise_scale_frame=0.5

# 3. SLC2A4 (GLUT4): Exofacial Vestibule (EL1 G65-W81 + EL3 G191-L201 + EL4 Y309-P322)
run_inference.py \
  inference.input_pdb=structures/AF-P14672-F1-model_v6.pdb \
  inference.output_prefix=outputs/slc2a4_vestibule/slc2a4_binder \
  'contigmap.contigs=[A1-509/0 55-90]' \
  'ppi.hotspot_res=[A66,A77,A193,A195,A311,A314,A317]' \
  inference.num_designs=100 \
  denoiser.noise_scale_ca=0.5 \
  denoiser.noise_scale_frame=0.5
