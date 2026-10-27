# Simulated Ischemia Top-40 Surfaceome Deep-Dive & LLM Triage Report

![Simulated ischemia expression shift heatmap](/Users/mightenyip/.gemini/antigravity/brain/e0187f1f-85e4-4d34-a0eb-4283c615cada/.user_uploaded/media_1792654844107_6d56ed39.png)

> [!IMPORTANT]
> **Half of the 40 genes in the ischemia heatmap are false-positive surface annotations (`20/40` rejected at LLM triage, plus `3/20` deep-dived genes graded `no`).**
> Upstream pathway/ligand-receptor databases frequently misclassify intracellular kinases (`CAMK2B`, `STK40`, `ATR`), ER/mitochondrial sensors (`ERN1`/IRE1$\alpha$, `BCL2L1`, `NT5M`, `PLEKHN1`), nuclear factors (`SIRT7`, `RBMX`, `PRMT5`), cytosolic enzymes (`EGLN1`/PHD2, `ACACA`, `CALM2`, `PSMD2`, `NEB`, `TNFAIP2`), inner-leaflet proteins (`THEM4`), and soluble hormones/enzymes (`NPPB`/BNP, `PTGDS`, `IL4I1`, `TGM4`) as `"receptor"`, `"plasma_membrane"`, or `"ecm"`. Deliverome's Claude Sonnet triage + deep-dive pipeline isolates the **15 true/likely cell-surface targets** (`7 canonical`, `8 likely`) and **2 low-evidence candidates** (`LSMEM1`, `FSTL4`).

```mermaid
flowchart LR
    A["40 Simulated-Ischemia Shift Genes"] --> B["Claude Sonnet Triage (/v1/triage)"]
    B -->|20 rejected: 'no'| C["20 Non-Surface Proteins (Cytoplasmic, Nuclear, Mito, ER, Inner-Leaflet, Secreted)"]
    B -->|20 advanced to Deep Dive| D["Deliverome Deep-Dive (/v1/genes)"]
    D --> E["7 Canonical Surface Targets (PIK3IP1, MOSMO, BCAM, CD46, DDR1, KCNQ4, IFNGR2)"]
    D --> F["8 Likely / Context-Gated Targets (SLC2A4, CACNA1S, PDIA3, CTSD, CATSPER2, ADAM32, FREM2, TMEM131L)"]
    D --> G["2 Low-Evidence Candidates (LSMEM1, FSTL4)"]
    D --> H["3 Downgraded to 'no' (PTGDS, IL4I1, PLEKHN1)"]
```

---

## 1. Validated & Likely Cell-Surface Targets (`canonical` and `likely` tiers — 15 genes)

These 15 genes pass both LLM triage and per-gene deep-dive validation (`SurfaceomeRecord` v2.14). Notice that several top ischemia-shift genes (`SLC2A4`/GLUT4, `PDIA3`/ERp57, `CTSD`, `PIK3IP1`, `CACNA1S`, `KCNQ4`) have direct cardiovascular, skeletal/vascular muscle, or stress/ischemia-regulated surface biology:

| Symbol (Plot Input) | Plot Label | DB Votes | Triage | Deep-Dive Tier | Access. | Call Reason | State Dep. | Topology (TM / ECD aa) | Mouse ECD %ID | Intern. | Bind / Tag Sites | Headline Risks |
|---|---|---:|---|---|---|---|---|---|---:|---|---|---|
| **PIK3IP1** (`Pik3ip1`) | `cell_surface` | 2/5 | `yes` | **canonical** | `high` | `cell_state_induced` | `high` | 1 TM / 151 aa | 76.3% | `moderate` | 1 / 11 | `shed_form`, `epitope_masked` |
| **SLC2A4** (`Slc2a4`) | `receptor` | 4/5 | `contextual` | **likely** | `high` | `cell_state_induced` | `high` | 12 TM / 71 aa | 85.1% | `moderate` | 0 / 4 | `epitope_masked` |
| **CACNA1S** (`Cacna1s`) | `receptor` | 2/5 | `yes` | **likely** | `high` | `tissue_restricted_surface` | `moderate` | 24 TM / 461 aa | 94.1% | `moderate` | 0 / 0 | `epitope_masked` |
| **MOSMO** (`Mosmo`) | `cell_surface` | 3/5 | `yes` | **canonical** | `moderate` | `dual_localization` | `moderate` | 4 TM / 64 aa | 100.0% | `moderate` | 0 / 2 | `epitope_masked` |
| **PDIA3** (`Pdia3`) | `receptor` | 1/5 | `contextual` | **likely** | `high` | `cell_state_induced` | `high` | 0 TM / 481 aa | 93.8% | `low` | 0 / 0 | `secreted_form`, `epitope_masked` |
| **BCAM** (`Bcam`) | `receptor` | 4/5 | `yes` | **canonical** | `high` | `classical_surface_receptor` | `moderate` | 1 TM / 520 aa | 72.7% | `moderate` | 2 / 9 | `shed_form`, `secreted_form` |
| **CTSD** (`Ctsd`) | `receptor` | 1/5 | `contextual` | **likely** | `moderate` | `cell_state_induced` | `high` | 0 TM / 392 aa | 83.7% | `very_low` | 0 / 0 | `secreted_form` |
| **CD46** (`Cd46`) | `receptor` | 5/5 | `yes` | **canonical** | `high` | `classical_surface_receptor` | `moderate` | 1 TM / 311 aa | — | `moderate` | 3 / 6 | `shed_form`, `secreted_form` |
| **CATSPER2** (`Catsper2`) | `receptor` | 2/5 | `contextual` | **likely** | `high` | `tissue_restricted_surface` | `high` | 6 TM / 75 aa | 70.9% | `moderate` | 0 / 0 | `co_receptor`, `epitope_masked` |
| **ADAM32** (`Adam32`) | `receptor` | 2/5 | `contextual` | **likely** | `moderate` | `tissue_restricted_surface` | `high` | 1 TM / 666 aa | 60.9% | `moderate` | 3 / 32 | none |
| **DDR1** (`Ddr1`) | `receptor` | 3/5 | `yes` | **canonical** | `high` | `classical_surface_receptor` | `moderate` | 1 TM / 397 aa | 92.7% | `high` | 3 / 18 | `shed_form` |
| **FREM2** (`Frem2`) | `cell_surface` | 3/5 | `no` | **likely** | `moderate` | `cell_state_induced` | `high` | 1 TM / 3064 aa | 89.1% | `low` | 6 / 41 | `shed_form`, `epitope_masked` |
| **KCNQ4** (`Kcnq4`) | `receptor` | 2/5 | `yes` | **canonical** | `high` | `multipass_with_exposed_loops` | `moderate` | 6 TM / 56 aa | 100.0% | `moderate` | 0 / 3 | `epitope_masked` |
| **IFNGR2** (`Ifngr2`) | `receptor` | 2/5 | `yes` | **canonical** | `high` | `classical_surface_receptor` | `moderate` | 1 TM / 227 aa | 55.1% | `moderate` | 3 / 7 | `epitope_masked` |
| **TMEM131L** (`Tmem131l`) | `cell_surface` | 3/5 | `contextual` | **likely** | `moderate` | `tissue_restricted_surface` | `high` | 1 TM / 834 aa | 82.9% | `moderate` | 5 / 25 | none |

### Per-Gene Deep-Dive Dossiers (Top Ischemia-Relevant & Canonical Receptors)

- **PIK3IP1 (`Pik3ip1`, `canonical`, high accessibility, Bin 4 terminal up-shift):**
  - *Topology & Conservation:* Single-pass Type I TM protein (`Q96FE7`), 151 aa kringle-containing ECD, 72 aa cytoplasmic PI3K-inhibitory tail. Mouse ECD identity: **76.3%** (cyno: **94.0%**). 1 SURFACE-Bind site, 11 tag sites, `moderate` internalization.
  - *Biology & Risks:* Constitutively surface-accessible on resting T cells and expressed in cardiac/hepatic tissues; undergoes rapid ADAM-protease ectodomain shedding upon activation (`shed_form` moderate risk; `epitope_masked` via ECD oligomerization).

- **SLC2A4 / GLUT4 (`Slc2a4`, `likely`, high accessibility, state-gated translocation):**
  - *Topology & Conservation:* 12-TM glucose transporter (`P14672`), 71 aa total exposed ECD loops (mean pLDDT 87.0). Mouse ECD identity: **85.1%** (FL: **95.3%**, cyno ECD: **100%**). 4 tag sites, `moderate` internalization.
  - *Biology & Risks:* Classic ischemia- and insulin-translocated transporter in cardiomyocytes, skeletal muscle, and adipocytes. Sequestered in intracellular GLUT4 storage vesicles (GSVs) at baseline, but acute ischemia/hypoxia, contraction, or insulin drives rapid sarcolemmal translocation. Main caveat is N57 glycan-mediated `epitope_masked` on the large exofacial loop; no shed/secreted decoy pool.

- **BCAM (`Bcam`, `canonical`, high accessibility, early Bin 0 high -> ischemia down-shift):**
  - *Topology & Conservation:* Single-pass Type I TM receptor (`P50895`) with 5 Ig-like domains (520 aa ECD, pLDDT 90.4). Mouse ECD identity: **72.7%** (cyno: **96.5%**). 2 SURFACE-Bind sites, 9 tag sites, `moderate` internalization.
  - *Biology & Risks:* Constitutively expressed on endothelial, epithelial, and erythroid cells (laminin $\alpha$5 receptor). Exhibits ADAM10/MMP-mediated shedding (`shed_form`) releasing soluble sBCAM.

- **CD46 (`Cd46`, `canonical`, high accessibility, 5/5 DB consensus):**
  - *Topology & Conservation:* Single-pass Type I complement regulatory receptor (`P15529`, 311 aa ECD, pLDDT 87.1). Cyno ECD identity: **84.2%** (rodents lack a 1:1 somatic CD46 surface ortholog, restricting mouse CD46 expression to spermatozoa/testis via *Crry* functional divergence). 3 SURFACE-Bind sites, 6 tag sites, `moderate` internalization.
  - *Biology & Risks:* Pan-nucleated human cell surface complement regulator; MMP-mediated ectodomain shedding (`shed_form`) is documented in acute coronary syndrome (ACS) and inflammation.

- **DDR1 (`Ddr1`, `canonical`, high accessibility, Bin 4 terminal strong up-shift):**
  - *Topology & Conservation:* Single-pass Type I collagen receptor tyrosine kinase (`Q08345`, 397 aa ECD, pLDDT 85.5). Mouse ECD identity: **92.7%** (cyno: **99.5%**). 3 SURFACE-Bind sites, 18 tag sites, **`high` internalization**.
  - *Biology & Risks:* Pan-epithelial and vascular smooth muscle collagen sensor; strong internalization grade (`high`) makes it well-suited for internalizing payloads, though constitutive MT1-MMP cleavage (`shed_form`) produces soluble sDDR1 ectodomain.

- **IFNGR2 (`Ifngr2`, `canonical`, high accessibility, Bin 4 terminal up-shift):**
  - *Topology & Conservation:* Single-pass Type I cytokine receptor chain (`P38484`, 227 aa ECD, pLDDT 92.9). Mouse ECD identity: **55.1%** (cyno: **94.7%** — note low mouse-human cross-reactivity). 3 SURFACE-Bind sites, 7 tag sites, `moderate` internalization.
  - *Biology & Risks:* Constitutively surface-accessible across endothelial, epithelial, and myeloid cells; upregulated on M1 macrophages and during inflammatory stress, with rapid LI-motif endocytosis.

- **KCNQ4 (`Kcnq4`, `canonical`, high accessibility) & CACNA1S (`Cacna1s`, `likely`, high accessibility):**
  - *Topology & Conservation:* Multi-pass ion channels (`KCNQ4`: 6 TM, 56 aa ECD, **100%** mouse ECD identity; `CACNA1S`: 24 TM, 461 aa ECD, **94.1%** mouse ECD identity).
  - *Biology & Risks:* `KCNQ4` is surface-accessible on vascular smooth muscle cells (regulating vascular tone) and cochlear hair cells, though in cardiomyocytes it localizes to mitochondria. `CACNA1S` (Cav1.1) is restricted to skeletal muscle T-tubule/sarcolemmal membranes. Both have compact extracellular loops with `epitope_masked` caveats.

- **PDIA3 / ERp57 (`Pdia3`, `likely`, high stress-induced accessibility) & CTSD (`Ctsd`, `likely`, moderate stress-induced accessibility):**
  - *Topology & Conservation:* Non-TM luminal proteins (`PDIA3`: 481 aa, **93.8%** mouse identity; `CTSD`: 392 aa, **83.7%** mouse identity) that translocate to the outer plasma membrane during cellular stress, immunogenic cell death, platelet activation (`PDIA3`), or lysosomal exocytosis / M6P-receptor retention (`CTSD`).
  - *Biology & Risks:* Relevant to late-bin ischemic stress (`PDIA3` peaks in Bin 4), but both carry `secreted_form` decoy risk and `low`/`very_low` internalization.

- **MOSMO (`Mosmo`, `canonical`, 4 TM, 100% mouse ECD identity), FREM2 (`Frem2`, `likely`, 1 TM / 3064 aa ECD, 89.1% mouse ECD identity), TMEM131L (`Tmem131l`, `likely`, 1 TM / 834 aa ECD, 82.9% mouse ECD identity), CATSPER2 (`Catsper2`, sperm-restricted), ADAM32 (`Adam32`, sperm-restricted):**
  - `MOSMO` is a 4-TM plasma-membrane/ciliary Hedgehog modulator (complexed with MEGF8). `CATSPER2` and `ADAM32` are strictly testis/spermatozoa-restricted in normal physiology.

---

## 2. Deep-Dived Candidates Graded `low` or `no` (5 genes)

Five genes from the plot entered the Deliverome deep-dive pipeline because of borderline triage calls (`contextual`) or 1 DB vote, but the literature ledger downgraded them to `low` (no protein-level surface validation yet) or `no` (soluble secreted or mitochondrial):

| Symbol (Plot Input) | Plot Label | DB Votes | Triage Call | Deep-Dive Tier | Access. | Call Reason | Deep-Dive Verdict & Rationale |
|---|---|---:|---|---|---|---|---|
| **LSMEM1** (`Lsmem1`) | `cell_surface` | 1/5 | `contextual` (`dual_localization`) | **low** | `low` | `cell_state_induced` | 1-TM predicted topology (48 aa ECD, 85.4% mouse ECD identity), upregulated in CKD kidney tubules, but **zero direct experimental plasma-membrane readouts** exist in literature. |
| **FSTL4** (`Fstl4`) | `ecm` | 0/5 | `contextual` (`lysosomal_exocytosis`) | **low** | `low` | `cell_state_induced` | 0-TM follistatin-like protein (820 aa); only RNA-level upregulation in lung adenocarcinoma is reported, with no direct surface protein evidence. |
| **PTGDS** (`Ptgds`) | `receptor` | 0/5 | `contextual` (`stable_surface_attachment`) | **no** | `no` | `secreted_only` | L-PGDS (prostaglandin D2 synthase) is a **freely secreted lipocalin** (0 TM; 2nd most abundant protein in CSF) whose *product* ($\text{PGD}_2$) binds receptors, not PTGDS itself. |
| **IL4I1** (`Il4i1`) | `ecm` | 1/5 | `contextual` (`dual_localization`) | **no** | `no` | `secreted_only` | Secreted L-amino-acid oxidase (0 TM) released into the synaptic cleft and extracellular space with no stable membrane anchor. |
| **PLEKHN1** (`Plekhn1`) | `plasma_membrane` | 1/5 | `contextual` (`dual_localization`) | **no** | `no` | `mitochondrial_internal` | N-myristoylated PH-domain protein anchored to the **cytoplasmic face of the mitochondrial outer membrane** and nucleus (0 TM, 0 aa ECD); SIM imaging explicitly refutes plasma-membrane presence. |

---

## 3. Triage-Rejected Candidates (`20/40` Genes Filtered Out by Claude Sonnet)

These 20 genes were labeled `receptor`, `plasma_membrane`, `ecm`, or `ligand` in the source plot, but Deliverome's Claude Sonnet triage (`claude-sonnet-4-6`) rejected all 20 as non-surface (`verdict: "no"`):

| Symbol (Plot Input) | Plot Label | DB Votes | Claude Verdict | Triage Reason | Mechanistic Summary from `/v1/triage/{SYMBOL}` |
|---|---|---:|---|---|---|
| **THEM4** (`Them4`) | `plasma_membrane` | 0/5 | `no` | `inner_leaflet_anchored` | Thioesterase (CTMP) bound to the **cytoplasmic inner leaflet** of the PM and mitochondria where it regulates Akt; no TM domain or ectodomain. |
| **CALM2** (`Calm2`) | `plasma_membrane` | 0/5 | `no` | `cytoplasmic` | Calmodulin-2: soluble intracellular $\text{Ca}^{2+}$ sensor in the cytoplasm; binds cytoplasmic tails of channels but never crosses the bilayer. |
| **RBMX** (`Rbmx`) | `ecm` | 0/5 | `no` | `nuclear` | Heterogeneous nuclear ribonucleoprotein G (hnRNP-G); strictly nuclear pre-mRNA splicing factor. |
| **NT5M** (`Nt5m`) | `receptor` | 0/5 | `no` | `mitochondrial_internal` | Mitochondrial matrix $5',3'$-nucleotidase; resides inside the mitochondrial matrix. |
| **ERN1** (`Ern1`) | `receptor` | 0/5 | `no` | `endomembrane_resident` | IRE1$\alpha$: canonical **ER-resident** unfolded-protein-response transmembrane sensor; luminal domain is inside the ER lumen, not the plasma membrane. |
| **CAMK2B** (`Camk2b`) | `receptor` | 0/5 | `no` | `cytoplasmic` | Cytosolic/synaptic $\text{Ca}^{2+}$/calmodulin-dependent serine/threonine kinase II $\beta$. |
| **PSMD2** (`Psmd2`) | `ecm` | 0/5 | `no` | `cytoplasmic` | 26S proteasome non-ATPase regulatory subunit 2 (Rpn1); cytosolic/nuclear proteasome lid scaffold. |
| **SIRT7** (`Sirt7`) | `receptor` | 0/5 | `no` | `nuclear` | Nucleolar/chromatin $\text{NAD}^+$-dependent deacylase (class IV sirtuin). |
| **NEB** (`Neb`) | `ecm` | 0/5 | `no` | `cytoplasmic` | Nebulin: giant (~800 kDa) intrasarcomeric actin-binding scaffold inside skeletal muscle cytoplasm. |
| **ACACA** (`Acaca`) | `receptor` | 0/5 | `no` | `cytoplasmic` | Acetyl-CoA carboxylase 1 (ACC1): cytosolic rate-limiting enzyme in fatty acid synthesis. |
| **TGM4** (`Tgm4`) | `ecm` | 0/5 | `no` | `secreted_only` | Prostate transglutaminase: soluble secreted luminal enzyme with no TM or GPI anchor. |
| **BCL2L1** (`Bcl2l1`) | `receptor` | 0/5 | `no` | `cytoplasmic` | Bcl-xL: outer mitochondrial and ER membrane anti-apoptotic protein facing the cytosol. |
| **EGLN1** (`Egln1`) | `receptor` | 0/5 | `no` | `cytoplasmic` | PHD2 (prolyl hydroxylase domain-containing protein 2): cytosolic oxygen sensor that hydroxylates HIF-$\alpha$ (tagged `"receptor"` only because it is an intracellular oxygen sensor). |
| **PRMT5** (`Prmt5`) | `receptor` | 0/5 | `no` | `cytoplasmic` | Cytoplasmic/nuclear protein arginine methyltransferase 5. |
| **OSBPL2** (`Osbpl2`) | `plasma_membrane` | 1/5 | `no` | `endomembrane_resident` | Oxysterol-binding protein-related protein 2 (ORP2): cytosolic/Golgi/ER-membrane-contact-site lipid transfer protein facing the cytoplasm. |
| **TNFAIP2** (`Tnfaip2`) | `ecm` | 0/5 | `no` | `cytoplasmic` | Intracellular exocyst-associated protein (EXOC3L3) localized to cytosol and tunneling nanotubes. |
| **STK40** (`Stk40`) | `receptor` | 0/5 | `no` | `cytoplasmic` | Cytosolic/nuclear Tribbles-family pseudokinase. |
| **ATR** (`Atr`) | `receptor` | 0/5 | `no` | `nuclear` | Nuclear DNA-damage sensor kinase (ataxia telangiectasia and Rad3-related). |
| **TXLNA** (`Txlna`) | `ligand` | 0/5 | `no` | `cytoplasmic` | $\alpha$-Taxilin: cytoplasmic syntaxin-binding vesicle traffic regulator (formerly misnamed IL-14). |
| **NPPB** (`Nppb`) | `ligand` | 0/5 | `no` | `secreted_only` | Brain Natriuretic Peptide (BNP): soluble circulating cardiac hormone cleaved and secreted into plasma. |

---

## 4. Claude Sonnet / Opus Triage & SurfaceBench Evaluation Framework

We added two new modules to `surfaceome_mining` so you can both inspect Deliverome's benchmark/triage data and run your own Claude Sonnet/Opus triage passes:

- [deliverome.py](file:///Users/mightenyip/Documents/GitHub/surfaceomeMining/src/surfaceome_mining/sources/deliverome.py): Added `benchmark()`, `benchmark_gene()`, `benchmark_matrix()`, `triage()`, `triage_cell()`, `genes_index()`, automatic uppercase HGNC normalization (`Pik3ip1` $\rightarrow$ `PIK3IP1`), and disk caching for all per-gene routes.
- [deep_dive.py](file:///Users/mightenyip/Documents/GitHub/surfaceomeMining/src/surfaceome_mining/deep_dive.py): Added [`DeepDiveProfile`](file:///Users/mightenyip/Documents/GitHub/surfaceomeMining/src/surfaceome_mining/deep_dive.py#L31-L86), [`summarize_gene_bundle`](file:///Users/mightenyip/Documents/GitHub/surfaceomeMining/src/surfaceome_mining/deep_dive.py#L89-L214), [`deep_dive_genes`](file:///Users/mightenyip/Documents/GitHub/surfaceomeMining/src/surfaceome_mining/deep_dive.py#L217-L254), and [`format_deep_dive_markdown`](file:///Users/mightenyip/Documents/GitHub/surfaceomeMining/src/surfaceome_mining/deep_dive.py#L266-L335).
- [triage.py](file:///Users/mightenyip/Documents/GitHub/surfaceomeMining/src/surfaceome_mining/triage.py): Added the full SurfaceBench taxonomy (`yes` / `contextual` / `no` and all 19 canonical reasons), replicate majority aggregation ([`compute_majority`](file:///Users/mightenyip/Documents/GitHub/surfaceomeMining/src/surfaceome_mining/triage.py#L170-L217)), SurfaceBench scoring ([`score_benchmark_matrix`](file:///Users/mightenyip/Documents/GitHub/surfaceomeMining/src/surfaceome_mining/triage.py#L220-L280) and [`evaluate_predictions_against_benchmark`](file:///Users/mightenyip/Documents/GitHub/surfaceomeMining/src/surfaceome_mining/triage.py#L283-L334)), prompt construction ([`build_triage_prompt`](file:///Users/mightenyip/Documents/GitHub/surfaceomeMining/src/surfaceome_mining/triage.py#L337-L382)), and [`ClaudeTriageRunner`](file:///Users/mightenyip/Documents/GitHub/surfaceomeMining/src/surfaceome_mining/triage.py#L385-L493).

### SurfaceBench 147-Gene Accuracy Across Claude Models (`/v1/benchmark/matrix`)

| Model | Prompt Variant | Evaluated | Correct (`yes ≡ contextual`) | Refused | Accuracy | False Positives | False Negatives |
|---|---|---:|---:|---:|---:|---|---|
| `claude-sonnet-5-5` | `ncbi` | 147 | 145 | 0 | **98.6%** | ITPR3 | SRC |
| `claude-opus-4-8` | `ncbi` | 147 | 144 | 0 | **98.0%** | IZUMO4 | SRC, TMED10 |
| `claude-sonnet-4-6` | `ncbi` | 147 | 144 | 0 | **98.0%** | — | KLK2, SRC, TGOLN2 |
| `claude-sonnet-4-6` | `web_ncbi` | 147 | 143 | 0 | **97.3%** | ATG9A | ATP5F1B, C3, SRC |
| `claude-sonnet-4-6` | `pubmed_ncbi` | 147 | 143 | 0 | **97.3%** | ABCB9, ATG9A, PRAME | TGOLN2 |
| `claude-sonnet-5-5` | `naive` | 147 | 143 | 0 | **97.3%** | ATG9A, ITPR3, IZUMO4 | SRC |
| `claude-opus-5-5` | `ncbi_effort_low` | 145 | 141 | 2 | **97.2%** | ABCB9, ATG9A | FN1, SRC |
| `claude-opus-4-8` | `naive` | 147 | 142 | 0 | **96.6%** | ABCB9, IZUMO4 | ATP5F1B, KLK2, SRC |
| `claude-sonnet-4-6` | `naive` | 147 | 141 | 0 | **95.9%** | F2 | B4GALT1, SRC, TGOLN2, TMED10, VDAC1 |
| `claude-haiku-5-5` | `ncbi` | 147 | 141 | 0 | **95.9%** | ABCB9, ATG9A, IZUMO4 | KLK2, SRC, STIM1 |

### CLI & Python Usage Examples

```bash
# 1. Generate a deep-dive markdown or JSON report for any gene list
surfaceome-mine deep-dive \
  --symbols "Pik3ip1 (cell_surface),Slc2a4 (receptor),Cacna1s (receptor),Them4 (plasma_membrane)" \
  --cache-dir .scratch

# 2. View the SurfaceBench leaderboard across Claude models and prompt variants
surfaceome-mine benchmark --cache-dir .scratch

# 3. Inspect Deliverome's recorded Claude triage runs & reasoning for specific genes
surfaceome-mine triage PIK3IP1 SLC2A4 THEM4 --cache-dir .scratch

# 4. Run live Claude Sonnet / Opus triage (with 3-replicate majority voting & optional disease context)
ANTHROPIC_API_KEY="sk-ant-..." surfaceome-mine triage SLC2A4 PDIA3 THEM4 \
  --run-llm \
  --model claude-sonnet-4-6 \
  --variant ncbi \
  --replicates 3 \
  --context "simulated myocardial/skeletal ischemia"
```

---

## 5. *In Silico* Stress-Test & Pugh-Chart Scorecard for MI Reperfusion Border-Zone Anchoring

### Definitions: ECD & pLDDT
- **ECD (Extracellular Domain):** The portion(s) of a membrane protein that protrude outside the plasma membrane into the extracellular space (measured in amino acid residues, `aa`). For a single-pass transmembrane protein (like `BCAM` or `DDR1`), this is one continuous ectodomain; for a multi-pass channel/transporter (like `SLC2A4` or `KCNQ4`), it is the sum of the short loops connecting transmembrane helices on the outer face. *In silico* protein/antibody design requires sufficient exposed ECD surface area (typically $>150\text{ aa}$ for standard antibodies, or well-ordered loops for nanobodies/peptides).
- **pLDDT (predicted Local Distance Difference Test):** AlphaFold2's per-residue confidence score ($0\text{--}100$) averaged across the ECD (`ecd_mean_plddt`), which strongly correlates with structural rigidity vs. intrinsic disorder:
  - **$\text{pLDDT} \ge 90$:** Very high confidence; rigid, well-folded 3D backbone ideal for structure-based *in silico* docking and *de novo* binder design (RFdiffusion / BindCraft).
  - **$70 \le \text{pLDDT} < 90$:** Confident backbone fold; suitable for binder design.
  - **$\text{pLDDT} < 70$ (or high `Disord%`):** Flexible or intrinsically disordered region; high entropic penalty makes rigid *in silico* docking unreliable.

### Raw Stress-Test Metrics Across the 15 Surface Candidates ([scripts/stress_test_15.py](file:///Users/mightenyip/Documents/GitHub/surfaceomeMining/scripts/stress_test_15.py))

```text
Gene     | Tier      | TM | ECD(aa) | pLDDT | Disord% | SBind | PDBs | Shed/Sec  | Intern.  | Paralog% | MouseECD%
------------------------------------------------------------------------------------------------------------------
BCAM     | canonical |  1 |     520 |  90.4 |    3.5% |     2 |    2 | Shed+Sec  | moderate |    31.9% |     72.7%
SLC2A4   | likely    | 12 |      71 |  87.0 |    0.0% |     0 |    0 | Clean     | moderate |    59.2% |     85.1%
DDR1     | canonical |  1 |     397 |  85.5 |   15.1% |     3 |   33 | Shed+Sec  | high     |    50.1% |     92.7%
CD46     | canonical |  1 |     311 |  87.1 |   17.7% |     3 |    6 | Shed+Sec  | moderate |    31.5% |     None 
KCNQ4    | canonical |  6 |      56 |  91.4 |    5.4% |     0 |    0 | Clean     | moderate |    78.6% |    100.0%
IFNGR2   | canonical |  1 |     227 |  92.9 |    4.8% |     3 |    3 | Clean     | moderate |    30.8% |     55.1%
PIK3IP1  | canonical |  1 |     151 |  64.3 |   55.6% |     1 |    0 | Shed      | moderate |    14.0% |     76.3%
CACNA1S  | likely    | 24 |     461 |  85.8 |    8.9% |     0 |    0 | Clean     | moderate |    67.0% |     94.1%
MOSMO    | canonical |  4 |      64 |  85.1 |   12.5% |     0 |    0 | Clean     | moderate |    0.0%  |    100.0%
PDIA3    | likely    |  0 |     481 |  93.4 |    2.7% |     0 |    0 | Secr      | low      |    35.0% |     93.8%
CTSD     | likely    |  0 |     392 |  86.8 |   14.5% |     0 |    0 | Secr      | very_low |    47.5% |     83.7%
FREM2    | likely    |  1 |    3064 |  72.8 |   28.1% |     6 |    0 | Shed      | low      |    40.6% |     89.1%
TMEM131L | likely    |  1 |     834 |  70.3 |   39.2% |     5 |    0 | Clean     | moderate |    26.9% |     82.9%
ADAM32   | likely    |  1 |     666 |  82.4 |   17.9% |     3 |    0 | Clean     | moderate |    43.8% |     60.9%
CATSPER2 | likely    |  6 |      75 |  85.0 |    5.3% |     0 |    0 | Clean     | moderate |    28.0% |     70.9%
```

### Pugh-Chart Decision Matrix: Anchor to Healthy / Viable Peri-Infarct Cells During MI Reperfusion

**Scoring Key:**Each criterion is scored relative to a baseline reference (`0` = neutral/acceptable, `+1` = favorable, `+2` = ideal, `-1` = unfavorable, `-2` = disqualifying bottleneck), weighted by importance for **anchoring a cytoprotective therapeutic to healthy/viable peri-infarct myocardium or coronary microvasculature surrounding infarcted cells**:
1. **Healthy / Peri-Infarct Border-Zone Trajectory (Weight 3$\times$):** High in Bin 0 (healthy baseline) and/or Bins 2–3 (viable ischemic border zone), and *downregulated* in Bin 4 (dead infarct core) so the drug anchors to salvageable cells rather than necrotic debris.
2. **Cardiovascular / Myocardial Expression (Weight 3$\times$):** Expressed on cardiomyocytes, coronary microvascular endothelium, or arterial smooth muscle.
3. ***In Silico* Binder Druggability (Weight 2$\times$):** Large, rigid ECD ($\text{pLDDT} > 85$, low disorder), experimental PDB structures, and validated SURFACE-Bind pockets.
4. **Low Soluble Sink / Shedding Risk (Weight 2$\times$):** Absence of circulating decoy protein, or presence of a non-shed membrane-proximal stalk domain.
5. **Surface Retention vs. Rapid Clearance (Weight 1$\times$):** Low-to-moderate internalization keeps the protective payload anchored at the cell membrane or pericellular space.
6. **Paralog Selectivity & Preclinical Mouse Conservation (Weight 1$\times$):** Low paralog ECD identity ($<50\%$) to avoid off-target binding, and high mouse ECD identity ($>70\%$) for *in vivo* mouse LAD ischemia-reperfusion validation.

| Gene | 1. Healthy / Border-Zone Trajectory (3$\times$) | 2. Cardiac / Vascular Expression (3$\times$) | 3. *In Silico* Druggability: ECD, pLDDT, PDBs (2$\times$) | 4. Low Soluble Decoy Risk (2$\times$) | 5. Surface Retention (1$\times$) | 6. Selectivity & Mouse Ortho (1$\times$) | **Weighted Total** | Rank & Verdict for MI Reperfusion Anchor |
|---|:---:|:---:|:---:|:---:|:---:|:---:|---:|---|
| **BCAM** | **+2** (Peak Bin 0 healthy; lost in infarct) | **+1** (Coronary/vascular endothelium) | **+2** (520 aa, pLDDT 90.4, 2 PDBs, 2 SBind) | **-1** (ADAM10 shed; target Ig4–5 stalk) | **+1** (`moderate` intern.) | **+1** (31.9% paralog, 72.7% mouse) | **+13** | **#1 Best Overall for *De Novo* Protein/mAb Design** (Coronary microvascular endothelial anchor) |
| **SLC2A4** (GLUT4) | **+2** (Bin 0 + strong Bins 2–3 viable ischemic border zone) | **+2** (Cardiomyocyte sarcolemma) | **-1** (71 aa across 12-TM loops, N57 glycan) | **+2** (Zero shedding/secretion) | **+1** (Regulated cycling) | **+1** (59.2% SLC2A1, 85.1% mouse) | **+13** | **#1 (Tied) Best Cardiomyocyte-Specific Anchor** (Requires loop-specific nanobody/peptide) |
| **KCNQ4** | **+1** (Peak Bin 2 mid-ischemia border zone; low Bin 4) | **+1** (Coronary/arterial smooth muscle & endothelium; *mito in CM*) | **-1** (56 aa across 6-TM loops, pLDDT 91.4) | **+2** (Zero shedding) | **+1** (`moderate` intern.) | **0** (78.6% KCNQ5 paralog, 100% mouse) | **+8** | **#3 Arterial / Microvascular Tone Anchor** (Small loop target; watch KCNQ5 cross-reactivity) |
| **CD46** | **+2** (Peak Bin 0 healthy; shed in infarct) | **+1** (Endothelium + all nucleated cells) | **+2** (311 aa, pLDDT 87.1, 6 PDBs, 3 SBind) | **-1** (MMP-1 sheds sCD46 in ACS plasma) | **+1** (`moderate` intern.) | **-2** (Pan-tissue sink; **no mouse ortholog**) | **+7** | **#4 Biologically Ideal (Complement Regulator) but Preclinically Blocked** (No mouse ortholog; systemic sink) |
| **DDR1** | **-1** (Low Bin 0; spikes in Bin 4 terminal infarct) | **+1** (Vascular smooth muscle & endothelium) | **+2** (397 aa, pLDDT 85.5, 33 PDBs, 3 SBind) | **-1** (MT1-MMP shedding) | **-1** (`high` intern. clears surface) | **+1** (50.1% DDR2, 92.7% mouse) | **+2** | **#5 Best for Late-Infarct Intracellular Delivery (ADC/Degrader)**, not healthy-cell surface anchoring |
| **IFNGR2** | **-1** (Low Bin 0; spikes in Bin 4 terminal infarct) | **+1** (Endothelium & myeloid cells) | **+2** (227 aa, pLDDT 92.9, 3 PDBs, 3 SBind) | **+2** (No shedding) | **0** (Rapid LI-motif endocytosis) | **0** (30.8% paralog, only 55.1% mouse) | **+4** | **#6 Late Inflammatory Marker** (Peaks in terminal infarct rather than healthy cells) |
| **PIK3IP1** | **-1** (Low Bin 0; spikes in Bin 4 terminal infarct) | **0** (Mainly T cells / liver; low cardiac surface proof) | **-1** (151 aa, low pLDDT 64.3, 55.6% disordered) | **-1** (ADAM shedding) | **+1** (`moderate` intern.) | **+2** (14.0% paralog, 76.3% mouse) | **-2** | **#7 Poor Structural Rigidity** (55.6% disordered ECD; terminal Bin 4 peak) |
| **CACNA1S** | **-1** (Low Bin 0–2; spikes in Bin 3–4) | **-1** (Skeletal muscle Cav1.1; cardiac uses CACNA1C) | **0** (461 aa across 24-TM loops, pLDDT 85.8) | **+2** (Clean) | **+1** (`moderate` intern.) | **0** (67% CACNA1D, 94.1% mouse) | **-1** | **#8 Wrong Isoform for Heart** (Skeletal muscle restricted; heart uses Cav1.2) |
| **PDIA3** | **-1** (Spikes in Bin 4 dying/stressed cells) | **0** (ER resident; surface DAMP on injured cells/platelets) | **0** (481 aa, pLDDT 93.4, 0 TM anchor) | **-1** (Secreted under stress) | **+1** (`low` intern.) | **+1** (35.0% P4HB, 93.8% mouse) | **-3** | **#9 Injured-Cell DAMP Marker** (Targets dying Bin 4 cells, not healthy cells) |
| **MOSMO** | **-1** (Peaks in Bin 2 & Bin 4) | **-1** (Developmental/ciliary; no cardiac enrichment) | **-1** (64 aa across 4-TM loops, MEGF8-masked) | **+2** (Clean) | **+1** (`moderate` intern.) | **+2** (No paralogs, 100% mouse) | **-3** | **#10 Ciliary / MEGF8-Masked** |
| **CTSD** | **0** (Mixed Bin 1 & Bin 4) | **-1** (Lysosomal protease) | **-1** (0 TM; receptor-tethered pro-enzyme) | **-2** (High free plasma pro-CTSD decoy) | **+1** (`very_low` intern.) | **+1** (47.5% NAPSA, 83.7% mouse) | **-6** | **#11 Disqualified by Soluble Plasma Sink** |
| **FREM2** | **0** (Bin 1 & Bin 3) | **-2** (Basement membrane ECM; GSC surface only) | **0** (3064 aa, pLDDT 72.8, 6 SBind) | **-1** (Furin shed; ECM-buried) | **+1** (`low` intern.) | **+1** (40.6% FREM3, 89.1% mouse) | **-6** | **#12 Inaccessible Basement Membrane ECM** |
| **TMEM131L**| **-1** (Spikes in Bin 4) | **-2** (Immature thymocyte restricted) | **0** (834 aa, pLDDT 70.3, 39.2% disordered) | **+2** (Clean) | **+1** (`moderate` intern.) | **+1** (26.9% TMEM131, 82.9% mouse) | **-5** | **#13 Thymocyte-Restricted** |
| **ADAM32** | **-1** (Spikes in Bin 4) | **-2** (Sperm/testis restricted) | **+1** (666 aa, pLDDT 82.4, 3 SBind) | **+2** (Clean) | **+1** (`moderate` intern.) | **0** (43.8% ADAM18, 60.9% mouse) | **-5** | **#14 Sperm-Restricted** |
| **CATSPER2**| **0** (Bin 2–3) | **-2** (Sperm flagellum restricted) | **-2** (75 aa, masked by CATSPER$\varepsilon$ canopy) | **+2** (Clean) | **+1** (`moderate` intern.) | **+1** (28% CATSPER3, 70.9% mouse) | **-4** | **#15 Sperm-Restricted & Canopy-Masked** |

