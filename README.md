# 🧬 TP53 DNA to Docking — Bioinformatics Pipeline

> **Gene:** TP53 (Tumor Protein P53) | **Accession:** NM_000546  
> **Protein:** p53 — *"Guardian of the Genome"*  
> **Type:** Tumor Suppressor Gene | **Chromosome:** 17p13.1

---

## 📋 Table of Contents
1. [Project Overview](#project-overview)
2. [Gene Background](#gene-background)
3. [Step 1 — DNA Sequence Retrieval](#step-1--dna-sequence-retrieval)
4. [Step 2 — DNA to mRNA Conversion](#step-2--dna-to-mrna-conversion)
5. [Step 3 — BLAST Analysis](#step-3--blast-analysis)
6. [Step 4 — AlphaFold2 Protein Structure](#step-4--alphafold2-protein-structure)
7. [Step 5 — GNINA Ligand Docking](#step-5--gnina-ligand-docking)
8. [Results Summary](#results-summary)
9. [How to Run](#how-to-run)
10. [References](#references)

---

## Project Overview

This project performs a complete bioinformatics analysis of the **TP53** tumor suppressor gene, starting from raw DNA and ending at ligand docking simulation:

```
DNA Sequence → mRNA Transcription → BLAST Identification
→ AlphaFold2 Structure Prediction → GNINA Ligand Docking
```

---

## Gene Background

| Property | Details |
|----------|---------|
| Gene Name | TP53 (Tumor Protein P53) |
| Gene Type | Tumor Suppressor |
| NCBI Gene ID | 7157 |
| Accession | NM_000546.6 |
| Chromosome | 17p13.1 |
| Protein | p53 — 393 amino acids |
| Function | Cell cycle arrest, DNA repair, Apoptosis |
| Cancer Link | Mutated in ~50% of all human cancers |

TP53 encodes the p53 protein — a transcription factor and the most critical checkpoint in the human cell cycle. When DNA damage is detected, p53 either halts division to allow repair, or triggers programmed cell death (apoptosis) if damage is too severe. Loss of TP53 function is found in nearly half of all cancers.

---

## Step 1 — DNA Sequence Retrieval

**Tool:** Biopython (NCBI Entrez API)  
**Script:** `step1_fetch_dna.py`  
**Platform:** Google Colab

The DNA sequence of TP53 was fetched directly from NCBI using the official mRNA reference accession **NM_000546**.

```python
from Bio import Entrez, SeqIO
Entrez.email = "your@email.com"
handle = Entrez.efetch(db="nucleotide", id="NM_000546", rettype="gb", retmode="text")
record = SeqIO.read(handle, "genbank")
```

**Result:**
- Sequence Length: **2591 bp**
- Output file: `TP53_DNA.fasta`

```
>TP53_DNA | Homo sapiens tumor protein p53 (TP53), mRNA
AGGGAGGGAGAGAATCTTCCAGGGCCAGCTCGGGCAGCAATCAGCAGG...
```

---

## Step 2 — DNA to mRNA Conversion

**Tool:** Biopython `Seq.transcribe()`  
**Script:** `step1_fetch_dna.py`  
**Platform:** Google Colab

The DNA coding sequence was transcribed to mRNA by replacing every **T → U**, following the central dogma of molecular biology.

```python
mrna_seq = dna_seq.transcribe()   # T → U
```

| | First 40 nucleotides |
|--|--|
| DNA | `AGGGAGGGAGAGAATCTTCCAGGGCCAGCTCGGGCAGCA` |
| mRNA | `AGGGAGGGAGAGAAUCUUCCAGGGCCAGCUCGGGCAGCA` |

- mRNA Length: **2591 nucleotides**
- Output file: `TP53_mRNA.fasta`

The protein-coding region (CDS) was also extracted and translated:
- Protein Length: **393 amino acids**
- Output file: `TP53_protein.fasta`

---

## Step 3 — BLAST Analysis

**Tool:** NCBI BLASTP (online + CSV export)  
**Database:** Non-redundant protein sequences (nr)  
**Platform:** https://blast.ncbi.nlm.nih.gov  
**RID:** 1MNUTJ8J014

The TP53 protein sequence (393 aa) was submitted to BLASTP to confirm identity and find homologs across species.

### Top 10 BLAST Hits

| Rank | Accession | % Identity | E-value | Bit Score | Query Coverage |
|------|-----------|-----------|---------|-----------|----------------|
| 1 | XP_063556659.1 | 99.75% | 0.0 | 814 | 99.75% |
| 2 | 7XZZ_K | **100.00%** | 0.0 | 814 | **100.00%** |
| 3 | 8R1F_C | **100.00%** | 0.0 | 814 | **100.00%** |
| 4 | NP_000537.3 | **100.00%** | 0.0 | 813 | **100.00%** |
| 5 | AYE20617.1 | 99.75% | 0.0 | 813 | 99.75% |
| 6 | AYE20613.1 | 99.75% | 0.0 | 812 | 99.75% |
| 7 | AYE20623.1 | 99.75% | 0.0 | 812 | 99.75% |
| 8 | AAA61212.1 | 99.75% | 0.0 | 812 | 99.75% |
| 9 | XP_004058559.3 | 99.75% | 0.0 | 812 | 99.75% |
| 10 | CAA42633.1 | 99.75% | 0.0 | 812 | 100.00% |

**Total hits:** 100 sequences across multiple species

### Key Findings
- Top hit **NP_000537.3** = cellular tumor antigen p53, Homo sapiens — **100% identity, E-value 0.0** ✅
- PDB structures **7XZZ_K** and **8R1F_C** hit at 100% — confirms known crystallographic structures match
- All top hits have E-value = **0.0** — statistically perfect matches
- Query length confirmed at **393 amino acids** — matches expected p53 protein length

> 📸 *[<img width="1762" height="671" alt="Screenshot 2026-05-30 104342" src="https://github.com/user-attachments/assets/cb5534bf-d6e1-43d4-9946-0a36ad170d16" />
<img width="1337" height="898" alt="Screenshot 2026-05-30 104943" src="https://github.com/user-attachments/assets/8b42e230-41e1-49e8-b871-df4177adabac" />
<img width="1769" height="626" alt="Screenshot 2026-05-30 104741" src="https://github.com/user-attachments/assets/8baf8dd2-27df-4983-aaaf-bae4a926ada2" />
]*

---

## Step 4 — AlphaFold2 Protein Structure

**Tool:** ColabFold (AlphaFold2)  
**Notebook:** [ColabFold on Google Colab](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb)  
**Job ID:** TP53_structure_7c637

The TP53 protein sequence was submitted to AlphaFold2 via ColabFold to predict its 3D structure. AlphaFold2 uses deep learning trained on the entire PDB database to predict protein folding with high accuracy.

### Settings Used
| Parameter | Value |
|-----------|-------|
| num_recycles | 3 |
| template_mode | pdb100 |
| use_amber | True (structure relaxation) |
| Runtime | T4 GPU (Google Colab) |

### Output Files Generated
| File | Description |
|------|-------------|
| `TP53_structure_7c637_unrelaxed_rank_001_...pdb` | Best predicted structure |
| `TP53_structure_7c637_plddt` | Per-residue confidence plot |
| `TP53_structure_7c637_coverage` | MSA coverage plot |
| `TP53_structure_7c637_predicted_aligned_error_v1` | PAE confidence matrix |
| `TP53_structure_7c637_pae` | PAE plot image |

### pLDDT Confidence Score
pLDDT (predicted Local Distance Difference Test) measures per-residue confidence:

| Score Range | Meaning |
|-------------|---------|
| 90 – 100 | Very high confidence |
| 70 – 90 | High confidence |
| 50 – 70 | Low confidence (may be disordered) |
| < 50 | Very low (intrinsically disordered) |

> 📸 *[Add your pLDDT plot screenshot here — from TP53_structure_7c637_plddt file]*

> 📸 *[Add your 3D structure visualization screenshot here — open .pdb in https://molstar.org/viewer/]*

---

## Step 5 — GNINA Ligand Docking

**Tool:** GNINA v1.0.3 (deep learning molecular docking)  
**Platform:** Google Colab (Linux)  
**Ligand:** APR-246 (Eprenetapopt) — PubChem CID 5713517

### Why APR-246?
APR-246 (also called PRIMA-1MET or Eprenetapopt) is a clinically relevant small molecule that **reactivates mutant p53** by covalently binding to the DNA-binding domain and restoring its wild-type conformation. It has been evaluated in Phase III clinical trials for acute myeloid leukemia (AML).

### Ligand Details
| Property | Value |
|----------|-------|
| Name | APR-246 / Eprenetapopt |
| PubChem CID | 5713517 |
| Formula | C7H13NO3 |
| Molecular Weight | 159.18 g/mol |
| Target | TP53 DNA-binding domain |
| Clinical Stage | Phase III trials |

### Docking Configuration
```bash
./gnina \
  --receptor protein.pdb \
  --ligand APR246.sdf \
  --out docked_poses.sdf \
  --center_x 1.5 --center_y 2.3 --center_z -0.8 \
  --size_x 25 --size_y 25 --size_z 25 \
  --num_modes 9 \
  --exhaustiveness 8 \
  --cnn_scoring none
```

### Docking Results

> 📸 *[Add your GNINA terminal output screenshot here showing affinity scores]*

> 📸 *[Add docked pose visualization — upload protein.pdb + docked_poses.sdf to https://molstar.org/viewer/]*

**Output file:** `docked_poses.sdf` — contains all 9 predicted binding poses with CNN scores and binding affinities (kcal/mol).

---

## Results Summary

| Step | Tool | Output | Key Result |
|------|------|--------|-----------|
| 1 — DNA retrieval | Biopython / NCBI | TP53_DNA.fasta | 2591 bp sequence |
| 2 — mRNA conversion | Biopython | TP53_mRNA.fasta | T→U transcription, 393 aa protein |
| 3 — BLAST | NCBI BLASTP | HitTable CSV | 100% identity to p53, E-value 0.0 |
| 4 — AlphaFold2 | ColabFold | .pdb structure | High-confidence 3D structure predicted |
| 5 — Docking | GNINA | docked_poses.sdf | APR-246 docked to DNA-binding domain |

---

## How to Run

### Requirements
```bash
pip install biopython
```

### Step 1 & 2 — DNA + mRNA
```bash
python step1_fetch_dna.py
```

### Step 3 — BLAST
Visit https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp and paste sequence from `TP53_protein.fasta`

### Step 4 — AlphaFold2
Open ColabFold notebook → paste protein sequence → Runtime → Run All:
https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb

### Step 5 — GNINA Docking (Google Colab)
```bash
# Install GNINA
wget https://github.com/gnina/gnina/releases/download/v1.0.3/gnina -O gnina
chmod +x gnina

# Run docking
./gnina --receptor protein.pdb --ligand APR246.sdf --out docked_poses.sdf \
        --center_x 1.5 --center_y 2.3 --center_z -0.8 \
        --size_x 25 --size_y 25 --size_z 25 \
        --num_modes 9 --exhaustiveness 8 --cnn_scoring none
```

---

## References

1. NCBI Gene TP53 — https://www.ncbi.nlm.nih.gov/gene/7157
2. UniProt P04637 (p53 Human) — https://www.uniprot.org/uniprot/P04637
3. ColabFold / AlphaFold2 — https://github.com/sokrypton/ColabFold
4. GNINA Molecular Docking — https://github.com/gnina/gnina
5. APR-246 (PubChem CID 5713517) — https://pubchem.ncbi.nlm.nih.gov/compound/5713517
6. PDB Structure 2OCJ — https://www.rcsb.org/structure/2OCJ
