# Universal Plasmid Construction Pipeline

## 📌 Overview

This repository implements a **computational plasmid construction pipeline** that generates a plasmid DNA sequence from:

* An **input genome sequence** (`.fa`)
* A **design specification file** (`Design.txt`)
* A library of **biological markers** (antibiotic genes, screening genes)

The pipeline:

* Automatically detects the **origin of replication (ORI)** from the input genome using GC-skew.
* Assembles a plasmid using the detected ORI and user-specified components.
* Removes **all known restriction enzyme sites** from the final plasmid.

---

## 🧬 Biological Background

In real molecular cloning:

* A plasmid requires an **origin of replication (ORI)** to replicate in a host.
* Selection and screening are done using:

  * Antibiotic resistance genes (AmpR, KanR, etc.)
  * Reporter genes (lacZα, GFP, etc.)
* Restriction enzymes define cloning sites.

This project simulates these principles **computationally**.

---

## 📁 Repository Structure

```
Assignment-1/
├── main.py                # Entry point
├── ori_finder.py          # Multi-scale ORI detection (GC-skew)
├── plasmid_builder.py     # Plasmid construction logic
├── restriction_sites.py   # Known restriction enzyme motifs
├── markers/               # Marker gene FASTA files
│   ├── Ampicillin.fa
│   ├── Kanamycin.fa
│   └── Chloramphenicol.fa
├── markers.tab            # Supported markers list
├── Design_pUC19.txt       # Example design file
├── pUC19.fa               # Test genome (E. coli plasmid)
└── tests/
    └── test_pUC19.py      # Automated test
```

---

## 🧾 Input Files

### 1️⃣ Genome Input (`input.fa`)

A FASTA file of any bacterial genome or plasmid.

Example:

```fasta
>Genome
ATGCGTAGCTAGCTAG...
```

---

### 2️⃣ Design File (`Design.txt`)

Each line specifies a plasmid component:

```
Label, Value
```

Example:

```
BamHI_site, BamHI
HindIII_site, HindIII
AmpR_gene, Ampicillin
lacZ_alpha, Blue_White_Selection
```

Rules:

* Restriction sites → must exist in `restriction_sites.py`
* Marker genes → must exist as FASTA in `markers/`
* Unknown markers are **skipped with warning**
* Unknown enzymes are **skipped with warning**

---

## 🧠 ORI Detection Algorithm

The ORI is detected using:

* GC-skew
* Sliding window
* Multi-scale consensus

Window sizes:

```
150 bp, 200 bp, 250 bp, 300 bp
```

Step size:

```
10 bp
```

Final ORI is selected using the **median start position** across scales.

This provides:

* Noise resistance
* No arbitrary ORI length assumption
* Robust inference

---

## ⚙️ How to Run

### Install dependency

```bash
pip install biopython
```

---

### Run full pipeline

```bash
python main.py input.fa Design.txt
```

Output:

```
Output.fa
```

This contains the final plasmid sequence.

---

### Run automated test

```bash
python -m tests.test_pUC19
```

Expected output:

```
Test passed: EcoRI successfully removed.
```

---

## 🔬 What the Pipeline Does

1. Reads genome FASTA
2. Detects ORI using GC-skew
3. Appends marker genes
4. Appends restriction motifs
5. Removes all known restriction sites globally
6. Outputs final plasmid

---

## ⚠️ Known Limitations

* Chromosomal ORI may not function as plasmid ORI in real biology.
* GC-skew only provides statistical ORI estimate.
* No transcriptional regulation included.
* No plasmid circularization simulation.

---

## 🧪 Testing Strategy

Test case uses:

* Input genome: `pUC19.fa`
* Design: `Design_pUC19.txt`

Expected behavior:

* EcoRI removed from final sequence
* ORI detected automatically
* Missing markers skipped safely
