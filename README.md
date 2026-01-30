# Plasmid Builder - Universal Plasmid Constructor

A computational tool for automated plasmid design and construction using bioinformatics approaches.

**Author:** Tamme Mokshagna (2023BB11055)  
**Course:** BBL434-Assignment 1

---

## 📋 Overview

This Python-based tool automatically constructs plasmid sequences by:

1. **Identifying the origin of replication (ORI)** from an input genome using GC-skew analysis
2. **Assembling plasmid components** including resistance genes, reporter genes, and restriction sites
3. **Optimizing the sequence** by removing unwanted restriction sites while preserving the multiple cloning site (MCS)

The tool uses computational biology techniques including Position Weight Matrices (PWM), GC-skew analysis, and motif discovery to create functional plasmid designs.

---

## 🧬 Biological Background

### What is a Plasmid?

A plasmid is a small, circular DNA molecule that replicates independently of chromosomal DNA. Key components include:

- **Origin of Replication (ORI):** Enables autonomous replication in host cells
- **Selection Markers:** Antibiotic resistance genes (e.g., Ampicillin, Kanamycin)
- **Reporter Genes:** Genes for screening (e.g., lacZα for blue-white selection)
- **Multiple Cloning Site (MCS):** Contains restriction enzyme recognition sites for gene insertion

### How This Tool Works

The tool mimics real molecular cloning workflows by:
- Extracting a functional ORI from any bacterial genome
- Adding user-specified genetic markers
- Creating a clean backbone free of interfering restriction sites
- Preserving essential cloning sites for downstream applications

---

## 🛠️ Features

### Core Capabilities

- **Automatic ORI Detection:** Uses GC-skew analysis to locate replication origins
- **Motif Discovery:** Implements Position Weight Matrix (PWM) for identifying conserved motifs
- **Flexible Design:** Accepts custom plasmid designs via simple text files
- **Sequence Optimization:** Removes unwanted restriction sites from functional regions
- **Modular Architecture:** Extensible class-based design

### Supported Restriction Enzymes

The tool recognizes 14 common restriction enzymes:

| Enzyme | Recognition Site |
|--------|------------------|
| EcoRI | GAATTC |
| BamHI | GGATCC |
| HindIII | AAGCTT |
| PstI | CTGCAG |
| KpnI | GGTACC |
| SacI | GAGCTC |
| SalI | GTCGAC |
| SphI | GCATGC |
| XbaI | TCTAGA |
| NotI | GCGGCCGC |
| SmaI | CCCGGG |
| BsaI | GGTCTC |
| BbsI | GAAGAC |
| BsmBI | CGTCTC |

---



## 🚀 Usage

### Basic Command

```bash
python plasmid_builder.py <genome.fa> <design.txt>
```




## 🔬 Algorithm Details

### 1. ORI Identification

The tool uses **GC-skew analysis** to locate the origin of replication:

```
GC-skew = (G - C) / (G + C)
```

**Steps:**
1. Calculate cumulative GC-skew across the genome
2. Find the minimum point (indicates ORI location)
3. Extract a 500 bp window centered on this position

**Why GC-Skew?**
- Leading strand accumulates G nucleotides during replication
- Lagging strand accumulates C nucleotides
- This creates a characteristic skew pattern with minimum at the ORI

### 2. Motif Discovery

Uses Position Weight Matrices (PWM) to identify conserved motifs:

**Process:**
1. Scan sequence with k-mers (lengths 6-15 bp)
2. Filter low-complexity sequences (must have >2 unique bases)
3. Build PWM from motif instances with pseudocounts
4. Calculate log-likelihood scores against background model
5. Select optimal k-mer length based on highest score

**PWM Formula:**
```
Score = Σ log(P(base|motif) / P(base|background))
```

### 3. Sequence Sanitization

**Goal:** Remove unwanted restriction sites from functional regions while preserving MCS

**Strategy:**
1. Assemble backbone (ORI + resistance genes + reporter genes)
2. Remove ALL known restriction sites from backbone
3. Add MCS restriction sites at the end (preserved)

This ensures:
- No accidental cutting during cloning
- MCS sites remain functional for gene insertion

---

## 📊 Example Workflow

### Input Setup

**Genome:** `pUC19.fa` (E. coli plasmid)

**Design:** `Design_pUC19.txt`
```
EcoRI_site, EcoRI
BamHI_site, BamHI
HindIII_site, HindIII
AmpR_gene, Ampicillin
lacZ_alpha, Blue_White_Selection
```

### Execution

```bash
python plasmid_builder.py pUC19.fa Design_pUC19.txt
```

### Console Output

```
═══ ORIGIN IDENTIFICATION ═══
→ Central position: 1247
→ Extracted region: 997 to 1497
→ Optimal k-mer: 9

✓ Plasmid sequence saved to Output.fa
✓ Total length: 3842 bp
```

### Result

A plasmid containing:
- 500 bp ORI region
- Ampicillin resistance gene
- lacZα reporter gene
- MCS with EcoRI, BamHI, and HindIII sites
- All other restriction sites removed

---

## 🏗️ Code Architecture

### Class Structure

```
RestrictionEnzymeDatabase
├── ENZYMES (dict)
├── get_site(enzyme_name)
└── all_enzymes()

SequenceAnalyzer
├── calculate_nucleotide_frequencies(sequence)
├── assess_sequence_complexity(kmer)
├── calculate_cumulative_gc_skew(sequence)
└── locate_replication_origin_position(sequence)

MotifPWM
├── __init__(motif_instances, motif_length)
├── _construct_matrix(instances)
└── calculate_sequence_score(sequence, background_probs)

OriginFinder
├── discover_motifs(sequence, kmer_size)
├── select_optimal_kmer_length(sequence)
└── identify_origin_region(fasta_path, window_size)

PlasmidAssembler
├── load_genetic_marker(marker_identifier)
├── parse_plasmid_design(design_filepath)
├── sanitize_backbone_sequence(sequence, enzyme_list)
└── construct_plasmid(input_genome, design_spec)
```



## ⚙️ Configuration

### Adjustable Parameters

#### ORI Window Size
Default: 500 bp

Modify in `identify_origin_region()`:
```python
origin_data = ori_finder.identify_origin_region(input_genome, window_size=500)
```

#### K-mer Range
Default: 6-15 bp

Modify in `OriginFinder.__init__()`:
```python
def __init__(self, k_range=(6, 15)):
```

#### Pseudocount Values
Default: 1 for each nucleotide

Modify in `MotifPWM._construct_matrix()`:
```python
nucleotide_counts = {"A": 1, "C": 1, "G": 1, "T": 1}
```

---

## 🧪 Testing

### Test Case: pUC19 Reconstruction

The repository includes `test_pUC19.py` for validation:

```bash
python test_pUC19.py
```

**Expected Output:**
```
Test passed: EcoRI successfully removed from backbone.
```

### Manual Verification

Check for restriction sites:
```python
from Bio import SeqIO

record = SeqIO.read("Output.fa", "fasta")
sequence = str(record.seq)

# Should return empty list (except MCS sites at end)
print(sequence.count("GAATTC"))  # EcoRI
```

---

## ⚠️ Limitations

### Biological Limitations

1. **ORI Functionality:** Computationally identified ORIs may not function in vivo
   - GC-skew provides statistical estimate, not experimental validation
   - Chromosomal ORIs may not work as plasmid ORIs

2. **No Regulatory Elements:** Tool does not include:
   - Promoters for marker gene expression
   - Terminators for transcription control
   - Ribosome binding sites (RBS)

3. **No Circularization:** Output is linear sequence
   - Real plasmids are circular
   - Requires additional processing for true plasmid function

### Technical Limitations

1. **Simple Sanitization:** Removes restriction sites by deletion
   - May disrupt gene function
   - Better approach: silent mutations

2. **No Size Optimization:** Does not consider:
   - Transformation efficiency (smaller plasmids transform better)
   - Replication stability

3. **Limited Enzyme Support:** Only 14 restriction enzymes
   - Many Type IIS enzymes not included
   - No support for methylation-sensitive enzymes

---

## 📝 File Structure

```
BBL434-Assignment-1/
├── plasmid_builder.py      # Main program
├── test_pUC19.py           # Test suite
├── markers/                # Genetic marker sequences
│   ├── Ampicillin.fa
│   ├── Kanamycin.fa
│   └── Chloramphenicol.fa
├── pUC19.fa               # Example genome
├── Design_pUC19.txt       # Example design
├── Output.fa              # Generated plasmid (after running)
├── README.md              # This file
└── LICENSE                # MIT License
```

---


##  Acknowledgments

- Course instructor for providing marker sequences and specifications
- Biopython developers for sequence analysis tools
- New England Biolabs for restriction enzyme data

---

