# 🧬 ORF Detection & Annotation Pipeline (Prodigal + NCBI)

Welcome! This repository contains a **student-friendly Google Colab pipeline** for:

1. Validating bacterial genome FASTA input,
2. Predicting ORFs with **Prodigal**,
3. Annotating translated proteins against **NCBI Protein (Entrez)**,
4. Exporting a clean `annotated_orfs.csv` for downstream analysis.

---

## 🎯 Learning Goals

By running this notebook script, students will practice:

- Working with FASTA and sequence quality checks,
- Understanding ORF prediction in prokaryotic genomes,
- Reading GFF/protein/CDS outputs from Prodigal,
- Performing practical homology-based annotation using NCBI,
- Building a reproducible mini bioinformatics workflow.

---

## 🧪 What the Pipeline Does

The script `colab_orf_pipeline.py` runs these stages:

| Stage | Action | Output |
|---|---|---|
| 1. Input | Upload FASTA in Colab | local uploaded file |
| 2. QC | Validate IUPAC bases + calculate GC/AT/N/ambiguity metrics | on-screen QC report |
| 3. ORF calling | Run Prodigal (`-p single`) | `.gff`, `.faa`, `.fna` |
| 4. Annotation | Query NCBI Entrez Protein per ORF sequence | accession/gene/organism/product |
| 5. Merge/export | Join ORF coordinates + annotation | `annotated_orfs.csv` |

---

## 🚀 Quick Start (Google Colab)

1. Open a new Colab notebook.
2. Install dependencies:

```bash
!apt-get update -y
!apt-get install -y prodigal
!pip install -q biopython pandas requests
```

3. Upload `colab_orf_pipeline.py` into Colab (left sidebar → Files → Upload).
4. Run the script:

```bash
!python colab_orf_pipeline.py
```

5. Provide:
   - a FASTA file when prompted,
   - your email for NCBI Entrez,
   - optional NCBI API key.

6. Download generated outputs from Colab.

---

## 📦 Output Files

After a successful run, students receive:

- `annotated_orfs.csv` — ORF table with coordinates and annotation,
- `prodigal.faa` — predicted protein sequences,
- `prodigal.fna` — nucleotide CDS sequences,
- `prodigal.gff` — genomic feature coordinates.

---

## 🧭 Suggested Classroom Flow

- **Part A (10–15 min):** Explain ORFs, codons, and prokaryotic gene prediction.
- **Part B (20–30 min):** Run the pipeline on a small bacterial genome.
- **Part C (20 min):** Inspect `annotated_orfs.csv` and discuss hypothetical proteins vs known products.
- **Part D (optional):** Compare annotation yield with/without API key and different genomes.

---

## 💡 Tips for Students

- Use high-quality FASTA inputs to reduce "No hit" annotations.
- NCBI requests are rate-limited; an API key generally improves throughput.
- Very short proteins may not yield reliable hits.
- Treat annotation as evidence-based inference, not absolute truth.

---

## 🧑‍🔬 Citation / Credit

If you use this in teaching materials, please cite:

- Prodigal: Hyatt et al., 2010.
- NCBI Entrez Programming Utilities (E-utilities).
- This repository/workflow.

---

## 🤝 Contributions

Issues and pull requests are welcome—especially improvements that make this workflow clearer for students and instructors.
