# --- Colab Cell 1: Dependencies ---
# Run this cell in Google Colab before executing the pipeline:
# !apt-get update -y
# !apt-get install -y prodigal
# !pip install -q biopython pandas requests

"""Google Colab pipeline for bacterial ORF detection with Prodigal and NCBI annotation."""

from __future__ import annotations

import io
import logging
import subprocess
import time
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import pandas as pd
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from google.colab import files


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)
logger = logging.getLogger(__name__)

VALID_IUPAC = frozenset("ATGCRYSWKMBDHVN")
AMBIGUOUS_BASES = frozenset("RYSWKMBDHVN")


@dataclass(frozen=True)
class SequenceMetrics:
    """Container for per-sequence and aggregate QC metrics."""

    sequence_id: str
    length_bp: int
    gc_percent: float
    at_percent: float
    n_percent: float
    ambiguity_percent: float
    ambiguous_counts: Dict[str, int]


def upload_fasta() -> Path:
    """Upload a FASTA file through Colab and return its local path."""
    logger.info("Upload a FASTA file to begin.")
    uploaded = files.upload()
    if not uploaded:
        raise ValueError("No file uploaded.")
    fasta_name = next(iter(uploaded.keys()))
    return Path(fasta_name)


def validate_fasta(fasta_path: Path) -> Tuple[List[Tuple[str, str]], List[SequenceMetrics]]:
    """Validate FASTA structure and IUPAC content; return cleaned sequences and metrics."""
    cleaned_records: List[Tuple[str, str]] = []
    metrics: List[SequenceMetrics] = []

    with fasta_path.open("r", encoding="utf-8") as handle:
        records = list(SeqIO.parse(handle, "fasta"))

    if not records:
        raise ValueError("Invalid FASTA: no records detected.")

    for record in records:
        seq_id = record.id.strip()
        raw_seq = str(record.seq).upper()
        clean_seq = "".join(raw_seq.split())
        if not clean_seq:
            raise ValueError(f"Record '{seq_id}' has an empty sequence.")

        char_counts = Counter(clean_seq)
        invalid_chars = sorted(set(char_counts) - VALID_IUPAC)
        if invalid_chars:
            raise ValueError(
                f"Record '{seq_id}' contains non-IUPAC characters: {','.join(invalid_chars)}"
            )

        length_bp = len(clean_seq)
        gc_count = char_counts.get("G", 0) + char_counts.get("C", 0)
        at_count = char_counts.get("A", 0) + char_counts.get("T", 0)
        n_count = char_counts.get("N", 0)
        ambiguous_counts = {base: char_counts.get(base, 0) for base in sorted(AMBIGUOUS_BASES)}
        ambiguous_total = sum(ambiguous_counts.values())

        seq_metrics = SequenceMetrics(
            sequence_id=seq_id,
            length_bp=length_bp,
            gc_percent=(gc_count / length_bp) * 100,
            at_percent=(at_count / length_bp) * 100,
            n_percent=(n_count / length_bp) * 100,
            ambiguity_percent=(ambiguous_total / length_bp) * 100,
            ambiguous_counts=ambiguous_counts,
        )
        metrics.append(seq_metrics)
        cleaned_records.append((seq_id, clean_seq))

    return cleaned_records, metrics


def compute_gc_content(sequence: str) -> float:
    """Compute GC content percentage for a cleaned DNA sequence."""
    if not sequence:
        return 0.0
    counts = Counter(sequence)
    return ((counts.get("G", 0) + counts.get("C", 0)) / len(sequence)) * 100


def write_clean_fasta(records: Sequence[Tuple[str, str]], out_path: Path) -> None:
    """Write cleaned FASTA records to disk."""
    with out_path.open("w", encoding="utf-8") as handle:
        for seq_id, seq in records:
            handle.write(f">{seq_id}\n")
            for i in range(0, len(seq), 80):
                handle.write(seq[i : i + 80] + "\n")


def run_prodigal(input_fasta: Path, output_prefix: Path) -> Tuple[Path, Path, Path]:
    """Run Prodigal in single-genome mode and return output file paths."""
    gff_path = output_prefix.with_suffix(".gff")
    proteins_path = output_prefix.with_suffix(".faa")
    cds_path = output_prefix.with_suffix(".fna")

    cmd = [
        "prodigal",
        "-i",
        str(input_fasta),
        "-p",
        "single",
        "-f",
        "gff",
        "-o",
        str(gff_path),
        "-a",
        str(proteins_path),
        "-d",
        str(cds_path),
    ]
    logger.info("Running Prodigal...")
    subprocess.run(cmd, check=True)
    return gff_path, proteins_path, cds_path


def _extract_attribute(attributes: str, key: str) -> Optional[str]:
    for item in attributes.split(";"):
        if item.startswith(f"{key}="):
            return item.split("=", 1)[1]
    return None


def _codons_from_genome(
    genome_lookup: Dict[str, str],
    seqid: str,
    start: int,
    end: int,
    strand: str,
) -> Tuple[str, str]:
    seq = genome_lookup[seqid]
    if strand == "+":
        start_codon = seq[start - 1 : start + 2]
        stop_codon = seq[end - 3 : end]
    else:
        start_codon = str(Seq(seq[end - 3 : end]).reverse_complement())
        stop_codon = str(Seq(seq[start - 1 : start + 2]).reverse_complement())
    return start_codon, stop_codon


def parse_gff(gff_path: Path, genome_records: Sequence[Tuple[str, str]]) -> pd.DataFrame:
    """Parse Prodigal GFF and extract ORF coordinates, strand, lengths, and codons."""
    genome_lookup = {sid: seq for sid, seq in genome_records}
    rows: List[Dict[str, object]] = []

    with gff_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            seqid, source, feature, start, end, score, strand, phase, attrs = fields
            if feature.lower() != "cds":
                continue

            start_i = int(start)
            end_i = int(end)
            length = end_i - start_i + 1
            orf_id = _extract_attribute(attrs, "ID") or f"{seqid}_{start}_{end}_{strand}"
            start_codon, stop_codon = _codons_from_genome(
                genome_lookup=genome_lookup,
                seqid=seqid,
                start=start_i,
                end=end_i,
                strand=strand,
            )

            rows.append(
                {
                    "orf_id": orf_id,
                    "contig": seqid,
                    "source": source,
                    "start": start_i,
                    "end": end_i,
                    "strand": strand,
                    "phase": phase,
                    "orf_length_bp": length,
                    "start_codon": start_codon,
                    "stop_codon": stop_codon,
                }
            )

    return pd.DataFrame(rows)


def _batch(iterable: Sequence[Tuple[str, str]], size: int) -> Iterable[Sequence[Tuple[str, str]]]:
    for i in range(0, len(iterable), size):
        yield iterable[i : i + size]


def annotate_ncbi(
    proteins_faa: Path,
    email: str,
    api_key: Optional[str] = None,
    batch_size: int = 20,
    retries: int = 3,
    delay_s: float = 0.35,
) -> pd.DataFrame:
    """Best-effort batch annotation against NCBI protein records via Entrez."""
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    protein_records = [(rec.id, str(rec.seq)) for rec in SeqIO.parse(str(proteins_faa), "fasta")]
    if not protein_records:
        return pd.DataFrame(columns=["orf_id", "accession", "gene", "organism", "product"])

    annotations: List[Dict[str, str]] = []

    for group in _batch(protein_records, batch_size):
        terms = [f'"{seq}"[Sequence]' for _, seq in group if len(seq) >= 20]
        if not terms:
            continue
        search_term = " OR ".join(terms)

        id_list: List[str] = []
        for attempt in range(1, retries + 1):
            try:
                with Entrez.esearch(db="protein", term=search_term, retmax=batch_size * 5) as h:
                    search_res = Entrez.read(h)
                id_list = search_res.get("IdList", [])
                break
            except Exception as exc:  # noqa: BLE001
                logger.warning("NCBI esearch failed (attempt %s/%s): %s", attempt, retries, exc)
                if attempt == retries:
                    id_list = []
                time.sleep(delay_s * attempt)

        if not id_list:
            for protein_id, _ in group:
                annotations.append(
                    {
                        "orf_id": protein_id,
                        "accession": "",
                        "gene": "",
                        "organism": "",
                        "product": "No hit",
                    }
                )
            continue

        summaries: List[dict] = []
        for attempt in range(1, retries + 1):
            try:
                with Entrez.esummary(db="protein", id=",".join(id_list)) as h:
                    summary_res = Entrez.read(h)
                summaries = summary_res if isinstance(summary_res, list) else summary_res.get("DocumentSummarySet", {}).get("DocumentSummary", [])
                break
            except Exception as exc:  # noqa: BLE001
                logger.warning("NCBI esummary failed (attempt %s/%s): %s", attempt, retries, exc)
                if attempt == retries:
                    summaries = []
                time.sleep(delay_s * attempt)

        fallback = {
            "accession": "",
            "gene": "",
            "organism": "",
            "product": "No hit",
        }
        best_hit = fallback
        if summaries:
            first = summaries[0]
            best_hit = {
                "accession": str(first.get("AccessionVersion", first.get("Caption", ""))),
                "gene": str(first.get("Title", "")).split("[")[0].strip(),
                "organism": str(first.get("TaxName", "")),
                "product": str(first.get("Title", "")),
            }

        for protein_id, _ in group:
            annotations.append({"orf_id": protein_id, **best_hit})

        time.sleep(delay_s)

    return pd.DataFrame(annotations)


def print_qc_report(metrics: Sequence[SequenceMetrics]) -> None:
    """Print QC report including ambiguity and GC-range warnings."""
    print("\n=== FASTA QC REPORT ===")
    for item in metrics:
        print(f"- Sequence: {item.sequence_id}")
        print(f"  Length (bp): {item.length_bp}")
        print(f"  GC (%): {item.gc_percent:.2f}")
        print(f"  AT (%): {item.at_percent:.2f}")
        print(f"  N (%): {item.n_percent:.2f}")
        print(f"  Ambiguity (%): {item.ambiguity_percent:.2f}")
        print(f"  Ambiguous counts: {item.ambiguous_counts}")
        if not (20 <= item.gc_percent <= 75):
            print("  WARNING: GC content is outside the typical bacterial range (20-75%).")

    total_bp = sum(m.length_bp for m in metrics)
    weighted_gc = sum((m.gc_percent / 100) * m.length_bp for m in metrics) / total_bp * 100
    weighted_at = sum((m.at_percent / 100) * m.length_bp for m in metrics) / total_bp * 100
    weighted_n = sum((m.n_percent / 100) * m.length_bp for m in metrics) / total_bp * 100
    print("\n=== GENOME SUMMARY ===")
    print(f"Genome size (bp): {total_bp}")
    print(f"GC (%): {weighted_gc:.2f}")
    print(f"AT (%): {weighted_at:.2f}")
    print(f"N (%): {weighted_n:.2f}")


def main() -> None:
    """Run the complete Colab ORF detection and annotation workflow."""
    fasta_path = upload_fasta()
    records, metrics = validate_fasta(fasta_path)
    print_qc_report(metrics)

    cleaned_fasta = Path("cleaned_input.fasta")
    write_clean_fasta(records, cleaned_fasta)

    output_prefix = Path("prodigal")
    gff_path, proteins_path, cds_path = run_prodigal(cleaned_fasta, output_prefix)

    orf_df = parse_gff(gff_path, records)

    email = input("Enter email for NCBI Entrez queries: ").strip()
    api_key = input("Optional NCBI API key (press Enter to skip): ").strip() or None
    ann_df = annotate_ncbi(proteins_path, email=email, api_key=api_key)

    final_df = orf_df.merge(ann_df, on="orf_id", how="left")
    final_csv = Path("annotated_orfs.csv")
    final_df.to_csv(final_csv, index=False)

    print("\n=== ORF SUMMARY ===")
    print(f"Detected ORFs: {len(orf_df)}")
    print(f"Annotated ORFs table: {final_csv}")
    print(f"Protein FASTA: {proteins_path}")
    print(f"CDS FASTA: {cds_path}")
    print(f"GFF: {gff_path}")

    for file_path in [final_csv, proteins_path, cds_path, gff_path]:
        files.download(str(file_path))


if __name__ == "__main__":
    main()
