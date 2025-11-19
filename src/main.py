#!/usr/bin/env python3
"""
LORF Finder

Reads nucleotide sequences from a FASTA file, finds Long Open Reading Frames (LORFs)
using a user-selected translation table, and writes amino acid sequences to a FASTA file.
"""

import argparse
from typing import List, Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Find Long Open Reading Frames (LORFs) and output amino acid FASTA."
    )
    parser.add_argument(
        "--input",
        "-i",
        required=True,
        help="Input FASTA file with DNA sequences.",
    )
    parser.add_argument(
        "--output",
        "-o",
        default="lorfs.faa",
        help="Output FASTA file for protein LORFs (default: lorfs.faa).",
    )
    parser.add_argument(
        "--table",
        "-t",
        type=int,
        default=1,
        help=(
            "NCBI translation table ID (default: 1 = Standard Code). "
            "Examples: 1 (Standard), 2 (Vertebrate Mitochondrial), 11 (Bacterial)."
        ),
    )
    parser.add_argument(
        "--min-aa",
        type=int,
        default=100,
        help="Minimum amino acid length to keep an ORF as a LORF (default: 100).",
    )
    return parser.parse_args()


def find_lorfs_in_record(
    record: SeqIO.SeqRecord, table_id: int, min_aa: int
) -> List[Tuple[str, Seq]]:
    """Find Long Open Reading Frames (LORFs) in a DNA sequence."""
    dna: Seq = record.seq.upper()

    try:
        table = CodonTable.unambiguous_dna_by_id[table_id]
    except KeyError:
        raise ValueError(f"Unknown translation table ID: {table_id}")

    start_codons = set(table.start_codons)
    stop_codons = set(table.stop_codons)

    lorfs: List[Tuple[str, Seq]] = []
    seq_len = len(dna)

    for frame in range(3):
        in_orf = False
        nt_start = None

        i = frame
        while i + 3 <= seq_len:
            codon = str(dna[i:i+3])

            if not in_orf:
                if codon in start_codons:
                    in_orf = True
                    nt_start = i
            else:
                if codon in stop_codons:
                    nt_end = i + 3
                    orf_dna = dna[nt_start:nt_end]

                    aa_seq = orf_dna.translate(table=table_id, to_stop=True)

                    if len(aa_seq) >= min_aa:
                        lorf_id = (
                            f"{record.id}|frame={frame}|"
                            f"nt={nt_start+1}-{nt_end}|aa_len={len(aa_seq)}"
                        )
                        lorfs.append((lorf_id, aa_seq))

                    in_orf = False
                    nt_start = None

            i += 3

    return lorfs


def main() -> None:
    args = parse_args()

    input_path = args.input
    output_path = args.output
    table_id = args.table
    min_aa = args.min_aa

    all_lorfs: List[Tuple[str, Seq]] = []

    for record in SeqIO.parse(input_path, "fasta"):
        lorfs = find_lorfs_in_record(record, table_id, min_aa)
        all_lorfs.extend(lorfs)

    with open(output_path, "w") as out:
        for lorf_id, aa_seq in all_lorfs:
            out.write(f">{lorf_id}\n")
            out.write(str(aa_seq) + "\n")


if __name__ == "__main__":
    main()
