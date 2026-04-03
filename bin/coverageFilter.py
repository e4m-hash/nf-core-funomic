#!/usr/bin/env python3

import argparse
import re
import pandas as pd


SAM_MANDATORY_COLUMNS = [
    "QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR",
    "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL"
]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Get read names that satisfy minimum alignment coverage and read length."
    )
    parser.add_argument("-i", "--input", dest="input_path", required=True, help="Input SAM file")
    parser.add_argument("-o", "--output", dest="output_path", required=True, help="Output file")
    parser.add_argument("-c", "--coverage", type=float, default=80.0, help="Min coverage %% (default: 80.0)")
    parser.add_argument("-m", "--min_length", type=int, default=60, help="Min read length in bp (default: 60)")
    return parser.parse_args()


def cigar_to_metrics(cigar):
    if pd.isna(cigar) or str(cigar).strip() in ("*", ""):
        return 0.0, 0

    cigar_str = str(cigar).strip()
    tokens = re.findall(r"(\d+)([MIDNSHP=X])", cigar_str)
    
    matched_bases = 0
    query_length = 0
    
    for length_str, op in tokens:
        length = int(length_str)
        if op in ("M", "="):
            matched_bases += length
        if op in ("M", "I", "S", "=", "X"):
            query_length += length
    
    coverage = (matched_bases / query_length * 100) if query_length > 0 else 0.0
    return coverage, query_length


def main():
    args = parse_args()
    
    q30sam = pd.read_csv(
        args.input_path,
        sep="\t",
        comment="@",
        header=None,
        names=SAM_MANDATORY_COLUMNS,
        usecols=range(11),
        on_bad_lines="skip",
        dtype=str,
    )
    
    metrics = q30sam["CIGAR"].apply(cigar_to_metrics)
    q30sam["coverage"] = metrics.str[0]
    q30sam["RLEN"] = metrics.str[1]
    
    filtered = q30sam[
        (q30sam["coverage"] >= args.coverage) &
        (q30sam["RLEN"] >= args.min_length)
    ]
    
    # paired-end reads는 R1/R2가 동일 QNAME을 공유한다.
    # QNAME 중복을 제거하여 id_list에 각 read pair가 한 번만 등장하도록 보장한다.
    # samtools view -N은 QNAME 단위로 매칭하므로, 중복이 있어도 동일하게 동작하지만
    # 명시적 dedup으로 리스트 크기를 최소화하고 의도를 명확히 한다.
    filtered["QNAME"].drop_duplicates().to_csv(
        args.output_path, sep="\t", index=False, header=False
    )

    print(f"Filtered {len(filtered)} alignments ({filtered['QNAME'].nunique()} unique reads) out of {len(q30sam)}")
    print(f"Saved to {args.output_path}")


if __name__ == "__main__":
    main()