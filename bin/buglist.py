#!/usr/bin/env python3
import pandas as pd
import argparse

#### Import from arguments
parser = argparse.ArgumentParser(description="process number of reads to buglist")

parser.add_argument(
    "-i",
    "--input",
    dest="input_path",
    action="store",
    default=None,
    help="Sample to import in this script \n",
)

parser.add_argument(
    "-o",
    "--output",
    dest="output_path",
    action="store",
    default=None,
    help="Sample to output in this script \n",
)

parser.add_argument(
    "-t",
    "--taxa",
    dest="taxonomy_path",
    action="store",
    default=None,
    help="Path of the taxanomy file \n",
)

option = parser.parse_args()

sample_to_process = option.input_path
output_path = option.output_path
taxa_path = option.taxonomy_path

bt_counts = pd.read_csv(sample_to_process, sep="\t", names=["Id", "Abundance"])

if len(bt_counts) >= 1:
    # 동일 prefix를 가진 어셈블리가 복수인 경우 median()은 total reads를 과소 계산한다.
    # 생물학적으로 올바른 집계는 매핑된 reads의 합산이다.
    raw_counts = (
        bt_counts.groupby(bt_counts.columns[0])[bt_counts.columns[-1]]
        .sum()
        .reset_index()
    )
    taxa = pd.read_csv(taxa_path, sep="\t", names=["Id", "taxa"], index_col=False)

    raw_counts.columns = ["Id", "Abundance"]
    taxa.columns = ["Id", "taxa"]
    buglist = pd.merge(taxa, raw_counts)

    column_names = [
        "kingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
        "strain",
    ]

    # expand=True 파라미터를 사용하면 리스트로 변환할 필요 없이 곧바로 데이터프레임이 생성됩니다.
    stratified_counts = buglist["taxa"].str.split("|", expand=True)
    # 간혹 계통 단계가 8단계(strain까지) 미만인 경우 컬럼 개수가 맞지 않을 수 있으므로,
    # 존재하는 열의 수만큼만 column_names를 매핑합니다.
    stratified_counts.columns = column_names[: stratified_counts.shape[1]]

    results = pd.concat([buglist, stratified_counts], axis=1)

    # 최신 Pandas 에러 방지를 위한 groupby 최적화: numeric_only=True 추가
    kingdom_counts = results.groupby(["kingdom"], as_index=False)["Abundance"].sum()

    phylum_counts = results.groupby(["kingdom", "phylum"], as_index=False)[
        "Abundance"
    ].sum()
    phylum_counts["combined"] = (
        phylum_counts["kingdom"].astype(str) + "|" + phylum_counts["phylum"]
    )

    class_counts = results.groupby(["kingdom", "phylum", "class"], as_index=False)[
        "Abundance"
    ].sum()
    class_counts["combined"] = (
        class_counts["kingdom"].astype(str)
        + "|"
        + class_counts["phylum"]
        + "|"
        + class_counts["class"]
    )

    order_counts = results.groupby(
        ["kingdom", "phylum", "class", "order"], as_index=False
    )["Abundance"].sum()
    order_counts["combined"] = (
        order_counts["kingdom"].astype(str)
        + "|"
        + order_counts["phylum"]
        + "|"
        + order_counts["class"]
        + "|"
        + order_counts["order"]
    )

    family_counts = results.groupby(
        ["kingdom", "phylum", "class", "order", "family"], as_index=False
    )["Abundance"].sum()
    family_counts["combined"] = (
        family_counts["kingdom"].astype(str)
        + "|"
        + family_counts["phylum"]
        + "|"
        + family_counts["class"]
        + "|"
        + family_counts["order"]
        + "|"
        + family_counts["family"]
    )

    genus_counts = results.groupby(
        ["kingdom", "phylum", "class", "order", "family", "genus"], as_index=False
    )["Abundance"].sum()
    genus_counts["combined"] = (
        genus_counts["kingdom"].astype(str)
        + "|"
        + genus_counts["phylum"]
        + "|"
        + genus_counts["class"]
        + "|"
        + genus_counts["order"]
        + "|"
        + genus_counts["family"]
        + "|"
        + genus_counts["genus"]
    )

    species_counts = results.groupby(
        ["kingdom", "phylum", "class", "order", "family", "genus", "species"],
        as_index=False,
    )["Abundance"].sum()
    species_counts["combined"] = (
        species_counts["kingdom"].astype(str)
        + "|"
        + species_counts["phylum"]
        + "|"
        + species_counts["class"]
        + "|"
        + species_counts["order"]
        + "|"
        + species_counts["family"]
        + "|"
        + species_counts["genus"]
        + "|"
        + species_counts["species"]
    )

    # strain 데이터가 존재할 경우에만 처리 (taxa 문자열 포맷에 따라 예외 방지)
    if "strain" in results.columns:
        strain_counts = results.groupby(
            [
                "kingdom",
                "phylum",
                "class",
                "order",
                "family",
                "genus",
                "species",
                "strain",
            ],
            as_index=False,
        )["Abundance"].sum()
        strain_counts["combined"] = (
            strain_counts["kingdom"].astype(str)
            + "|"
            + strain_counts["phylum"]
            + "|"
            + strain_counts["class"]
            + "|"
            + strain_counts["order"]
            + "|"
            + strain_counts["family"]
            + "|"
            + strain_counts["genus"]
            + "|"
            + strain_counts["species"]
            + "|"
            + strain_counts["strain"]
        )

    # taxa 리스트와 abundance 리스트를 구성한다.
    # strain 데이터가 있으면 포함하여 완전한 8단계 계층을 출력한다.
    taxa_parts = [
        kingdom_counts["kingdom"],
        phylum_counts["combined"],
        class_counts["combined"],
        order_counts["combined"],
        family_counts["combined"],
        genus_counts["combined"],
        species_counts["combined"],
    ]
    abundance_parts = [
        kingdom_counts["Abundance"],
        phylum_counts["Abundance"],
        class_counts["Abundance"],
        order_counts["Abundance"],
        family_counts["Abundance"],
        genus_counts["Abundance"],
        species_counts["Abundance"],
    ]

    if "strain" in results.columns:
        taxa_parts.append(strain_counts["combined"])
        abundance_parts.append(strain_counts["Abundance"])

    taxa_combined = pd.concat(taxa_parts, ignore_index=True)
    abundance_combined = pd.concat(abundance_parts, ignore_index=True)

    new = pd.concat([taxa_combined, abundance_combined], axis=1)
    new.columns = ["taxa", "Abundance"]

    # 상대 풍부도(RelAbundance) 계산: kingdom 레벨 합계를 전체 할당 reads로 사용한다.
    # 원본 스크립트는 .median 을 사용하였다.
    total_reads = kingdom_counts["Abundance"].sum()
    if total_reads > 0:
        new["RelAbundance"] = (new["Abundance"] / total_reads * 100).round(6)
    else:
        new["RelAbundance"] = 0.0

    new.to_csv(output_path, sep="\t", index=False)

else:
    print("no fungi found\n")
    empty_df = pd.DataFrame(columns=["taxa", "Abundance"])
    empty_df.to_csv(output_path, sep="\t", index=False)