from __future__ import annotations

from typing import Final, Mapping, Any

CHROMOSOME_MAP: Final[Mapping[str, int]] = {"X": 23, "Y": 24, "MT": 25}
VALID_ALLELES: Final[frozenset[str]] = frozenset(["A", "C", "T", "G"])


def chromosome_to_integer(chromosome: Any) -> int:
    """Remap chromosomes to integers"""
    chrom_string = str(chromosome).strip()
    try:
        chrom = int(chrom_string)
    except ValueError:
        try:
            chrom = CHROMOSOME_MAP[chrom_string]
        except KeyError as bad_remap:
            raise ValueError(f"Invalid chromosome {chromosome}") from bad_remap

    return chrom


def is_valid_sequence(allele: str) -> str:
    if bad_allele := set(allele) - VALID_ALLELES:
        raise ValueError(f"Invalid allele: {bad_allele}")
    return allele


def coerce_na_to_none(x: Any) -> Any:
    """R's default missing value is NA"""
    if x == "NA" or x == "#NA":
        return None
    return x
