from __future__ import annotations

from typing import Annotated, Optional

from pydantic import Field, BeforeValidator, AfterValidator

from ._type_funcs import chromosome_to_integer, is_valid_sequence, coerce_na_to_none
from .enums import RefAlleleState

# https://www.ebi.ac.uk/gwas/docs/summary-statistics-format

Chromosome = Annotated[
    int,
    Field(description="Chromosome where the variant is located", ge=1, le=26),
    BeforeValidator(chromosome_to_integer),
]

BasePairLocation = Annotated[
    int,
    Field(
        description="The first position of the variant in the reference, counting on the bases, from 1",
        gt=0,
    ),
]

EffectAllele = Annotated[
    str,
    Field(description="The allele associated with the effect", min_length=1),
    AfterValidator(is_valid_sequence),
]

OtherAllele = Annotated[
    str,
    Field(description="The non-effect allele", min_length=1),
    AfterValidator(is_valid_sequence),
]

# effect size measurements

Beta = Annotated[
    float,
    Field(
        description="Effect size of numeric traits",
        default=None,
    ),
    BeforeValidator(coerce_na_to_none),
]

OddsRatio = Annotated[
    float,
    Field(
        description="Effect measured as odds ratio",
        ge=0,
    ),
]

HazardRatio = Annotated[
    float, Field(description="Effect measured as hazard ratio", ge=0)
]

StandardError = Annotated[float, Field(description="Standard error of the effect")]

EffectAlleleFrequency = Annotated[
    float,
    Field(
        description="Frequency of the effect allele in the control population",
        # TODO: are 0/1 included or excluded?
        ge=0,
        le=0,
    ),
]

PValue = Annotated[
    float,
    Field(
        description="P-value of the association statistic",
        ge=0,
        le=1,
    ),
]

NegLog10PValue = Annotated[
    float,
    Field(
        description="Negative log10 p-value of the association statistic",
        ge=0,  # zero or more
    ),
]


def validate_variant_id(variant_id: str) -> str:
    parts = variant_id.split("_", 3)

    if len(parts) != 4:
        raise ValueError(f"Invalid variant_id: {variant_id} (bad delimiters _)")
    else:
        chrom, pos, ref, alt = parts

    # reuse existing validators
    # TODO: are X / Y / MT OK in VariantId? what about patches?
    _ = Chromosome(chrom)
    _ = BasePairLocation(pos)
    # just checking they're valid sequences
    _ = EffectAllele(ref)
    _ = EffectAllele(alt)

    return variant_id


VariantId = Annotated[
    str,
    Field(
        description="An internal variant identifier in the form of <chromosome>_<base_pair_location>_ <reference_allele>_<alternate_allele>",
        pattern=r"^([0-9]{1,2}|X|Y|MT)_[0-9]+_[ACGT]+_[ACGT]+$",
    ),
    AfterValidator(validate_variant_id),
]

RsID = Annotated[
    str,
    Field(description="The rsID of the variant", pattern=r"^rs[0-9]+$"),
]

Info = Annotated[float, Field(description="Imputation information metric", ge=0, le=1)]

# TODO: is this _just_ the odds ratio?
CI_Upper = Annotated[
    float, Field(description="Upper confidence interval for the odds ratio", ge=0)
]

CI_Lower = Annotated[
    float, Field(description="Lower confidence interval for the odds ratio", ge=0)
]


RefAllele = Annotated[
    Optional[RefAlleleState],
    Field(description="State which of the alleles is the reference allele"),
    BeforeValidator(coerce_na_to_none),
]

N = Annotated[int, Field(description="Sample size per variant", gt=0)]
