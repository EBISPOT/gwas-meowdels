from __future__ import annotations

from enum import Enum


class RefAlleleState(Enum, str):
    EFFECT_ALLELE = "EA"
    OTHER_ALLELE = "OA"
