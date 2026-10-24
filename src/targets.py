"""Which pattern words the decoy scripts test -- one definition for all of them.

Used by count_genome.py, count_chromosome.py and periodic_decoy.py, so the three
always test the same words with the same names and the same decoys.

Selection
---------
select_targets(targets, families, min_len, max_len) -> {tag: PatternInfo}

  targets = None or ['all']   every word of `families` with min_len <= length <= max_len
                              (default 8 families, 13-100 bits -> 36 words);
                              tag = 'family/variant/label', e.g. 'tribonacci/isolate_0/T6'
  targets = ['readme']        the five words of the original analysis (README section 6.3),
                              with their original tags (trib0_T6, trib0_T7, peri_n5,
                              fibo_F7, lucas_L6)
  targets = [names...]        any mix of README tags and 'family/variant/label' names

Why 100 bits by default: the longest word observed anywhere in T2T is 81 bits
(Tribonacci T7, README 6.1). A word that never occurs has 0 real loci, every
decoy ties with it, and its empirical p-value is exactly 1 -- the test cannot say
anything about it. Use max_len=10000 to include every word anyway.

Decoy seeds
-----------
target_rng(seed, tag) gives each target its own random stream, derived from the
run seed and the target's name. A target's decoys therefore do not change when
other targets are added or removed (in the original scripts one generator was
shared by all targets in turn).
"""
from __future__ import annotations

import zlib
from typing import Dict, Iterable, Optional

import numpy as np

import pattern_genome as PG

README_TARGETS = {
    'trib0_T6': ('tribonacci', 'isolate_0', 'T6'),
    'trib0_T7': ('tribonacci', 'isolate_0', 'T7'),
    'peri_n5':  ('period_doubling', 'standard', 'n5'),
    'fibo_F7':  ('fibonacci', 'standard', 'F7'),
    'lucas_L6': ('lucas', 'standard', 'L6'),
}


def tag_of(p: PG.PatternInfo) -> str:
    return f'{p.family}/{p.variant}/{p.label}'


def select_targets(targets: Optional[Iterable[str]] = None,
                   families: Iterable[str] = PG.DEFAULT_FAMILIES,
                   min_len: int = 13, max_len: int = 100) -> Dict[str, PG.PatternInfo]:
    """{tag: PatternInfo} for the requested targets (see module docstring)."""
    targets = list(targets) if targets else ['all']
    families = list(families)
    if targets == ['all']:
        return {tag_of(p): p for p in PG.build_patterns(families, min_len, max_len, [])}
    every = {(p.family, p.variant, p.label): p
             for p in PG.build_patterns(PG.DEFAULT_FAMILIES, 1, 10000, [])}
    by_tag = {tag_of(p): p for p in every.values()}
    out: Dict[str, PG.PatternInfo] = {}
    for name in targets:
        if name == 'readme':
            for tag, key in README_TARGETS.items():
                out[tag] = every[key]
        elif name in README_TARGETS:
            out[name] = every[README_TARGETS[name]]
        elif name in by_tag:
            out[name] = by_tag[name]
        else:
            raise ValueError(f'unknown target {name!r}: use all, readme, a README tag '
                             f'({", ".join(README_TARGETS)}) or family/variant/label')
    return out


def target_rng(seed: int, tag: str) -> np.random.Generator:
    """Random generator for one target's decoys: depends only on (seed, tag)."""
    return np.random.default_rng([seed, zlib.crc32(tag.encode())])
