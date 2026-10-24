#!/usr/bin/env python3
"""
Scan a DNA FASTA for mathematical/automatic binary words after DNA projection.

DNA -> binary projections (defined once, in projections.py):
  RY: purine/pyrimidine  A,G -> 1 ; C,T -> 0
  WS: strong/weak H-bond C,G -> 1 ; A,T -> 0
  MK: amino/keto        A,C -> 1 ; G,T -> 0
Non-ACGT characters become 'N', so no motif can cross a gap or ambiguity code.

Forms searched for each word (--strands)
----------------------------------------
  both (default)  word, bit complement, reverse, reversed complement. The
                  complement covers the arbitrary 0/1 labelling; the reverse
                  forms cover the opposite DNA strand (projections.py explains
                  why the set is the same in RY, WS and MK). Forms that coincide
                  (e.g. Thue-Morse prefixes, whose reverse equals the word or
                  its complement) are searched once, so no hit is counted twice.
  forward         word and complement only (the original behaviour).

Built-in pattern families
-------------------------
  fibonacci        F0=0, F1=01, Fn=F(n-1)F(n-2); lengths 1,2,3,5,8,13,21,34,...
  lucas            L0=0, L1=011, Ln=L(n-1)L(n-2); lengths 1,3,4,7,11,18,29,47,...
  tribonacci       ternary fixed point 0->01, 1->02, 2->0, tested under all three
                   one-vs-two binary codings (isolate_0, isolate_1, isolate_2)
  thue_morse       fixed point of 0->01, 1->10; prefix lengths 1,2,4,8,...
  period_doubling  fixed point of 0->01, 1->00; prefix lengths 1,2,4,8,...
  rudin_shapiro    parity of overlapping '11' pairs in binary(n); lengths 2^k
  paperfolding     regular paperfolding sequence; lengths 2^k
  baum_sweet       1 iff binary(n) has no odd-length run of 0s (b(0)=1); lengths 2^k
Custom exact motifs: --pattern NAME=BITS. Each built-in family is a chain of
nested prefixes; build_patterns() checks this, because the scan relies on it:
once a word and all its forms are absent from a contig, every longer word of the
chain is absent too (this also holds for the reverse forms, since the reverse
of a longer prefix ends with the reverse of the shorter one).

Background models, per contig, summed over contigs
--------------------------------------------------
  expected_iid        i.i.d. bits with the contig's 0/1 frequencies
  expected_markov{k}  order-k Markov chain on the bits (--markov-order, default 1,
                      which reproduces the original expected_markov1 exactly)
Both are multiplied by the number of valid start positions (windows not
crossing an N). They cannot produce tandem repeats, so repeat-shaped words look
enormously enriched against them: treat the O/E and p-values here as a
screening statistic. The decoy-word tests are the control for the rule.

Observed counts, clumps and p-values
------------------------------------
observed_total counts every (overlapping) occurrence of every form. Inside a
tandem array a self-overlapping word hits again and again, so raw hits are far
from independent and a Poisson test on them is much too optimistic.
observed_clumps merges overlapping hits (of any form) into one clump.
expected_clumps_* is the exact expected number of clump starts under the same
null: an occurrence starts a clump unless another occurrence begins fewer than
L positions earlier; partitioning by the NEAREST preceding occurrence gives

  E[clumps] = E[hits] - sum over (v, d, u) of starts(L+d) * P(v[:d] + u)

for every pair of forms v, u and shift d < L at which v at i-d and u at i are
compatible and the merged word contains no other occurrence in between.
(tests/test_pattern_genome.py checks this against simulation.)

p-values (Poisson, computed in log space, no scipy needed):
  p_enrich_hits_*     P(X >= observed_total | E[hits])     optimistic, for reference
  p_enrich_clumps_*   P(X >= observed_clumps | E[clumps])  the recommended one
  p_deplete_clumps_*  P(X <= observed_clumps | E[clumps])
  log10_p_enrich_clumps_*  same as p_enrich_clumps, without underflow
  q_enrich_clumps_*   Benjamini-Hochberg q-value over all rows of the run
The original p_ge_1_*_rare_event_approx columns (= 1 - exp(-E), the probability
of at least one hit) are kept for comparison; they are not p-values for the
observed count.

Example loci (--loci, --loci-sampling)
--------------------------------------
Up to --loci (default 20) hit positions per mapping x word are kept in the
example_loci column, for inspection with loci_context.py. They do not enter any
count. --loci-sampling random (default) draws them uniformly from ALL hits in
the genome (seeded, --seed); 'first' keeps the first ones in file order (the
original behaviour, biased toward the first chromosomes).

Reproducing the original numbers
--------------------------------
  --strands forward --loci-sampling first   (with the default --markov-order 1)
gives the original columns unchanged; the new columns are added.

Output
------
--csv FILE writes the table and FILE's stem + '_run.json' with the settings.
"""

from __future__ import annotations

import argparse
import csv
import datetime
import gzip
import json
import math
import platform
import re
import sys
from collections import Counter
from dataclasses import dataclass, field
from typing import Dict, Iterator, List, Optional, Sequence, Tuple

import numpy as np

import projections
import seqio

VERSION = '2.0'

MAPPINGS = projections.MAPPINGS
TRANSLATE_TABLES = projections.TRANSLATE_TABLES
make_translate_table = projections.make_translate_table
VALID_RUN_RE = re.compile(rb"[01]+")

DEFAULT_FAMILIES = [
    "fibonacci",
    "lucas",
    "tribonacci",
    "thue_morse",
    "period_doubling",
    "rudin_shapiro",
    "paperfolding",
    "baum_sweet",
]

FAMILY_DESCRIPTIONS = {
    "fibonacci": "F0=0, F1=01, Fn=F(n-1)F(n-2)",
    "lucas": "L0=0, L1=011, Ln=L(n-1)L(n-2)",
    "tribonacci": "ternary fixed point 0->01,1->02,2->0, tested under all 3 one-vs-two binary codings",
    "thue_morse": "binary fixed point 0->01,1->10",
    "period_doubling": "binary fixed point 0->01,1->00",
    "rudin_shapiro": "parity of overlapping 11-pairs in binary(n)",
    "paperfolding": "regular paperfolding sequence",
    "baum_sweet": "1 iff binary(n) has no odd-length zero run; b(0)=1",
}

FORM_LABELS = ("word", "complement", "reverse", "reverse_complement")


@dataclass(frozen=True)
class PatternInfo:
    family: str
    variant: str
    order: int
    label: str
    word: bytes
    comp: bytes
    length: int
    c0: int
    c1: int
    t00: int
    t01: int
    t10: int
    t11: int
    nested_prefix: bool = True

    @property
    def chain_key(self) -> Tuple[str, str]:
        return self.family, self.variant


@dataclass
class Result:
    mapping: str
    family: str
    variant: str
    order: int
    label: str
    length: int
    word: str
    complement: str
    observed_word: int = 0
    observed_complement: int = 0
    observed_reverse: int = 0
    observed_reverse_complement: int = 0
    observed_clumps: int = 0
    n_forms: int = 2
    expected_iid: float = 0.0
    expected_markov1: float = 0.0          # expected count under the Markov order of the run
    expected_clumps_iid: float = 0.0
    expected_clumps_markov: float = 0.0
    loci: List[Tuple[str, int, str]] = field(default_factory=list)
    q_enrich_clumps: float = float("nan")
    _sample: list = field(default_factory=list, repr=False)   # (key, record_idx, name, pos, form)

    @property
    def observed(self) -> int:
        return (self.observed_word + self.observed_complement
                + self.observed_reverse + self.observed_reverse_complement)


# ----------------------------- FASTA / projection -----------------------------

def open_maybe_gzip(path: str):
    """Kept for compatibility; reading goes through seqio (gzip detected by content)."""
    if path.endswith(".gz"):
        return gzip.open(path, "rb")
    return open(path, "rb")


def read_fasta(path: str) -> Iterator[Tuple[str, str, bytes]]:
    """Yield (record_name, full_header, sequence_bytes) -- via seqio.read_fasta."""
    for rec in seqio.read_fasta(path):
        yield rec.name, rec.header.decode("utf-8", errors="replace").strip(), rec.seq


def bit_complement(word: bytes) -> bytes:
    return word.translate(bytes.maketrans(b"01", b"10"))


def search_forms(word: bytes, strands: str = "both") -> List[Tuple[str, bytes]]:
    """Distinct forms to search: [(label, word), ...]; same set as projections.search_set."""
    c = bit_complement(word)
    cands = [("word", word), ("complement", c)]
    if strands == "both":
        r = word[::-1]
        cands += [("reverse", r), ("reverse_complement", bit_complement(r))]
    elif strands != "forward":
        raise ValueError("strands must be 'both' or 'forward'")
    seen, out = set(), []
    for lab, w in cands:
        if w not in seen:
            seen.add(w)
            out.append((lab, w))
    return out


# ------------------------------- pattern stats --------------------------------

def transition_counts(word: bytes) -> Tuple[int, int, int, int]:
    t00 = t01 = t10 = t11 = 0
    for a, b in zip(word, word[1:]):
        if a == 48 and b == 48:
            t00 += 1
        elif a == 48 and b == 49:
            t01 += 1
        elif a == 49 and b == 48:
            t10 += 1
        else:
            t11 += 1
    return t00, t01, t10, t11


def make_pattern(family: str, variant: str, order: int, label: str, word: bytes,
                 nested_prefix: bool = True) -> PatternInfo:
    if not word or any(ch not in (48, 49) for ch in word):
        raise ValueError(f"Pattern {family}/{variant}/{label} is not a nonempty binary word")
    c0 = word.count(b"0")
    c1 = len(word) - c0
    t00, t01, t10, t11 = transition_counts(word)
    return PatternInfo(family=family, variant=variant, order=order, label=label, word=word,
                       comp=bit_complement(word), length=len(word), c0=c0, c1=c1,
                       t00=t00, t01=t01, t10=t10, t11=t11, nested_prefix=nested_prefix)


def maybe_add(out: List[PatternInfo], min_len: int, max_len: int, family: str, variant: str,
              order: int, label: str, word: bytes, nested_prefix: bool = True) -> None:
    if min_len <= len(word) <= max_len:
        out.append(make_pattern(family, variant, order, label, word, nested_prefix))


# ---------------------------- built-in word families --------------------------

def fibonacci_patterns(min_len: int, max_len: int) -> List[PatternInfo]:
    out: List[PatternInfo] = []
    a, b = b"0", b"01"
    order = 0
    maybe_add(out, min_len, max_len, "fibonacci", "standard", 0, "F0", a)
    maybe_add(out, min_len, max_len, "fibonacci", "standard", 1, "F1", b)
    while True:
        c = b + a
        order += 1
        n = order + 1
        if len(c) > max_len:
            break
        maybe_add(out, min_len, max_len, "fibonacci", "standard", n, f"F{n}", c)
        a, b = b, c
    return out


def lucas_patterns(min_len: int, max_len: int) -> List[PatternInfo]:
    out: List[PatternInfo] = []
    a, b = b"0", b"011"
    maybe_add(out, min_len, max_len, "lucas", "standard", 0, "L0", a)
    maybe_add(out, min_len, max_len, "lucas", "standard", 1, "L1", b)
    n = 2
    while True:
        c = b + a
        if len(c) > max_len:
            break
        maybe_add(out, min_len, max_len, "lucas", "standard", n, f"L{n}", c)
        a, b = b, c
        n += 1
    return out


def apply_morphism(word: bytes, morphism: Dict[int, bytes]) -> bytes:
    return b"".join(morphism[ch] for ch in word)


def binary_morphic_patterns(family: str, morphism: Dict[int, bytes], min_len: int,
                            max_len: int) -> List[PatternInfo]:
    out: List[PatternInfo] = []
    word = b"0"
    order = 0
    while len(word) <= max_len:
        maybe_add(out, min_len, max_len, family, "standard", order, f"n{order}", word)
        nxt = apply_morphism(word, morphism)
        if len(nxt) <= len(word):
            raise ValueError(f"Morphism for {family} did not grow")
        word = nxt
        order += 1
    return out


def tribonacci_patterns(min_len: int, max_len: int) -> List[PatternInfo]:
    """Standard ternary Tribonacci word, tested under all 3 binary partitions."""
    out: List[PatternInfo] = []
    morphism = {ord("0"): b"01", ord("1"): b"02", ord("2"): b"0"}
    word = b"0"
    order = 0
    while len(word) <= max_len:
        for isolated in (48, 49, 50):
            variant = f"isolate_{chr(isolated)}"
            binary = bytes(49 if ch == isolated else 48 for ch in word)
            maybe_add(out, min_len, max_len, "tribonacci", variant, order, f"T{order}", binary)
        word = apply_morphism(word, morphism)
        order += 1
    return out


def rudin_shapiro_bit(n: int) -> int:
    return (n & (n >> 1)).bit_count() & 1


def paperfolding_bit(n: int) -> int:
    m = n + 1
    while (m & 1) == 0:
        m >>= 1
    return 1 if (m & 3) == 1 else 0


def baum_sweet_bit(n: int) -> int:
    if n == 0:
        return 1
    s = bin(n)[2:]
    i = 0
    while i < len(s):
        if s[i] == "1":
            i += 1
            continue
        j = i
        while j < len(s) and s[j] == "0":
            j += 1
        if (j - i) & 1:
            return 0
        i = j
    return 1


def power2_prefix_patterns(family: str, bit_func, min_len: int, max_len: int) -> List[PatternInfo]:
    out: List[PatternInfo] = []
    length = 1
    order = 0
    buf = bytearray()
    generated = 0
    while length <= max_len:
        while generated < length:
            buf.append(49 if bit_func(generated) else 48)
            generated += 1
        maybe_add(out, min_len, max_len, family, "standard", order, f"n{order}", bytes(buf))
        length *= 2
        order += 1
    return out


def custom_patterns(specs: Sequence[str], min_len: int, max_len: int) -> List[PatternInfo]:
    """Parse repeated --pattern NAME=BITS specifications."""
    out: List[PatternInfo] = []
    seen = set()
    for spec in specs:
        if "=" not in spec:
            raise SystemExit(f"Invalid --pattern {spec!r}; expected NAME=BITS")
        name, bits = spec.split("=", 1)
        name = name.strip()
        bits = "".join(bits.split())
        if not name:
            raise SystemExit("Custom pattern name cannot be empty")
        if name in seen:
            raise SystemExit(f"Duplicate custom pattern name: {name}")
        seen.add(name)
        if not bits or set(bits) - {"0", "1"}:
            raise SystemExit(f"Custom pattern {name!r} must contain only 0 and 1")
        word = bits.encode("ascii")
        if min_len <= len(word) <= max_len:
            out.append(make_pattern("custom", name, 0, name, word, nested_prefix=False))
    return out


def check_prefix_chains(patterns: Sequence[PatternInfo]) -> None:
    """Raise ValueError unless every nested_prefix chain is really prefix-nested."""
    chains: Dict[Tuple[str, str], List[PatternInfo]] = {}
    for p in patterns:
        if p.nested_prefix:
            chains.setdefault(p.chain_key, []).append(p)
    for key, pats in chains.items():
        pats = sorted(pats, key=lambda p: (p.length, p.order))
        for a, b in zip(pats, pats[1:]):
            if not b.word.startswith(a.word):
                raise ValueError(f"chain {key}: {b.label} does not start with {a.label}; "
                                 "the early stop in scan() would be wrong")


def build_patterns(families: Sequence[str], min_len: int, max_len: int,
                   custom_specs: Sequence[str]) -> List[PatternInfo]:
    if min_len < 1 or max_len < min_len:
        raise ValueError("Require 1 <= min_word_len <= max_word_len")

    generators = {
        "fibonacci": fibonacci_patterns,
        "lucas": lucas_patterns,
        "tribonacci": tribonacci_patterns,
        "thue_morse": lambda a, b: binary_morphic_patterns("thue_morse", {48: b"01", 49: b"10"}, a, b),
        "period_doubling": lambda a, b: binary_morphic_patterns("period_doubling", {48: b"01", 49: b"00"}, a, b),
        "rudin_shapiro": lambda a, b: power2_prefix_patterns("rudin_shapiro", rudin_shapiro_bit, a, b),
        "paperfolding": lambda a, b: power2_prefix_patterns("paperfolding", paperfolding_bit, a, b),
        "baum_sweet": lambda a, b: power2_prefix_patterns("baum_sweet", baum_sweet_bit, a, b),
    }

    out: List[PatternInfo] = []
    for family in families:
        out.extend(generators[family](min_len, max_len))
    out.extend(custom_patterns(custom_specs, min_len, max_len))
    check_prefix_chains(out)
    return sorted(out, key=lambda p: (p.family, p.variant, p.length, p.order))


# ----------------------------- background models ------------------------------

def run_lengths(binary: bytes) -> List[int]:
    return [m.end() - m.start() for m in VALID_RUN_RE.finditer(binary)]


def binary_stats(binary: bytes) -> Tuple[int, int, int, int, int, int]:
    """Return n0,n1,n00,n01,n10,n11 for valid runs only."""
    n0 = binary.count(b"0")
    n1 = binary.count(b"1")
    if n0 + n1 == 0:
        return 0, 0, 0, 0, 0, 0
    n01 = binary.count(b"01")
    n10 = binary.count(b"10")
    end0 = binary.count(b"0N") + (1 if binary.endswith(b"0") else 0)
    end1 = binary.count(b"1N") + (1 if binary.endswith(b"1") else 0)
    out0 = n0 - end0
    out1 = n1 - end1
    n00 = out0 - n01
    n11 = out1 - n10
    return n0, n1, n00, n01, n10, n11


def valid_start_count(runs: List[int], motif_len: int) -> int:
    return sum((r - motif_len + 1) for r in runs if r >= motif_len)


def safe_pow(p: float, k: int) -> float:
    if k == 0:
        return 1.0
    if p <= 0.0:
        return 0.0
    return math.exp(k * math.log(p))


def iid_prob(c0: int, c1: int, p0: float, p1: float) -> float:
    return safe_pow(p0, c0) * safe_pow(p1, c1)


def markov_prob(first: int, t00: int, t01: int, t10: int, t11: int, pi0: float, pi1: float,
                p00: float, p01: float, p10: float, p11: float) -> float:
    p = pi0 if first == 0 else pi1
    if p <= 0.0:
        return 0.0
    terms = ((p00, t00), (p01, t01), (p10, t10), (p11, t11))
    logp = math.log(p)
    for q, k in terms:
        if k == 0:
            continue
        if q <= 0.0:
            return 0.0
        logp += k * math.log(q)
    return math.exp(logp)


def expected_for_contig(stats, starts: int, pat: PatternInfo) -> Tuple[float, float]:
    """Original helper (word + complement, i.i.d. and Markov-1); kept for compatibility."""
    n0, n1, n00, n01, n10, n11 = stats
    nvalid = n0 + n1
    if nvalid == 0 or starts == 0:
        return 0.0, 0.0
    p0, p1 = n0 / nvalid, n1 / nvalid
    exp_iid = starts * (iid_prob(pat.c0, pat.c1, p0, p1) + iid_prob(pat.c1, pat.c0, p0, p1))
    out0, out1 = n00 + n01, n10 + n11
    p00 = n00 / out0 if out0 else 0.0
    p01 = n01 / out0 if out0 else 0.0
    p10 = n10 / out1 if out1 else 0.0
    p11 = n11 / out1 if out1 else 0.0
    first = 0 if pat.word[0] == 48 else 1
    pw = markov_prob(first, pat.t00, pat.t01, pat.t10, pat.t11, p0, p1, p00, p01, p10, p11)
    pc = markov_prob(1 - first, pat.t11, pat.t10, pat.t01, pat.t00, p0, p1, p00, p01, p10, p11)
    return exp_iid, starts * (pw + pc)


def kgram_counts(binary: bytes, k: int, chunk: int = 1 << 22) -> np.ndarray:
    """Counts of all 2^k bit k-grams lying entirely inside valid (non-N) runs."""
    a = np.frombuffer(binary, np.uint8)
    n = len(a)
    out = np.zeros(1 << k, np.int64)
    for s in range(0, n, chunk):
        seg = a[s:min(n, s + chunk + k - 1)]
        m = len(seg) - k + 1
        if m <= 0:
            break
        bit = (seg == 49).astype(np.uint32)
        bad = ((seg != 48) & (seg != 49)).astype(np.int64)
        cb = np.concatenate(([0], np.cumsum(bad)))
        ok = (cb[k:k + m] - cb[:m]) == 0
        ids = np.zeros(m, np.uint32)
        for j in range(k):
            ids = (ids << np.uint32(1)) | bit[j:j + m]
        out += np.bincount(ids[ok], minlength=1 << k)
    return out


@dataclass(frozen=True)
class WordStats:
    """What the background models need to know about one binary string."""
    length: int
    c0: int
    c1: int
    first: int
    t: Tuple[int, int, int, int]                 # t00, t01, t10, t11
    kfirst: int = 0                              # id of the first k bits (order k >= 2)
    kgrams: Tuple[Tuple[int, int], ...] = ()     # ((k+1)-gram id, count) (order k >= 2)


def word_stats(w: bytes, order: int) -> WordStats:
    c0 = w.count(b"0")
    kfirst, kgrams = 0, ()
    if order >= 2:
        bits = [ch - 48 for ch in w]
        kfirst = int(w[:order], 2)
        cnt: Counter = Counter()
        g = int(w[:order + 1], 2)
        mask = (1 << (order + 1)) - 1
        cnt[g] += 1
        for b in bits[order + 1:]:
            g = ((g << 1) | b) & mask
            cnt[g] += 1
        kgrams = tuple(sorted(cnt.items()))
    return WordStats(len(w), c0, len(w) - c0, w[0] - 48, transition_counts(w), kfirst, kgrams)


class ContigModel:
    """i.i.d. and order-k Markov models fitted to one projected contig."""

    def __init__(self, binary: bytes, order: int):
        self.order = order
        self.runs = run_lengths(binary)
        self._starts: Dict[int, int] = {}
        n0, n1, n00, n01, n10, n11 = binary_stats(binary)
        nvalid = n0 + n1
        self.nvalid = nvalid
        self.p0 = n0 / nvalid if nvalid else 0.0
        self.p1 = n1 / nvalid if nvalid else 0.0
        out0, out1 = n00 + n01, n10 + n11
        self.p00 = n00 / out0 if out0 else 0.0
        self.p01 = n01 / out0 if out0 else 0.0
        self.p10 = n10 / out1 if out1 else 0.0
        self.p11 = n11 / out1 if out1 else 0.0
        if order >= 2 and nvalid:
            ck = kgram_counts(binary, order).astype(float)
            ck1 = kgram_counts(binary, order + 1).astype(float)
            with np.errstate(divide="ignore", invalid="ignore"):
                self.logP0 = np.log(ck / ck.sum()) if ck.sum() else np.full(len(ck), -np.inf)
                ctx = ck1.reshape(-1, 2).sum(1)
                self.logT = np.log(ck1 / np.repeat(ctx, 2))
            self.logT[~np.isfinite(self.logT)] = -np.inf

    def starts(self, L: int) -> int:
        s = self._starts.get(L)
        if s is None:
            s = self._starts[L] = valid_start_count(self.runs, L)
        return s

    def prob_iid(self, ws: WordStats) -> float:
        return iid_prob(ws.c0, ws.c1, self.p0, self.p1)

    def prob_markov(self, ws: WordStats) -> float:
        if self.order == 0:
            return self.prob_iid(ws)
        if self.order == 1:
            return markov_prob(ws.first, *ws.t, self.p0, self.p1,
                               self.p00, self.p01, self.p10, self.p11)
        lp = self.logP0[ws.kfirst]
        if lp == -np.inf:
            return 0.0
        for g, c in ws.kgrams:
            q = self.logT[g]
            if q == -np.inf:
                return 0.0
            lp += c * q
        return math.exp(lp)


# ------------------------------------------------------------------- clumps

def _z_array(s: bytes) -> List[int]:
    n = len(s)
    z = [0] * n
    z[0] = n
    l = r = 0
    for i in range(1, n):
        if i < r:
            z[i] = min(r - i, z[i - l])
        while i + z[i] < n and s[z[i]] == s[i + z[i]]:
            z[i] += 1
        if i + z[i] > r:
            l, r = i, i + z[i]
    return z


def compatible_shifts(v: bytes, u: bytes) -> List[int]:
    """Shifts d in 1..L-1 such that v at i-d and u at i agree on their overlap
    (v[d:] == u[:L-d])."""
    L = len(u)
    z = _z_array(u + b"#" + v)
    return [d for d in range(1, L) if z[L + 1 + d] >= L - d]


def clump_terms(forms: Sequence[bytes]) -> List[bytes]:
    """Merged words v[:d] + u whose probability is subtracted from E[hits] to get
    E[clumps]: v at i-d is the NEAREST occurrence before u at i (no other form
    starts strictly between them)."""
    L = len(forms[0])
    terms = []
    for u in forms:
        for v in forms:
            for d in compatible_shifts(v, u):
                merged = v[:d] + u
                inner = any(0 <= merged.find(f, 1, d - 1 + L) for f in forms) if d > 1 else False
                if not inner:
                    terms.append(merged)
    return terms


# ----------------------------------------------------------- Poisson tails

def _log_gamma_p_q(a: float, x: float) -> Tuple[float, float]:
    """(log P(a,x), log Q(a,x)) for the regularized incomplete gamma functions."""
    if x <= 0:
        return -math.inf, 0.0
    gln = math.lgamma(a)
    pre = -x + a * math.log(x) - gln
    if x < a + 1.0:                       # series for P
        ap, s, d = a, 1.0 / a, 1.0 / a
        for _ in range(10_000_000):
            ap += 1.0
            d *= x / ap
            s += d
            if abs(d) < abs(s) * 1e-15:
                break
        lp = math.log(s) + pre
        lq = math.log1p(-math.exp(lp)) if lp < -1e-300 else -math.inf
        return lp, lq
    tiny = 1e-300                          # continued fraction for Q (Lentz)
    b = x + 1.0 - a
    c = 1.0 / tiny
    d = 1.0 / b
    h = d
    for i in range(1, 10_000_000):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < tiny:
            d = tiny
        c = b + an / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        de = d * c
        h *= de
        if abs(de - 1.0) < 1e-15:
            break
    lq = math.log(h) + pre
    lp = math.log1p(-math.exp(lq)) if lq < -1e-300 else -math.inf
    return lp, lq


def poisson_log_sf(k: int, lam: float) -> float:
    """log P(X >= k) for X ~ Poisson(lam)."""
    if k <= 0:
        return 0.0
    if lam <= 0:
        return -math.inf
    return _log_gamma_p_q(float(k), lam)[0]


def poisson_log_cdf(k: int, lam: float) -> float:
    """log P(X <= k) for X ~ Poisson(lam)."""
    if k < 0:
        return -math.inf
    if lam <= 0:
        return 0.0
    return _log_gamma_p_q(float(k + 1), lam)[1]


def bh_qvalues(p: Sequence[float]) -> List[float]:
    """Benjamini-Hochberg q-values; NaN entries are ignored and returned as NaN."""
    idx = [i for i, x in enumerate(p) if not math.isnan(x)]
    m = len(idx)
    q = [float("nan")] * len(p)
    order = sorted(idx, key=lambda i: p[i])
    running = 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        running = min(running, p[i] * m / rank)
        q[i] = running
    return q


# ---------------------------------- scanner ------------------------------------

def find_all(binary: bytes, pattern: bytes, cap: int) -> Tuple[int, List[int]]:
    """Count overlapping occurrences; return up to cap 0-based starts."""
    count = 0
    kept: List[int] = []
    pos = binary.find(pattern)
    while pos != -1:
        count += 1
        if len(kept) < cap:
            kept.append(pos)
        pos = binary.find(pattern, pos + 1)
    return count, kept


def find_positions(binary: bytes, pattern: bytes) -> np.ndarray:
    """All overlapping 0-based start positions of pattern."""
    out = []
    pos = binary.find(pattern)
    while pos != -1:
        out.append(pos)
        pos = binary.find(pattern, pos + 1)
    return np.asarray(out, dtype=np.int64)


def fmt(x: float) -> str:
    if x == 0:
        return "0"
    if math.isnan(x):
        return "nan"
    if x < 1e-3 or x >= 1e6:
        return f"{x:.4e}"
    return f"{x:.4f}"


@dataclass
class _Prepared:
    pattern: PatternInfo
    forms: List[Tuple[str, bytes]]
    form_stats: List[WordStats]
    term_stats: List[WordStats]


def _prepare(p: PatternInfo, strands: str, order: int) -> _Prepared:
    forms = search_forms(p.word, strands)
    words = [w for _, w in forms]
    return _Prepared(p, forms, [word_stats(w, order) for w in words],
                     [word_stats(t, order) for t in clump_terms(words)])


def _sample_loci(r: Result, positions: Dict[str, np.ndarray], name: str, rec_idx: int,
                 cap: int, mode: str, rng: np.random.Generator) -> None:
    if cap <= 0:
        return
    if mode == "first":
        for lab, pos in positions.items():
            remaining = cap - len(r.loci)
            if remaining <= 0:
                break
            r.loci.extend((name, int(x) + 1, lab) for x in pos[:remaining])
        return
    for lab, pos in positions.items():
        if len(pos) == 0:
            continue
        keys = rng.random(len(pos))
        if len(pos) > cap:
            sel = np.argpartition(keys, cap - 1)[:cap]
            keys, pos = keys[sel], pos[sel]
        r._sample.extend(zip(keys.tolist(), [rec_idx] * len(pos), [name] * len(pos),
                             (pos + 1).tolist(), [lab] * len(pos)))
        if len(r._sample) > cap:
            r._sample.sort(key=lambda t: t[0])
            del r._sample[cap:]


def scan(args) -> List[Result]:
    strands = getattr(args, "strands", "forward")
    order = getattr(args, "markov_order", 1)
    loci_mode = getattr(args, "loci_sampling", "first")
    rng = np.random.default_rng(getattr(args, "seed", 1))

    patterns = build_patterns(args.families, args.min_word_len, args.max_word_len, args.pattern)
    if not patterns:
        raise SystemExit("No requested patterns fall inside the requested length range.")
    if order < 0 or order > 16:
        raise SystemExit("--markov-order must be between 0 and 16")
    if order >= min(p.length for p in patterns):
        raise SystemExit("--markov-order must be smaller than the shortest word length")

    prepared = {p: _prepare(p, strands, order) for p in patterns}
    chains: Dict[Tuple[str, str], List[PatternInfo]] = {}
    for p in patterns:
        chains.setdefault(p.chain_key, []).append(p)
    for pats in chains.values():
        pats.sort(key=lambda p: (p.length, p.order))

    results: Dict[Tuple[str, str, str, int], Result] = {}
    for mapping in args.mappings:
        for p in patterns:
            results[(mapping, p.family, p.variant, p.order)] = Result(
                mapping=mapping, family=p.family, variant=p.variant, order=p.order, label=p.label,
                length=p.length, word=p.word.decode(), complement=p.comp.decode(),
                n_forms=len(prepared[p].forms))

    include_re = re.compile(args.include_regex) if args.include_regex else None
    exclude_re = re.compile(args.exclude_regex) if args.exclude_regex else None

    record_count = 0
    base_count = 0
    for name, header, seq in read_fasta(args.fasta):
        if include_re and not include_re.search(header):
            continue
        if exclude_re and exclude_re.search(header):
            continue
        record_count += 1
        base_count += len(seq)
        if args.progress:
            print(f"[record {record_count}] {name}  {len(seq):,} bp", file=sys.stderr)

        for mapping in args.mappings:
            binary = seq.translate(TRANSLATE_TABLES[mapping])
            model = ContigModel(binary, order)

            for chain_key in sorted(chains):
                alive = True
                for p in chains[chain_key]:
                    prep = prepared[p]
                    r = results[(mapping, p.family, p.variant, p.order)]
                    L = p.length
                    starts = model.starts(L)
                    if model.nvalid and starts:
                        pi = [model.prob_iid(s) for s in prep.form_stats]
                        pm = [model.prob_markov(s) for s in prep.form_stats]
                        e_iid = starts * sum(pi)
                        e_mk = starts * sum(pm)
                        r.expected_iid += e_iid
                        r.expected_markov1 += e_mk
                        c_iid = e_iid - (sum(model.starts(t.length) * model.prob_iid(t)
                                             for t in prep.term_stats) if e_iid > 0 else 0.0)
                        c_mk = e_mk - (sum(model.starts(t.length) * model.prob_markov(t)
                                           for t in prep.term_stats) if e_mk > 0 else 0.0)
                        r.expected_clumps_iid += max(c_iid, 0.0)
                        r.expected_clumps_markov += max(c_mk, 0.0)

                    if p.nested_prefix and not alive:
                        continue

                    positions = {lab: find_positions(binary, w) for lab, w in prep.forms}
                    total = 0
                    for lab, pos in positions.items():
                        setattr(r, "observed_" + lab, getattr(r, "observed_" + lab) + len(pos))
                        total += len(pos)
                    if total:
                        allpos = np.sort(np.concatenate(list(positions.values())))
                        r.observed_clumps += 1 + int((np.diff(allpos) >= L).sum())
                    _sample_loci(r, positions, name, record_count, args.loci, loci_mode, rng)

                    if p.nested_prefix and total == 0:
                        alive = False

    if loci_mode != "first":
        for r in results.values():
            r._sample.sort(key=lambda t: (t[1], t[3]))
            r.loci = [(nm, pos, lab) for _, _, nm, pos, lab in r._sample]
            r._sample = []

    rows = sorted(results.values(), key=lambda r: (r.mapping, r.family, r.variant, r.length, r.order))
    q = bh_qvalues([math.exp(poisson_log_sf(r.observed_clumps, r.expected_clumps_markov)) for r in rows])
    for r, qq in zip(rows, q):
        r.q_enrich_clumps = qq

    args._run_stats = dict(records=record_count, fasta_characters=base_count)
    if args.progress:
        print(f"Scanned {record_count} records, {base_count:,} FASTA characters.", file=sys.stderr)
    return rows


# ---------------------------------- output -------------------------------------

def _ratio(o: float, e: float) -> float:
    return o / e if e > 0 else math.inf


def row_values(r: Result, k: int) -> Dict[str, object]:
    """All output columns for one result (column names use the Markov order k)."""
    mk = f"markov{k}"
    lp_hits = poisson_log_sf(r.observed, r.expected_markov1)
    lp_clump = poisson_log_sf(r.observed_clumps, r.expected_clumps_markov)
    lp_clump_iid = poisson_log_sf(r.observed_clumps, r.expected_clumps_iid)
    lp_dep = poisson_log_cdf(r.observed_clumps, r.expected_clumps_markov)
    g = lambda x: f"{x:.12g}"
    return {
        "mapping": r.mapping, "family": r.family, "variant": r.variant, "order": r.order,
        "label": r.label, "length": r.length, "word": r.word, "complement": r.complement,
        "observed_word": r.observed_word, "observed_complement": r.observed_complement,
        "observed_total": r.observed,
        "expected_iid": g(r.expected_iid), f"expected_{mk}": g(r.expected_markov1),
        "obs_over_iid": g(_ratio(r.observed, r.expected_iid)),
        f"obs_over_{mk}": g(_ratio(r.observed, r.expected_markov1)),
        "p_ge_1_iid_rare_event_approx": g(-math.expm1(-r.expected_iid)),
        f"p_ge_1_{mk}_rare_event_approx": g(-math.expm1(-r.expected_markov1)),
        "n_forms": r.n_forms,
        "observed_reverse": r.observed_reverse,
        "observed_reverse_complement": r.observed_reverse_complement,
        "observed_clumps": r.observed_clumps,
        "expected_clumps_iid": g(r.expected_clumps_iid),
        f"expected_clumps_{mk}": g(r.expected_clumps_markov),
        f"clumps_obs_over_{mk}": g(_ratio(r.observed_clumps, r.expected_clumps_markov)),
        f"p_enrich_hits_{mk}": g(math.exp(lp_hits)),
        "p_enrich_clumps_iid": g(math.exp(lp_clump_iid)),
        f"p_enrich_clumps_{mk}": g(math.exp(lp_clump)),
        f"p_deplete_clumps_{mk}": g(math.exp(lp_dep)),
        f"log10_p_enrich_clumps_{mk}": g(lp_clump / math.log(10)),
        f"q_enrich_clumps_{mk}": g(r.q_enrich_clumps),
        "example_loci": ";".join(f"{c}:{p}:{ori}" for c, p, ori in r.loci),
    }


def write_csv(path: str, rows: List[Result], markov_order: int = 1) -> None:
    with open(path, "w", newline="") as fh:
        w = None
        for r in rows:
            vals = row_values(r, markov_order)
            if w is None:
                w = csv.DictWriter(fh, fieldnames=list(vals))
                w.writeheader()
            w.writerow(vals)


def write_run_info(path: str, args) -> None:
    info = dict(
        script="pattern_genome.py", version=VERSION,
        date=datetime.datetime.now().isoformat(timespec="seconds"),
        fasta=args.fasta, mappings=args.mappings, families=args.families, custom_patterns=args.pattern,
        min_word_len=args.min_word_len, max_word_len=args.max_word_len,
        strands=args.strands, markov_order=args.markov_order,
        loci=args.loci, loci_sampling=args.loci_sampling, seed=args.seed,
        include_regex=args.include_regex, exclude_regex=args.exclude_regex,
        python=platform.python_version(), numpy=np.__version__,
        **getattr(args, "_run_stats", {}),
    )
    with open(path, "w") as fh:
        json.dump(info, fh, indent=2)


def print_report(rows: List[Result], markov_order: int = 1) -> None:
    groups: Dict[Tuple[str, str, str], List[Result]] = {}
    for r in rows:
        groups.setdefault((r.mapping, r.family, r.variant), []).append(r)
    mk = f"m{markov_order}"
    for (mapping, family, variant), rr in groups.items():
        print(f"\n=== {mapping} | {family} | {variant} ===")
        print(f"label\tlen\tobserved\tE_iid\tE_markov{markov_order}\tO/E_iid\tO/E_{mk}"
              f"\tclumps\tE_clumps_{mk}\tp_clumps_{mk}")
        for r in rr:
            p = math.exp(poisson_log_sf(r.observed_clumps, r.expected_clumps_markov))
            print(f"{r.label}\t{r.length}\t{r.observed}\t{fmt(r.expected_iid)}\t"
                  f"{fmt(r.expected_markov1)}\t{fmt(_ratio(r.observed, r.expected_iid))}\t"
                  f"{fmt(_ratio(r.observed, r.expected_markov1))}\t{r.observed_clumps}\t"
                  f"{fmt(r.expected_clumps_markov)}\t{fmt(p)}")
        observed = [r for r in rr if r.observed > 0]
        if observed:
            longest = max(observed, key=lambda r: (r.length, r.order))
            print(f"Longest observed: {longest.label}, length {longest.length}, "
                  f"total hits {longest.observed}, clumps {longest.observed_clumps}")
            print(f"  word:       {longest.word}")
            print(f"  complement: {longest.complement}")
            if longest.loci:
                print("  example loci (1-based start):")
                for chrom, pos, ori in longest.loci:
                    print(f"    {chrom}:{pos}  {ori}")
        else:
            print("No requested word in this family/variant was observed.")

    observed_all = [r for r in rows if r.observed > 0]
    if observed_all:
        max_len = max(r.length for r in observed_all)
        print(f"\n=== OVERALL LONGEST OBSERVED: length {max_len} ===")
        for r in (r for r in observed_all if r.length == max_len):
            print(f"{r.mapping}\t{r.family}\t{r.variant}\t{r.label}\thits={r.observed}\t"
                  f"clumps={r.observed_clumps}\tO/E_{mk}={fmt(_ratio(r.observed, r.expected_markov1))}")


def print_family_list() -> None:
    print("Built-in families:")
    for name in DEFAULT_FAMILIES:
        print(f"  {name:18s} {FAMILY_DESCRIPTIONS[name]}")
    print("\nCustom exact binary motifs can also be supplied with --pattern NAME=BITS.")


# ------------------------------------ CLI --------------------------------------

def add_scan_arguments(ap: argparse.ArgumentParser) -> None:
    """Options shared with fibword_genome.py."""
    ap.add_argument("--mappings", nargs="+", choices=sorted(MAPPINGS), default=["RY", "WS"],
                    help="Binary DNA projections to scan (default: RY WS; add MK for all three)")
    ap.add_argument("--min-word-len", type=int, default=13, help="Smallest pattern length to count (default: 13)")
    ap.add_argument("--max-word-len", type=int, default=10000, help="Largest pattern length to consider (default: 10000)")
    ap.add_argument("--strands", choices=["both", "forward"], default="both",
                    help="both (default): also search the opposite strand; forward: word + complement only (original)")
    ap.add_argument("--markov-order", type=int, default=1,
                    help="Order of the binary Markov background model (default: 1, the original)")
    ap.add_argument("--loci", type=int, default=20,
                    help="Keep up to this many example loci per mapping/word (default: 20)")
    ap.add_argument("--loci-sampling", choices=["random", "first"], default="random",
                    help="random (default): uniform over all hits, seeded; first: first hits in file order (original)")
    ap.add_argument("--seed", type=int, default=1, help="Seed for --loci-sampling random (default: 1)")
    ap.add_argument("--include-regex", help="Only FASTA records whose full header matches this regex")
    ap.add_argument("--exclude-regex", help="Skip FASTA records whose full header matches this regex")
    ap.add_argument("--csv", help="Write detailed results to CSV (and the settings to <stem>_run.json)")
    ap.add_argument("--progress", action="store_true", help="Print per-record progress to stderr")


def parse_args(argv=None):
    ap = argparse.ArgumentParser(
        description=("Find mathematical/automatic binary word patterns in DNA under RY/WS/MK "
                     "projections and compare observed counts with IID/Markov expectations."))
    ap.add_argument("fasta", nargs="?", help="DNA FASTA, optionally gzip-compressed")
    ap.add_argument("--families", nargs="+", choices=DEFAULT_FAMILIES + ["all"], default=DEFAULT_FAMILIES,
                    help="Pattern families to scan (default: all built-ins)")
    ap.add_argument("--pattern", action="append", default=[], metavar="NAME=BITS",
                    help="Add an exact custom binary motif; may be repeated, e.g. --pattern myword=01001101")
    add_scan_arguments(ap)
    ap.add_argument("--list-families", action="store_true", help="Show built-in family definitions and exit")
    args = ap.parse_args(argv)
    if "all" in args.families:
        args.families = list(DEFAULT_FAMILIES)
    args.families = list(dict.fromkeys(args.families))
    return args


def main(argv=None):
    args = parse_args(argv)
    if args.list_families:
        print_family_list()
        return
    if not args.fasta:
        raise SystemExit("A FASTA path is required unless --list-families is used.")
    rows = scan(args)
    print_report(rows, args.markov_order)
    if args.csv:
        write_csv(args.csv, rows, args.markov_order)
        info = args.csv.rsplit(".", 1)[0] + "_run.json"
        write_run_info(info, args)
        print(f"\nWrote {args.csv} and {info}")


if __name__ == "__main__":
    main()
