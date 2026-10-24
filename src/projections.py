"""Binary projections of DNA (RY, WS, MK), the opposite-strand rule, and gaps.

This is the single definition of the binary encodings for the project. Both the
array form (codes 0..3 from seqio) and the byte form (for fast bytes.find
scanning in pattern_genome and the count scripts) live here, and
verify_identity() checks that they agree.

Encodings
---------
Codes A=0, C=1, G=2, T=3 (as produced by seqio.encode).

  RY  purine?        A,G -> 1   C,T -> 0
  WS  strong pair?   C,G -> 1   A,T -> 0
  MK  amino?         A,C -> 1   G,T -> 0

The three are the non-trivial characters of the Klein four-group on {A,C,G,T}:
RY xor WS xor MK = 0 at every base, and any two of them determine the base
exactly (unproject). One projection alone keeps half the information.

Opposite strand
---------------
Reading the other strand means reverse-complementing the DNA (A<->T, C<->G,
and reverse the order). In binary:

  WS  A<->T and C<->G stay in the same class -> the word is REVERSED
  RY  A<->T and C<->G swap purine/pyrimidine -> REVERSED and BIT-COMPLEMENTED
  MK  A<->T and C<->G swap amino/keto        -> REVERSED and BIT-COMPLEMENTED

Scanners already search a word together with its bit complement (which member
of a class is called 1 is arbitrary). Adding the opposite strand therefore means
adding the reversed word, in every projection. The full set to search is

    {w, ~w, reverse(w), ~reverse(w)}     with duplicates removed

search_set() returns exactly that. Removing duplicates matters: a word equal to
its own reverse (or to the complement of its reverse) would otherwise have
every hit counted twice. reverse_strand() and strand_variants() say which
strand and labelling each member corresponds to.

Gaps
----
project() and project_bytes() never turn a gap into a bit:
  * project() raises ValueError on any code outside 0..3 (seqio.GAP = 4,
    seqio.INVALID = 255). Project ACGT segments, e.g. with project_segments().
  * project_bytes() maps every non-ACGT character to b'N', so bytes.find can
    never match across it.

Functions
---------
project(v, name)               codes 0..3 -> uint8 bits
project_segments(x, name, L)   yield (start, bits) per ACGT stretch >= L
project_bytes(seq, name)       raw sequence bytes -> b'0'/b'1'/b'N' bytes
unproject(a, b, na, nb)        two projections -> codes 0..3
reverse_strand(word, name)     how a binary word reads on the opposite strand
strand_variants(word, name)    [(strand, labelling, word), ...] without duplicates
search_set(word, name)         the distinct words to search (both strands, both labellings)
bit_complement(word)           swap 0 and 1
verify_identity(v)             check XOR identity, invertibility, strand rule and
                               byte/array agreement on real data; raises ValueError
"""
from __future__ import annotations

from typing import Dict, Iterator, List, Tuple, Union

import numpy as np

import seqio

NAMES = ('RY', 'WS', 'MK')

# Which bases are 1 in each projection (everything else among ACGT is 0).
ONES: Dict[str, str] = {'RY': 'AG', 'WS': 'CG', 'MK': 'AC'}

# Whether the opposite strand also bit-complements the word (see docstring).
STRAND_COMPLEMENTS: Dict[str, bool] = {'RY': True, 'WS': False, 'MK': True}

Word = Union[str, bytes, np.ndarray]


# ---------------------------------------------------------------- array form

def ry(v): return ((v == 0) | (v == 2)).astype(np.uint8)   # A,G purine
def ws(v): return ((v == 1) | (v == 2)).astype(np.uint8)   # C,G strong
def mk(v): return ((v == 0) | (v == 1)).astype(np.uint8)   # A,C amino

PROJ = {'RY': ry, 'WS': ws, 'MK': mk}


def _check_name(name: str) -> None:
    if name not in PROJ:
        raise ValueError(f'unknown projection {name!r}; expected one of {NAMES}')


def project(v: np.ndarray, name: str) -> np.ndarray:
    """Project ACGT codes (0..3) to a uint8 0/1 array.

    Raises ValueError if v contains any other code, so a gap (seqio.GAP or
    seqio.INVALID) can never silently become a 0 bit.
    """
    _check_name(name)
    v = np.asarray(v)
    if v.size and (v.min() < 0 or v.max() > 3):
        raise ValueError(
            f'project({name}): input contains codes outside 0..3 (gaps/ambiguous bases). '
            'Project ACGT segments instead, e.g. with project_segments().')
    return PROJ[name](v)


def project_segments(x, name: str, min_len: int = 1) -> Iterator[Tuple[int, np.ndarray]]:
    """Yield (start, bits) for every ACGT stretch of length >= min_len.

    x is raw sequence bytes or codes from seqio.encode(). start is the 0-based
    position of bits[0] in the record, so hits map back to genome coordinates.
    """
    _check_name(name)
    for start, seg in seqio.iter_segments(x, min_len):
        yield start, PROJ[name](seg)


def unproject(a: np.ndarray, b: np.ndarray, na: str, nb: str) -> np.ndarray:
    """Rebuild ACGT codes from two different projections of the same sequence."""
    _check_name(na); _check_name(nb)
    if na == nb:
        raise ValueError('unproject needs two different projections')
    a = np.asarray(a); b = np.asarray(b)
    if a.shape != b.shape:
        raise ValueError(f'unproject: length mismatch {a.shape} vs {b.shape}')
    for arr, nm in ((a, na), (b, nb)):
        if arr.size and not np.isin(arr, (0, 1)).all():
            raise ValueError(f'unproject: {nm} input contains values other than 0/1')
    table = {}
    for code in range(4):
        vv = np.array([code], np.uint8)
        table[(int(PROJ[na](vv)[0]), int(PROJ[nb](vv)[0]))] = code
    out = np.empty(len(a), np.uint8)
    for (x, y), code in table.items():
        out[(a == x) & (b == y)] = code
    return out


# ----------------------------------------------------------------- byte form

# Byte maps used by the bytes.find scanners: base -> b'1' / b'0'.
MAPPINGS: Dict[str, Dict[int, int]] = {
    name: {ord(base): ord('1') if base in ONES[name] else ord('0') for base in 'ACGT'}
    for name in NAMES
}


def make_translate_table(mapping: Dict[int, int]) -> bytes:
    """256-byte table for bytes.translate. Upper and lower case ACGT map to
    b'0'/b'1'; every other byte maps to b'N' so no match can cross it."""
    table = bytearray(b'N' * 256)
    for base, bit in mapping.items():
        table[base] = bit
        if 65 <= base <= 90:
            table[base + 32] = bit
    return bytes(table)


TRANSLATE_TABLES: Dict[str, bytes] = {name: make_translate_table(MAPPINGS[name]) for name in NAMES}


def project_bytes(seq: bytes, name: str) -> bytes:
    """Raw sequence bytes -> b'0'/b'1' bytes, with b'N' at every non-ACGT position."""
    _check_name(name)
    return seq.translate(TRANSLATE_TABLES[name])


# ------------------------------------------------------- words and strands

def bit_complement(word: Word) -> Word:
    """Swap 0 and 1. Works on str, bytes or a 0/1 numpy array."""
    if isinstance(word, str):
        return word.translate(str.maketrans('01', '10'))
    if isinstance(word, (bytes, bytearray)):
        return bytes(word).translate(bytes.maketrans(b'01', b'10'))
    return (1 - np.asarray(word)).astype(np.uint8)


def _reverse(word: Word) -> Word:
    if isinstance(word, (str, bytes, bytearray)):
        return word[::-1]
    return np.asarray(word)[::-1].copy()


def reverse_strand(word: Word, name: str) -> Word:
    """How a binary word reads on the opposite DNA strand in projection `name`.

    WS: reversed. RY, MK: reversed and bit-complemented.
    """
    _check_name(name)
    r = _reverse(word)
    return bit_complement(r) if STRAND_COMPLEMENTS[name] else r


def _key(word: Word):
    return word.tobytes() if isinstance(word, np.ndarray) else word


def strand_variants(word: Word, name: str) -> List[Tuple[str, str, Word]]:
    """All distinct forms of `word` to search, labelled by origin.

    Returns [(strand, labelling, word), ...] where strand is '+' or '-' and
    labelling is 'as_is' or 'complement'. Forms that coincide are listed once
    (the first label wins), so no hit is counted twice.
    """
    _check_name(name)
    rs = reverse_strand(word, name)
    candidates = [('+', 'as_is', word), ('+', 'complement', bit_complement(word)),
                  ('-', 'as_is', rs), ('-', 'complement', bit_complement(rs))]
    seen, out = set(), []
    for strand, lab, w in candidates:
        k = _key(w)
        if k not in seen:
            seen.add(k); out.append((strand, lab, w))
    return out


def search_set(word: Word, name: str) -> List[Word]:
    """The distinct words to search so both strands and both 0/1 labellings are
    covered: {w, ~w, reverse(w), ~reverse(w)} without duplicates."""
    return [w for _, _, w in strand_variants(word, name)]


# -------------------------------------------------------------- self-check

def verify_identity(v: np.ndarray) -> bool:
    """Check the encoding rules on real data (codes 0..3). Raises ValueError.

    1. RY xor WS xor MK == 0 at every position
    2. every pair of projections inverts back to the original codes
    3. projecting the reverse complement equals reverse_strand() of the projection
    4. the byte form (project_bytes) agrees with the array form (project)
    Uses explicit checks, not assert, so it still runs under python -O.
    """
    v = np.asarray(v)
    p = {name: project(v, name) for name in NAMES}
    if not np.all((p['RY'] ^ p['WS'] ^ p['MK']) == 0):
        raise ValueError('Klein identity RY^WS^MK == 0 broken')
    for na, nb in (('RY', 'WS'), ('RY', 'MK'), ('WS', 'MK')):
        if not np.array_equal(unproject(p[na], p[nb], na, nb), v):
            raise ValueError(f'unproject({na},{nb}) does not recover the sequence')
    revcomp = (3 - v[::-1]).astype(v.dtype)          # A0<->T3, C1<->G2, reversed
    raw = np.frombuffer(b'ACGT', np.uint8)[v].tobytes()
    for name in NAMES:
        if not np.array_equal(project(revcomp, name), reverse_strand(p[name], name)):
            raise ValueError(f'opposite-strand rule broken for {name}')
        b = np.frombuffer(project_bytes(raw, name), np.uint8) - ord('0')
        if not np.array_equal(b, p[name]):
            raise ValueError(f'byte and array forms of {name} disagree')
    return True
