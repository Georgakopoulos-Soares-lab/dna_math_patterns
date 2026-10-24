"""Decoy words: controls that keep a pattern's SHAPE but destroy its rule.

Why decoys
----------
The enriched hits sit in tandem repeats, so a word looks "enriched" against an
i.i.d. or Markov null whenever it is repeat-shaped. To ask whether the
mathematical RULE matters, compare each word with decoys that look like it but
were not produced by that rule. If the real word scores like its decoys, the
enrichment is explained by its shape, not by the rule.

Methods
-------
  run (default, strict)
      Split the word into runs (maximal blocks of equal bits). Shuffle the
      lengths of the 0-runs among the 0-run positions and the lengths of the
      1-runs among the 1-run positions. Keeps EXACTLY: length, 0/1 composition,
      number of runs, and the multiset of 0-run and 1-run lengths (so also the
      longest run). Destroys only their order. For the Fibonacci word this is
      every arrangement of the same blocks '01' and '001'.

  segment (original method)
      Cut the word at every "slip" (two equal adjacent bits), permute the
      alternating segments, re-join so that each joint is again a slip. Keeps
      length, number of slips and the multiset of segment lengths. It does NOT
      by itself keep the composition or the run lengths (a segment of length 1
      placed next to another can create a longer run), hence the filters below
      are on by default. Few segments -> few decoys (Tribonacci T6: 4 segments,
      at most 5 distinct decoys).

  slip (weaker)
      Same length and number of slips, slip positions drawn at random. Also
      destroys near-periodicity.

Filters (all methods; defaults shown)
-------------------------------------
  match_composition=True   number of 1s equals the real word's number of 1s or
                           of 0s (the complement is always searched too)
  match_runs='max'         'max': no run longer than the real word's longest;
                           'exact': same multiset of run lengths (up to 0/1
                           swap); None: no run filter
For method='run' both always hold by construction.

Equivalent decoys (both strands, both labellings)
-------------------------------------------------
The scanners search a word together with its bit complement, and -- once the
opposite strand is searched -- its reverse and reversed complement (see
projections.search_set; the set {w, ~w, rev w, ~rev w} is the same in RY, WS
and MK). Decoys in the same group give identical counts, so only one per group
is kept, and the real word's group is excluded. canonical(w) names the group.
strands=False uses {w, ~w} only (complement, no reversal).

How many decoys exist
---------------------
If the number of arrangements a method can produce is at most
`enumerate_limit` (default 100,000), all of them are listed, so the number of
distinct decoys is EXACT, and the n returned are drawn without replacement.
Otherwise decoys are sampled at random until n distinct ones are found (or
`max_tries` is reached). Measured for the five genome-wide targets (distinct,
both strands):

  target               run (default)   segment + filters
  Tribonacci T6              669                5
  Tribonacci T7          >19,000              209
  period-doubling n5         461              159
  Fibonacci F7             1,286               69
  Lucas L6                    65 (all)         46

The slip method's space is too large to enumerate for most targets, and with
the filters on few random slip patterns pass (e.g. 19 for F7 in 100,000
draws); raise max_tries or turn the filters off if more are needed.

When fewer than n decoys exist, all of them are returned and a warning states
the count and the smallest achievable empirical p-value, 1/(N+1). Lucas L6 has
exactly 65 words of its shape, so its p-value floor is 1/66 = 0.015.

Interface
---------
make_decoys(w, n, rng, method='run', ...)  -> DecoySet (words + counts)
run_permutation_decoys(w, n, rng, ...)     -> list of words
segment_permutation_decoys(w, n, rng, ...) -> list of words (original signature)
random_slip_decoys(w, n, rng, ...)         -> list of words (original signature)
canonical(w, strands=True)                 group representative
empirical_p(real, decoys, tail='greater')  (1 + #{decoy >= real}) / (N + 1)
segments(w), rejoin(segs), n_slips(w), runs(w)

Tests: tests/test_decoy_words.py
"""
from __future__ import annotations

import itertools
import math
import warnings
from collections import Counter
from typing import Iterable, Iterator, List, NamedTuple, Optional, Sequence, Tuple

import numpy as np

FLIP = str.maketrans('01', '10')
METHODS = ('run', 'segment', 'slip')


# ------------------------------------------------------------- word helpers

def flip(w: str) -> str:
    return w.translate(FLIP)


def canonical(w: str, strands: bool = True) -> str:
    """Representative of w's group: min of {w, ~w, rev w, ~rev w} (strands=True)
    or of {w, ~w} (strands=False)."""
    forms = [w, flip(w)]
    if strands:
        r = w[::-1]
        forms += [r, flip(r)]
    return min(forms)


def runs(w: str) -> List[Tuple[str, int]]:
    """[(bit, length), ...] of the maximal runs of equal bits."""
    return [(k, len(list(g))) for k, g in itertools.groupby(w)]


def max_run(w: str) -> int:
    return max(n for _, n in runs(w))


def run_signature(w: str) -> Tuple[Tuple[int, ...], Tuple[int, ...]]:
    """(sorted 0-run lengths, sorted 1-run lengths)."""
    r = runs(w)
    return (tuple(sorted(n for b, n in r if b == '0')), tuple(sorted(n for b, n in r if b == '1')))


def segments(w: str) -> List[str]:
    """Cut w at every slip (two equal adjacent bits); each piece alternates."""
    cuts = [i for i in range(1, len(w)) if w[i] == w[i - 1]]
    b = [0] + cuts + [len(w)]
    return [w[b[i]:b[i + 1]] for i in range(len(b) - 1)]


def rejoin(segs: Sequence[str]) -> str:
    """Concatenate segments, complementing each one if needed so that every
    joint is a slip."""
    out = segs[0]
    for s in segs[1:]:
        if s[0] != out[-1]:
            s = flip(s)
        out += s
    return out


def n_slips(w: str) -> int:
    return sum(w[i] == w[i + 1] for i in range(len(w) - 1))


# ------------------------------------------------------------------ filters

def _composition_ok(d: str, ones: int, L: int) -> bool:
    return d.count('1') in (ones, L - ones)


def _runs_ok(d: str, w: str, mode: Optional[str], w_max: int, w_sig) -> bool:
    if mode is None:
        return True
    if mode == 'max':
        return max_run(d) <= w_max
    sig = run_signature(d)
    return sig == w_sig or (sig[1], sig[0]) == w_sig


# ------------------------------------------------------------ combinatorics

def _n_distinct_perms(items: Sequence) -> int:
    c = Counter(items)
    out = math.factorial(len(items))
    for v in c.values():
        out //= math.factorial(v)
    return out


def _distinct_perms(items: Sequence) -> Iterator[Tuple]:
    """All distinct permutations of a multiset, in lexicographic order."""
    a = sorted(items)
    n = len(a)
    while True:
        yield tuple(a)
        i = n - 2
        while i >= 0 and a[i] >= a[i + 1]:
            i -= 1
        if i < 0:
            return
        j = n - 1
        while a[j] <= a[i]:
            j -= 1
        a[i], a[j] = a[j], a[i]
        a[i + 1:] = reversed(a[i + 1:])


# --------------------------------------------------------- method backends
# Each backend returns (space_size, enumerate_all(), draw_one(rng)).

def _run_backend(w: str):
    r = runs(w)
    first = r[0][0]
    other = flip(first)
    a = [n for i, (_, n) in enumerate(r) if i % 2 == 0]     # runs of bit `first`
    b = [n for i, (_, n) in enumerate(r) if i % 2 == 1]     # runs of the other bit

    def build(pa, pb):
        out = []
        for i in range(len(r)):
            out.append(first * pa[i // 2] if i % 2 == 0 else other * pb[i // 2])
        return ''.join(out)

    def enum():
        for pa in _distinct_perms(a):
            for pb in _distinct_perms(b):
                yield build(pa, pb)

    def draw(rng):
        return build([a[i] for i in rng.permutation(len(a))], [b[i] for i in rng.permutation(len(b))])

    return _n_distinct_perms(a) * _n_distinct_perms(b), enum, draw


def _segment_backend(w: str):
    segs = segments(w)

    def enum():
        for p in _distinct_perms(segs):
            yield rejoin(p)

    def draw(rng):
        return rejoin([segs[i] for i in rng.permutation(len(segs))])

    return _n_distinct_perms(segs), enum, draw


def _slip_backend(w: str):
    L, k = len(w), n_slips(w)
    start = int(w[0])

    def build(slips):
        b = [start]
        for i in range(1, L):
            b.append(b[-1] if i in slips else 1 - b[-1])
        return ''.join(map(str, b))

    def enum():
        for c in itertools.combinations(range(1, L), k):
            yield build(set(c))

    def draw(rng):
        return build(set(rng.choice(np.arange(1, L), size=k, replace=False).tolist()))

    return math.comb(L - 1, k), enum, draw


_BACKENDS = {'run': _run_backend, 'segment': _segment_backend, 'slip': _slip_backend}


# ------------------------------------------------------------------- driver

class DecoySet(NamedTuple):
    words: List[str]            # the decoys, one per equivalence group
    method: str
    n_requested: int
    n_possible: Optional[int]   # exact number of distinct decoys if enumerated, else None
    exhaustive: bool            # True if every arrangement was examined
    p_floor: float              # smallest achievable empirical p-value, 1/(len(words)+1)

    @property
    def n_returned(self) -> int:
        return len(self.words)


def make_decoys(w: str, n: int, rng: np.random.Generator, method: str = 'run',
                match_composition: bool = True, match_runs: Optional[str] = 'max',
                strands: bool = True, enumerate_limit: int = 100_000,
                max_tries: Optional[int] = None, min_decoys: int = 19,
                warn: bool = True) -> DecoySet:
    """Up to n distinct decoys of binary word w (a str of '0'/'1').

    method        'run' (default), 'segment' or 'slip' -- see module docstring
    match_*       filters (module docstring); ignored where they hold by construction
    strands       treat reverse and reversed complement as equivalent too
    enumerate_limit  list every arrangement when there are at most this many
    max_tries     random draws allowed when sampling (default 1000 * n)
    min_decoys    warn if fewer decoys than this exist (19 -> p < 0.05 possible)
    """
    if method not in METHODS:
        raise ValueError(f'method must be one of {METHODS}')
    if match_runs not in (None, 'max', 'exact'):
        raise ValueError("match_runs must be None, 'max' or 'exact'")
    if not w or set(w) - {'0', '1'}:
        raise ValueError('w must be a non-empty string of 0/1')
    L, ones = len(w), w.count('1')
    w_max, w_sig = max_run(w), run_signature(w)
    real_key = canonical(w, strands)

    def keep(d):
        return ((not match_composition or _composition_ok(d, ones, L))
                and _runs_ok(d, w, match_runs, w_max, w_sig))

    space, enum, draw = _BACKENDS[method](w)
    found = {}
    if space <= enumerate_limit:
        for d in enum():
            key = canonical(d, strands)
            if key != real_key and key not in found and keep(d):
                found[key] = d
        pool = [found[k] for k in sorted(found)]
        n_possible, exhaustive = len(pool), True
        if len(pool) > n:
            pool = [pool[i] for i in sorted(rng.choice(len(pool), size=n, replace=False).tolist())]
        words = pool
    else:
        tries = max_tries if max_tries is not None else 1000 * n
        for _ in range(tries):
            d = draw(rng)
            key = canonical(d, strands)
            if key != real_key and key not in found and keep(d):
                found[key] = d
                if len(found) >= n:
                    break
        words = list(found.values())
        n_possible, exhaustive = None, False

    ds = DecoySet(words, method, n, n_possible, exhaustive, 1.0 / (len(words) + 1))
    if warn and len(words) < n:
        what = (f'only {n_possible} distinct decoys exist' if exhaustive
                else f'only {len(words)} distinct decoys found in the allowed random draws')
        msg = (f'{method} decoys for a {L}-bit word: {what} (requested {n}); '
               f'smallest achievable empirical p = 1/{len(words) + 1} = {ds.p_floor:.3f}')
        if len(words) < min_decoys:
            msg += ' -- p < 0.05 is impossible with this many decoys'
        warnings.warn(msg, RuntimeWarning, stacklevel=2)
    return ds


# ------------------------------------------------- convenience wrappers

def run_permutation_decoys(w: str, n: int, rng: np.random.Generator, **kw) -> List[str]:
    """Strict control: same length, composition and run lengths, order shuffled."""
    return make_decoys(w, n, rng, method='run', **kw).words


def segment_permutation_decoys(w: str, n: int, rng: np.random.Generator,
                               match_composition: bool = True, match_runs: Optional[str] = 'max',
                               **kw) -> List[str]:
    """Original method, now without equivalent duplicates and, by default,
    with the composition and longest-run filters."""
    return make_decoys(w, n, rng, method='segment', match_composition=match_composition,
                       match_runs=match_runs, **kw).words


def random_slip_decoys(w: str, n: int, rng: np.random.Generator,
                       match_composition: bool = True, match_runs: Optional[str] = 'max',
                       **kw) -> List[str]:
    """Weaker control: same length and number of slips, slip positions random;
    distinct, real word excluded, filters as above."""
    return make_decoys(w, n, rng, method='slip', match_composition=match_composition,
                       match_runs=match_runs, **kw).words


# --------------------------------------------------------------- statistics

def empirical_p(real: float, decoys: Iterable[float], tail: str = 'greater') -> float:
    """Empirical p-value with the +1 correction: (1 + #{decoy >= real}) / (N + 1).

    tail='greater' tests enrichment; tail='less' uses #{decoy <= real} (depletion).
    """
    d = np.asarray(list(decoys), dtype=float)
    if tail == 'greater':
        k = int((d >= real).sum())
    elif tail == 'less':
        k = int((d <= real).sum())
    else:
        raise ValueError("tail must be 'greater' or 'less'")
    return (1 + k) / (len(d) + 1)
