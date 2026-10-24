"""Is a word's excess of loci explained by it being NEARLY PERIODIC?

For each target word (targets.py; default: all 36 words, 13-100 bits) find its
best periodic approximation: the period p and unit u (majority vote per residue
class) minimising the Hamming distance d between the word and the tiling
u u u ... Then count, genome-wide and exactly as count_genome.py does, the real
word, its tiling, and two kinds of decoys:

  random_unit  same p and d, a DIFFERENT random unit. Asks: is this repeat more
               common than other repeats of the same period with as many defects?
  same_unit    the SAME unit, the d defects moved to other positions. Asks: does
               the word's particular defect pattern -- where its rule shows --
               beat other defective copies of its own repeat?
Only same_unit tests the rule against the word's own periodic approximation;
random_unit tests the repeat unit against other units. Both run by default.

Nearest period (profile)
------------------------
Periods p = 1 .. min(--pmax, L // --min-repeats) are tried (default pmax 16,
min-repeats 3). Requiring >= 3 repeats matters: any word "repeats" with a period
close to its own length (only the overlap is compared), so with no such
requirement 17 of the 36 words looked perfectly periodic (e.g. Fibonacci F5,
13 bits: "period 8, 0 defects"). Ties in d go to the smaller p; a residue class
with as many 0s as 1s votes 0 (as in the original script).
The profile reports d/L; a word with d/L > --max-defect-rate (default 0.10) is
flagged near_periodic = False -- the test is not meaningful for it, but it is
still run.

Decoy rules (both modes)
------------------------
  * same length L, period p and number of defects d as the real word
  * same 0/1 composition: number of 1s equals the real word's 1s or 0s
    (the complement is always searched too)
  * no run longer than the real word's longest run -- dropped with --no-run-rule
    (composition is still kept). Use it where the strict rules leave too few
    decoys: period-doubling n5 has 6 random-unit and 3 same-unit decoys with
    all rules, 48 and 11 without the run rule (smallest p 1/7 -> 1/49). Report
    which rule set was used (summary column run_rule).
  * random_unit: primitive units only (01010101 is period 2, not 8), never the
    real unit or any rotation, reversal or complement of it (those give the same
    tandem repeat)
  * one decoy per group {w, ~w, reverse, reversed complement}, never the real
    word's group (decoy_words.canonical)
Valid defect placements are counted and drawn by a small dynamic program, so
every draw satisfies all rules (uniform over valid (unit, placement) pairs).
When there are at most --enumerate-limit of them they are all listed and the
number of distinct decoys is exact; otherwise they are sampled. Repeat units are scarce at small periods: up to rotation, reversal
and complement there are only 14 distinct primitive 8-bit units, so the random
-unit decoys of period-doubling n5 (p = 8) come from 13 other repeats. The
summary reports n_units_used and n_units_available.

Counting and statistics are those of count_genome.py (one-pass search, both
strands by default, loci = hits chained within --locus-gap bp):
  frac_decoys_loci_ge_real (original column), p_emp_loci = (1 + #{decoys >= real})
  / (N + 1), q_emp_loci (Benjamini-Hochberg over all rows of the run), p_floor,
  testable (N >= 19), and tiling_loci -- the loci of the exact tiling u u u ...

usage
-----
  python3 src/periodic_decoy.py <chrom_dir | fasta> <n_decoys> <out.csv> [options]
  options (defaults in brackets)
     --targets all|readme|NAME ... [all]   --min-len [13] --max-len [100] --families
     --modes random_unit same_unit [both]  --pmax [16]  --min-repeats [3]
     --max-defect-rate [0.10]  --mappings [RY WS MK]  --strands [both]
     --locus-gap [1000]  --seed [23]  --workers [12]  --chroms ...  --exclude-chrom REGEX
     --enumerate-limit [200000]  --no-run-rule  --progress

Output
------
  <out.csv>           per word x mapping: wid, mapping, hits, loci, pattern, kind
                      (REAL / TILING / DECOY), word, mode, ...
  <out>_profile.csv   per target: length, best_period, unit, hamming_defects, ones
                      (original columns), repeats, defect_rate, near_periodic
  <out>_summary.csv   per target x mode x mapping: real vs decoys
  <out>_run.json      settings, decoy counts, versions
"""
from __future__ import annotations

import argparse
import datetime
import itertools
import json
import os
import platform
import sys
import time
from multiprocessing import Pool
from functools import lru_cache
from typing import Dict, List, NamedTuple, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import count_genome as CG     # noqa: E402
import decoy_words as DW      # noqa: E402
import pattern_genome as PG   # noqa: E402
import targets as TG          # noqa: E402

VERSION = '2.0'
MODES = ('random_unit', 'same_unit')
FLIP = str.maketrans('01', '10')


# ------------------------------------------------------------- periodicity

class PeriodFit(NamedTuple):
    d: int          # Hamming distance to the tiling
    p: int          # period
    unit: str       # repeat unit (length p)


def best_period(w: str, pmax: int = 16, min_repeats: int = 3) -> PeriodFit:
    """Nearest periodic word: smallest d over p = 1 .. min(pmax, L // min_repeats)."""
    L = len(w)
    a = np.array([int(c) for c in w])
    best = None
    for p in range(1, max(1, min(pmax, L // max(min_repeats, 1))) + 1):
        u = [int(round(a[r::p].mean())) for r in range(p)]
        d = int(sum(a[i] != u[i % p] for i in range(L)))
        if best is None or d < best[0]:
            best = PeriodFit(d, p, ''.join(map(str, u)))
    return best


def tile(u: str, L: int) -> str:
    return (u * (L // len(u) + 1))[:L]


def primitive_root(u: str) -> str:
    p = len(u)
    for k in range(1, p + 1):
        if p % k == 0 and u[:k] * (p // k) == u:
            return u[:k]
    return u


def canonical_unit(u: str) -> str:
    """One name per tandem repeat: min over rotations of u, its reverse and complements."""
    u = primitive_root(u)
    forms = []
    for x in (u, u[::-1]):
        for y in (x, x.translate(FLIP)):
            forms += [y[i:] + y[:i] for i in range(len(y))]
    return min(forms)


@lru_cache(maxsize=None)
def _units_of_period(p: int) -> Tuple[Tuple[str, str], ...]:
    """((canonical name, representative), ...) of all primitive non-constant p-bit units."""
    seen = {}
    for bits in itertools.product('01', repeat=p):
        u = ''.join(bits)
        if len(set(u)) < 2 or primitive_root(u) != u:
            continue
        seen.setdefault(canonical_unit(u), u)
    return tuple(sorted(seen.items()))


def available_units(p: int, real_unit: str) -> List[str]:
    """Primitive, non-constant p-bit units, one per tandem repeat, excluding the real one."""
    real = canonical_unit(real_unit)
    return [u for c, u in _units_of_period(p) if c != real]


# ---------------------------------------------------------------- decoys

class PeriodicDecoys(NamedTuple):
    words: List[str]
    units: List[str]              # unit of each decoy
    mode: str
    n_possible: Optional[int]     # distinct decoys, exact when enumerated
    exhaustive: bool
    n_units_available: int
    n_candidates: int = 0         # valid (unit, defect placement) pairs, before merging equivalents


def _flip(base: str, positions) -> str:
    b = list(base)
    for i in positions:
        b[i] = '1' if b[i] == '0' else '0'
    return ''.join(b)


class _Placements:
    """Defect placements on one tiling `base` that give a valid decoy.

    Exactly A flips 0->1 and B flips 1->0 (this fixes the composition) and no run
    longer than m. A dynamic program over positions counts the valid placements
    (`total`) and draws them uniformly (`sample`) or lists them (`enumerate`), so
    no random draw is ever rejected.
    State at position i: (a, b, previous bit (2 = none), current run length r).
    """

    def __init__(self, base: str, A: int, B: int, m: Optional[int]):
        # m = None: no limit on run length (the run-length state is then not needed)
        self.unlimited = m is None
        m = 1 if m is None else m
        self.base, self.A, self.B, self.m, self.L = base, A, B, m, len(base)
        self._full = None
        self.total = float(self._backward(keep=False)[0][0, 0, 2, 0])

    def _moves(self, i):
        ti = self.base[i]
        for x in (0, 1):
            flip = str(x) != ti
            yield x, int(flip and ti == '0'), int(flip and ti == '1')

    def _backward(self, keep: bool):
        A, B, m = self.A, self.B, self.m
        nxt = np.zeros((A + 1, B + 1, 3, m + 1))
        nxt[A, B] = 1.0
        out = [nxt] if keep else None
        for i in range(self.L - 1, -1, -1):
            cur = np.zeros_like(nxt)
            for x, da, db in self._moves(i):
                sh = np.zeros((A + 1, B + 1, m + 1))
                sh[:A + 1 - da, :B + 1 - db] = nxt[da:, db:, x]
                diff = sh[:, :, 1][:, :, None]                 # other bit / start: run = 1
                if self.unlimited:
                    cur[:, :, x] += diff                       # run length not tracked
                else:
                    cur[:, :, x, :m] += sh[:, :, 1:]           # same bit: run grows
                cur[:, :, 1 - x] += diff
                cur[:, :, 2] += diff
            nxt = cur
            if keep:
                out.append(cur)
        return out[::-1] if keep else [nxt]

    def _table(self):
        if self._full is None:
            self._full = self._backward(keep=True)
        return self._full

    def _options(self, N, i, a, b, bit, r):
        opts = []
        for x, da, db in self._moves(i):
            r2 = 1 if (self.unlimited or x != bit) else r + 1
            if a + da > self.A or b + db > self.B or r2 > self.m:
                continue
            wgt = N[i + 1][a + da, b + db, x, r2]
            if wgt > 0:
                opts.append((wgt, x, a + da, b + db, r2))
        return opts

    def sample(self, rng) -> str:
        if self.total <= 0:
            raise ValueError('no valid placement to sample from')
        N = self._table()
        a = b = r = 0
        bit, out = 2, []
        for i in range(self.L):
            opts = self._options(N, i, a, b, bit, r)
            tot = sum(o[0] for o in opts)
            u = rng.random() * tot
            for wgt, x, a2, b2, r2 in opts:
                u -= wgt
                if u <= 0:
                    break
            out.append(str(x)); a, b, bit, r = a2, b2, x, r2
        return ''.join(out)

    def enumerate(self):
        N = self._table()
        stack = [(0, 0, 0, 2, 0, '')]
        while stack:
            i, a, b, bit, r, pre = stack.pop()
            if i == self.L:
                yield pre
                continue
            for _, x, a2, b2, r2 in self._options(N, i, a, b, bit, r):
                stack.append((i + 1, a2, b2, x, r2, pre + str(x)))


def make_periodic_decoys(w: str, fit: PeriodFit, n: int, rng, mode: str,
                         enumerate_limit: int = 200_000, max_tries: Optional[int] = None,
                         run_rule: bool = True) -> PeriodicDecoys:
    """Up to n distinct decoys of w with the same period, number of defects,
    composition and (run_rule=True) no run longer than w's (see module docstring)."""
    if mode not in MODES:
        raise ValueError(f'mode must be one of {MODES}')
    L, d = len(w), fit.d
    ones = w.count('1')
    m = DW.max_run(w) if run_rule else None
    real_key = DW.canonical(w)
    units = [fit.unit] if mode == 'same_unit' else available_units(fit.p, fit.unit)

    strata = []                                   # (unit, placements)
    for u in units:
        base = tile(u, L)
        t = base.count('1')
        for f in sorted({ones, L - ones}):
            if (d + f - t) % 2:
                continue
            A = (d + f - t) // 2
            B = d - A
            if A < 0 or B < 0 or A > L - t or B > t:
                continue
            pl = _Placements(base, A, B, m)
            if pl.total > 0:
                strata.append((u, pl))
    total = sum(pl.total for _, pl in strata)
    found: Dict[str, Tuple[str, str]] = {}
    if total <= enumerate_limit:
        for u, pl in strata:
            for x in pl.enumerate():
                k = DW.canonical(x)
                if k != real_key and k not in found:
                    found[k] = (x, u)
        pool = [found[k] for k in sorted(found)]
        n_possible, exhaustive = len(pool), True
        if len(pool) > n:
            pool = [pool[i] for i in sorted(rng.choice(len(pool), size=n, replace=False).tolist())]
    else:
        weights = np.array([pl.total for _, pl in strata])
        weights = weights / weights.sum()
        tries = max_tries if max_tries is not None else 50 * n
        for _ in range(tries):
            u, pl = strata[int(rng.choice(len(strata), p=weights))]
            x = pl.sample(rng)
            k = DW.canonical(x)
            if k != real_key and k not in found:
                found[k] = (x, u)
                if len(found) >= n:
                    break
        pool = list(found.values())
        n_possible, exhaustive = None, False
    return PeriodicDecoys([x for x, _ in pool], [u for _, u in pool], mode, n_possible, exhaustive,
                          len(units), int(round(total)))


# ------------------------------------------------------------------ build

def build(tgts, n_decoys, modes, pmax, min_repeats, max_defect_rate, seed, enumerate_limit,
          run_rule=True):
    rows, prof, info = [], [], {}
    for tag, pinfo in tgts.items():
        w = pinfo.word.decode()
        fit = best_period(w, pmax, min_repeats)
        L = len(w)
        prof.append(dict(pattern=tag, length=L, best_period=fit.p, unit=fit.unit,
                         hamming_defects=fit.d, ones=w.count('1'),
                         repeats=round(L / fit.p, 2), defect_rate=round(fit.d / L, 4),
                         near_periodic=bool(fit.d / L <= max_defect_rate),
                         family=pinfo.family, variant=pinfo.variant, label=pinfo.label))
        base = dict(family=pinfo.family, variant=pinfo.variant, label=pinfo.label, length=L,
                    best_period=fit.p, defects=fit.d)
        rows.append(dict(wid=len(rows), pattern=tag, kind='REAL', mode='', unit=fit.unit, word=w, **base))
        rows.append(dict(wid=len(rows), pattern=tag, kind='TILING', mode='', unit=fit.unit,
                         word=tile(fit.unit, L), **base))
        for mode in modes:
            ds = make_periodic_decoys(w, fit, n_decoys, TG.target_rng(seed, f'{tag}|{mode}'), mode,
                                      enumerate_limit, run_rule=run_rule)
            info[(tag, mode)] = dict(n_decoys=len(ds.words), n_possible=ds.n_possible,
                                     n_candidates=ds.n_candidates, run_rule=run_rule,
                                     exhaustive=ds.exhaustive, n_units_available=ds.n_units_available,
                                     n_units_used=len({canonical_unit(u) for u in ds.units}))
            if len(ds.words) < n_decoys:
                what = f'{ds.n_possible} exist' if ds.exhaustive else f'{len(ds.words)} found'
                print(f'warning: {tag} {mode}: only {what} (requested {n_decoys}); '
                      f'p_emp >= 1/{len(ds.words) + 1}', file=sys.stderr)
            for x, u in zip(ds.words, ds.units):
                rows.append(dict(wid=len(rows), pattern=tag, kind='DECOY', mode=mode, unit=u, word=x, **base))
    return pd.DataFrame(rows), pd.DataFrame(prof), info


def summarise(g: pd.DataFrame, prof: pd.DataFrame, info) -> pd.DataFrame:
    pr = prof.set_index('pattern')
    rows = []
    for (tag, m), x in g.groupby(['pattern', 'mapping'], sort=False):
        r = x[x.kind == 'REAL'].iloc[0]
        tl = x[x.kind == 'TILING'].iloc[0]
        for mode in sorted(set(x.loc[x.kind == 'DECOY', 'mode']) | {k[1] for k in info if k[0] == tag}):
            d = x[(x.kind == 'DECOY') & (x['mode'] == mode)]
            n = len(d)
            i = info[(tag, mode)]
            rows.append(dict(
                pattern=tag, mapping=m, real_loci=int(r.loci), n_decoys=n,
                decoy_loci_median=d.loci.median() if n else np.nan,
                decoy_loci_max=d.loci.max() if n else np.nan,
                frac_decoys_loci_ge_real=round(float((d.loci >= r.loci).mean()), 3) if n else np.nan,
                p_emp_loci=DW.empirical_p(r.loci, d.loci) if n else np.nan,
                decoy_mode=mode, run_rule=i['run_rule'], tiling_loci=int(tl.loci),
                best_period=int(r.best_period), defects=int(r.defects),
                defect_rate=pr.loc[tag, 'defect_rate'], near_periodic=bool(pr.loc[tag, 'near_periodic']),
                n_decoys_possible=i['n_possible'], n_candidates=i['n_candidates'],
                n_units_used=i['n_units_used'],
                n_units_available=i['n_units_available'],
                p_floor=1.0 / (n + 1), testable=bool(n >= 19),
                real_hits=int(r.hits), decoy_hits_median=d.hits.median() if n else np.nan,
                p_emp_hits=DW.empirical_p(r.hits, d.hits) if n else np.nan,
                family=r.family, variant=r.variant, label=r.label, length=int(r.length),
            ))
    s = pd.DataFrame(rows)
    s.insert(s.columns.get_loc('p_emp_loci') + 1, 'q_emp_loci', PG.bh_qvalues(s.p_emp_loci.tolist()))
    return s


# -------------------------------------------------------------------- CLI

def parse_args(argv=None):
    ap = argparse.ArgumentParser(description='Real vs near-periodic decoys, genome-wide locus counts.')
    ap.add_argument('chrom_dir', help='folder with one FASTA per chromosome, or a single FASTA file')
    ap.add_argument('n_decoys', type=int)
    ap.add_argument('out')
    ap.add_argument('--targets', nargs='+', default=['all'])
    ap.add_argument('--families', nargs='+', default=list(PG.DEFAULT_FAMILIES), choices=PG.DEFAULT_FAMILIES)
    ap.add_argument('--min-len', type=int, default=13)
    ap.add_argument('--max-len', type=int, default=100)
    ap.add_argument('--modes', nargs='+', default=list(MODES), choices=MODES)
    ap.add_argument('--pmax', type=int, default=16)
    ap.add_argument('--min-repeats', type=int, default=3)
    ap.add_argument('--max-defect-rate', type=float, default=0.10)
    ap.add_argument('--mappings', nargs='+', default=['RY', 'WS', 'MK'], choices=sorted(PG.MAPPINGS))
    ap.add_argument('--strands', choices=['both', 'forward'], default='both')
    ap.add_argument('--locus-gap', type=int, default=1000)
    ap.add_argument('--seed', type=int, default=23)
    ap.add_argument('--workers', type=int, default=12)
    ap.add_argument('--chroms', nargs='+')
    ap.add_argument('--exclude-chrom')
    ap.add_argument('--enumerate-limit', type=int, default=200_000)
    ap.add_argument('--no-run-rule', action='store_true',
                    help='allow decoys with runs longer than the real word\'s (composition is still kept)')
    ap.add_argument('--progress', action='store_true')
    a = ap.parse_args(argv)
    if not 1 <= a.pmax <= 20:
        ap.error('--pmax must be between 1 and 20')
    return a


def main(argv=None):
    a = parse_args(argv)
    t0 = time.time()
    tgts = TG.select_targets(a.targets, a.families, a.min_len, a.max_len)
    words, prof, info = build(tgts, a.n_decoys, a.modes, a.pmax, a.min_repeats, a.max_defect_rate,
                              a.seed, a.enumerate_limit, run_rule=not a.no_run_rule)
    stem = a.out[:-4] if a.out.endswith('.csv') else a.out
    prof.to_csv(stem + '_profile.csv', index=False)
    for r in prof.itertuples():
        print(f'{r.pattern:32s} L={r.length:3d} period {r.best_period:2d} unit {r.unit:<16s} '
              f'defects {r.hamming_defects:2d} ({r.defect_rate:.2f}){"" if r.near_periodic else "  not nearly periodic"}',
              file=sys.stderr)

    specs = [(int(r.wid), tuple(x for _, x in PG.search_forms(r.word.encode(), a.strands)), int(r.length))
             for r in words.itertuples()]
    chroms = CG.chrom_files(a.chrom_dir, a.exclude_chrom, a.chroms)
    if not chroms:
        raise SystemExit(f'no FASTA files in {a.chrom_dir}')
    tasks = [(c, m) for c in chroms for m in a.mappings]
    rows = []
    with Pool(a.workers, initializer=CG._init, initargs=(specs, a.locus_gap)) as pool:
        for i, part in enumerate(pool.imap_unordered(CG._work, tasks, chunksize=1), 1):
            rows.extend(part)
            if a.progress:
                print(f'  [{i}/{len(tasks)}] {part[0][1]} {part[0][2]}  {time.time() - t0:.0f}s',
                      file=sys.stderr, flush=True)
    per = pd.DataFrame(rows, columns=['wid', 'chrom', 'mapping', 'hits', 'loci', 'locus_bp', 'max_locus_bp'])
    g = (per.groupby(['wid', 'mapping'])
            .agg(hits=('hits', 'sum'), loci=('loci', 'sum'), locus_bp=('locus_bp', 'sum'),
                 max_locus_bp=('max_locus_bp', 'max'))
            .reset_index().merge(words, on='wid'))
    cols = ['wid', 'mapping', 'hits', 'loci', 'pattern', 'kind', 'word', 'mode', 'unit',
            'locus_bp', 'max_locus_bp', 'best_period', 'defects', 'family', 'variant', 'label', 'length']
    g = g[cols].sort_values(['wid', 'mapping']).reset_index(drop=True)
    g.to_csv(a.out, index=False)
    s = summarise(g, prof, info)
    s.to_csv(stem + '_summary.csv', index=False)
    run = dict(script='periodic_decoy.py', version=VERSION,
               date=datetime.datetime.now().isoformat(timespec='seconds'), args=vars(a),
               n_targets=len(tgts), n_words=len(words), chromosomes=[CG.chrom_name(c) for c in chroms],
               decoys={f'{t}|{m}': v for (t, m), v in info.items()}, runtime_s=round(time.time() - t0, 1),
               python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__)
    with open(stem + '_run.json', 'w') as fh:
        json.dump(run, fh, indent=2, default=str)
    with pd.option_context('display.width', 220, 'display.max_columns', 14):
        print(s[['pattern', 'decoy_mode', 'mapping', 'real_loci', 'tiling_loci', 'n_decoys',
                 'decoy_loci_median', 'decoy_loci_max', 'p_emp_loci', 'q_emp_loci', 'near_periodic']]
              .to_string(index=False))
    print(f'\nWrote {a.out}, {stem}_profile.csv, {stem}_summary.csv, {stem}_run.json '
          f'in {time.time() - t0:.0f}s', file=sys.stderr)


if __name__ == '__main__':
    main()
