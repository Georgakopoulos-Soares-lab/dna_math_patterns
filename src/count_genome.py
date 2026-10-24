"""Genome-wide decoy test: does a pattern word occur at more LOCI than words of
the same shape?

For every target word (targets.py; default: all 36 words of the 8 families with
13-100 bits) the script builds decoys (decoy_words.py; default: 100
run-permutation decoys -- same length, composition and run lengths, order
shuffled), counts every word in every chromosome file, and compares the real
word's number of loci with its decoys'.

Counting
--------
  forms   every distinct form of the word is searched (pattern_genome.search_forms):
          --strands both (default) = word, complement, reverse, reversed complement;
          --strands forward = word and complement only (the original)
  hits    all overlapping occurrences of all forms
  loci    hits whose start positions are within --locus-gap bp (default 1000) of
          the previous hit on the same record are merged into one locus (chained).
          Hits inside one tandem array are not independent, and counting them
          separately would inflate a word's frequency by the array's length.
  locus_bp, max_locus_bp
          total and largest span (first hit start to last hit end) covered by the
          loci -- shows when a word's count is driven by one huge array.
Each FASTA record is translated and clustered separately (records are never
joined), and no match crosses an N.

Statistics (summary file)
-------------------------
  frac_decoys_loci_ge_real   #{decoys with loci >= real} / N   (original column)
  decoy_hits_median, decoy_hits_max, frac_decoys_hits_ge_real, p_emp_hits
                             the same comparison on raw hits (the measure of the
                             original count_chromosome.py)
  p_emp_loci                 (1 + #{decoys >= real}) / (N + 1)
  q_emp_loci                 Benjamini-Hochberg q-value over all rows of the run
  p_floor                    smallest achievable p_emp, 1 / (N + 1)
  testable                   True if p_emp < 0.05 is achievable (N >= 19 decoys).
                             Short words have few rearrangements with the same
                             runs: e.g. Tribonacci isolate_2 T4 (13 bits) has only
                             2 run decoys, Fibonacci F5 9; 14 of the 36 default
                             targets have fewer than 100 (all are 13-29 bits).
  n_decoys, n_decoys_possible
                             decoys used, and how many distinct ones exist for this
                             word (exact when the method could list them all)
A word with 0 real loci has p_emp = 1 by construction (every decoy ties).

usage
-----
  python3 src/count_genome.py <chrom_dir | fasta> <n_decoys> <out.csv> [options]
     chrom_dir   folder with one FASTA per chromosome (see src/split_t2t.sh);
                 *.fa, *.fasta, *.fa.gz, *.fasta.gz
     fasta       or a single FASTA file, e.g. data/t2t_chroms/chr1.fa, to run the
                 test on one chromosome (this replaces count_chromosome.py)
  options (defaults in brackets)
     --targets all|readme|NAME ...  [all]        --min-len [13]  --max-len [100]
     --families ...                 [all 8]      --mappings [RY WS MK]
     --strands both|forward         [both]       --decoy-method run|segment|slip [run]
     --no-decoy-filters  turn off the composition/longest-run filters of segment/slip
     --locus-gap [1000]  --seed [11]  --workers [12]  --per-chrom  --progress
     --chroms NAME ...      only these chromosomes (file names without extension)
     --exclude-chrom REGEX  skip chromosomes whose name matches, e.g. chrM

  Approximate original run:
     ... --targets readme --strands forward --decoy-method segment --no-decoy-filters
  (not exact: the decoy duplicates and composition fixes in decoy_words.py change
  the decoy set).

Output
------
  <out.csv>              per word x mapping: hits, loci, locus_bp ... (original columns first)
  <out>_summary.csv      per target x mapping: real vs decoys, p-values
  <out>_per_chrom.csv    per word x mapping x chromosome (with --per-chrom)
  <out>_run.json         settings, decoy counts per target, versions

job() keeps its original signature and counting (word + complement, 1 kb
chaining) so that older scripts importing it still work.
"""
from __future__ import annotations

import argparse
import datetime
import glob
import json
import os
import platform
import re
import sys
import time
import warnings
from multiprocessing import Pool
from typing import Dict, List, Sequence, Tuple

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import decoy_words as DW      # noqa: E402
import pattern_genome as PG   # noqa: E402
import seqio                  # noqa: E402
import targets as TG          # noqa: E402

VERSION = '2.0'
TARGETS = TG.README_TARGETS          # kept for scripts that import it from here
FASTA_GLOBS = ('*.fa', '*.fasta', '*.fa.gz', '*.fasta.gz')


# ------------------------------------------------------------------ counting

def chrom_name(path: str) -> str:
    base = os.path.basename(path)
    for ext in ('.fa.gz', '.fasta.gz', '.fa', '.fasta'):
        if base.endswith(ext):
            return base[:-len(ext)]
    return base


def chrom_files(cdir: str, exclude: str = None, only: Sequence[str] = None) -> List[str]:
    """FASTA files to scan: `cdir` itself if it is a file, else the FASTA files in
    the folder; optionally only the chromosomes named in `only` and without those
    matching the regex `exclude`. Largest first."""
    if os.path.isfile(cdir):
        files = [cdir]
    else:
        files = sorted({f for g in FASTA_GLOBS for f in glob.glob(os.path.join(cdir, g))})
    if only:
        missing = set(only) - {chrom_name(f) for f in files}
        if missing:
            raise SystemExit(f'--chroms: not found in {cdir}: {", ".join(sorted(missing))}')
        files = [f for f in files if chrom_name(f) in set(only)]
    if exclude:
        rx = re.compile(exclude)
        files = [f for f in files if not rx.search(chrom_name(f))]
    return sorted(files, key=os.path.getsize, reverse=True)


def load_records(path: str, mapping: str) -> List[bytes]:
    """Projected binary of every record in the file (b'0'/b'1'/b'N'), not joined."""
    table = PG.TRANSLATE_TABLES[mapping]
    return [rec.seq.translate(table) for rec in seqio.read_fasta(path)]


def count_forms(binaries: Sequence[bytes], forms: Sequence[bytes], L: int,
                gap: int) -> Tuple[int, int, int, int]:
    """(hits, loci, locus_bp, max_locus_bp) of all forms over all records."""
    hits = loci = bp = mx = 0
    for b in binaries:
        pos = [PG.find_positions(b, f) for f in forms]
        p = np.sort(np.concatenate(pos)) if pos else np.empty(0, np.int64)
        if len(p) == 0:
            continue
        hits += len(p)
        new = np.flatnonzero(np.diff(p) > gap) + 1          # first hit of each later locus
        starts = np.concatenate(([0], new))
        ends = np.concatenate((new - 1, [len(p) - 1]))
        spans = p[ends] + L - p[starts]
        loci += len(starts)
        bp += int(spans.sum())
        mx = max(mx, int(spans.max()))
    return hits, loci, bp, mx


def _cluster(p: np.ndarray, L: int, gap: int) -> Tuple[int, int, int, int]:
    """(hits, loci, locus_bp, max_locus_bp) for sorted start positions p on one record."""
    new = np.flatnonzero(np.diff(p) > gap) + 1
    starts = np.concatenate(([0], new))
    ends = np.concatenate((new - 1, [len(p) - 1]))
    spans = p[ends] + L - p[starts]
    return len(p), len(starts), int(spans.sum()), int(spans.max())


class WordIndex:
    """All words to count, grouped by length, for a one-pass multi-word search.

    For every window of K = min(L, 64) bits the packed 64-bit value is looked up
    among the keys of all words of length L at once (np.searchsorted). Words
    longer than 64 bits are matched on their first 64 bits and then verified in
    full. Every overlapping occurrence of every form is found, exactly as with
    bytes.find (tests/test_count_genome.py checks the two against each other).
    """

    def __init__(self, specs: Sequence[Tuple[int, Sequence[bytes], int]]):
        by_len: Dict[int, Dict[int, List[Tuple[int, bytes]]]] = {}
        for wid, forms, L in specs:
            K = min(L, 64)
            for f in forms:
                by_len.setdefault(L, {}).setdefault(int(f[:K], 2), []).append((wid, f))
        self.lengths = sorted(by_len)
        self.max_len = max(self.lengths)
        self.keys = {L: np.array(sorted(by_len[L]), dtype=np.uint64) for L in self.lengths}
        self.owners = by_len
        self.length_of = {wid: L for wid, _, L in specs}

    def find(self, binary: bytes, chunk: int = 1 << 21) -> Dict[int, np.ndarray]:
        """{wid: sorted start positions} for one record (b'0'/b'1'/b'N' bytes)."""
        a = np.frombuffer(binary, np.uint8)
        n = len(a)
        found: Dict[int, List[np.ndarray]] = {}
        powers = [p for p in (1, 2, 4, 8, 16, 32, 64)]
        for s in range(0, n, chunk):
            seg = a[s:min(n, s + chunk + self.max_len - 1)]
            bits = (seg == 49).astype(np.uint64)
            bad = np.concatenate(([0], np.cumsum((seg != 48) & (seg != 49), dtype=np.int64)))
            w = {1: bits}
            for p in powers[1:]:
                h = p // 2
                prev = w[h]
                if len(prev) <= h:
                    break
                w[p] = (prev[:-h] << np.uint64(h)) | prev[h:]
            for L in self.lengths:
                m = min(chunk, n - s - L + 1)
                if m <= 0:
                    continue
                K = min(L, 64)
                acc, off = None, 0
                for p in reversed(powers):
                    if K & p:
                        piece = w[p][off:off + m]
                        acc = piece if acc is None else (acc << np.uint64(p)) | piece
                        off += p
                keys = self.keys[L]
                idx = np.searchsorted(keys, acc)
                idx[idx >= len(keys)] = len(keys) - 1
                hit = (keys[idx] == acc) & ((bad[L:L + m] - bad[:m]) == 0)
                cand = np.flatnonzero(hit)
                if not len(cand):
                    continue
                owners = self.owners[L]
                for key in np.unique(acc[cand]).tolist():
                    pos = cand[acc[cand] == key] + s
                    for wid, f in owners[key]:
                        if L > 64:
                            pos_ok = np.array([q for q in pos.tolist() if binary[q:q + L] == f], dtype=np.int64)
                        else:
                            pos_ok = pos
                        if len(pos_ok):
                            found.setdefault(wid, []).append(pos_ok)
        return {wid: np.sort(np.concatenate(v)) for wid, v in found.items()}

    def count(self, binaries: Sequence[bytes], gap: int) -> Dict[int, Tuple[int, int, int, int]]:
        """{wid: (hits, loci, locus_bp, max_locus_bp)} summed over records."""
        tot: Dict[int, List[int]] = {}
        for b in binaries:
            for wid, p in self.find(b).items():
                h, l, bp, mx = _cluster(p, self.length_of[wid], gap)
                t = tot.setdefault(wid, [0, 0, 0, 0])
                t[0] += h; t[1] += l; t[2] += bp; t[3] = max(t[3], mx)
        return {wid: tuple(v) for wid, v in tot.items()}


_INDEX: WordIndex = None
_WIDS: List[int] = []
_GAP = 1000


def _init(specs, gap):
    global _INDEX, _WIDS, _GAP
    _INDEX, _WIDS, _GAP = WordIndex(specs), [wid for wid, _, _ in specs], gap


def _work(task):
    path, mapping = task
    recs = load_records(path, mapping)
    c = chrom_name(path)
    res = _INDEX.count(recs, _GAP)
    return [(wid, c, mapping) + res.get(wid, (0, 0, 0, 0)) for wid in _WIDS]


def job(j):
    """Original interface: j = (path, mapping, [(wid, word_str), ...]) ->
    [(wid, chrom, mapping, hits, loci), ...], counting word + complement with
    1 kb chaining (records are no longer joined)."""
    path, m, words = j
    recs = load_records(path, m)
    res = []
    for wid, w in words:
        wb = w.encode()
        h, l, _, _ = count_forms(recs, (wb, PG.bit_complement(wb)), len(wb), 1000)
        res.append((wid, chrom_name(path), m, h, l))
    return res


# ----------------------------------------------------------------- building

def build_words(tgts: Dict[str, PG.PatternInfo], n_decoys: int, method: str, filters: bool,
                seed: int) -> Tuple[pd.DataFrame, Dict[str, dict]]:
    """Rows (wid, pattern, kind, word, family, variant, label, length) and decoy info per target."""
    rows, info = [], {}
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        for tag, p in tgts.items():
            w = p.word.decode()
            ds = DW.make_decoys(w, n_decoys, TG.target_rng(seed, tag), method=method,
                                match_composition=filters, match_runs='max' if filters else None)
            info[tag] = dict(n_decoys=ds.n_returned, n_decoys_possible=ds.n_possible,
                             exhaustive=ds.exhaustive, p_floor=ds.p_floor)
            base = dict(family=p.family, variant=p.variant, label=p.label, length=p.length)
            rows.append(dict(wid=len(rows), pattern=tag, kind='REAL', word=w, **base))
            for d in ds.words:
                rows.append(dict(wid=len(rows), pattern=tag, kind='DECOY', word=d, **base))
    for c in caught:
        print('warning:', c.message, file=sys.stderr)
    return pd.DataFrame(rows), info


def summarise(g: pd.DataFrame, info: Dict[str, dict], method: str) -> pd.DataFrame:
    rows = []
    for (tag, m), x in g.groupby(['pattern', 'mapping'], sort=False):
        r = x[x.kind == 'REAL'].iloc[0]
        d = x[x.kind == 'DECOY']
        n = len(d)
        rows.append(dict(
            pattern=tag, mapping=m, real_hits=int(r.hits), real_loci=int(r.loci), n_decoys=n,
            decoy_loci_median=d.loci.median() if n else np.nan,
            decoy_loci_max=d.loci.max() if n else np.nan,
            frac_decoys_loci_ge_real=round(float((d.loci >= r.loci).mean()), 3) if n else np.nan,
            p_emp_loci=DW.empirical_p(r.loci, d.loci) if n else np.nan,
            family=r.family, variant=r.variant, label=r.label, length=int(r.length),
            decoy_method=method, n_decoys_possible=info[tag]['n_decoys_possible'],
            p_floor=1.0 / (n + 1), testable=bool(n >= 19),
            real_locus_bp=int(r.locus_bp), decoy_locus_bp_median=d.locus_bp.median() if n else np.nan,
            real_max_locus_bp=int(r.max_locus_bp),
            decoy_hits_median=d.hits.median() if n else np.nan,
            decoy_hits_max=d.hits.max() if n else np.nan,
            frac_decoys_hits_ge_real=round(float((d.hits >= r.hits).mean()), 3) if n else np.nan,
            p_emp_hits=DW.empirical_p(r.hits, d.hits) if n else np.nan,
        ))
    s = pd.DataFrame(rows)
    s.insert(s.columns.get_loc('p_emp_loci') + 1, 'q_emp_loci', PG.bh_qvalues(s.p_emp_loci.tolist()))
    return s


# ---------------------------------------------------------------------- CLI

def parse_args(argv=None):
    ap = argparse.ArgumentParser(description='Genome-wide real-vs-decoy locus counts for pattern words.')
    ap.add_argument('chrom_dir', help='folder with one FASTA per chromosome (src/split_t2t.sh), '
                                      'or a single FASTA file')
    ap.add_argument('n_decoys', type=int, help='decoys per target (e.g. 100)')
    ap.add_argument('out', help='output CSV')
    ap.add_argument('--targets', nargs='+', default=['all'],
                    help="all (default), readme (the original 5), README tags or family/variant/label")
    ap.add_argument('--families', nargs='+', default=list(PG.DEFAULT_FAMILIES), choices=PG.DEFAULT_FAMILIES)
    ap.add_argument('--min-len', type=int, default=13)
    ap.add_argument('--max-len', type=int, default=100)
    ap.add_argument('--mappings', nargs='+', default=['RY', 'WS', 'MK'], choices=sorted(PG.MAPPINGS))
    ap.add_argument('--strands', choices=['both', 'forward'], default='both')
    ap.add_argument('--decoy-method', choices=list(DW.METHODS), default='run')
    ap.add_argument('--no-decoy-filters', action='store_true',
                    help='segment/slip decoys without the composition and longest-run filters')
    ap.add_argument('--locus-gap', type=int, default=1000)
    ap.add_argument('--seed', type=int, default=11)
    ap.add_argument('--workers', type=int, default=12)
    ap.add_argument('--chroms', nargs='+', help='only these chromosomes (file names without extension)')
    ap.add_argument('--exclude-chrom', help='regex of chromosome file names to skip, e.g. chrM')
    ap.add_argument('--per-chrom', action='store_true', help='also write <out>_per_chrom.csv')
    ap.add_argument('--progress', action='store_true')
    return ap.parse_args(argv)


def main(argv=None):
    a = parse_args(argv)
    t0 = time.time()
    tgts = TG.select_targets(a.targets, a.families, a.min_len, a.max_len)
    if not tgts:
        raise SystemExit('no targets selected')
    words, info = build_words(tgts, a.n_decoys, a.decoy_method, not a.no_decoy_filters, a.seed)
    specs = [(int(r.wid), tuple(w for _, w in PG.search_forms(r.word.encode(), a.strands)), int(r.length))
             for r in words.itertuples()]
    chroms = chrom_files(a.chrom_dir, a.exclude_chrom, a.chroms)
    if not chroms:
        raise SystemExit(f'no FASTA files in {a.chrom_dir}')
    tasks = [(c, m) for c in chroms for m in a.mappings]
    print(f'{len(tgts)} targets, {len(words)} words (real + decoys), {len(chroms)} chromosome files, '
          f'{len(tasks)} jobs on {a.workers} workers', file=sys.stderr, flush=True)

    rows = []
    with Pool(a.workers, initializer=_init, initargs=(specs, a.locus_gap)) as pool:
        for i, part in enumerate(pool.imap_unordered(_work, tasks, chunksize=1), 1):
            rows.extend(part)
            if a.progress:
                print(f'  [{i}/{len(tasks)}] {part[0][1]} {part[0][2]}  {time.time() - t0:.0f}s',
                      file=sys.stderr, flush=True)

    per = pd.DataFrame(rows, columns=['wid', 'chrom', 'mapping', 'hits', 'loci', 'locus_bp', 'max_locus_bp'])
    g = (per.groupby(['wid', 'mapping'])
            .agg(hits=('hits', 'sum'), loci=('loci', 'sum'), locus_bp=('locus_bp', 'sum'),
                 max_locus_bp=('max_locus_bp', 'max'))
            .reset_index().merge(words, on='wid'))
    cols = ['wid', 'mapping', 'hits', 'loci', 'pattern', 'kind', 'word',
            'locus_bp', 'max_locus_bp', 'family', 'variant', 'label', 'length']
    g = g[cols].sort_values(['wid', 'mapping']).reset_index(drop=True)
    g.to_csv(a.out, index=False)
    stem = a.out[:-4] if a.out.endswith('.csv') else a.out
    if a.per_chrom:
        per.merge(words[['wid', 'pattern', 'kind']], on='wid').to_csv(stem + '_per_chrom.csv', index=False)

    s = summarise(g, info, a.decoy_method)
    s.to_csv(stem + '_summary.csv', index=False)
    run = dict(script='count_genome.py', version=VERSION,
               date=datetime.datetime.now().isoformat(timespec='seconds'),
               args=vars(a), n_targets=len(tgts), n_words=len(words), chromosomes=[chrom_name(c) for c in chroms],
               decoys=info, runtime_s=round(time.time() - t0, 1),
               python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__)
    with open(stem + '_run.json', 'w') as fh:
        json.dump(run, fh, indent=2, default=str)
    with pd.option_context('display.width', 200, 'display.max_columns', 12):
        print(s[['pattern', 'mapping', 'real_loci', 'n_decoys', 'decoy_loci_median', 'decoy_loci_max',
                 'p_emp_loci', 'q_emp_loci']].to_string(index=False))
    print(f'\nWrote {a.out}, {stem}_summary.csv, {stem}_run.json in {time.time() - t0:.0f}s', file=sys.stderr)


if __name__ == '__main__':
    main()
