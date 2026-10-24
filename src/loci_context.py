"""What DNA lies under the enriched pattern hits? Tandem period of every locus.

Takes the strongly enriched rows of a pattern_genome.py results CSV, finds the
words in the same FASTA that was scanned, and reports for every locus the
DNA's best tandem period and how well that period explains it. This is the
evidence behind README section 6.2 ("every enriched hit lies in a tandem
repeat"), now for all loci instead of a few examples.

usage
-----
  python3 src/loci_context.py <fasta> <names.tsv | -> <scan.csv> [min_OE] [options]
     fasta      the FASTA that was scanned by pattern_genome.py (plain or gzip)
     names.tsv  optional, only to display chromosome names: records.tsv or
                headers.tsv from split_t2t.sh; '-' to show record IDs
     scan.csv   pattern_genome.py output
     min_OE     rows with O/E (vs the Markov model of the scan) above this,
                default 50; rows with O/E = inf (hits but expected 0) are
                always included and listed first
  options (defaults in brackets)
     --max-q Q             also require q_enrich_clumps <= Q
     --from-examples       analyse only the example_loci stored in the CSV (the
                           original behaviour) instead of searching all hits
     --strands both|forward [both]   forms searched (match the scan's --strands)
     --locus-gap [1000]    hits closer than this are one locus (as count_genome)
     --max-period [20]     tandem periods tried, 1..max-period bp
     --flank [30]          bp of context on each side
     --tandem-threshold [0.9]  self-match needed to call a locus a tandem repeat
     --loci-per-row [6]    examples printed per pattern (drawn at random, --seed [1])
     --out PREFIX          write PREFIX_loci.csv (one row per locus) and
                           PREFIX_summary.csv (one row per pattern x mapping)

What is measured
----------------
For each locus (hits chained within --locus-gap on the same record):
  hit_period, hit_self_match       on the first hit's DNA alone (the original
                                   measure: the fraction of positions i with
                                   x[i] == x[i+p], best p, ties to the smaller p)
  context_period, context_self_match
                                   on the whole locus plus --flank bp each side:
                                   is the hit embedded in a longer repeat?
  verified                         the DNA at the reported position, projected in
                                   the row's encoding, equals the reported form
                                   (word, complement, reverse, reverse_complement)
Summary per pattern x mapping: hits, loci, loci with hit_self_match >= threshold
(n_tandem_hit) and with context_self_match >= threshold (n_tandem_context), the
fractions, median self-matches and the most common period.

DNA is fetched by coordinates from the FASTA itself, one record at a time; no
assumption about line width, no header-line arithmetic, no temporary files.
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import sys
from collections import Counter, defaultdict
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import count_genome as CG     # noqa: E402
import pattern_genome as PG   # noqa: E402
import projections as P       # noqa: E402
import seqio                  # noqa: E402

csv.field_size_limit(sys.maxsize)


# ------------------------------------------------------------------- input

def read_names(path: Optional[str]) -> Tuple[Dict[str, str], List[str]]:
    """({record_id: chrom}, [chrom in file order]) from records.tsv or headers.tsv."""
    if not path or path == '-':
        return {}, []
    rows = [l.rstrip('\n').split('\t') for l in open(path) if l.strip()]
    if rows and rows[0][:2] == ['chrom', 'record_id']:          # records.tsv
        return {r[1]: r[0] for r in rows[1:]}, [r[0] for r in rows[1:]]
    return {}, [r[1] for r in rows]                               # headers.tsv: line, chrom


def oe_column(fields: Sequence[str]) -> str:
    cols = [c for c in fields if c.startswith('obs_over_markov')]
    if not cols:
        raise SystemExit('scan CSV has no obs_over_markov* column')
    return cols[0]


def select_rows(scan: str, min_oe: float, max_q: Optional[float]) -> List[dict]:
    rows = list(csv.DictReader(open(scan)))
    if not rows:
        return []
    oe = oe_column(rows[0].keys())
    qcol = next((c for c in rows[0] if c.startswith('q_enrich_clumps')), None)
    out = []
    for r in rows:
        if int(r['observed_total']) <= 0 or r[oe] in ('', 'nan'):
            continue
        v = float(r[oe])
        if not (v > min_oe or math.isinf(v)):
            continue
        if max_q is not None:
            if qcol is None:
                raise SystemExit('--max-q needs a scan CSV with q_enrich_clumps (pattern_genome.py 2.0)')
            if float(r[qcol]) > max_q:
                continue
        r['_oe'] = v
        out.append(r)
    out.sort(key=lambda r: -r['_oe'])
    return out


# ---------------------------------------------------------------- analysis

def best_period(x: np.ndarray, maxp: int) -> Tuple[int, float]:
    """(period, self-match) maximising the fraction of x[i] == x[i+p] over valid
    positions; ties go to the smaller period."""
    best = (1, 0.0)
    for p in range(1, maxp + 1):
        if len(x) <= p:
            break
        a, b = x[p:], x[:-p]
        ok = (a < 4) & (b < 4)
        if not ok.any():
            continue
        f = float((a[ok] == b[ok]).mean())
        if f > best[1] + 1e-9:
            best = (p, f)
    return best


def dna_string(seq: bytes, s: int, e: int, hs: int, he: int) -> str:
    """seq[s:e] with the hit seq[hs:he] in capitals and the flanks in lower case."""
    return (seq[s:hs].decode().lower() + seq[hs:he].decode().upper() + seq[he:e].decode().lower())


def analyse_locus(seq: bytes, codes: np.ndarray, hits: List[Tuple[int, str]], L: int, mapping: str,
                  forms: Dict[str, bytes], a) -> dict:
    p0, lab = hits[0]
    last = hits[-1][0]
    n = len(seq)
    hp, hf = best_period(codes[p0:p0 + L], a.max_period)
    cs, ce = max(0, p0 - a.flank), min(n, last + L + a.flank)
    cp, cf = best_period(codes[cs:ce], a.max_period)
    projected = P.project_bytes(seq[p0:p0 + L], mapping)
    return dict(start=p0 + 1, end=last + L, n_hits=len(hits), span_bp=last + L - p0, form=lab,
                verified=bool(forms.get(lab) == projected),
                hit_period=hp, hit_self_match=round(hf, 3), context_period=cp, context_self_match=round(cf, 3),
                dna=dna_string(seq, max(0, p0 - a.flank), min(n, p0 + L + a.flank), p0, p0 + L))


def cluster(hits: List[Tuple[int, str]], gap: int) -> List[List[Tuple[int, str]]]:
    hits = sorted(hits)
    loci, cur = [], []
    for h in hits:
        if cur and h[0] - cur[-1][0] > gap:
            loci.append(cur); cur = []
        cur.append(h)
    if cur:
        loci.append(cur)
    return loci


def scan_all(fasta: str, rows: List[dict], a) -> Dict[int, List[dict]]:
    """Search every selected word (all forms) in the FASTA; analyse every locus."""
    words = []                                     # (row index, mapping, L, {label: form})
    for i, r in enumerate(rows):
        forms = dict(PG.search_forms(r['word'].encode(), a.strands))
        words.append((i, r['mapping'], int(r['length']), forms))
    by_map = defaultdict(list)
    for w in words:
        by_map[w[1]].append(w)
    indexes, fid_info = {}, {}
    for m, ws in by_map.items():
        specs = []
        for i, _, L, forms in ws:
            for lab, f in forms.items():
                fid = len(fid_info)
                fid_info[fid] = (i, lab)
                specs.append((fid, (f,), L))
        indexes[m] = CG.WordIndex(specs)
    out: Dict[int, List[dict]] = defaultdict(list)
    for rec_i, rec in enumerate(seqio.read_fasta(fasta)):
        codes = seqio.encode(rec.seq)
        for m, idx in indexes.items():
            hits: Dict[int, List[Tuple[int, str]]] = defaultdict(list)
            for fid, pos in idx.find(P.project_bytes(rec.seq, m)).items():
                i, lab = fid_info[fid]
                hits[i].extend((int(x), lab) for x in pos)
            for i, hs in hits.items():
                _, _, L, forms = words[i]
                for locus in cluster(hs, a.locus_gap):
                    d = analyse_locus(rec.seq, codes, locus, L, m, forms, a)
                    d.update(record=rec.name, record_index=rec_i)
                    out[i].append(d)
        if a.progress:
            print(f'  {rec.name}: {sum(len(v) for v in out.values())} loci so far', file=sys.stderr, flush=True)
    return out


def from_examples(fasta: str, rows: List[dict], a) -> Dict[int, List[dict]]:
    """Original behaviour: the example_loci stored in the CSV, one locus each."""
    want = defaultdict(list)                       # record -> [(row, pos0, label)]
    for i, r in enumerate(rows):
        for loc in (r.get('example_loci') or '').split(';'):
            if loc:
                c, p, o = loc.rsplit(':', 2)
                want[c].append((i, int(p) - 1, o))
    out: Dict[int, List[dict]] = defaultdict(list)
    for rec_i, rec in enumerate(seqio.read_fasta(fasta)):
        if rec.name not in want:
            continue
        codes = seqio.encode(rec.seq)
        for i, p0, lab in want[rec.name]:
            r = rows[i]
            forms = dict(PG.search_forms(r['word'].encode(), 'both'))
            d = analyse_locus(rec.seq, codes, [(p0, lab)], int(r['length']), r['mapping'], forms, a)
            d.update(record=rec.name, record_index=rec_i)
            out[i].append(d)
    return out


def summarise(rows: List[dict], loci: Dict[int, List[dict]], thr: float) -> List[dict]:
    out = []
    for i, r in enumerate(rows):
        L = loci.get(i, [])
        n = len(L)
        th = sum(d['hit_self_match'] >= thr for d in L)
        tc = sum(d['context_self_match'] >= thr for d in L)
        per = Counter(d['context_period'] for d in L)
        out.append(dict(
            mapping=r['mapping'], pattern=f"{r['family']}/{r['variant']} {r['label']}", length=r['length'],
            obs_over_expected=r['_oe'], hits=sum(d['n_hits'] for d in L), loci=n,
            n_tandem_hit=th, frac_tandem_hit=round(th / n, 3) if n else float('nan'),
            n_tandem_context=tc, frac_tandem_context=round(tc / n, 3) if n else float('nan'),
            median_hit_self_match=float(np.median([d['hit_self_match'] for d in L])) if n else float('nan'),
            median_context_self_match=float(np.median([d['context_self_match'] for d in L])) if n else float('nan'),
            most_common_period=per.most_common(1)[0][0] if per else '',
            n_unverified=sum(not d['verified'] for d in L)))
    return out


# --------------------------------------------------------------------- CLI

def parse_args(argv=None):
    ap = argparse.ArgumentParser(description='DNA context and tandem period of enriched pattern loci.')
    ap.add_argument('fasta')
    ap.add_argument('names', help="records.tsv or headers.tsv from split_t2t.sh, or '-'")
    ap.add_argument('scan', help='pattern_genome.py results CSV')
    ap.add_argument('min_oe', nargs='?', type=float, default=50.0)
    ap.add_argument('--max-q', type=float)
    ap.add_argument('--from-examples', action='store_true')
    ap.add_argument('--strands', choices=['both', 'forward'], default='both')
    ap.add_argument('--locus-gap', type=int, default=1000)
    ap.add_argument('--max-period', type=int, default=20)
    ap.add_argument('--flank', type=int, default=30)
    ap.add_argument('--tandem-threshold', type=float, default=0.9)
    ap.add_argument('--loci-per-row', type=int, default=6)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--out')
    ap.add_argument('--progress', action='store_true')
    return ap.parse_args(argv)


def main(argv=None):
    a = parse_args(argv)
    id2chrom, order = read_names(a.names)
    rows = select_rows(a.scan, a.min_oe, a.max_q)
    if not rows:
        raise SystemExit(f'no rows with O/E > {a.min_oe} (or inf) and at least one hit')
    loci = from_examples(a.fasta, rows, a) if a.from_examples else scan_all(a.fasta, rows, a)
    for i, ds in loci.items():
        for d in ds:
            d['chrom'] = id2chrom.get(d['record']) or (order[d['record_index']] if d['record_index'] < len(order)
                                                       else d['record'])
    rng = np.random.default_rng(a.seed)
    print(f"{'proj':4s} {'pattern':28s} {'O/E':>9s} {'locus':22s} {'hits':>5s} {'period':>6s} {'self':>5s} "
          f"{'ctx p':>5s} {'ctx':>5s}  DNA (hit in CAPITALS, {a.flank} bp flanks)")
    for i, r in enumerate(rows):
        ds = loci.get(i, [])
        pick = sorted(rng.choice(len(ds), size=min(a.loci_per_row, len(ds)), replace=False).tolist()) if ds else []
        name = f"{r['family']}/{r['variant']} {r['label']}"
        for j in pick:
            d = ds[j]
            flag = '' if d['verified'] else '  [NOT VERIFIED]'
            print(f"{r['mapping']:4s} {name:28s} {r['_oe']:9.3g} {d['chrom'] + ':' + str(d['start']):22s} "
                  f"{d['n_hits']:5d} {d['hit_period']:6d} {d['hit_self_match']:5.2f} {d['context_period']:5d} "
                  f"{d['context_self_match']:5.2f}  {d['dna']}{flag}")
    summ = summarise(rows, loci, a.tandem_threshold)
    print('\nSummary per pattern (loci = hits within %d bp merged; tandem = self-match >= %.2f):'
          % (a.locus_gap, a.tandem_threshold))
    for s in summ:
        print(f"  {s['mapping']:3s} {s['pattern']:30s} loci {s['loci']:7d}  tandem (hit) {s['n_tandem_hit']:7d}"
              f" ({s['frac_tandem_hit']:.2f})  tandem (context) {s['n_tandem_context']:7d} ({s['frac_tandem_context']:.2f})"
              f"  period {s['most_common_period']}")
    if a.out:
        cols = ['mapping', 'pattern', 'label', 'length', 'obs_over_expected', 'chrom', 'record', 'start', 'end',
                'n_hits', 'span_bp', 'form', 'verified', 'hit_period', 'hit_self_match',
                'context_period', 'context_self_match', 'dna']
        with open(a.out + '_loci.csv', 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=cols, extrasaction='ignore')
            w.writeheader()
            for i, r in enumerate(rows):
                for d in loci.get(i, []):
                    w.writerow(dict(d, mapping=r['mapping'], pattern=f"{r['family']}/{r['variant']}",
                                    label=r['label'], length=r['length'], obs_over_expected=r['_oe']))
        with open(a.out + '_summary.csv', 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(summ[0]))
            w.writeheader(); w.writerows(summ)
        print(f'\nWrote {a.out}_loci.csv and {a.out}_summary.csv', file=sys.stderr)
    bad = sum(not d['verified'] for ds in loci.values() for d in ds)
    if bad:
        print(f'WARNING: {bad} loci whose DNA does not project to the reported form', file=sys.stderr)


if __name__ == '__main__':
    main()
