#!/usr/bin/env python3
"""Choose 100 reference genomes per group (bacteria, archaea, viruses) from NCBI RefSeq.

Rule (fixed seed, so the selection is reproducible and can be reported):
  * candidates: RefSeq reference genomes with assembly level "Complete Genome"
    (NCBI Datasets API, taxon 2 = Bacteria, 2157 = Archaea, 10239 = Viruses);
  * bacteria and archaea: one genome per genus (genus = first word of the organism
    name, skipping "Candidatus" and brackets), so no genus is counted twice; the
    genomes already analysed (the 10 + 10 of the first run) are always included,
    and the remaining places are filled by genera drawn at random;
  * viruses: drawn at random within four genome-length quartiles (25 each), so that
    large DNA viruses are represented as well as the many small genomes.

Writes data/sets/<group>100.tsv (name, accession, organism, length_bp, gc_percent,
assembly_level) and data/sets/<group>100_candidates.tsv (every candidate seen).

usage:  python3 select_genomes.py                 # all three groups
        python3 select_genomes.py bacteria        # one group
        N=100 SEED=1 OUT=data/sets python3 select_genomes.py
needs:  Python 3 only (standard library); internet access to api.ncbi.nlm.nih.gov
"""
import csv
import json
import os
import random
import sys
import time
import urllib.parse
import urllib.request

API = 'https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon'
ROOT = os.environ.get('ROOT', '/Users/kleiopapatheodorou/Documents/lab/fibonacci/patterns')
OUT = os.environ.get('OUT', os.path.join(ROOT, 'data', 'sets'))
N = int(os.environ.get('N', 100))
SEED = int(os.environ.get('SEED', 1))
TAXA = {'bacteria': 2, 'archaea': 2157, 'viruses': 10239}

# the genomes of the first run, always kept so the two analyses overlap
ALREADY = {
    'bacteria': ['GCF_000005845.2', 'GCF_000009045.1', 'GCF_000013425.1', 'GCF_000195955.2',
                 'GCF_000203835.1', 'GCF_000006765.1', 'GCF_000006745.1', 'GCF_000008525.1',
                 'GCF_000027325.1', 'GCF_000009725.1'],
    'archaea': ['GCF_000025685.1', 'GCF_000006805.1', 'GCF_000012285.1', 'GCF_000007005.1',
                'GCF_000009965.1', 'GCF_000007305.1', 'GCF_000011585.1', 'GCF_000007345.1',
                'GCF_000091665.1', 'GCF_000018465.1'],
    'viruses': [],
}


def get(url, tries=4):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return json.load(r)
        except Exception as e:  # network hiccup or rate limit: wait and retry
            if i == tries - 1:
                raise
            print(f'   retry after error: {e}', file=sys.stderr)
            time.sleep(3 * (i + 1))


def candidates(taxon, reference_only=True):
    params = {'filters.assembly_source': 'refseq', 'filters.assembly_level': 'complete_genome',
              'page_size': '1000'}
    if reference_only:
        params['filters.reference_only'] = 'true'
    out, token = [], None
    while True:
        p = dict(params, **({'page_token': token} if token else {}))
        d = get(f'{API}/{taxon}/dataset_report?' + urllib.parse.urlencode(p))
        for r in d.get('reports', []):
            st = r.get('assembly_stats', {})
            out.append({'accession': r['accession'],
                        'organism': r.get('organism', {}).get('organism_name', ''),
                        'length_bp': int(st.get('total_sequence_length') or 0),
                        'gc_percent': st.get('gc_percent', ''),
                        'assembly_level': r.get('assembly_info', {}).get('assembly_level', '')})
        token = d.get('next_page_token')
        print(f'   {len(out)} candidates so far', file=sys.stderr)
        if not token:
            return out
        time.sleep(0.4)


def by_accession(accs):
    """Reports for specific accessions (used for first-run genomes not flagged as reference)."""
    if not accs:
        return []
    d = get('https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/'
            + urllib.parse.quote(','.join(accs)) + '/dataset_report')
    out = []
    for r in d.get('reports', []):
        st = r.get('assembly_stats', {})
        out.append({'accession': r['accession'],
                    'organism': r.get('organism', {}).get('organism_name', ''),
                    'length_bp': int(st.get('total_sequence_length') or 0),
                    'gc_percent': st.get('gc_percent', ''),
                    'assembly_level': r.get('assembly_info', {}).get('assembly_level', '')})
    return out


def genus(name):
    w = name.replace('[', '').replace(']', '').split()
    if w and w[0] == 'Candidatus':
        w = w[1:]
    return w[0] if w else name


def choose(group, cands, rng):
    if group == 'viruses':
        cands = sorted(cands, key=lambda c: (c['length_bp'], c['accession']))
        q = [cands[i * len(cands) // 4:(i + 1) * len(cands) // 4] for i in range(4)]
        per = N // 4
        pick = []
        for i, part in enumerate(q):
            k = per + (N - 4 * per if i == 3 else 0)
            pick += rng.sample(part, min(k, len(part)))
        return pick
    by_genus = {}
    for c in sorted(cands, key=lambda c: c['accession']):
        by_genus.setdefault(genus(c['organism']), []).append(c)
    acc = {c['accession']: c for c in cands}
    keep = [acc[a] for a in ALREADY[group] if a in acc]
    missing = [a for a in ALREADY[group] if a not in acc]
    if missing:
        print(f'   note: {len(missing)} genomes of the first run are not in the candidate list '
              f'(not "reference" or not complete): {" ".join(missing)}', file=sys.stderr)
    used = {genus(c['organism']) for c in keep}
    rest = sorted(g for g in by_genus if g not in used)
    for g in rng.sample(rest, max(0, min(N - len(keep), len(rest)))):
        keep.append(rng.choice(by_genus[g]))
    if len(keep) < N:
        print(f'   note: only {len(keep)} genera available for {group}', file=sys.stderr)
    return keep


def write(path, rows):
    cols = ['name', 'accession', 'organism', 'length_bp', 'gc_percent', 'assembly_level']
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, cols, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        for r in rows:
            w.writerow(dict(r, name=r['accession']))


def main():
    groups = sys.argv[1:] or list(TAXA)
    os.makedirs(OUT, exist_ok=True)
    for group in groups:
        print(f'{group}: reading candidates from NCBI', file=sys.stderr)
        cands = candidates(TAXA[group])
        if len(cands) < N:   # e.g. viral assemblies may not carry the "reference" flag
            print(f'   only {len(cands)} reference genomes; using all complete RefSeq assemblies',
                  file=sys.stderr)
            cands = candidates(TAXA[group], reference_only=False)
        have = {c['accession'] for c in cands}
        extra = by_accession([a for a in ALREADY[group] if a not in have])
        if extra:
            print(f'   added {len(extra)} first-run genomes by accession', file=sys.stderr)
            cands += extra
        write(os.path.join(OUT, f'{group}{N}_candidates.tsv'), cands)
        pick = sorted(choose(group, cands, random.Random(SEED)), key=lambda c: c['organism'])
        write(os.path.join(OUT, f'{group}{N}.tsv'), pick)
        mb = sum(c['length_bp'] for c in pick) / 1e6
        print(f'{group}: {len(pick)} genomes chosen from {len(cands)} candidates, {mb:.1f} Mb in total '
              f'-> {os.path.join(OUT, f"{group}{N}.tsv")}')


if __name__ == '__main__':
    main()
