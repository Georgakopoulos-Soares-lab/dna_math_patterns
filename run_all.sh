#!/bin/bash
# Run the complete analysis (code version 2.0).
#
# usage:  bash run_all.sh                    # every step, in order
#         bash run_all.sh decoys periodic    # only the named steps
#         bash run_all.sh --list             # list the steps
#
# Settings below can also be given as environment variables, e.g.
#         WORKERS=24 REPLICATES=20 bash run_all.sh longest
#         CORPORA="hs_chr20:data/hs_chr20.fa" bash run_all.sh fibmod4
#
# Each step logs to $RESULTS/logs/<step>.log (per corpus for longest/fibmod4)
# and stops the script on error. Run the steps in order the first time: context
# needs scan, and decoys/periodic need split. Times are estimates (12 cores).
set -euo pipefail

# ------------------------------------------------------------------ settings
ROOT=/Users/kleiopapatheodorou/Documents/lab/fibonacci/patterns
HERE=$ROOT/src                            # the .py files and split_t2t.sh
DATA=$ROOT/data
T2T=$DATA/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna
CHROMS=$DATA/t2t_chroms                   # per-chromosome T2T files, written by step split
RESULTS=$ROOT/results_new
WORKERS=${WORKERS:-12}
N_DECOYS=${N_DECOYS:-100}               # decoys per target word
REPLICATES=${REPLICATES:-100}           # null replicates (longest, fibmod4)
LM_BASES=${LM_BASES:-8000000}           # longest: ACGT bases per corpus (same for all corpora)
# corpora for longest and fibmod4: name:path, space separated
CORPORA=${CORPORA:-"hs_chr20:$DATA/hs_chr20.fa hs_chr21:$DATA/hs_chr21.fa mm_chr19:$DATA/mm_chr19.fa wheat_1A:$DATA/wheat_1A.fa t2t_chr1_cen:$DATA/t2t_chr1_cen.fa"}
PY=/opt/anaconda3/bin/python
SRC=$HERE                                 # folder with the .py files and split_t2t.sh

STEPS=(tests split scan context decoys periodic longest fibmod4)

describe() {
  case $1 in
    tests)    echo "unit tests of all scripts (~5 min)";;
    split)    echo "T2T -> one FASTA per chromosome + headers.tsv + records.tsv (~5-10 min)";;
    scan)     echo "pattern_genome: 8 families, 13-10000 bits, RY WS MK, both strands, whole T2T (hours)";;
    context)  echo "loci_context: DNA and tandem period of every locus of the enriched rows (needs scan)";;
    decoys)   echo "count_genome: 36 words x $N_DECOYS run decoys, genome-wide loci (needs split)";;
    periodic) echo "periodic_decoy: nearest period + periodic decoys, strict and --no-run-rule (needs split)";;
    longest)  echo "longest_match: 26 generators vs $LM_BASES bases of each corpus, 3 nulls x $REPLICATES (~15 min/corpus)";;
    fibmod4)  echo "fibmod4: mod-4 recurrence on each corpus, 3 nulls x $REPLICATES (~20-40 min/corpus, ~3 GB/worker)";;
  esac
}

log_run() {                               # log_run <log name> <command...>
  local name=$1; shift
  mkdir -p "$RESULTS/logs"
  local t0=$SECONDS
  if "$@" > "$RESULTS/logs/$name.log" 2>&1; then
    echo "   $name: done in $((SECONDS - t0)) s (log: $RESULTS/logs/$name.log)"
  else
    echo "   $name: FAILED after $((SECONDS - t0)) s -- see $RESULTS/logs/$name.log" >&2
    tail -n 20 "$RESULTS/logs/$name.log" >&2
    exit 1
  fi
}

need() { for f in "$@"; do [ -e "$f" ] || { echo "missing $f -- run the step that makes it first (or set the path)" >&2; exit 1; }; done; }

corpus_paths() { for c in $CORPORA; do echo "${c#*:}"; done; }

combine() {                               # combine <output.csv> <glob of <corpus>_summary.csv> -- adds a 'corpus' column
  $PY - "$@" <<'PY'
import sys, glob, os, pandas as pd
out, pattern = sys.argv[1], sys.argv[2]
frames = []
for f in sorted(glob.glob(pattern)):
    corpus = os.path.basename(f)[:-len('_summary.csv')]
    d = pd.read_csv(f); d.insert(0, 'corpus', corpus); frames.append(d)
pd.concat(frames, ignore_index=True).to_csv(out, index=False)
print(f'   combined {len(frames)} corpora -> {out}')
PY
}

# --------------------------------------------------------------------- steps
step_tests() {
  log_run tests bash -c "cd '$ROOT' && for t in tests/test_*.py; do echo \"### \$t\"; PYTHONPATH='$SRC' $PY \"\$t\" || exit 1; done"
}
step_split() { need "$T2T"; log_run split bash "$SRC/split_t2t.sh" "$T2T" "$CHROMS"; }
step_scan() {
  need "$T2T"
  log_run scan $PY "$SRC/pattern_genome.py" "$T2T" --mappings RY WS MK --families all \
      --min-word-len 13 --max-word-len 10000 --progress \
      --csv "$RESULTS/pattern_scan_T2T.csv"
}
step_context() {
  need "$T2T" "$RESULTS/pattern_scan_T2T.csv"
  local names=-; [ -e "$CHROMS/records.tsv" ] && names=$CHROMS/records.tsv
  log_run context $PY "$SRC/loci_context.py" "$T2T" "$names" "$RESULTS/pattern_scan_T2T.csv" 50 \
      --out "$RESULTS/enriched_loci_context" --progress
}
step_decoys() {
  need "$CHROMS"
  log_run decoys $PY "$SRC/count_genome.py" "$CHROMS" "$N_DECOYS" "$RESULTS/decoy_test_genome.csv" \
      --workers "$WORKERS" --per-chrom --progress
}
step_periodic() {
  need "$CHROMS"
  log_run periodic $PY "$SRC/periodic_decoy.py" "$CHROMS" "$N_DECOYS" "$RESULTS/periodic_test_genome.csv" \
      --workers "$WORKERS" --progress
  log_run periodic_norunrule $PY "$SRC/periodic_decoy.py" "$CHROMS" "$N_DECOYS" \
      "$RESULTS/periodic_test_genome_norunrule.csv" --no-run-rule --workers "$WORKERS" --progress
}
step_longest() {
  need $(corpus_paths)
  mkdir -p "$RESULTS/longest_match"
  for c in $CORPORA; do
    local name=${c%%:*} path=${c#*:}
    log_run "longest__$name" $PY "$SRC/longest_match.py" "$path" "$LM_BASES" \
        "$RESULTS/longest_match/$name.csv" \
        --replicates "$REPLICATES" --workers "$WORKERS" --progress
  done
  combine "$RESULTS/longest_match_summary_all_corpora.csv" "$RESULTS/longest_match/*_summary.csv"
}
step_fibmod4() {
  need $(corpus_paths)
  mkdir -p "$RESULTS/fibmod4"
  for c in $CORPORA; do
    local name=${c%%:*} path=${c#*:}
    log_run "fibmod4__$name" $PY "$SRC/fibmod4.py" "$path" --out "$RESULTS/fibmod4/$name" \
        --replicates "$REPLICATES" --workers "$WORKERS" --progress
  done
  combine "$RESULTS/fibmod4_summary_all_corpora.csv" "$RESULTS/fibmod4/*_summary.csv"
}

# ---------------------------------------------------------------------- main
if [ "${1:-}" = "--list" ]; then
  for s in "${STEPS[@]}"; do printf '  %-9s %s\n' "$s" "$(describe "$s")"; done; exit 0
fi
selected=("$@")
[ ${#selected[@]} -eq 0 ] && selected=("${STEPS[@]}")
for s in "${selected[@]}"; do
  if ! printf '%s\n' "${STEPS[@]}" | grep -qx "$s"; then
    echo "unknown step '$s'; steps: ${STEPS[*]}" >&2; exit 1
  fi
done
mkdir -p "$RESULTS"
echo "run_all.sh: ${selected[*]}  (T2T=$T2T, workers=$WORKERS, decoys=$N_DECOYS, replicates=$REPLICATES)"
for s in "${selected[@]}"; do
  echo "== $s: $(describe "$s")"
  "step_$s"
done
echo "finished; results in $RESULTS/"
