#!/bin/bash
# Run the complete analysis (code version 2.1) on T2T-CHM13 and the six T2T apes.
#
# usage:  bash run_all.sh                    # every step, every genome, in order
#         bash run_all.sh decoys periodic    # only the named steps
#         bash run_all.sh --list             # list the steps and genomes
#
# Which genomes (GENOMES, default all): all | human | apes | names, e.g.
#         GENOMES=apes bash run_all.sh split scan context decoys periodic combine
#         GENOMES="chimpanzee gorilla" bash run_all.sh scan
#         GENOMES=human bash run_all.sh decoys
# The per-genome steps (split scan context decoys periodic) run genome by genome:
# all selected steps for the first genome, then the next. tests runs first;
# longest, fibmod4 (human corpora only) and combine run last, once.
#
# Settings below can also be given as environment variables, e.g.
#         WORKERS=24 REPLICATES=20 bash run_all.sh longest
#
# Each step logs to $RESULTS/logs/<step>__<genome>.log (per corpus for
# longest/fibmod4) and stops the script on error. Run the steps in order the
# first time: context needs scan, decoys/periodic need split (and scan, for the
# word range, see DECOY_MAX_LEN), and the apes' shared context needs the human
# scan. Times are estimates (12 cores).
set -euo pipefail

# ------------------------------------------------------------------ settings
ROOT=${ROOT:-/Users/kleiopapatheodorou/Documents/lab/fibonacci/patterns}
HERE=$ROOT/src                            # the .py files and split_t2t.sh
DATA=${DATA:-$ROOT/data}
T2T=${T2T:-$DATA/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna}
CHROMS=${CHROMS:-$DATA/t2t_chroms}        # per-chromosome T2T files, written by step split
APES=${APES:-$DATA/apes}                  # the six ape FASTA files (NCBI primary assemblies v2.1)
RESULTS=${RESULTS:-$ROOT/results_new}
WORKERS=${WORKERS:-12}
N_DECOYS=${N_DECOYS:-100}               # decoys per target word
REPLICATES=${REPLICATES:-100}           # null replicates (longest, fibmod4)
LM_BASES=${LM_BASES:-8000000}           # longest: ACGT bases per corpus (same for all corpora)
# corpora for longest and fibmod4: name:path, space separated
CORPORA=${CORPORA:-"hs_chr20:$DATA/hs_chr20.fa hs_chr21:$DATA/hs_chr21.fa mm_chr19:$DATA/mm_chr19.fa wheat_1A:$DATA/wheat_1A.fa t2t_chr1_cen:$DATA/t2t_chr1_cen.fa"}
PY=${PY:-/opt/anaconda3/bin/python}
SRC=$HERE                                 # folder with the .py files and split_t2t.sh

# Records left out of every genome (split, scan, context): unplaced and
# unlocalized scaffolds and the mitochondrion, so that only the assembled
# chromosomes are analysed. Matched against the FASTA header line (lower case
# in NCBI headers). EXCLUDE= (empty) keeps every record.
EXCLUDE=${EXCLUDE-unlocalized|unplaced|mitochondri}

# Words tested by decoys and periodic: the 8 families from 13 bits up to
# DECOY_MAX_LEN bits. auto (default) = 100, or the longest word observed in that
# genome's scan if longer, so every word present in a genome is tested. The 36
# words of 13-100 bits are tested in every genome (the shared set).
DECOY_MAX_LEN=${DECOY_MAX_LEN:-auto}

# name:FASTA of every genome; chm13 is the human reference (results in $RESULTS,
# as before), the apes go to $RESULTS/apes/<name>/.
GENOME_LIST="chm13:$T2T
chimpanzee:$APES/GCF_028858775.2_NHGRI_mPanTro3-v2.1_pri_genomic.fna
bonobo:$APES/GCF_029289425.2_NHGRI_mPanPan1-v2.1_pri_genomic.fna
gorilla:$APES/GCF_029281585.2_NHGRI_mGorGor1-v2.1_pri_genomic.fna
sumatran_orangutan:$APES/GCF_028885655.2_NHGRI_mPonAbe1-v2.1_pri_genomic.fna
bornean_orangutan:$APES/GCF_028885625.2_NHGRI_mPonPyg2-v2.1_pri_genomic.fna
siamang:$APES/GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri_genomic.fna"
GENOMES=${GENOMES:-all}

STEPS=(tests split scan context decoys periodic longest fibmod4 combine)
GENOME_STEPS=" split scan context decoys periodic "

describe() {
  case $1 in
    tests)    echo "unit tests of all scripts (~5 min)";;
    split)    echo "genome -> one FASTA per chromosome + headers.tsv + records.tsv + skipped.tsv (~5-10 min)";;
    scan)     echo "pattern_genome: 8 families, 13-10000 bits, RY WS MK, both strands, whole genome (hours)";;
    context)  echo "loci_context: DNA and tandem period of every locus of the enriched rows (needs scan; apes also the human rows)";;
    decoys)   echo "count_genome: words of 13-100+ bits x $N_DECOYS run decoys, genome-wide loci (needs split, scan)";;
    periodic) echo "periodic_decoy: nearest period + periodic decoys, strict and --no-run-rule (needs split, scan)";;
    longest)  echo "longest_match: 26 generators vs $LM_BASES bases of each corpus, 3 nulls x $REPLICATES (~15 min/corpus)";;
    fibmod4)  echo "fibmod4: mod-4 recurrence on each corpus, 3 nulls x $REPLICATES (~20-40 min/corpus, ~3 GB/worker)";;
    combine)  echo "one table per analysis with all genomes, in $RESULTS/all_genomes/ (seconds)";;
  esac
}

# ---------------------------------------------------------------- genomes
all_genomes() { for e in $GENOME_LIST; do echo "${e%%:*}"; done; }
genome_fasta() { for e in $GENOME_LIST; do [ "${e%%:*}" = "$1" ] && { echo "${e#*:}"; return; }; done; }
gdir()    { if [ "$1" = chm13 ]; then echo "$RESULTS"; else echo "$RESULTS/apes/$1"; fi; }
gchroms() { if [ "$1" = chm13 ]; then echo "$CHROMS"; else echo "$APES/chroms/$1"; fi; }
gscan()   { if [ "$1" = chm13 ]; then echo "$RESULTS/pattern_scan_T2T.csv"; else echo "$(gdir "$1")/pattern_scan_$1.csv"; fi; }

expand_genomes() {                        # GENOMES -> list of names, checked
  local out="" g
  for g in $GENOMES; do
    case $g in
      all)   out="$out $(all_genomes | tr '\n' ' ')";;
      human) out="$out chm13";;
      apes)  out="$out $(all_genomes | grep -vx chm13 | tr '\n' ' ')";;
      *)     [ -n "$(genome_fasta "$g")" ] || { echo "unknown genome '$g'; genomes: $(all_genomes | tr '\n' ' ')(or all, human, apes)" >&2; exit 1; }
             out="$out $g";;
    esac
  done
  echo $out | tr ' ' '\n' | awk '!seen[$0]++' | tr '\n' ' '
}

# ----------------------------------------------------------------- helpers
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

exclude_args() { [ -n "$EXCLUDE" ] && printf '%s\n' --exclude-regex "$EXCLUDE"; return 0; }

decoy_max_len() {                         # decoy_max_len <genome> -> upper word length for decoys/periodic
  if [ "$DECOY_MAX_LEN" != auto ]; then echo "$DECOY_MAX_LEN"; return; fi
  local scan; scan=$(gscan "$1")
  if [ ! -e "$scan" ]; then
    echo "   note: no scan for $1 yet, so decoys/periodic use words up to 100 bits" >&2
    echo 100; return
  fi
  $PY - "$scan" <<'PY'
import sys, pandas as pd
d = pd.read_csv(sys.argv[1], usecols=['length', 'observed_total'])
seen = d.loc[d.observed_total > 0, 'length']
print(max(100, int(seen.max()) if len(seen) else 100))
PY
}

# --------------------------------------------------------------------- steps
step_tests() {
  log_run tests bash -c "cd '$ROOT' && for t in tests/test_*.py; do echo \"### \$t\"; PYTHONPATH='$SRC' $PY \"\$t\" || exit 1; done"
}

step_split() {                            # step_split <genome>
  local g=$1 fa; fa=$(genome_fasta "$g")
  need "$fa"
  log_run "split__$g" bash "$SRC/split_t2t.sh" "$fa" "$(gchroms "$g")" "$EXCLUDE"
}

step_scan() {
  local g=$1 fa; fa=$(genome_fasta "$g")
  need "$fa"; mkdir -p "$(gdir "$g")"
  local ex=(); while IFS= read -r x; do ex+=("$x"); done < <(exclude_args)
  log_run "scan__$g" $PY "$SRC/pattern_genome.py" "$fa" --mappings RY WS MK --families all \
      --min-word-len 13 --max-word-len 10000 ${ex[@]+"${ex[@]}"} --progress \
      --csv "$(gscan "$g")"
}

step_context() {
  local g=$1 fa; fa=$(genome_fasta "$g")
  local dir; dir=$(gdir "$g")
  need "$fa" "$(gscan "$g")"
  local names=-; [ -e "$(gchroms "$g")/records.tsv" ] && names=$(gchroms "$g")/records.tsv
  local ex=(); while IFS= read -r x; do ex+=("$x"); done < <(exclude_args)
  # own: the rows enriched in this genome's own scan
  log_run "context__$g" $PY "$SRC/loci_context.py" "$fa" "$names" "$(gscan "$g")" 50 \
      --out "$dir/enriched_loci_context" ${ex[@]+"${ex[@]}"} --progress
  # shared (apes only): the rows enriched in the HUMAN scan, searched in this genome
  if [ "$g" != chm13 ]; then
    if [ -e "$(gscan chm13)" ]; then
      log_run "context_shared__$g" $PY "$SRC/loci_context.py" "$fa" "$names" "$(gscan chm13)" 50 \
          --out "$dir/enriched_loci_context_shared" ${ex[@]+"${ex[@]}"} --progress
    else
      echo "   context_shared__$g: skipped, no human scan yet ($(gscan chm13))" >&2
    fi
  fi
}

step_decoys() {
  local g=$1 chroms; chroms=$(gchroms "$g")
  need "$chroms/records.tsv"; mkdir -p "$(gdir "$g")"
  local maxlen; maxlen=$(decoy_max_len "$g")
  log_run "decoys__$g" $PY "$SRC/count_genome.py" "$chroms" "$N_DECOYS" "$(gdir "$g")/decoy_test_genome.csv" \
      --max-len "$maxlen" --workers "$WORKERS" --per-chrom --progress
}

step_periodic() {
  local g=$1 chroms; chroms=$(gchroms "$g")
  need "$chroms/records.tsv"; mkdir -p "$(gdir "$g")"
  local maxlen; maxlen=$(decoy_max_len "$g")
  log_run "periodic__$g" $PY "$SRC/periodic_decoy.py" "$chroms" "$N_DECOYS" "$(gdir "$g")/periodic_test_genome.csv" \
      --max-len "$maxlen" --workers "$WORKERS" --progress
  log_run "periodic_norunrule__$g" $PY "$SRC/periodic_decoy.py" "$chroms" "$N_DECOYS" \
      "$(gdir "$g")/periodic_test_genome_norunrule.csv" --max-len "$maxlen" --no-run-rule \
      --workers "$WORKERS" --progress
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

step_combine() {                          # all genomes with results, whatever GENOMES says
  local args=() g
  for g in $(all_genomes); do
    args+=(--genome "$g" "$(genome_fasta "$g")" "$(gdir "$g")" "$(gscan "$g")" "$(gchroms "$g")")
  done
  log_run combine $PY "$SRC/combine_genomes.py" --out "$RESULTS/all_genomes" "${args[@]}"
  cat "$RESULTS/logs/combine.log"
}

# ---------------------------------------------------------------------- main
if [ "${1:-}" = "--list" ]; then
  for s in "${STEPS[@]}"; do printf '  %-9s %s\n' "$s" "$(describe "$s")"; done
  echo "genomes (GENOMES=all | human | apes | names):"
  for e in $GENOME_LIST; do printf '  %-19s %s\n' "${e%%:*}" "${e#*:}"; done
  exit 0
fi
selected=("$@")
[ ${#selected[@]} -eq 0 ] && selected=("${STEPS[@]}")
for s in "${selected[@]}"; do
  if ! printf '%s\n' "${STEPS[@]}" | grep -qx "$s"; then
    echo "unknown step '$s'; steps: ${STEPS[*]}" >&2; exit 1
  fi
done
genomes=$(expand_genomes)
mkdir -p "$RESULTS"
echo "run_all.sh: ${selected[*]}"
echo "   genomes: $genomes"
echo "   workers=$WORKERS, decoys=$N_DECOYS, replicates=$REPLICATES, exclude='$EXCLUDE', decoy max length=$DECOY_MAX_LEN"

is_selected() { printf '%s\n' "${selected[@]}" | grep -qx "$1"; }

is_selected tests && { echo "== tests: $(describe tests)"; step_tests; }
for g in $genomes; do
  for s in "${selected[@]}"; do
    case "$GENOME_STEPS" in *" $s "*) echo "== $s [$g]: $(describe "$s")"; "step_$s" "$g";; esac
  done
done
for s in longest fibmod4 combine; do
  is_selected "$s" && { echo "== $s: $(describe "$s")"; "step_$s"; }
done
echo "finished; results in $RESULTS/"
