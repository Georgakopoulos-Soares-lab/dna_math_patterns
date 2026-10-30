#!/bin/bash
# Run the complete analysis (code version 2.1) on T2T-CHM13, the six T2T apes,
# the viral, bacterial and archaeal genomes.
#
# usage:  bash run_all.sh                    # every step, every primate genome, in order
#         bash run_all.sh decoys periodic    # only the named steps
#         bash run_all.sh --list             # list the steps and genomes
#
# Which genomes (GENOMES, default all): all | human | apes | viruses | bacteria | archaea | names, e.g.
#         GENOMES=apes bash run_all.sh split scan context decoys periodic combine
#         GENOMES="chimpanzee gorilla" bash run_all.sh scan
#         GENOMES=human bash run_all.sh decoys
#         GENOMES=viruses bash run_all.sh split scan context decoys periodic combine
#         GENOMES=bacteria bash run_all.sh split scan context decoys periodic combine
#         GENOMES=archaea bash run_all.sh split scan context decoys periodic combine
#         GENOMES="bacteria100 archaea100 viruses100" PARALLEL=12 bash run_all.sh split scan context decoys periodic combine
#         GENOMES="e_coli h_pylori" bash run_all.sh scan
# all = the seven primate genomes of the paper (human + apes), as before; the
# viruses, bacteria and archaea run only when named (the group, or their names).
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
# word range, see DECOY_MAX_LEN), and the shared context of apes and viruses
# needs the human scan. Times are estimates (12 cores).
set -euo pipefail

# ------------------------------------------------------------------ settings
ROOT=${ROOT:-/Users/kleiopapatheodorou/Documents/lab/fibonacci/patterns}
HERE=$ROOT/src                            # the .py files and split_t2t.sh
DATA=${DATA:-$ROOT/data}
T2T=${T2T:-$DATA/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna}
CHROMS=${CHROMS:-$DATA/t2t_chroms}        # per-chromosome T2T files, written by step split
APES=${APES:-$DATA/apes}                  # the six ape FASTA files (NCBI primary assemblies v2.1)
VIRUSES=${VIRUSES:-$DATA/viruses}         # one FASTA per virus, all segments (get_viruses.sh)
BACTERIA=${BACTERIA:-$DATA/bacteria}      # one FASTA per bacterium, whole assembly (get_bacteria.sh)
ARCHAEA=${ARCHAEA:-$DATA/archaea}         # one FASTA per archaeon, whole assembly (get_archaea.sh)
RESULTS=${RESULTS:-$ROOT/results_new}
SETS=${SETS:-$DATA/sets}                  # genome sets from select_genomes.py (<set>.tsv)
PARALLEL=${PARALLEL:-1}                   # genomes run at the same time (WORKERS are shared between them)
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
# in NCBI headers). EXCLUDE= (empty) keeps every record. No viral record
# matches it, so every viral segment is kept.
EXCLUDE=${EXCLUDE-unlocalized|unplaced|mitochondri}
# Bacteria: the same, plus plasmids, so that only the chromosomes are analysed
# (as for the primates). EXCLUDE_BACTERIA='unlocalized|unplaced' keeps plasmids.
EXCLUDE_BACTERIA=${EXCLUDE_BACTERIA-unlocalized|unplaced|plasmid}
EXCLUDE_ARCHAEA=${EXCLUDE_ARCHAEA-unlocalized|unplaced|plasmid}        # archaea: the same rule

# Words tested by decoys and periodic: the 8 families from 13 bits up to
# DECOY_MAX_LEN bits. auto (default) = 100, or the longest word observed in that
# genome's scan if longer, so every word present in a genome is tested. The 36
# words of 13-100 bits are tested in every genome (the shared set).
DECOY_MAX_LEN=${DECOY_MAX_LEN:-auto}

# name:FASTA of every primate genome; chm13 is the human reference (results in
# $RESULTS, as before), the apes go to $RESULTS/apes/<name>/.
GENOME_LIST="chm13:$T2T
chimpanzee:$APES/GCF_028858775.2_NHGRI_mPanTro3-v2.1_pri_genomic.fna
bonobo:$APES/GCF_029289425.2_NHGRI_mPanPan1-v2.1_pri_genomic.fna
gorilla:$APES/GCF_029281585.2_NHGRI_mGorGor1-v2.1_pri_genomic.fna
sumatran_orangutan:$APES/GCF_028885655.2_NHGRI_mPonAbe1-v2.1_pri_genomic.fna
bornean_orangutan:$APES/GCF_028885625.2_NHGRI_mPonPyg2-v2.1_pri_genomic.fna
siamang:$APES/GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri_genomic.fna"

# name:FASTA of every virus (RefSeq; segmented viruses with all their segments).
# Results go to $RESULTS/viruses/<name>/, per-segment files to $VIRUSES/chroms/<name>/.
VIRUS_LIST="sars_cov_2:$VIRUSES/sars_cov_2.fna
sars_cov:$VIRUSES/sars_cov.fna
dengue_2:$VIRUSES/dengue_2.fna
zika:$VIRUSES/zika.fna
vsv:$VIRUSES/vsv.fna
hiv_1:$VIRUSES/hiv_1.fna
influenza_a:$VIRUSES/influenza_a.fna
rsv_a:$VIRUSES/rsv_a.fna
hbv:$VIRUSES/hbv.fna
ebv:$VIRUSES/ebv.fna
rotavirus_a:$VIRUSES/rotavirus_a.fna
bluetongue:$VIRUSES/bluetongue.fna"

# name:FASTA of every bacterium (RefSeq reference assemblies, get_bacteria.sh).
# Results go to $RESULTS/bacteria/<name>/, per-chromosome files to $BACTERIA/chroms/<name>/.
BACTERIA_LIST="e_coli:$BACTERIA/e_coli.fna
b_subtilis:$BACTERIA/b_subtilis.fna
s_aureus:$BACTERIA/s_aureus.fna
m_tuberculosis:$BACTERIA/m_tuberculosis.fna
s_coelicolor:$BACTERIA/s_coelicolor.fna
p_aeruginosa:$BACTERIA/p_aeruginosa.fna
v_cholerae:$BACTERIA/v_cholerae.fna
h_pylori:$BACTERIA/h_pylori.fna
m_genitalium:$BACTERIA/m_genitalium.fna
synechocystis:$BACTERIA/synechocystis.fna"

# name:FASTA of every archaeon (RefSeq reference assemblies, get_archaea.sh).
# Results go to $RESULTS/archaea/<name>/, per-chromosome files to $ARCHAEA/chroms/<name>/.
ARCHAEA_LIST="h_volcanii:$ARCHAEA/h_volcanii.fna
h_salinarum:$ARCHAEA/h_salinarum.fna
s_acidocaldarius:$ARCHAEA/s_acidocaldarius.fna
s_solfataricus:$ARCHAEA/s_solfataricus.fna
t_kodakarensis:$ARCHAEA/t_kodakarensis.fna
p_furiosus:$ARCHAEA/p_furiosus.fna
m_maripaludis:$ARCHAEA/m_maripaludis.fna
m_acetivorans:$ARCHAEA/m_acetivorans.fna
m_jannaschii:$ARCHAEA/m_jannaschii.fna
n_maritimus:$ARCHAEA/n_maritimus.fna"
GENOMES=${GENOMES:-all}

STEPS=(tests split scan context decoys periodic longest fibmod4 combine)
GENOME_STEPS=" split scan context decoys periodic "

describe() {
  case $1 in
    tests)    echo "unit tests of all scripts (~5 min)";;
    split)    echo "genome -> one FASTA per chromosome + headers.tsv + records.tsv + skipped.tsv (~5-10 min)";;
    scan)     echo "pattern_genome: 8 families, 13-10000 bits, RY WS MK, both strands, whole genome (hours; minutes for a virus)";;
    context)  echo "loci_context: DNA and tandem period of every locus of the enriched rows (needs scan; apes and viruses also the human rows)";;
    decoys)   echo "count_genome: words of 13-100+ bits x $N_DECOYS run decoys, genome-wide loci (needs split, scan)";;
    periodic) echo "periodic_decoy: nearest period + periodic decoys, strict and --no-run-rule (needs split, scan)";;
    longest)  echo "longest_match: 26 generators vs $LM_BASES bases of each corpus, 3 nulls x $REPLICATES (~15 min/corpus)";;
    fibmod4)  echo "fibmod4: mod-4 recurrence on each corpus, 3 nulls x $REPLICATES (~20-40 min/corpus, ~3 GB/worker)";;
    combine)  echo "one table per analysis with all genomes, in $RESULTS/all_genomes/ (primates), $RESULTS/viruses/all_viruses/, $RESULTS/bacteria/all_bacteria/ and $RESULTS/archaea/all_archaea/ (seconds)";;
  esac
}

# ---------------------------------------------------------------- genomes
all_genomes()   { for e in $GENOME_LIST; do echo "${e%%:*}"; done; }   # the seven primates
virus_genomes() { for e in $VIRUS_LIST;  do echo "${e%%:*}"; done; }
bacteria_genomes() { for e in $BACTERIA_LIST; do echo "${e%%:*}"; done; }
archaea_genomes() { for e in $ARCHAEA_LIST; do echo "${e%%:*}"; done; }
is_virus()      { for e in $VIRUS_LIST;  do [ "${e%%:*}" = "$1" ] && return 0; done; return 1; }
is_bacterium()  { for e in $BACTERIA_LIST; do [ "${e%%:*}" = "$1" ] && return 0; done; return 1; }
is_archaeon()   { for e in $ARCHAEA_LIST; do [ "${e%%:*}" = "$1" ] && return 0; done; return 1; }

# Genome sets: every $SETS/<set>.tsv written by select_genomes.py is a group of
# genomes (first column = name = accession). FASTA in $DATA/<set>/<name>.fna,
# per-chromosome files in $DATA/<set>/chroms/<name>/, results in $RESULTS/<set>/<name>/.
set_list()    { local f b; for f in "$SETS"/*.tsv; do [ -e "$f" ] || continue; b=$(basename "$f" .tsv); case $b in *_candidates) ;; *) echo "$b";; esac; done; }
set_genomes() { tail -n +2 "$SETS/$1.tsv" | cut -f1; }
gset()        { local s; for s in $(set_list); do tail -n +2 "$SETS/$s.tsv" | cut -f1 | grep -xF "$1" >/dev/null && { echo "$s"; return; }; done; return 0; }
is_set()      { [ -n "$1" ] && set_list | grep -xF "$1" >/dev/null; }

genome_fasta()  {
  local e s
  for e in $GENOME_LIST $VIRUS_LIST $BACTERIA_LIST $ARCHAEA_LIST; do [ "${e%%:*}" = "$1" ] && { echo "${e#*:}"; return; }; done
  s=$(gset "$1"); [ -n "$s" ] && echo "$DATA/$s/$1.fna"; return 0
}
gdir() {
  if [ "$1" = chm13 ]; then echo "$RESULTS"
  elif is_virus "$1"; then echo "$RESULTS/viruses/$1"
  elif is_bacterium "$1"; then echo "$RESULTS/bacteria/$1"
  elif is_archaeon "$1"; then echo "$RESULTS/archaea/$1"
  elif [ -n "$(gset "$1")" ]; then echo "$RESULTS/$(gset "$1")/$1"
  else echo "$RESULTS/apes/$1"; fi
}
gchroms() {
  if [ "$1" = chm13 ]; then echo "$CHROMS"
  elif is_virus "$1"; then echo "$VIRUSES/chroms/$1"
  elif is_bacterium "$1"; then echo "$BACTERIA/chroms/$1"
  elif is_archaeon "$1"; then echo "$ARCHAEA/chroms/$1"
  elif [ -n "$(gset "$1")" ]; then echo "$DATA/$(gset "$1")/chroms/$1"
  else echo "$APES/chroms/$1"; fi
}
gexclude() {                              # records left out of <genome>
  if is_bacterium "$1"; then echo "$EXCLUDE_BACTERIA"
  elif is_archaeon "$1"; then echo "$EXCLUDE_ARCHAEA"
  else
    case $(gset "$1") in
      bacteria*) echo "$EXCLUDE_BACTERIA";;
      archaea*)  echo "$EXCLUDE_ARCHAEA";;
      *)         echo "$EXCLUDE";;
    esac
  fi
}
gscan()   { if [ "$1" = chm13 ]; then echo "$RESULTS/pattern_scan_T2T.csv"; else echo "$(gdir "$1")/pattern_scan_$1.csv"; fi; }

expand_genomes() {                        # GENOMES -> list of names, checked
  local out="" g
  for g in $GENOMES; do
    case $g in
      all)     out="$out $(all_genomes | tr '\n' ' ')";;
      human)   out="$out chm13";;
      apes)    out="$out $(all_genomes | grep -vx chm13 | tr '\n' ' ')";;
      viruses) out="$out $(virus_genomes | tr '\n' ' ')";;
      bacteria) out="$out $(bacteria_genomes | tr '\n' ' ')";;
      archaea) out="$out $(archaea_genomes | tr '\n' ' ')";;
      *)       if is_set "$g"; then out="$out $(set_genomes "$g" | tr '\n' ' ')"; continue; fi
               [ -n "$(genome_fasta "$g")" ] || { echo "unknown genome '$g'; genomes: $(all_genomes | tr '\n' ' ')$(virus_genomes | tr '\n' ' ')$(bacteria_genomes | tr '\n' ' ')$(archaea_genomes | tr '\n' ' ')(or all, human, apes, viruses, bacteria, archaea, or a set: $(set_list | tr '\n' ' '))" >&2; exit 1; }
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

log_run_context() {                       # like log_run, but "no enriched rows" is a result, not an error
  local name=$1; shift
  mkdir -p "$RESULTS/logs"
  local t0=$SECONDS
  if "$@" > "$RESULTS/logs/$name.log" 2>&1; then
    echo "   $name: done in $((SECONDS - t0)) s (log: $RESULTS/logs/$name.log)"
  elif grep -q '^no rows with O/E' "$RESULTS/logs/$name.log"; then
    echo "   $name: no enriched rows (O/E > 50 with a hit) in this genome, nothing to describe; continuing"
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

exclude_args() {                          # exclude_args <genome>
  local ex; ex=$(gexclude "$1")
  [ -n "$ex" ] && printf '%s\n' --exclude-regex "$ex"; return 0
}

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
  log_run "split__$g" bash "$SRC/split_t2t.sh" "$fa" "$(gchroms "$g")" "$(gexclude "$g")"
}

step_scan() {
  local g=$1 fa; fa=$(genome_fasta "$g")
  need "$fa"; mkdir -p "$(gdir "$g")"
  local ex=(); while IFS= read -r x; do ex+=("$x"); done < <(exclude_args "$g")
  log_run "scan__$g" $PY "$SRC/pattern_genome.py" "$fa" --mappings RY WS MK --families all \
      --min-word-len 13 --max-word-len 10000 ${ex[@]+"${ex[@]}"} --progress \
      --csv "$(gscan "$g")"
}

step_context() {
  local g=$1 fa; fa=$(genome_fasta "$g")
  local dir; dir=$(gdir "$g")
  need "$fa" "$(gscan "$g")"
  local names=-; [ -e "$(gchroms "$g")/records.tsv" ] && names=$(gchroms "$g")/records.tsv
  local ex=(); while IFS= read -r x; do ex+=("$x"); done < <(exclude_args "$g")
  # own: the rows enriched in this genome's own scan
  log_run_context "context__$g" $PY "$SRC/loci_context.py" "$fa" "$names" "$(gscan "$g")" 50 \
      --out "$dir/enriched_loci_context" ${ex[@]+"${ex[@]}"} --progress
  # shared (apes and viruses): the rows enriched in the HUMAN scan, searched in this genome
  if [ "$g" != chm13 ]; then
    if [ -e "$(gscan chm13)" ]; then
      log_run_context "context_shared__$g" $PY "$SRC/loci_context.py" "$fa" "$names" "$(gscan chm13)" 50 \
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

genome_args() {                           # genome_args <names...> -> --genome arguments for combine_genomes.py
  local g
  for g in "$@"; do
    printf '%s\n' --genome "$g" "$(genome_fasta "$g")" "$(gdir "$g")" "$(gscan "$g")" "$(gchroms "$g")"
  done
}

combine_group() {                         # combine_group <group> <names...>: every genome with a scan -> $RESULTS/<group>/all_<group>/
  local group=$1; shift
  local have=() g
  for g in "$@"; do [ -e "$(gscan "$g")" ] && have+=("$g"); done
  if [ ${#have[@]} -eq 0 ]; then
    echo "   combine_$group: skipped, no genome of this group has a scan yet" >&2; return 0
  fi
  local args=(); while IFS= read -r x; do args+=("$x"); done < <(genome_args "${have[@]}")
  log_run "combine_$group" $PY "$SRC/combine_genomes.py" --out "$RESULTS/$group/all_$group" "${args[@]}"
  cat "$RESULTS/logs/combine_$group.log"
}

step_combine() {
  # primates: all seven, whatever GENOMES says (unchanged), into $RESULTS/all_genomes/
  # viruses, bacteria: every genome of the group that has a scan, into $RESULTS/<group>/all_<group>/
  local sel=" $genomes " g
  local want_primates=0 want_viruses=0 want_bacteria=0 want_archaea=0
  for g in $(all_genomes);      do case $sel in *" $g "*) want_primates=1;; esac; done
  for g in $(virus_genomes);    do case $sel in *" $g "*) want_viruses=1;; esac; done
  for g in $(bacteria_genomes); do case $sel in *" $g "*) want_bacteria=1;; esac; done
  for g in $(archaea_genomes);  do case $sel in *" $g "*) want_archaea=1;; esac; done
  if [ $want_primates = 1 ]; then
    local args=(); while IFS= read -r x; do args+=("$x"); done < <(genome_args $(all_genomes))
    log_run combine $PY "$SRC/combine_genomes.py" --out "$RESULTS/all_genomes" "${args[@]}"
    cat "$RESULTS/logs/combine.log"
  fi
  [ $want_viruses = 1 ]  && combine_group viruses  $(virus_genomes)
  [ $want_bacteria = 1 ] && combine_group bacteria $(bacteria_genomes)
  [ $want_archaea = 1 ]  && combine_group archaea  $(archaea_genomes)
  # genome sets: each set with a selected genome, into $RESULTS/<set>/all_<set>/
  local s
  for s in $(for g in $genomes; do gset "$g"; done | awk 'NF && !seen[$0]++'); do
    combine_group "$s" $(set_genomes "$s")
  done
  return 0
}

# ---------------------------------------------------------------------- main
if [ "${1:-}" = "--list" ]; then
  for s in "${STEPS[@]}"; do printf '  %-9s %s\n' "$s" "$(describe "$s")"; done
  echo "genomes (GENOMES=all | human | apes | viruses | bacteria | archaea | names):"
  for e in $GENOME_LIST; do printf '  %-19s %s\n' "${e%%:*}" "${e#*:}"; done
  echo "viruses (GENOMES=viruses, not part of all):"
  for e in $VIRUS_LIST;  do printf '  %-19s %s\n' "${e%%:*}" "${e#*:}"; done
  echo "bacteria (GENOMES=bacteria, not part of all; exclude='$EXCLUDE_BACTERIA'):"
  for e in $BACTERIA_LIST; do printf '  %-19s %s\n' "${e%%:*}" "${e#*:}"; done
  echo "archaea (GENOMES=archaea, not part of all; exclude='$EXCLUDE_ARCHAEA'):"
  for e in $ARCHAEA_LIST; do printf '  %-19s %s\n' "${e%%:*}" "${e#*:}"; done
  echo "genome sets in $SETS (GENOMES=<set>):"
  for s in $(set_list); do printf '  %-19s %s genomes, FASTA in %s\n' "$s" "$(set_genomes "$s" | wc -l | tr -d ' ')" "$DATA/$s/"; done
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
echo "   workers=$WORKERS, decoys=$N_DECOYS, replicates=$REPLICATES, exclude='$EXCLUDE' (bacteria, archaea: '$EXCLUDE_BACTERIA', '$EXCLUDE_ARCHAEA'), decoy max length=$DECOY_MAX_LEN"

is_selected() { printf '%s\n' "${selected[@]}" | grep -qx "$1"; }

is_selected tests && { echo "== tests: $(describe tests)"; step_tests; }
per_genome=(); for s in "${selected[@]}"; do case "$GENOME_STEPS" in *" $s "*) per_genome+=("$s");; esac; done
n_genomes=$(echo $genomes | wc -w | tr -d ' ')
if [ "$PARALLEL" -gt 1 ] && [ ${#per_genome[@]} -gt 0 ] && [ "$n_genomes" -gt 1 ]; then
  # several genomes at once, each in its own run of this script with WORKERS/PARALLEL
  # workers; a genome that fails is reported and the others carry on
  w=$(( WORKERS / PARALLEL )); [ "$w" -lt 1 ] && w=1
  mkdir -p "$RESULTS/logs"
  echo "== ${per_genome[*]}: $n_genomes genomes, $PARALLEL at a time, $w worker(s) each (per-genome log: $RESULTS/logs/run__<genome>.log)"
  fails="$RESULTS/logs/parallel_failed.txt"; : > "$fails"; k=0
  for g in $genomes; do
    while [ "$(jobs -rp | wc -l | tr -d ' ')" -ge "$PARALLEL" ]; do sleep 1; done
    k=$((k + 1))
    ( if GENOMES="$g" PARALLEL=1 WORKERS="$w" bash "$0" "${per_genome[@]}" > "$RESULTS/logs/run__$g.log" 2>&1; then
        echo "   [$k/$n_genomes] $g: done"
      else
        echo "   [$k/$n_genomes] $g: FAILED -- see $RESULTS/logs/run__$g.log"; echo "$g" >> "$fails"
      fi ) &
  done
  wait
  nf=$(wc -l < "$fails" | tr -d ' ')
  [ "$nf" -gt 0 ] && echo "   $nf genome(s) failed, listed in $fails; rerun them with GENOMES=\"\$(cat $fails)\"" >&2
else
  for g in $genomes; do
    for s in "${selected[@]}"; do
      case "$GENOME_STEPS" in *" $s "*) echo "== $s [$g]: $(describe "$s")"; "step_$s" "$g";; esac
    done
  done
fi
for s in longest fibmod4 combine; do
  is_selected "$s" && { echo "== $s: $(describe "$s")"; "step_$s"; }
done
echo "finished; results in $RESULTS/"
