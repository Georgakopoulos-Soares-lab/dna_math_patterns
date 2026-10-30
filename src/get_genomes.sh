#!/usr/bin/env bash
# Download every genome of a set chosen by select_genomes.py, one FASTA per genome.
#
# usage: bash get_genomes.sh bacteria100          # reads data/sets/bacteria100.tsv
#        bash get_genomes.sh archaea100 viruses100
# Each genome is saved as data/<set>/<accession>.fna (whole assembly: all chromosomes,
# plasmids or segments; run_all.sh decides which records are analysed). Genomes
# already downloaded are skipped, so the script can simply be run again after an
# interruption.
# needs: curl, unzip
set -euo pipefail

ROOT=${ROOT:-/Users/kleiopapatheodorou/Documents/lab/fibonacci/patterns}
SETS=${SETS:-$ROOT/data/sets}
API=https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession
[ $# -ge 1 ] || { echo "usage: bash get_genomes.sh <set> [<set> ...]   (sets: $(ls "$SETS" 2>/dev/null | grep -v candidates | sed 's/\.tsv$//' | tr '\n' ' '))" >&2; exit 1; }
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT

for set in "$@"; do
  tsv=$SETS/$set.tsv
  [ -s "$tsv" ] || { echo "missing $tsv -- run select_genomes.py first" >&2; exit 1; }
  out=$ROOT/data/$set; mkdir -p "$out"
  total=$(($(wc -l < "$tsv") - 1)); i=0; failed=0
  # the list is read on file descriptor 3, so curl cannot consume it
  while IFS=$'\t' read -r name acc organism len rest <&3; do
    i=$((i + 1))
    f=$out/$name.fna
    if [ -s "$f" ]; then continue; fi
    if curl -sS --fail -L --retry 3 --retry-delay 3 -o "$TMP/g.zip" \
         "$API/$acc/download?include_annotation_type=GENOME_FASTA" \
       && unzip -p "$TMP/g.zip" 'ncbi_dataset/data/*/*.fna' > "$f.tmp" && [ -s "$f.tmp" ]; then
      mv "$f.tmp" "$f"
      printf '%s  %3d/%d  %s  %s\n' "$set" "$i" "$total" "$acc" "$organism"
    else
      rm -f "$f.tmp"; failed=$((failed + 1))
      echo "$set  $i/$total  $acc  FAILED -- run the script again to retry" >&2
    fi
    sleep 0.4      # NCBI allows a few requests per second without an API key
  done 3< <(tail -n +2 "$tsv")
  have=$(ls "$out"/*.fna 2>/dev/null | wc -l | tr -d ' ')
  echo "$set: $have of $total genomes in $out ($failed failed this run)"
done
