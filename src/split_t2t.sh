#!/bin/bash
# Split a genome FASTA (T2T-CHM13v2.0 or any other) into one FASTA per record.
# Needed by count_genome.py and periodic_decoy.py.
#
# usage: bash split_t2t.sh <genome fasta[.gz]> <outdir>
#   e.g. T2T=/scratch/10917/michalakis/compression_paper/data/fasta/work/T2T-CHM13v2.0.fa
#        bash src/split_t2t.sh "$T2T" data/t2t_chroms
#
# Output in <outdir>:
#   <chrom>.fa     each record exactly as in the input (header line included,
#                  line breaks unchanged), in a single pass over the input
#   headers.tsv    "<header line number>\t<chrom>" -- same format as before; the
#                  line number refers to the (decompressed) input file and is
#                  used by loci_context.py
#   records.tsv    chrom, record_id, header_line, lines -- links the record ID
#                  used in scan results (e.g. NC_060925.1) to the file name (chr1)
#
# File names (<chrom>), first rule that applies:
#   1. header contains "chromosome N"  (NCBI style)  -> chrN   (N = 1..22, X, Y, M)
#   2. record ID starts with "chr"     (UCSC style)  -> the record ID
#   3. header mentions "mitochondri"                 -> chrM
#   4. otherwise                                     -> the record ID
# Characters other than A-Z a-z 0-9 . _ - are replaced by "_". The script stops
# with an error if two records would get the same name, or if sequence appears
# before the first header. At the end it checks that every file has the
# expected number of lines and prints a summary.
#
# Input may be plain or gzip-compressed. Works with any POSIX awk (gawk, mawk,
# BSD awk). A final line without a newline is kept (the earlier version dropped
# it).
set -euo pipefail

usage="usage: bash split_t2t.sh <genome fasta[.gz]> <outdir>"
T=${1:?$usage}
OUT=${2:?$usage}
[ -r "$T" ] || { echo "error: cannot read $T" >&2; exit 1; }
mkdir -p "$OUT"

gzip -cdf -- "$T" | awk -v out="$OUT" '
function finish() {
    if (file != "") {
        close(file)
        print name "\t" id "\t" hline "\t" n > rec
    }
}
BEGIN {
    hdr = out "/headers.tsv"; rec = out "/records.tsv"; file = ""; bad = 0
    printf "" > hdr
    print "chrom\trecord_id\theader_line\tlines" > rec
}
/^>/ {
    finish()
    h = $0; sub(/\r$/, "", h)
    id = substr(h, 2); sub(/[ \t].*/, "", id)
    if (match(h, /chromosome [0-9XYM]+/)) name = "chr" substr(h, RSTART + 11, RLENGTH - 11)
    else if (id ~ /^chr/)                 name = id
    else if (tolower(h) ~ /mitochondri/)  name = "chrM"
    else                                  name = id
    gsub(/[^A-Za-z0-9._-]/, "_", name)
    if (name == "") name = "record_line" NR
    if (name in seen) {
        printf "error: the records at lines %d and %d would both be written to %s.fa\n", seen[name], NR, name > "/dev/stderr"
        bad = 1; exit 1
    }
    seen[name] = NR
    file = out "/" name ".fa"; hline = NR; n = 0
    print NR "\t" name > hdr
}
{
    if (file == "") {
        print "error: sequence found before the first FASTA header (line " NR ")" > "/dev/stderr"
        bad = 1; exit 1
    }
    print > file
    n++
}
END { if (!bad) finish() }
'

# check: every file has the number of lines recorded while splitting
tail -n +2 "$OUT/records.tsv" | while IFS=$'\t' read -r chrom id hl n; do
    got=$(wc -l < "$OUT/$chrom.fa"); got=${got//[[:space:]]/}
    if [ "$got" != "$n" ]; then
        echo "error: $OUT/$chrom.fa has $got lines, expected $n (disk full?)" >&2
        exit 1
    fi
done

nrec=$(( $(wc -l < "$OUT/records.tsv") - 1 ))
echo "Wrote $nrec FASTA files to $OUT (plus headers.tsv and records.tsv):"
awk -F'\t' 'NR > 1 { printf "  %-10s %-16s %12d lines\n", $1, $2, $4 }' "$OUT/records.tsv"
