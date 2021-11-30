#!/usr/bin/env bash
#
# Use MACS2 to build forward and reverse peaks files for given .bam file.
# Iterate peaks through set of criteria to determine UTR viability, before annotating in .gff file.
#
# Usage: annotate_3pUTR.sh [--max-distance=N] BAM_FILE GFF_FILE ...
#
# Arguments:
#   BAM_FILE            input bam file
#   GFF_FILE            input 'canonical' gff file.
#
# Options:
#   --max-distance=N    limit the number of line to display
#
# Examples:
#    ./annotate_3pUTR.sh --max-distance=2500 in.bam in.gff
#

source docopts.sh
help=$(docopt_get_help_string $0)
version='0.1'

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" || exit ; pwd -P )
export PATH=$parent_path:$PATH

parsed=$(docopts -G args -h "$help" -V $version : "$@")
eval "$parsed"

GFF_IN="$2"

PYTHON_CODE=$(cat <<END
import os.path

import gffutils

gff_db = os.path.splitext("$GFF_IN")[0] + '.db'
if not os.path.isfile(gff_db):
    gffutils.create_db("$GFF_IN", gff_db, force=True)
print(gff_db)
END
)

GFF_DB="$(python3 -c "$PYTHON_CODE")"

call_forward_peaks.sh "$1" &
call_reverse_peaks.sh "$1" &

wait

peaks_to_UTR.py "$GFF_DB" forward_peaks.broadPeak forward --max-distance "$3" &
peaks_to_UTR.py "$GFF_DB" reverse_peaks.broadPeak reverse --max-distance "$3" &

wait

cat forward_three_prime_UTRs.gff reverse_three_prime_UTRs.gff > three_prime_UTRs.gff

if [[ $GFF_IN == *.gff ]]; then
  cat "$GFF_IN" three_prime_UTRs.gff > full.gff

  gt gff3 -sort -retainids -tidy full.gff > full.sorted.gff
fi
