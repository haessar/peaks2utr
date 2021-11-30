#!/usr/bin/env bash
#
# Use MACS2 to build forward and reverse peaks files for given .bam file.
# Iterate peaks through set of criteria to determine UTR viability, before annotating in .gff file.
#
# Usage: annotate_3pUTR.sh [--max-distance=N] BAM_FILE GFF_FILE ...
#
# Arguments:
#   BAM_FILE            input reads file in bam format.
#   GFF_FILE            input 'canonical' annotations file in gff or gtf format.
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

re='^[0-9]+$'
if [[ ! (-z $args_max_distance) && ! ($args_max_distance =~ $re) ]] ; then
  echo "$help"
  echo "ERROR: max distance not a number" >&2; exit 1
fi

if [[ (! -f "$args_BAM_FILE") ||  ($args_BAM_FILE != *.bam)]]; then
  echo "$help"
  echo "ERROR: BAM_FILE doesn't exist or in incorrect format" >&2; exit 1
fi

if [[ (! -f "$args_GFF_FILE") || ! ($args_GFF_FILE = *.@(gff|gtf)) ]]; then
  echo "$help"
  echo "ERROR: GFF_FILE doesn't exist or in incorrect format" >&2; exit 1
fi

if [[ -z "$args_max_distance" ]] ; then
  echo "setting max distance to default 2500"
  args_max_distance=2500
fi

PYTHON_CODE=$(cat <<END
import os.path

import gffutils

gff_db = os.path.splitext("$args_GFF_FILE")[0] + '.db'
if not os.path.isfile(gff_db):
    gffutils.create_db("$args_GFF_FILE", gff_db, force=True)
print(gff_db)
END
)
GFF_DB="$(python3 -c "$PYTHON_CODE")"

BAM_BASENAME=$(basename $args_BAM_FILE .bam)
call_forward_peaks.sh "$BAM_BASENAME" &
call_reverse_peaks.sh "$BAM_BASENAME" &

wait

peaks_to_UTR.py "$GFF_DB" forward_peaks.broadPeak forward --max-distance "$args_max_distance" &
peaks_to_UTR.py "$GFF_DB" reverse_peaks.broadPeak reverse --max-distance "$args_max_distance" &

wait

cat forward_three_prime_UTRs.gff reverse_three_prime_UTRs.gff > three_prime_UTRs.gff

if [[ $args_GFF_FILE == *.gff ]]; then
  cat "$args_GFF_FILE" three_prime_UTRs.gff > full.gff

  gt gff3 -sort -retainids -tidy full.gff > full.sorted.gff
fi
