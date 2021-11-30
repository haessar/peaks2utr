#!/usr/bin/env bash
#
# cat -n all files
#
# Usage: cat-n_wrapper_example.sh [--count=N] FILE...
#
# Arguments:
#   FILE     input file, if FILLE equal - stdin is used instead.
#
# Options:
#   --count=N   limit the number of line to display
#
# Examples:
#    ./cat-n_wrapper.sh --count=3 cat-n_wrapper.sh  quick_example.sh
#
parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" || exit ; pwd -P )
export PATH=$parent_path:$PATH

source docopts.sh
help=$(docopt_get_help_string $0)
version='0.1'

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
