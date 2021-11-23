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

../../call_forward_peaks.sh "$1" &
../../call_reverse_peaks.sh "$1" &

wait

../../peaks_to_3pUTR/peaks_to_UTR.py "$GFF_DB" forward_peaks.broadPeak forward --max-distance "$3" &
../../peaks_to_3pUTR/peaks_to_UTR.py "$GFF_DB" reverse_peaks.broadPeak reverse --max-distance "$3" &

wait

cat forward_three_prime_UTRs.gff reverse_three_prime_UTRs.gff > three_prime_UTRs.gff

if [[ $GFF_IN == *.gff ]]; then
  cat "$GFF_IN" three_prime_UTRs.gff > full.gff

  gt gff3 -sort -retainids -tidy full.gff > full.sorted.gff
fi
