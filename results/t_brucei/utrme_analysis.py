#%%
import os.path
from pathlib import Path
import sys

import gffutils

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from peaks2utr.models import UTR

#%%
base_dir = os.path.dirname(__file__)
out_dir = os.path.join(base_dir, "output")
Path(out_dir).mkdir(parents=True, exist_ok=True)
new = os.path.join(base_dir, "input/utrme_output.gff3")

db_new = gffutils.create_db(new, os.path.join(out_dir, "utrme.db"), force=True, id_spec={"three_prime_UTR": "score", "gene": "ID", "mRNA": "ID"})


#%%
total_canonical_genes = 0
total_utrs = 0
total_matched = 0
total_extended = 0
total_reduced = 0
total_new = 0
total_existing = 0
total_missing = 0
for gene in db_new.all_features(featuretype="gene"):
    total_canonical_genes += 1
    gene_id = gene.id
    new_utrs = list(db_new.children(id=gene_id, featuretype="three_prime_UTR"))
    if new_utrs:
        total_utrs += 1
    if len(new_utrs) > 1:
        utrme_utrs = [f for f in new_utrs if f.source=="UTRme"]
        canonical_utrs = [f for f in new_utrs if f.source=="EuPathDB"]
        if any(utrme_utrs) and any(canonical_utrs):
            nUTR = UTR(start=utrme_utrs[0].start, end=utrme_utrs[0].end)
            cUTR = UTR(start=min(utr.start for utr in canonical_utrs), end=max(utr.end for utr in canonical_utrs))
            if gene.strand == "+":
                if nUTR.end > cUTR.end:
                    total_extended += 1
                elif nUTR.end < cUTR.end:
                    total_reduced += 1
                else:
                    total_matched += 1
            else:
                if nUTR.start < cUTR.start:
                    total_extended += 1
                elif nUTR.start > cUTR.start:
                    total_reduced += 1
                else:
                    total_matched += 1
        else:
            if any(utrme_utrs):
                total_new += 1
            elif any(canonical_utrs):
                total_existing +=  1
    elif len(new_utrs) == 1:
        if new_utrs[0].source == "UTRme":
            total_new += 1
        elif new_utrs[0].source == "EuPathDB":
            total_existing += 1
    else:
        total_missing += 1

total_altered = total_matched + total_extended + total_reduced

print("Total canonical genes: %d" % total_canonical_genes)
print("Total canonical UTRs: %d" % (total_altered + total_existing))
print("Total genes missing UTR: %d" % total_missing)
print("Total annotated UTRs: %d" % total_utrs)
print("Of those,")
print("... total new: %d" % total_new)
print("... total altered: %d" % total_altered)
print("... Of those,")
print("... ... total matched: %d" % total_matched)
print("... ... total extended: %d" % total_extended)
print("... ... total reduced: %d" % total_reduced)
print("... total pre-existing: %d" % total_existing)
assert total_utrs + total_missing == total_canonical_genes
assert total_new + total_altered + total_existing == total_utrs
assert total_altered + total_existing == total_utrs - total_new

# %%
