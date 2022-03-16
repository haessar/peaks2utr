import logging

import gffutils
import pysam


def matching_chr(db, args):
    db = gffutils.FeatureDB(db)
    gff_chrs = {f.chrom for f in db.all_features()}
    bam_chrs = set()

    samfile = pysam.AlignmentFile(args.BAM_IN, "rb", require_index=True)

    for chr in gff_chrs:
        try:
            samfile.fetch(chr)
        except ValueError:
            logging.warning("Chromosome {} from GFF_IN not found in BAM_IN.".format(chr))
        else:
            bam_chrs.add(str(chr))
    if len(bam_chrs) > 1:
        logging.warning(
            """
            Chromosomes {} are present in both GFF_IN and BAM_IN.
            Consider reducing both to a single chromosome to improve performance.
            """.format(', '.join(bam_chrs))
        )
    return bool(bam_chrs)
