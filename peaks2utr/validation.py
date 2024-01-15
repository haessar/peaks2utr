import logging

import pysam

from .utils import connect_db, index_bam_file


def matching_chr(db_path, args):
    """
    Check seqids in BAM and GFF input files to ensure at least one matches. Returns bool.
    """
    db = connect_db(db_path)
    gff_chrs = {f.seqid for f in db.all_features()}
    bam_chrs = set()

    index_bam_file(args.BAM_IN, args.processors)
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
