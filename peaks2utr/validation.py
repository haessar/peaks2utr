import logging
import os
import os.path

from gffutils.inspect import inspect
import pysam

from .constants import LOG_DIR
from .exceptions import EXCEPTIONS_MAP
from .utils import index_bam_file

LOG_FN = "validation.log"


def valid_bam(args):
    """
    Check that BAM file is valid / has valid header. Returns bool.
    """
    try:
        pysam.AlignmentFile(args.BAM_IN, "rb")
    except ValueError as e:
        with open(os.path.join(LOG_DIR, LOG_FN), 'w') as flog:
            flog.write(str(e))
        raise EXCEPTIONS_MAP.get(valid_bam.__name__, Exception)("Check %s." % LOG_FN)
    return True


def matching_chr(args):
    """
    Check seqids in BAM and GFF input files to ensure at least one matches. Returns bool.
    """
    gff_inspection = inspect(args.GFF_IN, look_for=["chrom"], verbose=False)
    gff_chrs = gff_inspection.get('chrom', {})
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
