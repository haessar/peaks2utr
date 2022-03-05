import asyncio
from collections import defaultdict
import logging
import math
import os.path
import pickle

from asgiref.sync import sync_to_async
import gffutils
import pysam

from .exceptions import EXCEPTIONS_MAP
from .models import SoftClippedRead
from .utils import cached, consume_lines, filter_nested_dict
from .constants import CACHE_DIR, LOG_DIR, STRAND_PYSAM_ARGS


def pysam_index(bam_file):
    if not os.path.isfile(bam_file + '.bai'):
        logging.info("Indexing BAM file.")
        pysam.index(bam_file)


def pysam_strand_split(strand, bam_basename, args):
    """
    Call 'samtools view' to split BAM_IN for each strand.
    """
    output_file = cached(bam_basename + '.%s.bam' % strand)
    if not os.path.isfile(output_file) :
        logging.info("Splitting %s strand from %s." % (strand, args.BAM_IN))
        processors = max(math.floor(args.processors / 2), 1)        
        try:
            pysam.view("--threads", str(processors), "-b", *STRAND_PYSAM_ARGS[strand], "-o", output_file, args.BAM_IN, catch_stdout=False)            
        except TypeError as e:
            logging.error("pysam returned an error: %s" % e)
            raise
        else:
            logging.info("Finished splitting %s strand." % strand)
            pysam_index(output_file)                        
            samfile = pysam.AlignmentFile(output_file, "rb")
            logging.info("Iterating over %s strand reads to determine pileups." % strand)
            unmapped = defaultdict(lambda: defaultdict(int))
            for seg in samfile.fetch():
                read = SoftClippedRead(chr=seg.reference_name, start=seg.reference_start, end=seg.reference_end, cigar=seg.cigarstring, seq=seg.query_sequence, strand="reverse" if seg.is_reverse else "forward")
                if read.poly_tail_exists():
                    unmapped[read.chr][read.extremity] += 1
            with open(cached(strand + "_unmapped.pickle"), "wb") as f:
                pickle.dump(filter_nested_dict(unmapped, 10), f)


async def create_db(gff_in):
    """
    Asynchronously create sqlite3 db for GFF_IN.
    """
    gff_db = cached(os.path.basename(os.path.splitext(gff_in)[0] + '.db'))
    if not os.path.isfile(gff_db):
        logging.info('Creating gff db.')
        await sync_to_async(gffutils.create_db)(gff_in, gff_db, force=True, verbose=True, disable_infer_genes=True, disable_infer_transcripts=True)
        logging.info('Finished creating gff db.')
    return gff_db


async def call_peaks(bam_basename, strand):
    """
    Call MACS2 asynchronously for stranded BAM file.
    """
    if not os.path.isfile(cached("%s_peaks.broadPeak" % strand)):
        logging.info("Calling peaks for %s strand with MACS2." % strand)
        process = await asyncio.create_subprocess_exec(
            "macs2", "callpeak", "-t", cached(bam_basename + '.%s.bam' % strand), "-n", strand, "--nomodel", "--extsize", "200", "--broad", "--outdir", CACHE_DIR,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.STDOUT
        )
        asyncio.create_task(consume_lines(process.stdout, os.path.join(LOG_DIR, "%s_macs2.log" % strand)))
        exit_code = await process.wait()
        if exit_code != 0:
            logging.error("MACS2 returned an error.")
            raise EXCEPTIONS_MAP.get("call_peaks", Exception)("Check %s_macs2.log." % strand)
        logging.info("Finished calling %s strand peaks." % strand)
    else:
        logging.info("Using cached %s strand peaks file." % strand)
