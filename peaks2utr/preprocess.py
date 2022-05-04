import asyncio
from collections import defaultdict
from glob import glob
import json
import logging
import os.path
import re

from asgiref.sync import sync_to_async
import gffutils
import pysam
from tqdm import tqdm

from .exceptions import EXCEPTIONS_MAP
from .models import SoftClippedRead
from .utils import cached, consume_lines, filter_nested_dict, sum_nested_dicts, multiprocess_over_dict
from .constants import CACHE_DIR, LOG_DIR, STRAND_PYSAM_ARGS


class BAMSplitter:
    def __init__(self, bam_basename, args):
        self.basename = bam_basename
        self.args = args
        self.pbar = None
    
    def split_strands(self):
        for strand in ["forward", "reverse"]:
            output_file = cached(self.basename + '.%s.bam' % strand)
            if not os.path.isfile(output_file):
                logging.info("Splitting %s strand from %s." % (strand, self.args.BAM_IN))
                try:
                    pysam.view("--threads", str(self.args.processors), "-b", *STRAND_PYSAM_ARGS[strand], "-o", output_file, self.args.BAM_IN, catch_stdout=False)            
                except TypeError as e:
                    logging.error("pysam returned an error: %s" % e)
                    raise
                else:
                    logging.info("Finished splitting %s strand." % strand)
            else:
                logging.info("Using cached %s strand BAM file." % strand)
    
    @staticmethod
    def num_read_groups(bam):
        header = pysam.view("-H", bam).split("\n")
        return len([h for h in header if h.startswith("@RG")])
    
    def split_read_groups(self):
        for strand in ["forward", "reverse"]:
            input_bam = cached(self.basename + '.%s.bam' % strand)
            if len(glob(cached(self.basename + ".%s_*.bam" % strand))) < self.num_read_groups(input_bam):
                logging.info("Splitting %s-stranded BAM file into read-groups." % strand)
                pysam.split("-@", str(self.args.processors), "-f", cached("%*_%#.%."), input_bam)
            
        self.read_group_bams = sorted(glob(cached(self.basename + ".forward_*.bam")) + glob(cached(self.basename + ".reverse_*.bam")),
                                      key = lambda x: os.stat(x).st_size,
                                      reverse=True)
        self.outputs = {bf: cached(re.search( r'%s.(.*).bam$' % self.basename, os.path.basename(bf)).group(1) + "_unmapped.json") for bf in self.read_group_bams}
        self.full_outputs = self.outputs.copy()
    
    def index_bam_file(self, bam_file):
        if not os.path.isfile(cached(bam_file + '.bai')):
            logging.info("Indexing %s." % bam_file)
            pysam.index("-@", str(self.args.processors), bam_file)
    
    def _get_max_reads_for_pbar(self):
        max_reads = 0
        for bf in self.read_group_bams:
            if not os.path.isfile(self.outputs[bf]):
                self.index_bam_file(bf)
                idxstats = pysam.idxstats(bf).split('\n')
                num_reads = sum([int(chr.split("\t")[2]) + int(chr.split("\t")[3]) for chr in idxstats[:-1]])
                if num_reads > max_reads:
                    max_reads = num_reads
                    self.max_bam = bf
            else:
                del self.outputs[bf]
        return max_reads        
    
    def pileup_soft_clipped_reads(self):
        if not os.path.isfile(cached("forward_unmapped.json")) or not os.path.isfile(cached("reverse_unmapped.json")):        
            max_reads = self._get_max_reads_for_pbar()
            if self.outputs and max_reads > 0:
                with tqdm(total=max_reads,
                        desc=f'{"INFO": <8} Iterating over reads to determine SPAT pileups',
                        bar_format='{l_bar}{bar}| [{elapsed}<{remaining}]') as self.pbar:
                    multiprocess_over_dict(self._count_unmapped_pileups, self.outputs)
            
            logging.info('Merging SPAT outputs.')
            for strand in ["forward", "reverse"]:
                strand_output = {}
                for output in self.full_outputs.values():
                    if strand in os.path.basename(output):
                        with open(output, 'r') as f:
                            strand_output = sum_nested_dicts(strand_output, json.load(f))
                with open(cached("%s_unmapped.json" % strand), "w") as f:
                    json.dump(filter_nested_dict(strand_output, self.args.min_pileups), f)
        else:
            logging.info("Using cached SPAT pileups.")
    
    def _count_unmapped_pileups(self, bam_file, output_file):
        samfile = pysam.AlignmentFile(bam_file, "rb")            
        unmapped = defaultdict(lambda: defaultdict(int))
        for seg in samfile.fetch(until_eof=True):
            read = SoftClippedRead(chr=seg.reference_name, start=seg.reference_start, end=seg.reference_end, cigar=seg.cigarstring, seq=seg.query_sequence, strand="reverse" if seg.is_reverse else "forward")
            if read.poly_tail_exists(self.args.min_poly_tail):
                unmapped[read.chr][read.extremity] += 1
            if bam_file == self.max_bam:
                self.pbar.update()
                    
        with open(output_file, "w") as f:
            json.dump(unmapped, f)


async def create_db(gff_in):
    """
    Asynchronously create sqlite3 db for GFF_IN.
    """
    gff_db = cached(os.path.basename(os.path.splitext(gff_in)[0] + '.db'))
    if not os.path.isfile(gff_db):
        logging.info('Creating gff db.')
        await sync_to_async(gffutils.create_db)(gff_in, gff_db, force=True, verbose=True, disable_infer_genes=True, disable_infer_transcripts=True)
        logging.info('Finished creating gff db.')
    else:        
        logging.info("Using cached gff db.")
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
