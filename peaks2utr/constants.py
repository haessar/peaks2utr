import os
import os.path


STRAND_MAP = {
    'forward': '+',
    'reverse': '-',
}

STRAND_CIGAR_SOFT_CLIP_REGEX = {
    "forward": r"([0-9]+)S$",
    "reverse": r"^([0-9]+)S"
}

STRAND_PYSAM_ARGS = {
    'forward': ["-F", "20"],
    'reverse': ["-f", "16"],
}

CACHE_DIR = os.path.join(os.getcwd(), '.cache')
LOG_DIR = os.path.join(os.getcwd(), '.log')

TMP_GFF_FN = "_tmp.gff"

PERC_ALLOCATED_VRAM= 75
