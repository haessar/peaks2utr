import os
import os.path


class AnnotationColour:
    artemis_colour_map = {
        "red": "4",
        "green": "3",
        "blue": "6",
    }
    Extended = artemis_colour_map["green"]
    ExtendedWithSPAT = artemis_colour_map["blue"]
    TruncatedZeroCoverage = artemis_colour_map["red"]


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

GFFUTILS_GTF_DIALECT = {
    'leading semicolon': False,
    'trailing semicolon': True,
    'quoted GFF2 values': True,
    'field separator': '; ',
    'keyval separator': ' ',
    'multival separator': ',',                
    'fmt': 'gtf',
    'repeated keys': False,
    'order': ['gene_id', 'transcript_id', 'ID', 'Name'],
}


CACHE_DIR = os.path.join(os.getcwd(), '.cache')
LOG_DIR = os.path.join(os.getcwd(), '.log')

TMP_GFF_FN = "_tmp.gff"

PERC_ALLOCATED_VRAM = 75
