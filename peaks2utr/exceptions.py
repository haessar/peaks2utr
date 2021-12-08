class PysamError(Exception):
    pass


class MACS2Error(Exception):
    pass


EXCEPTIONS_MAP = {
    "pysam_strand_split": PysamError,
    "call_peaks": MACS2Error,
}
