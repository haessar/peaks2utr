class AnnotationsError(Exception):
    pass


class PybedtoolsError(Exception):
    pass


class PysamError(Exception):
    pass


class MACS3Error(Exception):
    pass


EXCEPTIONS_MAP = {
    "_find_zero_coverage_intervals": PybedtoolsError,
    "_count_unmapped_pileups": PysamError,
    "call_peaks": MACS3Error,
}
