import collections
import json


class ZeroCoverageIntervalsDict(collections.UserDict):
    """
    Dictionary of zero coverage intervals per chromosome from parsed BED file.
    """
    class Interval:
        def __init__(self, start, end):
            self.start = int(start)
            self.end = int(end)

    def __init__(self, dict=None, bed_fn=None):
        super().__init__(dict)
        if bed_fn:
            with open(bed_fn, 'r') as f:
                for line in f.readlines():
                    chr, start, end = line.strip().split('\t')
                    if chr not in self.data:
                        self.data[chr] = []
                    self.data[chr].append(self.Interval(start, end))

    def filter(self, chr, base):
        """
        Filter intervals that contain base.

        This was motivated by https://github.com/haessar/peaks2utr/issues/9 - the native
        pybedtools.BedTool.filter method was causing "Too many files open" errors in some
        distributed systems.
        """
        if chr in self:
            return [i for i in self[chr] if i.start < base < i.end]
        return []


class SPATTruncationPointsDict(collections.UserDict):
    """
    Dictionary of SPAT "truncation points" per chromosome from json file.
    """
    def __init__(self, dict=None, json_fn=None):
        super().__init__(dict)
        if json_fn:
            with open(json_fn, 'r') as f:
                self.data.update(json.load(f) or {})


