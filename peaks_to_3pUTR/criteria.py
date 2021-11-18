class CriteriaFailure(Exception):
    pass


def assert_not_already_annotated(db, peak, gene):
    if any(db.children(gene, featuretype='three_prime_UTR')):
        raise CriteriaFailure("3' UTR already annotated for gene %s near peak %s" % (gene.id, peak.name))


def assert_not_a_subset(peak, gene):
    if peak.range.issubset(gene.range):
        raise CriteriaFailure("Peak %s wholly contained within gene %s" % (peak.name, gene.id))


def assert_3_prime_end_and_truncate(peak, gene, utr):
    if peak.strand == "+" and peak.end > gene.end:
        utr.start = gene.end
    elif peak.strand == "-" and peak.start < gene.start:
        utr.end = gene.start
    else:
        raise CriteriaFailure("Peak %s corresponds to 5'-end of gene %s" % (peak.name, gene.id))


def truncate_5_prime_end(peak, next_gene, utr):
    if utr.range.intersection(next_gene.range):
        print("Peak %s overlapping following gene %s: Truncating." % (peak.name, next_gene.id))
        if peak.strand == "+":
            utr.end = next_gene.start
        else:
            utr.start = next_gene.end


def belongs_to_next_gene(peak, next_gene):
    if (peak.strand == "+" and peak.start > next_gene.start) or (peak.strand == "-" and peak.end < next_gene.end):
        raise CriteriaFailure("Peak %s belongs entirely to following gene %s." % (peak.name, next_gene.id))
