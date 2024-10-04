import unittest
from unittest.mock import MagicMock

from peaks2utr import prepare_argparser
from peaks2utr.collections import AnnotationsDict
from peaks2utr.constants import FeatureTypes, GFFUTILS_GTF_DIALECT
from peaks2utr.models import Feature, UTR
from peaks2utr.utils import get_output_filename


class TestOutputFormatting(unittest.TestCase):

    def setUp(self):
        chr = "chr1"
        start = 1000
        end = 2000
        strand = "+"
        gene_id = "gene1"
        ncRNA_gene_id = "ncRNA_gene1"
        argparser = prepare_argparser()
        self.args = argparser.parse_args(["", ""])
        self.gene_gff = Feature(chr, id=gene_id, featuretype=FeatureTypes.Gene[0], start=start, end=end, strand=strand,
                                attributes={"ID": [gene_id]})
        self.transcript_gff = Feature(chr, id="gene1:mRNA", featuretype=FeatureTypes.GffTranscript[0], start=start, end=end,
                                      strand=strand, attributes={"ID": ["gene1:mRNA"], "Parent": [gene_id]})
        self.ncRNA_gene_gff = Feature(chr, id=ncRNA_gene_id, featuretype=FeatureTypes.NonCodingGene[0], start=start, end=end, strand=strand,
                                attributes={"ID": [ncRNA_gene_id]})
        self.ncRNA_feature_gff = Feature(chr, id=ncRNA_gene_id+":t1", featuretype=FeatureTypes.NonCodingTranscript[0], start=start, end=end,
                                         strand=strand, attributes={"ID": [ncRNA_gene_id+":t1"], "Parent": [ncRNA_gene_id]})
        self.gene_gtf = Feature(chr, id=gene_id, featuretype=FeatureTypes.Gene[0], source="gffutils_derived",
                                start=start, end=end, strand=strand, attributes={"gene_id": [gene_id]},
                                dialect=GFFUTILS_GTF_DIALECT)
        self.transcript_gtf = Feature(chr, id="gene1.1", featuretype=FeatureTypes.GtfTranscript[0], start=start, end=end,
                                      strand=strand, attributes={"transcript_id": ["gene1.1"], "gene_id": [gene_id]},
                                      dialect=GFFUTILS_GTF_DIALECT)
        self.ncRNA_gene_gtf = Feature(chr, id=ncRNA_gene_id, featuretype=FeatureTypes.NonCodingGene[0], source="gffutils_derived",
                                start=start, end=end, strand=strand, attributes={"gene_id": [ncRNA_gene_id]},
                                dialect=GFFUTILS_GTF_DIALECT)
        self.ncRNA_feature_gtf = Feature(chr, id=ncRNA_gene_id+".t1", featuretype=FeatureTypes.NonCodingTranscript[0], start=start, end=end,
                                         strand=strand, attributes={"transcript_id": [ncRNA_gene_id+".t1"], "gene_id": [ncRNA_gene_id]},
                                         dialect=GFFUTILS_GTF_DIALECT)
        self.utr = UTR(start=start, end=end)
        self.db = MagicMock()
        self.db.children = MagicMock(return_value=[Feature(id="utr_1", featuretype=FeatureTypes.FivePrimeUTR[0])])

    def test_gff_to_gff(self):
        self.args.gtf_in = False
        self.args.gtf_out = False
        expected_gene = ["chr1", ".", "gene", "1000", "2000", ".", "+", ".", "ID=gene1"]
        expected_transcript = ["chr1", ".", "mRNA", "1000", "2000", ".", "+", ".", "ID=gene1:mRNA;Parent=gene1"]
        expected_utr = ["chr1", "peaks2utr", "three_prime_UTR", "1000", "2000", ".", "+", ".",
                        "ID=utr_2;Parent=gene1:mRNA;colour=3"]
        self.utr.generate_feature(self.gene_gff, self.transcript_gff, self.db, gtf_in=self.args.gtf_in)
        annotations = AnnotationsDict(args=self.args)
        annotations.update({
            self.gene_gff.id: {"gene": self.gene_gff, "transcript": self.transcript_gff, "utr": self.utr.feature}
        })
        gene, transcript, utr = annotations.iter_feature_strings()
        self.assertListEqual(gene.strip().split("\t"), expected_gene)
        self.assertListEqual(transcript.strip().split("\t"), expected_transcript)
        self.assertListEqual(utr.strip().split("\t"), expected_utr)

    def test_gff_to_gtf(self):
        self.args.gtf_in = False
        self.args.gtf_out = True
        expected_transcript = ["chr1", ".", "transcript", "1000", "2000", ".", "+", ".",
                               'gene_id "gene1"; transcript_id "gene1:mRNA";']
        expected_utr = ["chr1", "peaks2utr", "three_prime_UTR", "1000", "2000", ".", "+", ".",
                        'gene_id "gene1"; transcript_id "gene1:mRNA"; colour "3";']
        expected_exon = expected_utr[:2] + ['exon'] + expected_utr[3:]
        self.utr.generate_feature(self.gene_gff, self.transcript_gff, self.db, gtf_in=self.args.gtf_in)
        annotations = AnnotationsDict(args=self.args)
        annotations.update({
            self.gene_gff.id: {"gene": self.gene_gff, "transcript": self.transcript_gff, "utr": self.utr.feature}
        })
        transcript, utr, exon = annotations.iter_feature_strings()
        self.assertListEqual(transcript.strip().split("\t"), expected_transcript)
        self.assertListEqual(utr.strip().split("\t"), expected_utr)
        self.assertListEqual(exon.strip().split("\t"), expected_exon)

    def test_gtf_to_gff(self):
        self.args.gtf_in = True
        self.args.gtf_out = False
        expected_gene = ["chr1", "gffutils_derived", "gene", "1000", "2000", ".", "+", ".", "ID=gene1"]
        expected_transcript = ["chr1", ".", "mRNA", "1000", "2000", ".", "+", ".", "ID=gene1.1;Parent=gene1"]
        expected_utr = ["chr1", "peaks2utr", "three_prime_UTR", "1000", "2000", ".", "+", ".",
                        "ID=utr_2;Parent=gene1.1;colour=3"]
        self.utr.generate_feature(self.gene_gtf, self.transcript_gtf, self.db, gtf_in=self.args.gtf_in)
        annotations = AnnotationsDict(args=self.args)
        annotations.update({
            self.gene_gtf.id: {"gene": self.gene_gtf, "transcript": self.transcript_gtf, "utr": self.utr.feature}
        })
        gene, transcript, utr = annotations.iter_feature_strings()
        self.assertListEqual(gene.strip().split("\t"), expected_gene)
        self.assertListEqual(transcript.strip().split("\t"), expected_transcript)
        self.assertListEqual(utr.strip().split("\t"), expected_utr)

    def test_gtf_to_gtf(self):
        self.args.gtf_in = True
        self.args.gtf_out = True
        expected_transcript = ["chr1", ".", "transcript", "1000", "2000", ".", "+", ".",
                               'gene_id "gene1"; transcript_id "gene1.1";']
        expected_utr = ["chr1", "peaks2utr", "three_prime_UTR", "1000", "2000", ".", "+", ".",
                        'gene_id "gene1"; transcript_id "gene1.1"; colour "3";']
        expected_exon = expected_utr[:2] + ['exon'] + expected_utr[3:]
        self.utr.generate_feature(self.gene_gtf, self.transcript_gtf, self.db, gtf_in=self.args.gtf_in)
        annotations = AnnotationsDict(args=self.args)
        annotations.update({
            self.gene_gtf.id: {"gene": self.gene_gtf, "transcript": self.transcript_gtf, "utr": self.utr.feature}
        })
        transcript, utr, exon = annotations.iter_feature_strings()
        self.assertListEqual(transcript.strip().split("\t"), expected_transcript)
        self.assertListEqual(utr.strip().split("\t"), expected_utr)
        self.assertListEqual(exon.strip().split("\t"), expected_exon)
    
    def test_gtf_output_without_gtf_out_flag(self):
        self.args.output = "test_output.gtf"
        self.assertFalse(self.args.gtf_out)
        output_fn = get_output_filename(self.args)
        self.assertEqual(output_fn, self.args.output)
        self.assertTrue(self.args.gtf_out)
    
    def test_gff_output_with_gtf_out_flag(self):
        self.args.gtf_out = True
        self.args.output = "test_output.gff"
        with self.assertLogs(level='WARNING') as cm:
            output_fn = get_output_filename(self.args)
        self.assertRegex(cm.output[0], "WARNING")
        self.assertEqual(output_fn, self.args.output)
    
    def test_gtf_to_gff_ncRNA_retention(self):
        self.args.gtf_in = True
        self.args.gtf_out = False
        expected_gene = ["chr1", "gffutils_derived", "ncRNA_gene", "1000", "2000", ".", "+", ".", "ID=ncRNA_gene1"]
        expected_feature_0 = ["chr1", ".", "ncRNA", "1000", "2000", ".", "+", ".", "ID=ncRNA_gene1.t1;Parent=ncRNA_gene1"]
        annotations = AnnotationsDict(args=self.args)
        annotations.update({
            self.ncRNA_gene_gtf.id: {"gene": self.ncRNA_gene_gtf, "feature_0": self.ncRNA_feature_gtf}
        })
        gene, feature_0 = annotations.iter_feature_strings()
        self.assertListEqual(gene.strip().split("\t"), expected_gene)
        self.assertListEqual(feature_0.strip().split("\t"), expected_feature_0)
    
    def test_gff_to_gtf_ncRNA_retention(self):
        self.args.gtf_in = False
        self.args.gtf_out = True
        expected_gene = ["chr1", ".", "ncRNA_gene", "1000", "2000", ".", "+", ".",
                               'gene_id "ncRNA_gene1";']
        expected_feature_0 = ["chr1", ".", "ncRNA", "1000", "2000", ".", "+", ".",
                               'gene_id "ncRNA_gene1"; transcript_id "ncRNA_gene1:t1";']
        annotations = AnnotationsDict(args=self.args)
        annotations.update({
            self.ncRNA_gene_gff.id: {"gene": self.ncRNA_gene_gff, "feature_0": self.ncRNA_feature_gff}
        })
        gene, feature_0 = annotations.iter_feature_strings()
        self.assertListEqual(gene.strip().split("\t"), expected_gene)
        self.assertListEqual(feature_0.strip().split("\t"), expected_feature_0)


if __name__ == '__main__':
    unittest.main()
