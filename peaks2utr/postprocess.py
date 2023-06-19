import logging
import os.path
import shutil
import subprocess

import sqlite3

from . import criteria
from .constants import FeatureTypes, LOG_DIR, TMP_GFF_FN
from .models import FeatureDB
from .utils import cached, features_dict_for_gene, format_stats_line


def write_summary_stats(annotations, pipeline):
    total_peaks = pipeline.total_peaks
    with open('summary_stats.txt', 'w') as fstats:
        logging.info("Writing summary statistics file.")
        fstats.write(format_stats_line("Total peaks", total_peaks))
        fstats.write(format_stats_line("\t...with no nearby features", total_peaks, int(pipeline.no_features_counter)))
        fstats.write(format_stats_line("\t...corresponding to an already annotated 3' UTR", total_peaks,
                                       int(criteria.assert_whether_utr_already_annotated.fails)))
        fstats.write(format_stats_line("\t...contained within a feature", total_peaks,
                                       int(criteria.assert_peak_not_a_subset_of_transcript.fails)))
        fstats.write(format_stats_line("\t...corresponding to 5'-end of a feature", total_peaks,
                                       int(criteria.assert_3_prime_end_and_truncate.fails)))
        fstats.write(format_stats_line("\t...corresponding to potential 3' UTR removed due to zero read coverage",
                                       total_peaks, int(pipeline.zero_coverage_removal_counter)))
        fstats.write(format_stats_line("Total 3' UTRs", len(annotations.filter(featuretype=FeatureTypes.ThreePrimeUTR))))
        fstats.write(format_stats_line("\t...annotated by {}".format(__package__),
                                       len(annotations.filter(featuretype=FeatureTypes.ThreePrimeUTR, source=__package__))))


def merge_annotations(db, annotations):
    """
    Update three_prime_UTR annotations dict with all features from GFF_IN file.
    """
    logging.info("Merging annotations with canonical gff file.")

    db = sqlite3.connect(db, check_same_thread=False)
    db = FeatureDB(db)
    for gene in db.all_features(featuretype=FeatureTypes.Gene):
        if gene.id not in annotations:
            features = features_dict_for_gene(db, gene)
            annotations[gene.id] = features


def gt_gff3_sort(annotations, new_gff_fn, force=False, gtf_out=False):
    """
    Use genometools (gt) binary to sort and tidy tmp file into new combined output gff3 file.
    """
    log_fn = "gt_gff3.log"
    with open(cached(TMP_GFF_FN), 'w') as fout:
        fout.writelines(annotations.iter_feature_strings())
    if not gtf_out:
        command = "gt gff3 -sort -retainids -tidy -o {} ".format(new_gff_fn)
        if force:
            command += "-force "
        with open(os.path.join(LOG_DIR, log_fn), 'w') as flog:
            try:
                output = subprocess.check_output(
                    command + cached(TMP_GFF_FN),
                    universal_newlines=True,
                    stderr=subprocess.STDOUT,
                    shell=True
                )
            except subprocess.CalledProcessError as e:
                flog.write(e.output)
                if e.returncode == 127:
                    logging.warning("Genometools binary can't be called. Please ensure it is installed.")
            except MemoryError:
                logging.warning("Process required too much memory. Aborting.")
            else:
                flog.write(output)
                if os.path.exists(new_gff_fn):
                    logging.info("Successfully formatted GFF3 output file %s using genometools." % new_gff_fn)
                    return
        logging.warning("Some issues were encountered when processing output file. Check %s." % log_fn)
    shutil.copy(cached(TMP_GFF_FN), new_gff_fn)
