import logging
import os.path
import shutil
import subprocess

import sqlite3
import gffutils

from . import criteria
from .constants import LOG_DIR, TMP_GFF_FN
from .utils import cached, format_stats_line, format_stats_line


def write_summary_stats(annotations, total_peaks):
    with open('summary_stats.txt', 'w') as fstats:
        logging.info("Writing summary statistics file.")
        fstats.write(format_stats_line("Total peaks", total_peaks))
        fstats.write(format_stats_line("Total 3' UTRs annotated", len([a for a in annotations if "utr" in a])))
        fstats.write(format_stats_line("Peaks with no nearby features", total_peaks, annotations.no_features_counter))
        fstats.write(format_stats_line("Peaks corresponding to an already annotated 3' UTR", total_peaks,
                                       criteria.assert_whether_utr_already_annotated.fails.value))
        fstats.write(format_stats_line("Peaks contained within a feature", total_peaks,
                                       criteria.assert_not_a_subset.fails.value))
        fstats.write(format_stats_line("Peaks corresponding to 5'-end of a feature", total_peaks,
                                       criteria.assert_3_prime_end_and_truncate.fails.value))


def merge_and_gt_gff3_sort(db, annotations, new_gff_fn, args):
    """
    Concatenate three_prime_UTRs GFF and original GFF_IN file in a tmp file.
    Use genometools (gt) binary to sort and tidy tmp file into new combined output gff3 file.
    """
    logging.info("Merging annotations with canonical gff file.")
    
    db = sqlite3.connect(db, check_same_thread=False)
    db = gffutils.FeatureDB(db)
    for gene in db.all_features(featuretype="gene"):
        if gene.id not in annotations:
            features = {"gene": gene}
            features.update({"feature_{}".format(idx): f for idx, f in enumerate(list(db.children(gene)))})
            annotations[gene.id] = features

    log_fn = "gt_gff3.log"
    with open(cached(TMP_GFF_FN), 'w') as fout:
        fout.writelines(annotations)
    
    command = "gt gff3 -sort -retainids -tidy -o {} ".format(new_gff_fn)
    if args.force:
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
                logging.info("Successfully merged three_prime_UTRs into canonical gff: %s." % new_gff_fn)
                return
    logging.warning("Some issues were encountered when processing output file. Check %s." % log_fn)
    shutil.copy(cached(TMP_GFF_FN), new_gff_fn)
