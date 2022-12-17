import asyncio
import os
import os.path


def prepare_argparser():
    import argparse
    import pkg_resources

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=r"""
        ____________________________________________________________________
                                                  __
                                   /            /    )
        ------__-----__-----__----/-__----__-----___/------------_/_----)__-
            /   )  /___)  /   )  /(      (_ `  /         /   /   /     /   )
        ___/___/__(___ __(___(__/___\___(__)__/____/____(___(___(_ ___/_____
          /
         /

        Use MACS3 to build forward and reverse peaks files for given .bam
        file.
        Iterate peaks through set of criteria to determine UTR viability,
        before annotating in .gff file.
        """
    )
    parser.add_argument('GFF_IN', help="input 'canonical' annotations file in gff or gtf format.")
    parser.add_argument('BAM_IN', help="input reads file in bam format.")
    parser.add_argument('--max-distance', type=int, default=200,
                        help='maximum distance in bases that UTR can be from a transcript.')
    parser.add_argument('--override-utr', action="store_true", help="ignore already annotated 3' UTRs in criteria.")
    parser.add_argument('--extend-utr', action="store_true",
                        help="extend previously existing 3' UTR annotations where possible.")
    parser.add_argument('--five-prime-ext', type=int, default=0,
                        help='a peak within this many bases of a gene\'s 5\'-end should be assumed to belong to it.')
    parser.add_argument('--skip-soft-clip', action="store_true",
                        help="skip the resource-intensive logic to pileup soft-clipped read edges.")
    parser.add_argument('--min-pileups', type=int, default=10, help='Minimum number of piled-up mapped reads for UTR cut-off.')
    parser.add_argument('--min-poly-tail', type=int, default=10,
                        help='Minimum length of poly-A/T tail considered in soft-clipped reads.')
    parser.add_argument('-p', '--processors', type=int, default=1, help="How many processor cores to use.")
    parser.add_argument('-f', '-force', '--force', action="store_true", help="Overwrite outputs if they exist.")
    parser.add_argument('-o', '--output', help="output filename.")
    parser.add_argument('--gtf', action="store_true", help="output in GTF format (rather than default GFF3).")
    parser.add_argument('--keep-cache', action="store_true", help="Keep cached files on run completion.")
    parser.add_argument('--version', action='version',
                        version='%(prog)s {version}'.format(version=pkg_resources.require(__package__)[0].version))
    return parser


def demo():
    """
    Entry-point for peaks2utr-demo
    """
    from glob import glob

    demo_dir = os.path.join(os.path.dirname(__file__), "demo")
    gff_in = glob(os.path.join(demo_dir, "*.gff"))[0]
    bam_in = glob(os.path.join(demo_dir, "*.bam"))[0]
    argparser = prepare_argparser()
    args = argparser.parse_args([gff_in, bam_in, "-f"])
    asyncio.run(_main(args))


def main():
    """
    Main entry-point
    """
    from psutil import virtual_memory

    from .constants import PERC_ALLOCATED_VRAM
    from .utils import limit_memory

    limit_memory(PERC_ALLOCATED_VRAM * virtual_memory().total / 100)
    argparser = prepare_argparser()
    args = argparser.parse_args()
    asyncio.run(_main(args))


async def _main(args):
    """
    The main function / pipeline for peaks2utr.
    """
    import logging
    import math
    import multiprocessing
    import shutil
    import sys

    from tqdm import tqdm

    from . import constants
    from .annotations import Annotations, NoNearbyFeatures, batch_annotate_strand
    from .collections import BroadPeaksList
    from .utils import cached, iter_batches, yield_from_process
    from .preprocess import BAMSplitter, call_peaks, create_db
    from .postprocess import merge_annotations, gt_gff3_sort, write_summary_stats

    try:
        ###################
        # Setup logging   #
        ###################

        # Change root logger level from WARNING (default) to NOTSET in order for all messages to be delegated.
        logging.getLogger().setLevel(logging.NOTSET)

        # Add stdout handler, with level INFO.
        console = logging.StreamHandler(sys.stdout)
        console.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        console.setFormatter(formatter)
        logging.getLogger().addHandler(console)

        if not os.path.exists(constants.LOG_DIR):
            logging.info("Make .log directory")
            os.mkdir(constants.LOG_DIR)

        # Add file handler, with level DEBUG.
        fileHandler = logging.FileHandler(
            filename=os.path.join(constants.LOG_DIR, '{}_debug.log'.format(__package__)), mode="w")
        fileHandler.setLevel(logging.DEBUG)
        fileHandler.setFormatter(formatter)
        logging.getLogger().addHandler(fileHandler)

        ###################
        # Define outputs  #
        ###################

        if not os.path.exists(constants.CACHE_DIR):
            logging.info("Make .cache directory")
            os.mkdir(constants.CACHE_DIR)

        gff_basename = os.path.basename(os.path.splitext(args.GFF_IN)[0])
        new_gff_fn = gff_basename + ".new.gff" if not args.output else args.output
        bam_basename = os.path.basename(os.path.splitext(args.BAM_IN)[0])
        if not args.output:
            new_gff_fn = gff_basename + ".new" + ".gtf" if args.gtf else ".gff3"
        else:
            new_gff_fn = args.output
        ###################
        # Perform checks  #
        ###################

        if os.path.exists(new_gff_fn) and not args.force:
            logging.error("%s already exists. Re-run with -f flag to force overwrite of output files. Aborting." % new_gff_fn)
            sys.exit(1)

        if args.override_utr and args.extend_utr:
            logging.error("Only one of --extend-utr and --override-utr can be used simultaneously. Aborting.")
            sys.exit(1)

        ###################
        # Pre-processing  #
        ###################

        BAMSplitter(bam_basename, args).process()

        db, _, _ = await asyncio.gather(
            create_db(args.GFF_IN),
            call_peaks(bam_basename, "forward"),
            call_peaks(bam_basename, "reverse")
        )
        peaks = \
            BroadPeaksList(broadpeak_fn=cached("forward_peaks.broadPeak"), strand="forward") + \
            BroadPeaksList(broadpeak_fn=cached("reverse_peaks.broadPeak"), strand="reverse")
        total_peaks = len(peaks)

        ###################
        # Process peaks   #
        ###################

        queue = multiprocessing.Queue()
        processes = [
            batch_annotate_strand(db, batch, queue, args)
            for batch in iter_batches(peaks, math.ceil(total_peaks/args.processors))
        ]
        for p in processes:
            p.start()

        annotations = Annotations()
        with tqdm(total=total_peaks, desc=f'{"INFO": <8} Iterating over peaks to annotate 3\' UTRs.') as pbar:
            for p in processes:
                for result in yield_from_process(queue, p, pbar):
                    if result:
                        if result is NoNearbyFeatures:
                            annotations.no_features_counter += 1
                        else:
                            annotations.update(result)

        ###################
        # Post-processing #
        ###################

        merge_annotations(db, annotations)
        gt_gff3_sort(annotations, new_gff_fn, args.force, args.gtf)
        write_summary_stats(annotations, total_peaks)

        logging.info("%s finished successfully." % __package__)
        await asyncio.sleep(1)
        sys.exit(0)
    except KeyboardInterrupt:
        logging.error("User interrupted processing. Aborting.")
        sys.exit(130)
    finally:
        try:
            if not args.keep_cache:
                logging.info("Clearing cache.")
                shutil.rmtree(constants.CACHE_DIR)
        except NameError:
            pass
