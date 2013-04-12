"""
example script for running a RNA-seq analysis

python rnaseq_pipeline.py rnaseq_pipeline.yaml

you will have to write a couple of functions to group the input
data in useful ways

"""
from bipy.cluster import start_cluster, stop_cluster
import sys
import yaml
from bipy.log import setup_logging, logger
from bcbio.utils import safe_makedir, file_exists
from bipy.utils import (combine_pairs, flatten, append_stem,
                        prepare_ref_file, replace_suffix)
from bipy.toolbox import (htseq_count, deseq, annotate, rseqc, sam)
from bcbio.broad import BroadRunner, picardrun
from bipy.toolbox.rseqc import RNASeqMetrics
from bipy.plugins import StageRepository
from cluster_helper.cluster import cluster_view

import glob
from itertools import product, repeat, islice
import sh
import os, fnmatch

def locate(pattern, root=os.curdir):
    '''Locate all files matching supplied filename pattern in and below
    supplied root directory.'''
    for path, dirs, files in os.walk(os.path.abspath(root)):
        for filename in fnmatch.filter(files, pattern):
            yield os.path.join(path, filename)

def main(config_file, view):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)

    # make the needed directories
    map(safe_makedir, config["dir"].values())

    # specific for project
    input_dir = config["input_dir"]
    logger.info("Loading files from %s" % (input_dir))
    input_files = list(locate("*.disambiguous*.sorted.bam", input_dir))
    logger.info("Input files: %s" % (input_files))

    results_dir = config["dir"]["results"]
    safe_makedir(results_dir)

    # make the stage repository
    repository = StageRepository(config)
    logger.info("Stages found: %s" % (repository.plugins))

    curr_files = input_files
    final_bamfiles = input_files
    view.map(sam.bamindex, final_bamfiles)
    logger.info("Running quantitation on %s." % (curr_files))

    for stage in config["run"]:

        if stage == "htseq-count":
            logger.info("Running htseq-count on %s." % (curr_files))
            htseq_args = zip(*product(curr_files, [config], [stage]))
            htseq_outputs = view.map(htseq_count.run_with_config,
                                     *htseq_args)
            htseq_count.combine_counts(htseq_outputs)

        if stage == "rnaseq_metrics":
            logger.info("Calculating RNASeq metrics on %s." % (curr_files))
            #coverage = repository[stage](config)
            curr_files = view.map(sam.bam2sam, curr_files)
            coverage = RNASeqMetrics(config)
            view.map(coverage, curr_files)

        if stage == "rseqc":
            logger.info("Running rseqc on %s." % (curr_files))
            rseq_args = zip(*product(curr_files, [config]))
            view.map(rseqc.bam_stat, *rseq_args)
            view.map(rseqc.genebody_coverage, *rseq_args)
            view.map(rseqc.junction_annotation, *rseq_args)
            view.map(rseqc.junction_saturation, *rseq_args)
            RPKM_args = zip(*product(final_bamfiles, [config]))
            RPKM_count_out = view.map(rseqc.RPKM_count, *RPKM_args)
            RPKM_count_fixed = view.map(rseqc.fix_RPKM_count_file,
                                        RPKM_count_out)
            """
                            annotate_args = zip(*product(RPKM_count_fixed,
                                         ["gene_id"],
                                         ["ensembl_gene_id"],
                                         ["human"]))
            view.map(annotate.annotate_table_with_biomart,
                     *annotate_args)
                     """
            view.map(rseqc.RPKM_saturation, *rseq_args)

    # end gracefully
    stop_cluster()


if __name__ == "__main__":
    # read in the config file and perform initial setup
    main_config_file = sys.argv[1]
    with open(main_config_file) as config_in_handle:
        startup_config = yaml.load(config_in_handle)
    setup_logging(startup_config)
    cluster_config = startup_config["cluster"]
    with cluster_view(cluster_config["scheduler"],
                      cluster_config["queue"],
                      cluster_config["cores"]) as view:
        main(main_config_file, view)
