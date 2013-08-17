"""
example script for running a RNA-seq analysis

python rnaseq_pipeline.py rnaseq_pipeline.yaml

you will have to write a couple of functions to group the input
data in useful ways

"""
from bipy.cluster import stop_cluster
import sys
import yaml
from itertools import product
import os, fnmatch

from bcbio.utils import safe_makedir
from bipy.toolbox import (htseq_count, rseqc, sam)
from bcbio.log import logger, setup_local_logging, create_base_logger
from bipy.toolbox.rseqc import RNASeqMetrics
from bipy.plugins import StageRepository
from cluster_helper.cluster import cluster_view


def locate(pattern, root=os.curdir):
    '''Locate all files matching supplied filename pattern in and below
    supplied root directory.'''
    for path, dirs, files in os.walk(os.path.abspath(root)):
        for filename in fnmatch.filter(files, pattern):
            yield os.path.join(path, filename)

def main(config, view):
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
    logger.info("Running quantitation on %s." % (curr_files))

    for stage in config["run"]:
        if stage == "htseq-count":
            logger.info("Running htseq-count on %s." % (input_files))
            name_sorted = view.map(sam.bam_name_sort, input_files)
            curr_files = view.map(sam.bam2sam, name_sorted)
            htseq_args = zip(*product(curr_files, [config], [stage]))
            htseq_outputs = view.map(htseq_count.run_with_config,
                                     *htseq_args)
            htseq_count.combine_counts(htseq_outputs)

        if stage == "rnaseq_metrics":
            logger.info("Calculating RNASeq metrics on %s." % (curr_files))
            #coverage = repository[stage](config)
            #curr_files = view.map(sam.bam2sam, curr_files)
            coverage = RNASeqMetrics(config)
            view.map(coverage, curr_files)

        if stage == "rseqc":
            logger.info("Running rseqc on %s." % (curr_files))
            rseq_args = zip(*product(curr_files, [config]))
            view.map(rseqc.bam_stat, *rseq_args)
            view.map(rseqc.genebody_coverage, *rseq_args)
            view.map(rseqc.junction_annotation, *rseq_args)
            view.map(sam.bamindex, curr_files)
            RPKM_count_out = view.map(rseqc.RPKM_count, *rseq_args)
            view.map(rseqc.fix_RPKM_count_file, RPKM_count_out)
            #view.map(rseqc.junction_saturation, *rseq_args)
            #RPKM_args = zip(*product(final_bamfiles, [config]))
            #RPKM_count_out = view.map(rseqc.RPKM_count, *RPKM_args)
            #RPKM_count_fixed = view.map(rseqc.fix_RPKM_count_file,
            #                            RPKM_count_out)
            #RPKM_count_fixed = view.map(rseqc.fix_RPKM_count_file,
            """
                            annotate_args = zip(*product(RPKM_count_fixed,
                                         ["gene_id"],
                                         ["ensembl_gene_id"],
                                         ["human"]))
            view.map(annotate.annotate_table_with_biomart,
                     *annotate_args)
                     """
                     #view.map(rseqc.RPKM_saturation, *rseq_args)

    # end gracefully
    stop_cluster()


if __name__ == "__main__":
    # read in the config file and perform initial setup
    main_config_file = sys.argv[1]
    with open(main_config_file) as config_in_handle:
        startup_config = yaml.load(config_in_handle)
    parallel = create_base_logger(startup_config, {"type": "ipython"})
    setup_local_logging(startup_config, parallel)
    startup_config["parallel"] = parallel

    if startup_config["cluster"].get("local", False):
        main(startup_config, DummyView())
    else:
        cluster_config = startup_config["cluster"]
        with cluster_view(cluster_config["scheduler"],
                          cluster_config["queue"],
                          cluster_config["cores"]) as view:
            main(startup_config, view)

class DummyView(object):

    def __init__(self):
        self.map = map
