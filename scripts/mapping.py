"""
example script for running a RNA-seq analysis

python rnaseq_pipeline.py rnaseq_pipeline.yaml

you will have to write a couple of functions to group the input
data in useful ways

"""
import sys
import yaml
from itertools import product, islice

from cluster_helper.cluster import cluster_view
from bcbio.utils import safe_makedir
from bipy.utils import (combine_pairs, append_stem)
from bipy.toolbox import (htseq_count, rseqc, sam)
from bipy.toolbox.trim import Cutadapt
from bipy.toolbox.fastqc import FastQC
from bipy.toolbox.fastq import HardClipper
from bipy.toolbox.tophat import Tophat
from bipy.toolbox.rseqc import RNASeqMetrics
from bipy.plugins import StageRepository

import os, fnmatch
from bcbio.log import create_base_logger, setup_local_logging, logger

def locate(pattern, root=os.curdir):
    '''Locate all files matching supplied filename pattern in and below
    supplied root directory.'''
    for path, dirs, files in os.walk(os.path.abspath(root)):
        for filename in fnmatch.filter(files, pattern):
            yield os.path.join(path, filename)


def make_test(in_file, config, lines=1000000):
    """
    take a small subset of the input files for testing. only makes sense for
    text files where lines gives an appopriate number of records, for example,
    FASTQ files should be a multiple of 4.

    """
    results_dir = config["dir"]["results"]
    out_dir = os.path.join(results_dir, "test", "data")
    safe_makedir(out_dir)
    out_file = os.path.join(out_dir,
                            append_stem(os.path.basename(in_file), "test"))
    with open(in_file) as in_handle, open(out_file, "w") as out_handle:
        for line in islice(in_handle, lines):
            out_handle.write(line)

    return out_file


def _get_stage_config(config, stage):
    return config["stage"][stage]


def _get_program(config, stage):
    return config["stage"][stage]["program"]


def _emit_stage_message(stage, curr_files):
    logger.info("Running %s on %s" % (stage, curr_files))


def main(config, view):

    # make the needed directories
    map(safe_makedir, config["dir"].values())

    # specific for project
    input_dir = config["dir"]["data"]
    logger.info("Loading files from %s" % (input_dir))
    input_files = list(locate("*.fq", input_dir))
    input_files += list(locate("*.fastq", input_dir))
    logger.info("Input files: %s" % (input_files))

    results_dir = config["dir"]["results"]
    safe_makedir(results_dir)

    # make the stage repository
    repository = StageRepository(config)
    logger.info("Stages found: %s" % (repository.plugins))

    if config.get("test_pipeline", False):
        logger.info("Running a test pipeline on a subset of the reads.")
        results_dir = os.path.join(results_dir, "test_pipeline")
        config["dir"]["results"] = results_dir
        safe_makedir(results_dir)
        curr_files = map(make_test, input_files, [config] * len(input_files))
        logger.info("Converted %s to %s. " % (input_files, curr_files))
    else:
        curr_files = input_files
        logger.info("Running RNASeq alignment pipeline on %s." % (curr_files))

    for stage in config["run"]:
        if stage == "fastqc":
            logger.info("Running fastqc on %s." % (curr_files))
            stage_runner = FastQC(config)
            view.map(stage_runner, curr_files)

        if stage == "cutadapt":
            curr_files = combine_pairs(curr_files)
            logger.info("Running cutadapt on %s." % (curr_files))
            stage_runner = Cutadapt(config)
            curr_files = view.map(stage_runner, curr_files)

        if stage == "tophat":
            logger.info("Running Tophat on %s." % (curr_files))
            #tophat = repository["tophat"](config)
            tophat = Tophat(config)
            tophat_outputs = view.map(tophat, curr_files)
            sortsam = view.map(sam.coordinate_sort_sam, tophat_outputs,
                               [config] * len(tophat_outputs))
            bamfiles = view.map(sam.sam2bam, sortsam)
            bamsort = view.map(sam.bamsort, bamfiles)
            view.map(sam.bamindex, bamsort)
            final_bamfiles = bamsort
            curr_files = tophat_outputs

        if stage == "disambiguate":
            logger.info("Disambiguating %s." % (curr_files))
            disambiguate = repository[stage](config)
            view.map(disambiguate, curr_files)

        if stage == "htseq-count":
            logger.info("Running htseq-count on %s." % (bamfiles))
            name_sorted = view.map(sam.bam_name_sort, bamfiles)
            curr_files = view.map(sam.bam2sam, name_sorted)
            htseq_args = zip(*product(curr_files, [config], [stage]))
            htseq_outputs = view.map(htseq_count.run_with_config,
                                     *htseq_args)
            htseq_count.combine_counts(htseq_outputs)

        if stage == "rnaseq_metrics":
            logger.info("Calculating RNASeq metrics on %s." % (curr_files))
            #coverage = repository[stage](config)
            coverage = RNASeqMetrics(config)
            view.map(coverage, curr_files)

        if stage == "hard_clip":
            logger.info("Trimming from the beginning of reads on %s." % (curr_files))
            hard_clipper = HardClipper(config)
            curr_files = view.map(hard_clipper, curr_files)

        if stage == "rseqc":
            logger.info("Running rseqc on %s." % (curr_files))
            curr_files = view.map(sam.sam2bam, curr_files)
            rseq_args = zip(*product(curr_files, [config]))
            view.map(rseqc.bam_stat, *rseq_args)
            view.map(rseqc.genebody_coverage, *rseq_args)
            view.map(rseqc.junction_annotation, *rseq_args)
            view.map(sam.bamindex, curr_files)
            RPKM_count_out = view.map(rseqc.RPKM_count, *rseq_args)
            view.map(rseqc.fix_RPKM_count_file, RPKM_count_out)
            """
                            annotate_args = zip(*product(RPKM_count_fixed,
                                         ["gene_id"],
                                         ["ensembl_gene_id"],
                                         ["human"]))
            view.map(annotate.annotate_table_with_biomart,
                     *annotate_args)
                     """
                     #view.map(rseqc.RPKM_saturation, *rseq_args)

if __name__ == "__main__":
    main_config_file = sys.argv[1]
    with open(main_config_file) as config_in_handle:
        startup_config = yaml.load(config_in_handle)
    parallel = create_base_logger(startup_config, {"type": "ipython"})
    setup_local_logging(startup_config, parallel)
    startup_config["parallel"] = parallel
    cluster_config = startup_config["cluster"]
    cores_per_job = cluster_config.get("cores_per_job", 1)
    if startup_config["cluster"].get("local", False):
        main(startup_config, DummyView())
    else:
        with cluster_view(cluster_config["scheduler"],
                          cluster_config["queue"],
                          cluster_config["cores"],
                          cores_per_job) as view:
            main(startup_config, view)

class DummyView(object):

    def __init__(self):
        self.map = map
