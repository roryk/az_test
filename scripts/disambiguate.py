"""
example script for running a RNA-seq analysis

python rnaseq_pipeline.py rnaseq_pipeline.yaml

you will have to write a couple of functions to group the input
data in useful ways

"""
import sys
import yaml
#from bipy.log import setup_logging, logger
from bcbio.log import create_base_logger, logger, setup_local_logging
from bipy.utils import (combine_pairs, flatten, append_stem,
                        prepare_ref_file, replace_suffix)
from bipy.toolbox import (htseq_count, deseq, annotate, rseqc, sam)
from bipy.toolbox.trim import Cutadapt
from bipy.toolbox.fastqc import FastQC
from bipy.toolbox.tophat import Tophat
from bipy.toolbox.rseqc import RNASeqMetrics
from bcbio.utils import safe_makedir
from bcbio.distributed.ipythontasks import _setup_logging
from az.plugins.disambiguate import Disambiguate
from cluster_helper.cluster import cluster_view

from itertools import product,  islice
import sh
import os
import fnmatch

def locate(pattern, root=os.curdir):
    '''Locate all files matching supplied filename pattern in and below
    supplied root directory.'''
    for path, dirs, files in os.walk(os.path.abspath(root)):
        for filename in fnmatch.filter(files, pattern):
            yield os.path.join(path, filename)

def _emit_stage_message(stage, curr_files):
    logger.info("Running %s on %s" % (stage, curr_files))


def find_sam_files(in_dir):
    logger.info("Finding SAM files in %s." % (in_dir))
    # the descriptive filenames are symbolic links so find all of them
    files = list(locate("*.sam", in_dir))
    logger.info("Found %s." % (files))
    links = filter(os.path.islink, files)
    links.sort()
    logger.info("Found %s." % (links))
    return links


def find_bam_files(in_dir):
    logger.info("Finding BAM files in %s." % (in_dir))
    files = list(locate("*.sorted.bam", in_dir))
    logger.info("Found %s." % (files))
    return files


def main(config, view):

    # make the needed directories
    map(safe_makedir, config["dir"].values())

    # specific for project
    human_input = find_sam_files(config["input_dir_human"])
    mouse_input = find_sam_files(config["input_dir_mouse"])
    if len(human_input) != len(mouse_input):
        logger.error("The length of the number of human SAM files does "
                     "not match the length of the number of mouse SAM "
                     "files, aborting.")
        sys.exit(1)
    input_files = zip(human_input, mouse_input)

    curr_files = input_files

    logger.info("Running disambiguation pipeline on %s." % (curr_files))

    for stage in config["run"]:
        if stage == "disambiguate":
            logger.info("Disambiguating %s." % (curr_files))
            disambiguate = Disambiguate(config)
            out_files = list(flatten(view.map(disambiguate, curr_files)))
            bam_files = view.map(sam.sam2bam, out_files)
            bam_sorted = view.map(sam.bamsort, bam_files)
            view.map(sam.bamindex, bam_sorted)

if __name__ == "__main__":
    # read in the config file and perform initial setup
    main_config_file = sys.argv[1]
    with open(main_config_file) as config_in_handle:
        startup_config = yaml.load(config_in_handle)
    parallel = create_base_logger(startup_config, {"type": "ipython"})
    setup_local_logging(startup_config, parallel)
    startup_config["parallel"] = parallel
         #setup_logging(startup_config)

    cluster_config = startup_config["cluster"]
    cores_per_job = cluster_config.get("cores_per_job", 1)
    with cluster_view(cluster_config["scheduler"],
                      cluster_config["queue"],
                      cluster_config["cores"],
                      cores_per_job) as view:
        main(startup_config, view)
