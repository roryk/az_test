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
from bipy.utils import (combine_pairs, flatten, append_stem,
                        prepare_ref_file, replace_suffix)
from bipy.toolbox import (htseq_count, deseq, annotate, rseqc, sam)
from bipy.toolbox.trim import Cutadapt
from bipy.toolbox.fastqc import FastQC
from bipy.toolbox.tophat import Tophat
from bipy.toolbox.rseqc import RNASeqMetrics
from bipy.toolbox.sam import Disambiguate
from bipy.plugins import StageRepository
from bcbio.utils import safe_makedir

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
    logger.info("Found %s." % (links))
    return links

def find_bam_files(in_dir):
    logger.info("Finding BAM files in %s." % (in_dir))
    files = list(locate("*.sorted.bam", in_dir))
    logger.info("Found %s." % (files))
    return files


def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)

    # make the needed directories
    map(safe_makedir, config["dir"].values())

    # specific for project
    #human_input = find_sam_files(config["input_dir_human"])
    #mouse_input = find_sam_files(config["input_dir_mouse"])
    human_input = find_bam_files(config["input_dir_human"])
    mouse_input = find_bam_files(config["input_dir_mouse"])
    input_files = zip(human_input, mouse_input)

    # make the stage repository
    repository = StageRepository(config)
    logger.info("Stages found: %s" % (repository.plugins))

    curr_files = input_files

    logger.info("Running disambiguation pipeline on %s." % (curr_files))

    for stage in config["run"]:
        # for now use my hack version
        if stage == "disambiguate":
            logger.info("Disambiguating %s." % (curr_files))
            disambiguate = Disambiguate(config)
            view.map(disambiguate, curr_files)
            #disambiguate = repository[stage](config)
            #view.map(disambiguate, curr_files)

    # end gracefully
    stop_cluster()


if __name__ == "__main__":
    # read in the config file and perform initial setup
    main_config_file = sys.argv[1]
    with open(main_config_file) as config_in_handle:
        startup_config = yaml.load(config_in_handle)
    setup_logging(startup_config)
    start_cluster(startup_config)
    from bipy.cluster import view
    main(main_config_file)
