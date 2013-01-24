"""
wrapper around Qi's diambiguate script

"""

from bipy.pipeline.stages import AbstractStage
from bipy.utils import flatten
from bcbio.utils import safe_makedir, file_exists
import os
import sh
from itertools import groupby, starmap
from bipy.log import logger
import shutil


class Disambiguate(AbstractStage):
    """
    wrapper class around Qi's disambiguate script. takes a tuple of SAM files,
    the results of mapping to the human (the first) and the mouse (the second)
    genomes and an output directory and disambiguates the reads mapping to
    each genome.

    requires disamb_byMapping2.pl or something that does a similar task to
    be specified in this section of the YAML file:

    stage:
        disambiguate:
            program: /path/to/disamb_byMapping2.pl

    example:
    stage_runner = Disambiguate(config)
    stage_runner(("human.sam", "mouse.sam") -> creates
    "human.disambiguousHuman.sam", "mouse.disambiguousMouse.sam",
    "human.ambiguousHuman.sam", "mouse.disambiguousMouse.sam"

    """

    stage = "disambiguate"
    organisms = ("Human", "Mouse")

    def __init__(self, config):
        # abstract class does some simple initialization for us
        super(Disambiguate, self).__init__(config)
        self.stage_config = config["stage"][self.stage]
        self.out_dir = os.path.join(config["dir"].get("results", "results"),
                                    self.stage)
        self.program = self.stage_config.get("program", "disamb_byMapping2.pl")
        safe_makedir(self.out_dir)

    def out_file(self, in_tuple):
        """
        returns the set of output filenames that will be made from
        running this stage on a set of input files
        """
        return map(self._disambiguate_to_bipy,
                   self._disambiguate_out(in_tuple))

    def _disambiguate_out(self, in_tuple):
        """
        returns the set of output filenames that will be made from
        running disambiguate on the tuple of input files
        """
        return list(flatten(map(self._organism_files,
                                in_tuple, self.organisms)))

    def _disambiguate_to_bipy(self, in_file):
        """
        convert the filenames into something easier to work with
        """
        base = os.path.basename(in_file)
        if "Human" in base:
            out_dir = os.path.join(self.config["dir"].get("results",
                                                          "results"),
                                   "human_mapping", self.stage)
            safe_makedir(out_dir)
            out_file = str.replace(base, "Human", ".human")
            return os.path.join(out_dir, out_file)
        elif "Mouse" in base:
            out_dir = os.path.join(self.config["dir"].get("results",
                                                          "results"),
                                   "mouse_mapping", self.stage)
            safe_makedir(out_dir)
            out_file = str.replace(base, "Mouse", ".mouse")
            return os.path.join(out_dir, out_file)
        else:
            logger.error("Could not find Human or Mouse in the filenames, "
                         "aborting.")
            exit(1)

    def _organism_files(self, in_file, organism):
        base, _ = os.path.splitext(os.path.basename(in_file))
        # remove organism name
        #sample_name, _ = os.path.splitext(base)
        sample_name = base
        disamb = sample_name + ".disambiguous" + organism + ".sam"
        ambig = sample_name + ".ambiguous" + organism + ".sam"
        return [os.path.join(self.out_dir, x) for x in
                (disamb, ambig)]

    def _disambiguate(self, org1_sam, org2_sam):
        run_disambiguate = sh.Command("perl")
        out_files = self.out_file((org1_sam, org2_sam))
        # if files exist already and are non-zero, skip this step
        if all(map(file_exists, out_files)):
            return out_files

        # disambiguate and return the output filenames
        run_disambiguate(self.program, org1_sam, org2_sam, self.out_dir)
        return out_files

    def __call__(self, in_files):
        self._start_message(in_files)
        # first is human, second is mouse
        out_files = self._disambiguate(in_files[0], in_files[1])
        dis_files = self._disambiguate_out(in_files)
        if all(map(file_exists, dis_files)):
            [shutil.move(x[0], x[1]) for x in zip(dis_files, out_files)]
        shutil.rmtree(self.out_dir)
        self._end_message(in_files)
        return out_files
