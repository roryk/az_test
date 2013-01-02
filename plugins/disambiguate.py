"""
wrapper around Qi's diambiguate script

"""

from bipy.pipeline.stages import AbstractStage
from bipy.utils import flatten
from bcbio.utils import safe_makedir
import os


class Disambiguate(AbstractStage):
    """
    wrapper class around Qi's disambiguate script. takes a tuple of SAM files,
    the results of mapping to the human (the first) and the mouse (the second)
    genomes and an output directory and disambiguates the reads mapping to
    each genome.

    requires disamb_byMapping2.pl to either be in your path or a program
    that does a similar task to be specified in this section of the YAML file:

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

    def __init__(self, config):
        self.config = config
        self.stage_config = config["stage"][self.stage]
        self.out_dir = os.path.join(config["dir"].get("results", "results"),
                                    self.stage)
        self.program = self.stage_config.get("program", "disamb_byMapping2.pl")
        safe_makedir(self.out_dir)

    def out_file(self, in_tuple):
        """
        returns the set of output filenames that will be made from
        running disambiguate on the tuple of input files
        """
        return list(flatten(map(self._organism_files,
                                in_tuple, ("Human", "Mouse"))))

    def _organism_files(self, in_file, organism):
        base, _ = os.path.splitext(os.path.basename(in_file))
        disamb = base + ".disambiguous" + organism + ".sam"
        ambig = base + ".ambiguous" + organism + ".sam"
        return [os.path.join(self.out_dir, x) for x in (disamb, ambig)]

    def __call__(self):
        pass
