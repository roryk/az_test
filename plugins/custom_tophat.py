from bipy.toolbox.tophat import Tophat, run_with_config
from bipy.utils import is_pair


class TophatMouse(Tophat):

    stage = "tophat_mouse"

    def __init__(self, config):
        super(TophatMouse, self).__init__(config)
        self.stage_config = self.config["stage"][self.stage]
        self.ref = self.stage_config.get("ref", self.ref)
        self.gtf = self.stage_config.get("gtf", self.config["gtf"])

    def __call__(self, in_file):
        self._start_message(in_file)
        if is_pair(in_file):
            out_file = run_with_config(in_file[0], in_file[1],
                                       self.ref, self.stage, self.config,
                                       gtf=self.gtf)
        else:
            out_file = run_with_config(in_file[0], None, self.ref,
                                       self.stage, self.config, gtf=self.gtf)


class TophatHuman(TophatMouse):

    stage = "tophat_human"
