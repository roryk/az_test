import yaml
import unittest
from bipy.plugins import StageRepository
from bipy.pipeline.stages import AbstractStage
from bcbio.utils import file_exists
from bipy.toolbox import sam
import os

STAGENAME = "disambiguate"


class TestDisambiguate(unittest.TestCase):

    def setUp(self):
        self.config_file = os.path.join("test", STAGENAME,
                                        "test_" + STAGENAME + ".yaml")
        with open(self.config_file) as in_handle:
            self.config = yaml.load(in_handle)

        self.input_file = self.config["input"]
        self.stage_config = self.config["stage"][STAGENAME]
        self.repository = StageRepository(self.config)

    def test_find_plugin(self):
        """
        test that the plugin finder can find the disambiguate plugin

        """
        stage_runner = self.repository[STAGENAME]
        self.assertTrue(issubclass(stage_runner, AbstractStage))

    def test_call(self):
        """
        test that running disambiguate works properly

        """
        stage = self.repository[STAGENAME]
        stage_runner = stage(self.config)
        out_files = stage_runner(self.input_file)
        self.assertTrue(all(map(self._is_valid_output, out_files)))
        map(os.remove, out_files)

    def _is_valid_output(self, in_file):
        return file_exists(in_file) and sam.is_sam(in_file)

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestDisambiguate)
    unittest.TextTestRunner(verbosity=2).run(suite)
