import yaml
import unittest
from bipy.plugins import StageRepository
from bipy.pipeline.stages import AbstractStage
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

    def test_out_file(self):
        """
        test that the disambiguate class makes the correct output
        filenames

        """
        stage = self.repository[STAGENAME]
        stage_runner = stage(self.config)
        print stage_runner.out_file(self.input_file)
        print stage_runner
        self.assertTrue(False)

    def test_call(self):
        """
        test that running disambiguate works properly

        """
        self.assertTrue(False)

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestDisambiguate)
    unittest.TextTestRunner(verbosity=2).run(suite)
