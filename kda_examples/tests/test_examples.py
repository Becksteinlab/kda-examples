"""
Check that the test model scripts are running without error.
"""

import os
import pytest

# TODO: cleanup imports
from kda_examples.test_model_3_state.script import main as tm_3_state
from kda_examples.test_model_4_state.script import main as tm_4_state
from kda_examples.test_model_4_state_leakage.script import main as tm_4_state_leakage
from kda_examples.test_model_5_state_leakage.script import main as tm_5_state_leakage
from kda_examples.test_model_6_state.script import main as tm_6_state
from kda_examples.test_model_6_state_leakage.script import main as tm_6_state_leakage


class TestModelTesting:

    def test_3_state(self, tmpdir):
        with tmpdir.as_cwd():
            os.mkdir("diagrams")
            tm_3_state()

    def test_4_state(self, tmpdir):
        with tmpdir.as_cwd():
            os.mkdir("diagrams")
            tm_4_state()

    def test_4_state_leakage(self, tmpdir):
        with tmpdir.as_cwd():
            os.mkdir("diagrams")
            tm_4_state_leakage()

    def test_5_state_leakage(self, tmpdir):
        with tmpdir.as_cwd():
            os.mkdir("diagrams")
            tm_5_state_leakage()

    def test_6_state(self, tmpdir):
        with tmpdir.as_cwd():
            os.mkdir("diagrams")
            tm_6_state()

    def test_6_state_leakage(self, tmpdir):
        with tmpdir.as_cwd():
            os.mkdir("diagrams")
            tm_6_state_leakage()
