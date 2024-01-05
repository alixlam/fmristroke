from pathlib import Path

import nibabel as nb
import numpy as np
import pytest
from nipype.pipeline.engine.utils import generate_expanded_graph

from .... import config, data
from ...tests import mock_config
from ..base import init_roi_preproc_wf


@pytest.fixture(scope="module", autouse=True)
def _quiet_logger():
    import logging

    logger = logging.getLogger("nipype.workflow")
    old_level = logger.getEffectiveLevel()
    logger.setLevel(logging.ERROR)
    yield
    logger.setLevel(old_level)


@pytest.mark.parametrize(
    "fmriprep_root", [data.load("tests/FCStroke_2009_deriv").absolute()]
)
@pytest.mark.parametrize("freesurfer", [False, True])
@pytest.mark.parametrize("level", ["minimal", "resampling", "full"])
def test_roi_preproc_wf(
    fmriprep_root: Path,
    tmp_path: Path,
    freesurfer: bool,
    level: str,
):
    """Test as many combinations of precomputed files and input
    configurations as possible."""
    output_dir = tmp_path / "output"
    output_dir.mkdir()

    with mock_config():
        config.workflow.level = level
        config.workflow.freesurfer = freesurfer
        wf = init_roi_preproc_wf()

    flatgraph = wf._create_flat_graph()
    generate_expanded_graph(flatgraph)
