from pathlib import Path

import nibabel as nb
import numpy as np
import pytest
from nipype.pipeline.engine.utils import generate_expanded_graph

from .... import config
from ...tests import mock_config
from ..base import init_lesion_preproc_wf

from .... import data


@pytest.fixture(scope="module", autouse=True)
def _quiet_logger():
    import logging

    logger = logging.getLogger("nipype.workflow")
    old_level = logger.getEffectiveLevel()
    logger.setLevel(logging.ERROR)
    yield
    logger.setLevel(old_level)

@pytest.mark.parametrize("fmriprep_root", [data.load('tests/FCStroke_2009_deriv').absolute()])
@pytest.mark.parametrize("session_level", [False, True])
@pytest.mark.parametrize("freesurfer", [False, True])
@pytest.mark.parametrize("level", ["minimal", "resampling", "full"])
def test_lesion_preproc_wf(
    fmriprep_root: Path,
    tmp_path: Path,
    freesurfer: bool,
    level: str,
    session_level: bool | None,
):
    """Test as many combinations of precomputed files and input
    configurations as possible."""
    output_dir = tmp_path / 'output'
    output_dir.mkdir()

    if session_level:
        bold_series = [
            str(fmriprep_root / 'sub-027' / 'ses-1week' / 'func'  / f'sub-027_ses-1week_task-rest_run-0{i}_desc-preproc_bold.nii.gz')
            for i in range(2, 4)
            ]
        boldref_file = [
            str(fmriprep_root / 'sub-027' / 'ses-1week' / 'func'  / f'sub-027_ses-1week_task-rest_run-0{i}_boldref.nii.gz')
            for i in range(2, 4)
        ]
        boldmask_file = [
            str(fmriprep_root / 'sub-027' / 'ses-1week' / 'func'  / f'sub-027_ses-1week_task-rest_run-0{i}_desc-brain_mask.nii.gz')
            for i in range(2, 4)
        ]
        confounds_file = [
            str(fmriprep_root / 'sub-027' / 'ses-1week' / 'func'  / f'sub-027_ses-1week_task-rest_run-0{i}_desc-confounds_timeseries.tsv')
            for i in range(2, 4)
        ]
        confounds_metadata = [
            str(fmriprep_root / 'sub-027' / 'ses-1week' / 'func'  / f'sub-027_ses-1week_task-rest_run-0{i}_desc-confounds_timeseries.json')
            for i in range(2, 4)
        ]
    
    else:
        bold_series = [
            str(fmriprep_root / 'sub-027' / 'ses-1week' / 'func'  / f'sub-027_ses-1week_task-rest_run-0{i}_desc-preproc_bold.nii.gz')
            for i in range(2, 3)
            ]
        boldref_file = [
            str(fmriprep_root / 'sub-027' / 'ses-1week' / 'func'  / f'sub-027_ses-1week_task-rest_run-0{i}_boldref.nii.gz')
            for i in range(2, 3)
        ]
        boldmask_file = [
            str(fmriprep_root / 'sub-027' / 'ses-1week' / 'func'  / f'sub-027_ses-1week_task-rest_run-0{i}_desc-brain_mask.nii.gz')
            for i in range(2, 3)
        ]
        confounds_file = [
            str(fmriprep_root / 'sub-027' / 'ses-1week' / 'func'  / f'sub-027_ses-1week_task-rest_run-0{i}_desc-confounds_timeseries.tsv')
            for i in range(2, 3)
        ]
        confounds_metadata = [
            str(fmriprep_root / 'sub-027' / 'ses-1week' / 'func'  / f'sub-027_ses-1week_task-rest_run-0{i}_desc-confounds_timeseries.json')
            for i in range(2, 3)
        ]
    

    with mock_config():
        config.workflow.level = level
        config.workflow.freesurfer = freesurfer
        config.execution.run_sessionlevel = session_level

        wf = init_lesion_preproc_wf(
            bold_file=bold_series,
            boldref_file=boldref_file,
            boldmask_file=boldmask_file,
            confounds_file=confounds_file,
            confounds_metadata=confounds_metadata
        )

    flatgraph = wf._create_flat_graph()
    generate_expanded_graph(flatgraph)