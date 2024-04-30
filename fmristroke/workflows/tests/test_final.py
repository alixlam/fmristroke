from copy import deepcopy
from pathlib import Path
from unittest.mock import patch

import bids
import nibabel as nb
import numpy as np
import pytest
from nipype.pipeline.engine.utils import generate_expanded_graph
from niworkflows.utils.testing import generate_bids_skeleton

from ... import config
from ..final import init_fmristroke_wf
from ..tests import mock_config


@pytest.fixture(scope="module", autouse=True)
def _quiet_logger():
    import logging

    logger = logging.getLogger("nipype.workflow")
    old_level = logger.getEffectiveLevel()
    logger.setLevel(logging.ERROR)
    yield
    logger.setLevel(old_level)


def _make_params(
    session_level = None,
    ncomp_method: str = "varexp",
    ica_method: str = "canica",
    maxlag: int = 10,
    croprun: int = 0,
    conn_measure: list[str] = ["correlation"],
    freesurfer: bool = True,
    bids_filters: dict = {},
):
    if bids_filters is None:
        bids_filters = {}
    return (
        session_level,
        ncomp_method,
        ica_method,
        maxlag,
        croprun,
        conn_measure,
        freesurfer,
        bids_filters,
    )


@pytest.mark.parametrize("level", ["minimal", "resampling", "full"])
@pytest.mark.parametrize(
    (
        "session_level",
        "ncomp_method",
        "ica_method",
        "maxlag",
        "croprun",
        "conn_measure",
        "freesurfer",
        "bids_filters",
    ),
    [
        _make_params(),
        _make_params(session_level=False),
        _make_params(session_level=True),
        _make_params(ncomp_method="aic"),
        _make_params(ncomp_method="kic"),
        _make_params(ncomp_method="mdl"),
        _make_params(croprun=4),
        _make_params(conn_measure=["correlation"]),
        _make_params(conn_measure=["correlation", "covariance"]),
        _make_params(bids_filters={"bold": {"task": "rest"}}),
        _make_params(freesurfer=False),
    ],
)
def test_init_fmristroke_wf(
    level: str,
    session_level: bool,
    ncomp_method: str,
    ica_method: str,
    maxlag: int,
    croprun: int,
    conn_measure: list[str],
    freesurfer: bool,
    bids_filters: dict,
):
    with mock_config():
        config.workflow.level = level
        config.workflow.freesurfer = freesurfer
        config.execution.run_sessionlevel = session_level
        config.workflow.ncomp_method = ncomp_method
        config.workflow.ica_method = ica_method
        config.workflow.maxlag = maxlag
        config.workflow.croprun = croprun
        config.workflow.conn_measure = conn_measure
        config.execution.bids_filters = bids_filters
        wf = init_fmristroke_wf()

    generate_expanded_graph(wf._create_flat_graph())
