"""Utilities and mocks for testing and documentation building."""
import os
import shutil
from contextlib import contextmanager
from pathlib import Path
from tempfile import mkdtemp

from toml import loads

from ... import data


@contextmanager
def mock_config(bids_dir=None, fmriprep_dir=None):
    """Create a mock config for documentation and testing purposes."""
    from ... import config

    settings = loads(data.load.readable('tests/mock_config.toml').read_text())
    for sectionname, configs in settings.items():
        if sectionname != 'environment':
            section = getattr(config, sectionname)
            section.load(configs, init=False)
    config.nipype.omp_nthreads = 1
    config.nipype.init()
    config.loggers.init()
    config.init_spaces()
    config.init_pipelines()
    config.init_atlases()

    bids_dir = bids_dir or data.load('tests/FCStroke_2009').absolute()
    fmriprep_dir = fmriprep_dir or data.load('tests/FCStroke_2009_deriv').absolute()

    config.execution.work_dir = Path(mkdtemp())
    config.execution.output_dir = Path(mkdtemp())
    config.execution.bids_dir = bids_dir
    config.execution.fmriprep_dir = fmriprep_dir
    config.execution.bids_database_dir = None
    config.execution._layout = None
    config.execution.init()

    yield

    shutil.rmtree(config.execution.work_dir)
    shutil.rmtree(config.execution.output_dir)
