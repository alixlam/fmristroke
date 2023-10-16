from pathlib import Path
from toml import loads
from tempfile import mkdtemp
import shutil
import os
from contextlib import contextmanager

@contextmanager
def mock_config():
	"""Create mock config for tests and documentation"""
	from .. import config

	config_file = Path(__file__).parent / "../../data/tests/mock_config.toml"
	settings = loads(config_file.read_text())

	for sectionname, configs in settings.items():
		if sectionname != 'environment':
			section = getattr(config, sectionname)
			section.load(configs, init=False)
	config.nipype.omp_nthreads = 1
	config.nipype.init()
	config.loggers.init()
	config.init_spaces()

	config.execution.work_dir = Path(mkdtemp())
	config.execution.bids_dir = Path(__file__.parent / "../../data/tests/FCStroke_2009").absolute()
	config.execution.fmriprep_dir = Path(mkdtemp())
	config.execution.init()

	yield

	shutil.rmtree(config.execution.work_dir)
	shutil.rmtree(config.execution.fmriprep_dir)