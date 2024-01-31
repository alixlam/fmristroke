"""
A Python module to maintain unique, run-wide *fMRIStroke* settings.

This module implements the memory structures to keep a consistent, singleton config.
Settings are passed across processes via filesystem, and a copy of the settings for
each run and subject is left under
``<fmristroke_dir>/sub-<participant_id>/log/<run_unique_id>/fmristroke.toml``.
Settings are stored using :abbr:`ToML (Tom's Markup Language)`.
The module has a :py:func:`~fmristroke.config.to_filename` function to allow writing out
the settings to hard disk in *ToML* format, which looks like:

This config file is used to pass the settings across processes,
using the :py:func:`~fmristroke.config.load` function.

Configuration sections
----------------------
.. autoclass:: environment
    :members:
.. autoclass:: execution
    :members:
.. autoclass:: workflow
    :members:
.. autoclass:: nipype
    :members:

Usage
-----
A config file is used to pass settings and collect information as the execution
graph is built across processes.


"""
import os
from multiprocessing import set_start_method

# Disable NiPype etelemetry always
_disable_et = bool(
    os.getenv("NO_ET") is not None or os.getenv("NIPYPE_NO_ET") is not None
)
os.environ["NIPYPE_NO_ET"] = "1"
os.environ["NO_ET"] = "1"

CONFIG_FILENAME = "fmristroke.toml"

try:
    set_start_method("forkserver")
except RuntimeError:
    pass  # context has been already set
finally:
    # Defer all custom import for after initializing the forkserver and
    # ignoring the most annoying warnings
    import logging
    import random
    import sys
    from pathlib import Path
    from time import strftime
    from uuid import uuid4

    from fmriprep import __version__ as _fmriprep_ver
    from nipype import __version__ as _nipype_ver
    from templateflow import __version__ as _tf_ver

    from . import __version__


if not hasattr(sys, "_is_pytest_session"):
    sys._is_pytest_session = False  # Trick to avoid sklearn's FutureWarnings


# Add a new level between INFO and WARNING
logging.addLevelName(25, "IMPORTANT")
logging.addLevelName(15, "VERBOSE")  # Add a new level between INFO and DEBUG

DEFAULT_MEMORY_MIN_GB = 0.01

# Ping NiPype eTelemetry once if env var was not set
# workers on the pool will have the env variable set from the master process
if not _disable_et:
    # Just get so analytics track one hit
    from contextlib import suppress

    from requests import ConnectionError, ReadTimeout
    from requests import get as _get_url

    with suppress((ConnectionError, ReadTimeout)):
        _get_url("https://rig.mit.edu/et/projects/nipy/nipype", timeout=0.05)

# Execution environment
_exec_env = os.name


_templateflow_home = Path(
    os.getenv(
        "TEMPLATEFLOW_HOME",
        os.path.join(os.getenv("HOME"), ".cache", "templateflow"),
    )
)

try:
    from psutil import virtual_memory

    _free_mem_at_start = round(virtual_memory().available / 1024**3, 1)
except Exception:
    _free_mem_at_start = None

_oc_limit = "n/a"
_oc_policy = "n/a"
try:
    # Memory policy may have a large effect on types of errors experienced
    _proc_oc_path = Path("/proc/sys/vm/overcommit_memory")
    if _proc_oc_path.exists():
        _oc_policy = {"0": "heuristic", "1": "always", "2": "never"}.get(
            _proc_oc_path.read_text().strip(), "unknown"
        )
        if _oc_policy != "never":
            _proc_oc_kbytes = Path("/proc/sys/vm/overcommit_kbytes")
            if _proc_oc_kbytes.exists():
                _oc_limit = _proc_oc_kbytes.read_text().strip()
            if (
                _oc_limit in ("0", "n/a")
                and Path("/proc/sys/vm/overcommit_ratio").exists()
            ):
                _oc_limit = "{}%".format(
                    Path("/proc/sys/vm/overcommit_ratio").read_text().strip()
                )
except Exception:
    pass


class _Config:
    """An abstract class forbidding instantiation."""

    _paths = tuple()

    def __init__(self):
        """Avert instantiation."""
        raise RuntimeError("Configuration type is not instantiable.")

    @classmethod
    def load(cls, settings, init=True, ignore=None):
        """Store settings from a dictionary."""
        ignore = ignore or {}
        for k, v in settings.items():
            if k in ignore or v is None:
                continue
            if k in cls._paths:
                setattr(cls, k, Path(v).absolute())
            elif hasattr(cls, k):
                setattr(cls, k, v)

        if init:
            try:
                cls.init()
            except AttributeError:
                pass

    @classmethod
    def get(cls):
        """Return defined settings."""
        from niworkflows.utils.spaces import Reference, SpatialReferences

        out = {}
        for k, v in cls.__dict__.items():
            if k.startswith("_") or v is None:
                continue
            if callable(getattr(cls, k)):
                continue
            if k in cls._paths:
                v = str(v)
            if isinstance(v, SpatialReferences):
                v = " ".join([str(s) for s in v.references]) or None
            if isinstance(v, Reference):
                v = str(v) or None
            out[k] = v
        return out


class environment(_Config):
    """
    Read-only options regarding the platform and environment.

    The ``environment`` section is not loaded in from file,
    only written out when settings are exported.
    This config section is useful when reporting issues.

    """

    cpu_count = os.cpu_count()
    """Number of available CPUs."""
    exec_env = _exec_env
    """A string representing the execution platform."""
    free_mem = _free_mem_at_start
    """Free memory at start."""
    overcommit_policy = _oc_policy
    """Linux's kernel virtual memory overcommit policy."""
    overcommit_limit = _oc_limit
    """Linux's kernel virtual memory overcommit limits."""
    nipype_version = _nipype_ver
    """Nipype's current version."""
    templateflow_version = _tf_ver
    """The TemplateFlow client version installed."""
    fmriprep_version = _fmriprep_ver
    """*fMRIPrep*'s version."""
    version = __version__
    """*fmriStroke's version."""


class nipype(_Config):
    """Nipype settings."""

    crashfile_format = "txt"
    """The file format for crashfiles, either text or pickle."""
    get_linked_libs = False
    """Run NiPype's tool to enlist linked libraries for every interface."""
    memory_gb = None
    """Estimation in GB of the RAM this workflow can allocate at any given time."""
    nprocs = os.cpu_count()
    """Number of processes (compute tasks) that can be run in parallel (multiprocessing only)."""
    omp_nthreads = None
    """Number of CPUs a single process can access for multithreaded execution."""
    plugin = "MultiProc"
    """NiPype's execution plugin."""
    plugin_args = {
        "maxtasksperchild": 1,
        "raise_insufficient": False,
    }
    """Settings for NiPype's execution plugin."""
    resource_monitor = False
    """Enable resource monitor."""
    stop_on_first_crash = True
    """Whether the workflow should stop or continue after the first error."""

    @classmethod
    def get_plugin(cls):
        """Format a dictionary for Nipype consumption."""
        out = {
            "plugin": cls.plugin,
            "plugin_args": cls.plugin_args,
        }
        if cls.plugin in ("MultiProc", "LegacyMultiProc"):
            out["plugin_args"]["n_procs"] = int(cls.nprocs)
            if cls.memory_gb:
                out["plugin_args"]["memory_gb"] = float(cls.memory_gb)
        return out

    @classmethod
    def init(cls):
        """Set NiPype configurations."""
        from nipype import config as ncfg

        # Configure resource_monitor
        if cls.resource_monitor:
            ncfg.update_config(
                {
                    "monitoring": {
                        "enabled": cls.resource_monitor,
                        "sample_frequency": "0.5",
                        "summary_append": True,
                    }
                }
            )
            ncfg.enable_resource_monitor()

        # Nipype config (logs and execution)
        ncfg.update_config(
            {
                "execution": {
                    "crashdump_dir": str(execution.log_dir),
                    "crashfile_format": cls.crashfile_format,
                    "get_linked_libs": cls.get_linked_libs,
                    "stop_on_first_crash": cls.stop_on_first_crash,
                    "check_version": False,  # disable future telemetry
                }
            }
        )

        if cls.omp_nthreads is None:
            cls.omp_nthreads = min(
                cls.nprocs - 1 if cls.nprocs > 1 else os.cpu_count(), 8
            )


class execution(_Config):
    """Configure run-level settings."""

    fmriprep_dir = None
    """A path where fmriprep derivatives are found."""
    bids_dir = None
    """An existing path to the dataset, which must be BIDS-compliant."""
    bids_filters = None
    """A dictionary of BIDS selection filters."""
    layout = None
    """A :py:class:`~bids.layout.BIDSLayout` object, see :py:func:`init`."""
    log_dir = None
    """The path to a directory that contains execution logs."""
    log_level = 25
    """Output verbosity."""
    low_mem = None
    """Utilize uncompressed NIfTIs and other tricks to minimize memory allocation."""
    output_dir = None
    """Folder where derivatives will be stored."""
    output_spaces = None
    """List of (non)standard spaces designated (with the ``--output-spaces`` flag of
    the command line) as spatial references for outputs."""
    output_pipelines = None
    """List of pipelines designated (with the ``--output-pipelines`` flag of the comand line) as pipelines for denoising"""
    output_atlases = None
    """List of pipelines designated (with the ``--output-atlases`` flag of the comand line) as atlases for connectivity"""
    run_sessionlevel = False
    """Concatenate runs from the same session and task"""
    reports_only = False
    """Only build the reports, based on the reportlets found in a cached working directory."""
    run_uuid = f"{strftime('%Y%m%d-%H%M%S')}_{uuid4()}"
    """Unique identifier of this particular run."""
    participant_label = None
    """List of participant identifiers that are to be preprocessed."""
    templateflow_home = _templateflow_home
    """The root folder of the TemplateFlow client."""
    work_dir = Path("work").absolute()
    """Path to a working directory where intermediate results will be available."""
    write_graph = False
    """Write out the computational graph corresponding to the planned preprocessing."""

    _layout = None

    _paths = (
        "bids_dir",
        "fmriprep_dir",
        "layout",
        "log_dir",
        "output_dir",
        "templateflow_home",
        "work_dir",
    )

    @classmethod
    def init(cls):
        """Create a new BIDS Layout accessible with :attr:`~execution.layout`."""
        if cls._layout is None:
            import re

            from bids.layout import BIDSLayout
            from bids.layout.index import BIDSLayoutIndexer

            _db_path = cls.work_dir / cls.run_uuid / "bids_db"
            _db_path.mkdir(exist_ok=True, parents=True)

            # Recommended after PyBIDS 12.1
            _indexer = BIDSLayoutIndexer(
                validate=False,
                ignore=(
                    "code",
                    "stimuli",
                    "sourcedata",
                    "models",
                    re.compile(r"^\."),
                    re.compile(
                        r"sub-[a-zA-Z0-9]+(/ses-[a-zA-Z0-9]+)?/(beh|dwi|eeg|ieeg|meg|perf)"
                    ),
                ),
            )
            cls._layout = BIDSLayout(
                str(cls.fmriprep_dir),
                database_path=_db_path,
                indexer=_indexer,
            )
            cls.bids_database_dir = _db_path
        cls.layout = cls._layout
        if cls.bids_filters:
            from bids.layout import Query

            # unserialize pybids Query enum values
            for acq, filters in cls.bids_filters.items():
                cls.bids_filters[acq] = {
                    k: getattr(Query, v[7:-4])
                    if not isinstance(v, Query) and "Query" in v
                    else v
                    for k, v in filters.items()
                }


# These variables are not necessary anymore
del _exec_env
del _nipype_ver
del _templateflow_home
del _tf_ver
del _free_mem_at_start
del _oc_limit
del _oc_policy


class workflow(_Config):
    """Configure the particular execution graph of this workflow."""

    freesurfer = None
    """Wether freesurfer must be used"""
    ncomp_method = None
    """Method to estimate components for lesion ICA"""
    ica_method = None
    """ICA method to use"""
    maxlag = None
    """Maxlag for hemodynamics computations"""
    spaces = None
    """Keeps the :py:class:`~niworkflows.utils.spaces.SpatialReferences`
    instance keeping standard and nonstandard spaces."""
    pipelines = None
    """Keeps the pipelines instances for denoising"""
    atlases = None
    """Keeps the atlases instances for connectivity"""
    croprun = None
    """Crops n first volumes in fMRI run before concatenating when session-level is true"""
    conn_measure = None
    """Measure to use for connectivity"""


class loggers:
    """Keep loggers easily accessible (see :py:func:`init`)."""

    _fmt = (
        "%(asctime)s,%(msecs)d %(name)-2s " "%(levelname)-2s:\n\t %(message)s"
    )
    _datefmt = "%y%m%d-%H:%M:%S"

    default = logging.getLogger()
    """The root logger."""
    cli = logging.getLogger("cli")
    """Command-line interface logging."""
    workflow = logging.getLogger("nipype.workflow")
    """NiPype's workflow logger."""
    interface = logging.getLogger("nipype.interface")
    """NiPype's interface logger."""
    utils = logging.getLogger("nipype.utils")
    """NiPype's utils logger."""

    @classmethod
    def init(cls):
        """
        Set the log level, initialize all loggers into :py:class:`loggers`.

            * Add new logger levels (25: IMPORTANT, and 15: VERBOSE).
            * Add a new sub-logger (``cli``).
            * Logger configuration.

        """
        from nipype import config as ncfg

        if not cls.cli.hasHandlers():
            _handler = logging.StreamHandler(stream=sys.stdout)
            _handler.setFormatter(
                logging.Formatter(fmt=cls._fmt, datefmt=cls._datefmt)
            )
            cls.cli.addHandler(_handler)
        cls.default.setLevel(execution.log_level)
        cls.cli.setLevel(execution.log_level)
        cls.interface.setLevel(execution.log_level)
        cls.workflow.setLevel(execution.log_level)
        cls.utils.setLevel(execution.log_level)
        ncfg.update_config(
            {
                "logging": {
                    "log_directory": str(execution.log_dir),
                    "log_to_file": True,
                }
            }
        )


class seeds(_Config):
    """Initialize the PRNG and track random seed assignments"""

    _random_seed = None
    master = None
    """Master random seed to initialize the Pseudorandom Number Generator (PRNG)"""
    numpy = None
    """Seed used by NumPy"""

    @classmethod
    def init(cls):
        if cls._random_seed is not None:
            cls.master = cls._random_seed
        if cls.master is None:
            cls.master = random.randint(1, 65536)
        random.seed(cls.master)  # initialize the PRNG
        # functions to set program specific seeds
        cls.numpy = _set_numpy_seed()


def _set_numpy_seed():
    """NumPy's random seed is independent from Python's `random` module"""
    import numpy as np

    val = random.randint(1, 65536)
    np.random.seed(val)
    return val


def from_dict(settings, init=True, ignore=None):
    """Read settings from a flat dictionary.

    Arguments
    ---------
    setting : dict
        Settings to apply to any configuration
    init : `bool` or :py:class:`~collections.abc.Container`
        Initialize all, none, or a subset of configurations.
    ignore : :py:class:`~collections.abc.Container`
        Collection of keys in ``setting`` to ignore
    """

    # Accept global True/False or container of configs to initialize
    def initialize(x):
        return init if init in (True, False) else x in init

    nipype.load(settings, init=initialize("nipype"), ignore=ignore)
    execution.load(settings, init=initialize("execution"), ignore=ignore)
    workflow.load(settings, init=initialize("workflow"), ignore=ignore)
    seeds.load(settings, init=initialize("seeds"), ignore=ignore)

    loggers.init()


def load(filename, skip=None, init=True):
    """Load settings from file.

    Arguments
    ---------
    filename : :py:class:`os.PathLike`
        TOML file containing fMRIStroke configuration.
    skip : dict or None
        Sets of values to ignore during load, keyed by section name
    init : `bool` or :py:class:`~collections.abc.Container`
        Initialize all, none, or a subset of configurations.
    """
    from toml import loads

    skip = skip or {}

    # Accept global True/False or container of configs to initialize
    def initialize(x):
        return init if init in (True, False) else x in init

    filename = Path(filename)
    settings = loads(filename.read_text())
    for sectionname, configs in settings.items():
        if sectionname != "environment":
            section = getattr(sys.modules[__name__], sectionname)
            ignore = skip.get(sectionname)
            section.load(configs, ignore=ignore, init=initialize(sectionname))
    init_spaces()
    init_pipelines()
    init_atlases()


def get(flat=False):
    """Get config as a dict."""
    settings = {
        "environment": environment.get(),
        "execution": execution.get(),
        "workflow": workflow.get(),
        "nipype": nipype.get(),
        "seeds": seeds.get(),
    }
    if not flat:
        return settings

    return {
        ".".join((section, k)): v
        for section, configs in settings.items()
        for k, v in configs.items()
    }


def dumps():
    """Format config into toml."""
    from toml import dumps

    return dumps(get())


def to_filename(filename):
    """Write settings to file."""
    filename = Path(filename)
    filename.write_text(dumps())


def init_spaces(checkpoint=True):
    """Initialize the :attr:`~workflow.spaces` setting."""
    from niworkflows.utils.spaces import Reference, SpatialReferences

    spaces = execution.output_spaces or SpatialReferences()
    if not isinstance(spaces, SpatialReferences):
        spaces = SpatialReferences(
            [
                ref
                for s in spaces.split(" ")
                for ref in Reference.from_string(s)
            ]
        )

    # Add the default standard space if not already present (required by
    # several sub-workflows)
    if "MNI152NLin2009cAsym" not in spaces.get_spaces(
        nonstandard=False, dim=(3,)
    ):
        spaces.add(Reference("MNI152NLin2009cAsym", {}))

    if checkpoint and not spaces.is_cached():
        spaces.checkpoint()

    # Make the SpatialReferences object available
    workflow.spaces = spaces


def init_pipelines():
    """Initialize the :attr:`~workflow.pipelines` setting."""
    from pkg_resources import resource_filename as pkgrf

    from .utils.pipelines import Pipeline, Pipelines

    pipelines = execution.output_pipelines or Pipelines()
    known_pipelines = ["SimpleGS", "ICLesionGS", "CompCorGS", "SimpleLesionGS", "CompCorLesionGS", "ICLesionCompCorGS"]
    if not isinstance(pipelines, Pipelines):
        pipelines_filenames = [
            pkgrf("fmristroke", f"data/denoising/{s}.json")
            if s in known_pipelines
            else s
            for s in pipelines
        ]
        pipelines = Pipelines(
            [Pipeline.from_json(f) for f in pipelines_filenames]
        )

    # Make the Pipelines object available
    workflow.pipelines = pipelines


def init_atlases():
    """Initialize the :attr:`~workflow.pipelines` setting."""
    from pkg_resources import resource_filename as pkgrf

    from .utils.atlas import Atlas, Atlases

    atlases = execution.output_atlases or Atlases()
    known_atlases = ["Scheafer"]
    if not isinstance(atlases, Atlases):
        atlases_filenames = [
            pkgrf("fmristroke", f"data/atlas/{s}.json")
            if s in known_atlases
            else s
            for s in atlases
        ]
        atlases = Atlases([Atlas.from_json(f) for f in atlases_filenames])

    # Make the Atlases object available
    workflow.atlases = atlases
