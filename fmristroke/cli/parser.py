"""Parser."""
import sys

from .. import config


def _build_parser(**kwargs):
    """Build parser object.

    ``kwargs`` are passed to ``argparse.ArgumentParser`` (mainly useful for debugging).
    """
    from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
    from functools import partial
    from pathlib import Path

    from niworkflows.utils.spaces import OutputReferencesAction, Reference
    from packaging.version import Version


    def _path_exists(path, parser):
        """Ensure a given path exists."""
        if path is None or not Path(path).exists():
            raise parser.error(f"Path does not exist: <{path}>.")
        return Path(path).absolute()

    def _is_file(path, parser):
        """Ensure a given path exists and it is a file."""
        path = _path_exists(path, parser)
        if not path.is_file():
            raise parser.error(f"Path should point to a file (or symlink of file): <{path}>.")
        return path

    def _min_one(value, parser):
        """Ensure an argument is not lower than 1."""
        value = int(value)
        if value < 1:
            raise parser.error("Argument can't be less than one.")
        return value

    def _to_gb(value):
        scale = {"G": 1, "T": 10**3, "M": 1e-3, "K": 1e-6, "B": 1e-9}
        digits = "".join([c for c in value if c.isdigit()])
        units = value[len(digits) :] or "M"
        return int(digits) * scale[units[0]]

    def _drop_sub(value):
        return value[4:] if value.startswith("sub-") else value

    def _filter_pybids_none_any(dct):
        import bids

        return {
            k: bids.layout.Query.NONE if v is None else (bids.layout.Query.ANY if v == "*" else v)
            for k, v in dct.items()
        }

    def _bids_filter(value, parser):
        from json import JSONDecodeError, loads

        if value:
            if Path(value).exists():
                try:
                    return loads(Path(value).read_text(), object_hook=_filter_pybids_none_any)
                except JSONDecodeError:
                    raise parser.error(f"JSON syntax error in: <{value}>.")
            else:
                raise parser.error(f"Path does not exist: <{value}>.")

    currentv = Version(config.environment.version)
    
    parser = ArgumentParser(
        description="fMRIStroke: fMRI Stroke workflows v{}".format(
            config.environment.version
        ),
        formatter_class=ArgumentDefaultsHelpFormatter,
        **kwargs,
    )
    PathExists = partial(_path_exists, parser=parser)
    IsFile = partial(_is_file, parser=parser)
    PositiveInt = partial(_min_one, parser=parser)
    BIDSFilter = partial(_bids_filter, parser=parser)

    # Arguments as specified by BIDS-Apps
    # required, positional arguments
    # IMPORTANT: they must go directly with the parser object
    parser.add_argument(
        "bids_dir",
        action="store",
        type=PathExists,
        help="The root folder of a BIDS valid dataset (sub-XXXXX folders should "
        "be found at the top level in this folder).",
    )
    parser.add_argument(
        "output_dir",
        action="store",
        type=Path,
        help="The output path for the outcomes of preprocessing and visual reports",
    )
    parser.add_argument(
        "analysis_level",
        choices=["participant"],
        help='Processing stage to be run, only "participant" in the case of '
        "fMRIStroke as fMRIprep (see BIDS-Apps specification).",
    )

    g_bids = parser.add_argument_group("Options for filtering BIDS queries")
    g_bids.add_argument(
        "--participant-label",
        "--participant_label",
        action="store",
        nargs="+",
        type=_drop_sub,
        help="A space delimited list of participant identifiers or a single "
        "identifier (the sub- prefix can be removed)",
    )
    g_bids.add_argument(
        "--bids-filter-file",
        dest="bids_filters",
        action="store",
        type=BIDSFilter,
        metavar="FILE",
        help="A JSON file describing custom BIDS input filters using PyBIDS. "
    )
    g_bids.add_argument(
        "--fmriprep_dir",
        action="store",
        metavar="PATH",
        type=PathExists,
        help="Derivatives from fMRIPrep run",
    )

    g_perfm = parser.add_argument_group("Options to handle performance")
    g_perfm.add_argument(
        "--nprocs",
        "--nthreads",
        "--n_cpus",
        "--n-cpus",
        dest='nprocs',
        action="store",
        type=PositiveInt,
        help="Maximum number of threads across all processes",
    )
    g_perfm.add_argument(
        "--omp-nthreads",
        action="store",
        type=PositiveInt,
        help="Maximum number of threads per-process",
    )
    g_perfm.add_argument(
        "--mem",
        "--mem_mb",
        "--mem-mb",
        dest="memory_gb",
        action="store",
        type=_to_gb,
        metavar="MEMORY_MB",
        help="Upper bound memory limit for fMRIPrep processes",
    )
    g_perfm.add_argument(
        "--low-mem",
        action="store_true",
        help="Attempt to reduce memory usage (will increase disk usage in working directory)",
    )
    g_perfm.add_argument(
        "--use-plugin",
        "--nipype-plugin-file",
        action="store",
        metavar="FILE",
        type=IsFile,
        help="Nipype plugin configuration file",
    )

    g_subset = parser.add_argument_group("Options for performing only a subset of the workflow")
    g_subset.add_argument(
        "--reports-only",
        action="store_true",
        default=False,
        help="Only generate reports, don't run workflows. This will only rerun report "
        "aggregation, not reportlet generation for specific nodes.",
    )

    g_conf = parser.add_argument_group("Workflow configuration")
    g_conf.add_argument(
        "--output-spaces",
        nargs="*",
        action=OutputReferencesAction,
        help="""\
Standard and non-standard spaces to resample anatomical and functional images to. \
Standard spaces may be specified by the form \
``<SPACE>[:cohort-<label>][:res-<resolution>][...]``, where ``<SPACE>`` is \
a keyword designating a spatial reference, and may be followed by optional, \
colon-separated parameters. \
Non-standard spaces imply specific orientations and sampling grids. \
Important to note, the ``res-*`` modifier does not define the resolution used for \
the spatial normalization. To generate no BOLD outputs, use this option without specifying \
any spatial references. For further details, please check out \
https://fmriprep.readthedocs.io/en/%s/spaces.html"""
        % ("latest"),
    )
    
    g_pipe = parser.add_argument_group("Denoising pipeline configuration")
    g_pipe.add_argument(
        "--output_pipelines",
        nargs="+",
        help="Pipelines to use for denoising, space separated pipeline name, accepts either json file describing pipeline or name of pipeline from package."
    )
    
    g_confounds = parser.add_argument_group("Options relating to confounds")
    g_confounds.add_argument(
        "--ncomp_method",
        dest="ncomp_method",
        choices=["varexp", "aic", "kic", "mdl"],
        required=False,
        default= "varexp",
        type=str,
        help="method to estimate number of components for ICA lesion confounds",
    )
    g_confounds.add_argument(
        "--ica_method",
        dest="ica_method",
        required=False,
        action="store",
        default="canica",
        choices=["canica", "dictlearning"],
        type=str,
        help="Method to run ICA lesion",
    )

    # Hemodynamics options
    g_lag = parser.add_argument_group("Specific options for hemodynmics analysis")
    g_lag.add_argument(
        "--maxlag",
        default=10,
        type=int,
        help="Max lag for hemodynamic analysis",
    )
    
    # FreeSurfer options
    g_fs = parser.add_argument_group("Specific options for FreeSurfer preprocessing")
    g_fs.add_argument(
        "--fs-license-file",
        metavar="FILE",
        type=IsFile,
        help="Path to FreeSurfer license key file. Get it (for free) by registering"
        " at https://surfer.nmr.mgh.harvard.edu/registration.html",
    )
    g_fs.add_argument(
        "--freesurfer",
        action="store_true",
        default=True,
        dest="freesurfer",
        help="Was freesurfer run",
    )

    g_other = parser.add_argument_group("Other options")
    g_other.add_argument(
        "-w",
        "--work-dir",
        action="store",
        type=Path,
        default=Path("work").absolute(),
        help="Path where intermediate results should be stored",
    )
    g_other.add_argument(
        "--clean-workdir",
        action="store_true",
        default=False,
        help="Clears working directory of contents. Use of this flag is not "
        "recommended when running concurrent processes of fMRIPrep.",
    )
    g_other.add_argument(
        "--resource-monitor",
        action="store_true",
        default=False,
        help="Enable Nipype's resource monitoring to keep track of memory and CPU usage",
    )
    g_other.add_argument(
        "--config-file",
        action="store",
        metavar="FILE",
        help="Use pre-generated configuration file. Values in file will be overridden "
        "by command-line arguments.",
    )
    g_other.add_argument(
        "--write-graph",
        action="store_true",
        default=False,
        help="Write workflow graph.",
    )
    g_other.add_argument(
        "-v",
        "--verbose",
        dest="verbose_count",
        action="count",
        default=0,
        help="Increases log verbosity for each occurrence, debug level is -vvv",
    )
    g_other.add_argument(
        "--stop-on-first-crash",
        action="store_true",
        default=False,
        help="Force stopping on first crash, even if a work directory was specified.",
    )
    return parser


def parse_args(args=None, namespace=None):
    """Parse args and run further checks on the command line."""
    import logging

    from niworkflows.utils.spaces import Reference, SpatialReferences

    parser = _build_parser()
    opts = parser.parse_args(args, namespace)

    if opts.config_file:
        skip = {} if opts.reports_only else {"execution": ("run_uuid",)}
        config.load(opts.config_file, skip=skip, init=False)
        config.loggers.cli.info(f"Loaded previous configuration file {opts.config_file}")

    config.execution.log_level = int(max(25 - 5 * opts.verbose_count, logging.DEBUG))
    config.from_dict(vars(opts), init=['nipype'])

    # Initialize --output-spaces if not defined
    if config.execution.output_spaces is None:
        config.execution.output_spaces = SpatialReferences(
            [Reference("MNI152NLin2009cAsym", {"res": "native"})]
        )
    
    # Initialize --output-pipelines if not defined
    if config.execution.output_pipelines is None:
        config.execution.output_pipelines = ["SimpleGS", "ICLesionGS", "CompCorGS", "SimpleLesionGS"]

    # Retrieve logging level
    build_log = config.loggers.cli

    # Load base plugin_settings from file if --use-plugin
    if opts.use_plugin is not None:
        import yaml

        with open(opts.use_plugin) as f:
            plugin_settings = yaml.load(f, Loader=yaml.FullLoader)
        _plugin = plugin_settings.get("plugin")
        if _plugin:
            config.nipype.plugin = _plugin
            config.nipype.plugin_args = plugin_settings.get("plugin_args", {})
            config.nipype.nprocs = opts.nprocs or config.nipype.plugin_args.get(
                "n_procs", config.nipype.nprocs
            )

    # Resource management options
    # Note that we're making strong assumptions about valid plugin args
    # This may need to be revisited if people try to use batch plugins
    if 1 < config.nipype.nprocs < config.nipype.omp_nthreads:
        build_log.warning(
            f"Per-process threads (--omp-nthreads={config.nipype.omp_nthreads}) exceed "
            f"total threads (--nthreads/--n_cpus={config.nipype.nprocs})"
        )

    bids_dir = config.execution.bids_dir
    output_dir = config.execution.output_dir
    work_dir = config.execution.work_dir
    version = config.environment.version

    # Wipe out existing work_dir
    if opts.clean_workdir and work_dir.exists():
        from niworkflows.utils.misc import clean_directory

        build_log.info(f"Clearing previous fMRIStroke working directory: {work_dir}")
        if not clean_directory(work_dir):
            build_log.warning(f"Could not clear all contents of working directory: {work_dir}")

    # Update the config with an empty dict to trigger initialization of all config
    # sections (we used `init=False` above).
    # This must be done after cleaning the work directory, or we could delete an
    # open SQLite database
    config.from_dict({})

    # Ensure input and output folders are not the same
    if output_dir == bids_dir:
        parser.error(
            "The selected output folder is the same as the input BIDS folder. "
            "Please modify the output path (suggestion: %s)."
            % bids_dir
            / "derivatives"
            / ("fmriprep-%s" % version.split("+")[0])
        )

    if bids_dir in work_dir.parents:
        parser.error(
            "The selected working directory is a subdirectory of the input BIDS folder. "
            "Please modify the output path."
        )

    # Setup directories
    config.execution.log_dir = config.execution.output_dir / "logs"
    # Check and create output and working directories
    config.execution.log_dir.mkdir(exist_ok=True, parents=True)
    work_dir.mkdir(exist_ok=True, parents=True)

    # Force initialization of the BIDSLayout
    config.execution.init()
    all_subjects = config.execution.layout.get_subjects()
    if config.execution.participant_label is None:
        config.execution.participant_label = all_subjects

    participant_label = set(config.execution.participant_label)
    missing_subjects = participant_label - set(all_subjects)
    if missing_subjects:
        parser.error(
            "One or more participant labels were not found in the BIDS directory: "
            "%s." % ", ".join(missing_subjects)
        )

    config.execution.participant_label = sorted(participant_label)
