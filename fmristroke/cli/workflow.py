"""
The workflow builder factory method.

All the checks and the construction of the workflow are done
inside this function that has pickleable inputs and output
dictionary (``retval``) to allow isolation using a
``multiprocessing.Process`` that allows fmristroke to enforce
a hard-limited memory-scope.

"""


def build_workflow(config_file, retval):
    """Create the Nipype Workflow that supports the whole execution graph."""
    from pathlib import Path

    from fmriprep.reports.core import generate_reports
    from fmriprep.utils.misc import check_deps
    from niworkflows.utils.bids import collect_participants
    from pkg_resources import resource_filename as pkgrf

    from .. import config
    from ..workflows.final import init_fmristroke_wf

    config.load(config_file)
    build_log = config.loggers.workflow

    fmriprep_dir = config.execution.fmriprep_dir
    version = config.environment.version

    retval["return_code"] = 1
    retval["workflow"] = None

    banner = [f"Running fMRIStroke version {version}"]

    build_log.log(25, f"\n{' ' * 9}".join(banner))
    # First check that bids_dir looks like a BIDS folder
    subject_list = collect_participants(
        config.execution.layout,
        participant_label=config.execution.participant_label,
    )

    # Called with reports only
    if config.execution.reports_only:
        build_log.log(
            25,
            "Running --reports-only on participants %s",
            ", ".join(subject_list),
        )
        retval["return_code"] = generate_reports(
            config.execution.participant_label,
            config.execution.fmriprep_dir,
            config.execution.run_uuid,
            config=pkgrf("fmristroke", "data/reports-specs.yml"),
            packagename="fmriprep",
        )
        return retval

    # Build main workflow
    init_msg = [
        "Building fMRIStroke's workflow:",
        f"BIDS dataset path: {config.execution.bids_dir}.",
        f"fMRIPrep path: {config.execution.fmriprep_dir}.",
        f"Participant list: {subject_list}.",
        f"Run identifier: {config.execution.run_uuid}.",
        f"Output spaces: {config.execution.output_spaces}.",
        f"Denoising pipelines: {config.execution.output_pipelines}.",
        f"Atlases: {config.execution.output_atlases}.",
    ]

    build_log.log(25, f"\n{' ' * 11}* ".join(init_msg))

    retval["workflow"] = init_fmristroke_wf()

    # Check workflow for missing commands
    missing = check_deps(retval["workflow"])
    if missing:
        build_log.critical(
            "Cannot run fMRIStroke. Missing dependencies:%s",
            "\n\t* ".join(
                [""]
                + [f"{cmd} (Interface: {iface})" for iface, cmd in missing]
            ),
        )
        retval["return_code"] = 127  # 127 == command not found.
        return retval

    config.to_filename(config_file)
    build_log.info(
        "fMRIStroke workflow graph with %d nodes built successfully.",
        len(retval["workflow"]._get_all_nodes()),
    )
    retval["return_code"] = 0
    return retval

def build_boilerplate(config_file, workflow):
    """Write boilerplate in an isolated process."""
    from fmristroke import config

    config.load(config_file)
    logs_path = config.execution.output_dir / "logsfmristroke"
    boilerplate = workflow.visit_desc()
    citation_files = {
        ext: logs_path / ("CITATION.%s" % ext) for ext in ("bib", "tex", "md", "html")
    }

    if boilerplate:
        # To please git-annex users and also to guarantee consistency
        # among different renderings of the same file, first remove any
        # existing one
        for citation_file in citation_files.values():
            try:
                citation_file.unlink()
            except FileNotFoundError:
                pass

    citation_files["md"].write_text(boilerplate)

    if not citation_files["md"].exists():
        from pathlib import Path
        from subprocess import CalledProcessError, TimeoutExpired, check_call

        from pkg_resources import resource_filename as pkgrf

        bib_text = Path(pkgrf("fmristroke", "data/boilerplate.bib")).read_text()
        citation_files["bib"].write_text(
            bib_text.replace("fMRIStroke <version>", f"fMRIStroke {config.environment.version}")
        )

        # Generate HTML file resolving citations
        cmd = [
            "pandoc",
            "-s",
            "--bibliography",
            str(citation_files["bib"]),
            "--citeproc",
            "--metadata",
            'pagetitle="fMRIPrep citation boilerplate"',
            str(citation_files["md"]),
            "-o",
            str(citation_files["html"]),
        ]

        config.loggers.cli.info("Generating an HTML version of the citation boilerplate...")
        try:
            check_call(cmd, timeout=10)
        except (FileNotFoundError, CalledProcessError, TimeoutExpired):
            config.loggers.cli.warning("Could not generate CITATION.html file:\n%s", " ".join(cmd))

        # Generate LaTex file resolving citations
        cmd = [
            "pandoc",
            "-s",
            "--bibliography",
            str(citation_files["bib"]),
            "--natbib",
            str(citation_files["md"]),
            "-o",
            str(citation_files["tex"]),
        ]
        config.loggers.cli.info("Generating a LaTeX version of the citation boilerplate...")
        try:
            check_call(cmd, timeout=10)
        except (FileNotFoundError, CalledProcessError, TimeoutExpired):
            config.loggers.cli.warning("Could not generate CITATION.tex file:\n%s", " ".join(cmd))