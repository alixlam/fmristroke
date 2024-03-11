"""fMRI Stroke workflow."""
from .. import config

EXITCODE: int = -1


def main():
    """Entry point."""
    import gc
    import sys
    import warnings
    from multiprocessing import Manager, Process
    from os import EX_SOFTWARE
    from pathlib import Path

    from fmriprep.utils.bids import (
        write_bidsignore,
        write_derivative_description,
    )

    from .parser import parse_args
    from .workflow import build_workflow

    parse_args()

    # CRITICAL Save the config to a file. This is necessary because the execution graph
    # is built as a separate process to keep the memory footprint low. The most
    # straightforward way to communicate with the child process is via the
    # filesystem.
    config_file = (
        config.execution.work_dir / config.execution.run_uuid / "config.toml"
    )
    config_file.parent.mkdir(exist_ok=True, parents=True)
    config.to_filename(config_file)

    # CRITICAL Call build_workflow(config_file, retval) in a subprocess.
    # Because Python on Linux does not ever free virtual memory (VM), running the
    # workflow construction jailed within a process preempts excessive VM
    # buildup.
    with Manager() as mgr:
        retval = mgr.dict()
        p = Process(target=build_workflow, args=(str(config_file), retval))
        p.start()
        p.join()
        retval = dict(retval.items())  # Convert to base dictionary

        if p.exitcode:
            retval["return_code"] = p.exitcode

    global EXITCODE
    EXITCODE = retval.get("return_code", 0)
    fmristroke_wf = retval.get("workflow", None)

    # CRITICAL Load the config from the file. This is necessary because the ``build_workflow``
    # function executed constrained in a process may change the config (and thus the global
    # state of fMRIStroke).
    config.load(config_file)

    if config.execution.reports_only:
        sys.exit(int(EXITCODE > 0))

    if fmristroke_wf and config.execution.write_graph:
        fmristroke_wf.write_graph(
            graph2use="colored", format="svg", simple_form=True
        )

    EXITCODE = EXITCODE or (fmristroke_wf is None) * EX_SOFTWARE
    if EXITCODE != 0:
        sys.exit(EXITCODE)

    # Clean up master process before running workflow, which may create forks
    gc.collect()

    # Sentry tracking
    config.loggers.workflow.log(
        15,
        "\n".join(
            ["fMRIStroke config:"]
            + ["\t\t%s" % s for s in config.dumps().splitlines()]
        ),
    )
    config.loggers.workflow.log(25, "fMRIStroke started!")
    errno = 1  # Default is error exit unless otherwise set
    try:
        fmristroke_wf.run(**config.nipype.get_plugin())
    except Exception as e:
        config.loggers.workflow.critical("fMRIStroke failed: %s", e)
        raise
    else:
        config.loggers.workflow.log(25, "fMRIStroke finished successfully!")

    finally:
        from fmriprep.reports.core import generate_reports
        from pkg_resources import resource_filename as pkgrf

        # Generate reports phase
        failed_reports = generate_reports(
            config.execution.participant_label,
            config.execution.output_dir,
            config.execution.run_uuid,
            config=pkgrf("fmristroke", "data/reports-specs.yml"),
            packagename="fmriprep",
        )
        write_derivative_description(
            config.execution.bids_dir, config.execution.fmriprep_dir
        )
        write_bidsignore(config.execution.fmriprep_dir)
        sys.exit(int((errno + failed_reports) > 0))


if __name__ == "__main__":
    main()
