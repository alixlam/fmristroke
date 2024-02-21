"""
fMRIStroke base processing workflows
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: init_fmristroke_wf
.. autofunction:: init_single_subject_wf

"""

import os
import sys
import warnings
from copy import deepcopy

from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from niworkflows.utils.connections import listify
from packaging.version import Version

from .. import config
from ..interfaces import DerivativesDataSink
from .anat.base import init_roi_preproc_wf
from .bold.base import init_lesion_connectivity_wf, init_lesion_preproc_wf


def init_fmristroke_wf():
    """
    Build *fMRIStroke*'s pipeline.

    This workflow organizes the execution of FMRISTROKE, with a sub-workflow for
    each subject.

    """
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow

    ver = Version(config.environment.version)

    fmriStroke_wf = Workflow(name=f"fmristroke_{ver.major}_{ver.minor}_wf")
    fmriStroke_wf.base_dir = config.execution.work_dir
    single_subject_wfs = []
    for subject_id in config.execution.participant_label:
        single_subject_wf = init_single_subject_wf(subject_id)

        single_subject_wf.config["execution"]["crashdump_dir"] = str(
            config.execution.fmriprep_dir
            / f"sub-{subject_id}"
            / "log"
            / config.execution.run_uuid
        )
        for node in single_subject_wf._get_all_nodes():
            node.config = deepcopy(single_subject_wf.config)
        fmriStroke_wf.add_nodes([single_subject_wf])

        single_subject_wfs.append(single_subject_wf)

        # Dump a copy of the config file into the log directory
        log_dir = (
            config.execution.fmriprep_dir
            / f"sub-{subject_id}"
            / "log"
            / config.execution.run_uuid
        )
        log_dir.mkdir(exist_ok=True, parents=True)
        config.to_filename(log_dir / "fmristroke.toml")

    return fmriStroke_wf


def init_single_subject_wf(subject_id: str):
    """
    Organize the preprocessing pipeline for a single subject.

    It collects and reports information about the subject, and prepares
    sub-workflows to perform pipelines..
    Functional is performed using a separate workflow for each
    individual BOLD series.

    Parameters
    ----------
    subject_id : :obj:`str`
        Subject label for this single-subject workflow.

    Inputs
    ------
    subjects_dir : :obj:`str`
        FreeSurfer's ``$SUBJECTS_DIR``.

    """
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow
    from niworkflows.interfaces.bids import BIDSInfo
    from niworkflows.utils.spaces import Reference
    from smriprep.utils.bids import collect_derivatives

    from ..interfaces.bids import BIDSDerivativeDataGrabber
    from ..utils.bids import (
        collect_bold_derivatives,
        collect_roi_mask,
        group_runs,
    )

    name = "single_subject_%s_wf" % subject_id
    fmriprep_dir = str(config.execution.fmriprep_dir)
    raw_data_dir = str(config.execution.bids_dir)
    spaces = config.workflow.spaces
    session_level = config.execution.run_sessionlevel
    if session_level:
        croprun = config.workflow.croprun

    bold_derivatives = collect_bold_derivatives(
        fmriprep_dir,
        subject_id,
        bids_filters=config.execution.bids_filters,
    )
    for imtype in ["bold_t1", "boldref_t1", "confounds_file"]:
        if not bold_derivatives[imtype]:
            config.loggers.workflow.warning(
                f"""\
        A ttempted to access pre-existing bold derivatives at \
        <{config.execution.fmriprep_dir}>, however not all expectations of fMRIPrep \
        were met (for participant <{subject_id}>."""
            )
    std_spaces = spaces.get_spaces(nonstandard=False, dim=(3,))

    anat_derivatives = collect_derivatives(
        fmriprep_dir,
        subject_id,
        std_spaces,
        config.workflow.freesurfer,
    )
    if anat_derivatives is None:
        config.loggers.workflow.warning(
            f"""\
        Attempted to access pre-existing anatomical derivatives at \
        <{config.execution.fmriprep_dir}>, however not all expectations of fMRIPrep \
        were met (for participant <{subject_id}>, spaces <{', '.join(std_spaces)}>, \
        >)."""
        )

    # Get roi mask, currently fmriprep does not preprocess the mask,
    # it is not in the derivatives folder
    roi = collect_roi_mask(bids_dir=raw_data_dir, subject_id=subject_id)

    workflow = Workflow(name=name)

    bidssrc = pe.Node(
        BIDSDerivativeDataGrabber(
            bold_derivatives=bold_derivatives,
            anat_derivatives=anat_derivatives,
            subject_id=subject_id,
        ),
        name="bidssrc",
    )

    # ROI resampling
    roi_anat_wf = init_roi_preproc_wf(name="roi_std_wf")
    roi_anat_wf.inputs.inputnode.roi = roi["roi"][0]

    # fmt:off
    workflow.connect([
            (bidssrc, roi_anat_wf, [
                    ("t1w_preproc", "inputnode.t1w_preproc"),
                    ("anat2std_xfm", "inputnode.anat2std_xfm"),
                    ("std2anat_xfm", "inputnode.std2anat_xfm"),
                    ("template", "inputnode.template"),
                ],
            )
    ])

    func_pre_desc = """
Functional data lesion specific data preprocessing

: For each of the {num_bold} BOLD runs found per subject (across all
tasks and sessions), the following lesion specific preprocessing was performed.
""".format(
        num_bold=len(bold_derivatives["bold_t1"])
    )

    func_preproc_wfs = []

    for (
        bold_file,
        boldref,
        boldmask,
        confounds_file,
        confounds_metadata,
    ) in zip(
        bold_derivatives["bold_t1"],
        bold_derivatives["boldref_t1"],
        bold_derivatives["boldmask_t1"],
        bold_derivatives["confounds_file"],
        bold_derivatives["confounds_metadata"],
    ):
        lesion_preproc_wf = init_lesion_preproc_wf(
            bold_file,
            boldref,
            boldmask,
            confounds_file,
            confounds_metadata,
        )
        if lesion_preproc_wf is None:
            continue
        lesion_preproc_wf.inputs.inputnode.roi = roi["roi"][0]
        lesion_preproc_wf.__desc__ = func_pre_desc + (
            lesion_preproc_wf.__desc__ or ""
        )

        # fmt:off
        workflow.connect([
            (bidssrc, lesion_preproc_wf, [
                ('t1w_preproc', 'inputnode.t1w_preproc'),
                ('t1w_mask', 'inputnode.t1w_mask'),
                ('t1w_dseg', 'inputnode.t1w_dseg'),
                ('t1w_tpms', 'inputnode.t1w_tpms'),
                ('std2anat_xfm', 'inputnode.std2anat_xfm'),
                ('anat2std_xfm', 'inputnode.anat2std_xfm'),
                ('template', 'inputnode.templates'),
                # Undefined if freesurfer was not run
                ('t1w_aseg', 'inputnode.t1w_aseg')

            ]),
        ])
        # fmt:on

        func_preproc_wfs.append(lesion_preproc_wf)

    node_names = [func_preproc_wf.name for func_preproc_wf in func_preproc_wfs]

    # Merge nodes
    merge_denoised_t1 = pe.Node(
        niu.Merge(
            len(func_preproc_wfs),
            no_flatten=True),
        name="merge_denoised_t1")
    merge_denoised_std = pe.Node(
        niu.Merge(
            len(func_preproc_wfs),
            no_flatten=True),
        name="merge_denoised_std")
    merge_boldmask = pe.Node(
        niu.Merge(
            len(func_preproc_wfs),
            no_flatten=True),
        name="merge_boldmask")

    for i, bold_wf in enumerate(func_preproc_wfs):
        # fmt:off
        workflow.connect([
            (bold_wf, merge_boldmask, [("outputnode.boldmask", f"in{i+1}")]),
            (bold_wf, merge_denoised_std, [("outputnode.bold_denoised_std", f"in{i+1}")]),
            (bold_wf, merge_denoised_t1, [("outputnode.bold_denoised_t1", f"in{i+1}")]),
        ])
        # fmt:on

    if session_level:
        bold_derivatives["bold_t1"] = group_runs(bold_derivatives["bold_t1"])

    grp_bold_denoised_t1 = pe.Node(
        niu.Function(
            function=_group_outputs),
        name="group_bold_denoised_t1")
    grp_bold_denoised_t1.inputs.node_names = node_names
    grp_bold_denoised_t1.inputs.session_level = session_level
    grp_bold_denoised_std = pe.Node(
        niu.Function(
            function=_group_outputs),
        name="group_bold_denoised_std")
    grp_bold_denoised_std.inputs.node_names = node_names
    grp_bold_denoised_std.inputs.session_level = session_level
    grp_bold_mask = pe.Node(
        niu.Function(
            function=_group_outputs),
        name="group_boldmask")
    grp_bold_mask.inputs.node_names = node_names
    grp_bold_mask.inputs.session_level = session_level

    # fmt:off
    workflow.connect([
        (merge_denoised_std, grp_bold_denoised_std, [("out", "in_outputs")]),
        (merge_denoised_t1, grp_bold_denoised_t1, [("out", "in_outputs")]),
        (merge_boldmask, grp_bold_mask, [("out", "in_outputs")]),
    ])
    # fmt:on

    for i, bold_file in enumerate(bold_derivatives["bold_t1"]):
        lesion_connectivity_wf = init_lesion_connectivity_wf(bold_file)

        if lesion_connectivity_wf is None:
            continue

        lesion_connectivity_wf.inputs.inputnode.roi = roi["roi"][0]

        select_denoised_t1 = pe.Node(
            niu.Select(index=i), name=f"select_denoised_t1_{i}"
        )
        select_denoised_std = pe.Node(
            niu.Select(index=i), name=f"select_denoised_std_{i}"
        )
        select_boldmask = pe.Node(
            niu.Select(index=i), name=f"select_boldmask_{i}"
        )

        # fmt:off
        workflow.connect([
            (bidssrc, lesion_connectivity_wf, [
                ('t1w_preproc', 'inputnode.t1w_preproc'),
                ('t1w_tpms', 'inputnode.t1w_tpms'),
                ('t1w_mask', 'inputnode.t1w_mask'),
                ('std2anat_xfm', 'inputnode.std2anat_xfm'),
                ('anat2std_xfm', 'inputnode.anat2std_xfm'),
                ('template', 'inputnode.templates'),
                # Undefined if freesurfer was not run
                ('t1w_aseg', 'inputnode.t1w_aseg')
            ]),
            (bold_wf, lesion_connectivity_wf, [
                ("outputnode.pipelines", "inputnode.pipelines")]),
            (grp_bold_denoised_std, select_denoised_std, [
                ("out", "inlist")]),
            (grp_bold_denoised_t1, select_denoised_t1, [
                ("out", "inlist")]),
            (grp_bold_mask, select_boldmask, [
                ("out", "inlist")]),
            (select_denoised_std, lesion_connectivity_wf, [
                ("out", "inputnode.bold_denoised_std")]),
            (select_denoised_t1, lesion_connectivity_wf, [
                ("out", "inputnode.bold_denoised_t1")]),
            (select_boldmask, lesion_connectivity_wf, [
                ("out", "inputnode.boldmask_t1")]),
        ])
        workflow.connect([
            (roi_anat_wf, lesion_connectivity_wf, [
                ("outputnode.roi_mask_std", "inputnode.roi_std")],)
        ])
        # fmt:on

    return workflow


def _group_outputs(in_outputs, node_names, session_level=False):
    import re
    from itertools import groupby

    nodes = [(name, out) for name, out in zip(node_names, in_outputs)]

    # Sort by node name
    nodes.sort(key=lambda x: x[0])

    def _grp_runs(x):
        if "_run_" not in x[0]:
            return x
        run = re.search("_run_\\d*", x[0]).group(0)
        return x[0].replace(run, "_run_?")

    grouped = []
    for _, nodes in groupby(nodes, key=_grp_runs):
        nodes = list(nodes)
        grouped.append(nodes)

    if session_level:
        grouped_outputs = [[run[1] for run in outputs] for outputs in grouped]
    else:
        grouped_outputs = [[run[1]] for outputs in grouped for run in outputs]
    return grouped_outputs
