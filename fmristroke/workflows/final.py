"""
fMRIStroke base processing workflows
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
from .bold.base import init_lesion_preproc_wf


def init_fmristroke_wf():
    """
    Build *fMRIStroke*'s pipeline.

    This workflow organizes the execution of FMRISTROKE, with a sub-workflow for
    each subject.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from fmriprep.workflows.tests import mock_config
            from fmriprep.workflows.base import init_fmriprep_wf
            with mock_config():
                wf = init_fmriprep_wf()

    """
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow

    ver = Version(config.environment.version)

    fmriStroke_wf = Workflow(name=f'fmristroke_{ver.major}_{ver.minor}_wf')
    fmriStroke_wf.base_dir = config.execution.work_dir

    for subject_id in config.execution.participant_label:
        single_subject_wf = init_single_subject_wf(subject_id)

        single_subject_wf.config['execution']['crashdump_dir'] = str(
            config.execution.fmriprep_dir / f"sub-{subject_id}" / "log" / config.execution.run_uuid
        )
        for node in single_subject_wf._get_all_nodes():
            node.config = deepcopy(single_subject_wf.config)
        fmriStroke_wf.add_nodes([single_subject_wf])

        # Dump a copy of the config file into the log directory
        log_dir = (
            config.execution.fmriprep_dir / f"sub-{subject_id}" / 'log' / config.execution.run_uuid
        )
        log_dir.mkdir(exist_ok=True, parents=True)
        config.to_filename(log_dir / 'fmristroke.toml')

    return fmriStroke_wf


def init_single_subject_wf(subject_id: str):
    """
    Organize the preprocessing pipeline for a single subject.

    It collects and reports information about the subject, and prepares
    sub-workflows to perform pipelines..
    Functional is performed using a separate workflow for each
    individual BOLD series.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from fmriprep.workflows.tests import mock_config
            from fmriprep.workflows.base import init_single_subject_wf
            with mock_config():
                wf = init_single_subject_wf('01')

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
    from ..utils.bids import collect_bold_derivatives, collect_roi_mask
    from ..interfaces.bids import BIDSDerivativeDataGrabber

    name = "single_subject_%s_wf" % subject_id
    fmriprep_dir = str(config.execution.fmriprep_dir)
    raw_data_dir = str(config.execution.bids_dir)
    spaces = config.workflow.spaces


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
    roi = collect_roi_mask(
        bids_dir=raw_data_dir,
        subject_id=subject_id
    )

    workflow = Workflow(name=name)
    #workflow.__desc__ = """
#Results included in this manuscript come from preprocessing
#performed using *fMRIPrep* {fmriprep_ver}
#(@fmriprep1; @fmriprep2; RRID:SCR_016216),
#which is based on *Nipype* {nipype_ver}
#(@nipype1; @nipype2; RRID:SCR_002502).

#""".format(
#        fmriprep_ver=config.environment.version, nipype_ver=config.environment.nipype_version
#    )
#    workflow.__postdesc__ = """

#Many internal operations of *fMRIPrep* use
#*Nilearn* {nilearn_ver} [@nilearn, RRID:SCR_001362],
#mostly within the functional processing workflow.
#For more details of the pipeline, see [the section corresponding
#to workflows in *fMRIPrep*'s documentation]\
#(https://fmriprep.readthedocs.io/en/latest/workflows.html \
#"FMRIPrep's documentation")."""


    bidssrc = pe.Node(
        BIDSDerivativeDataGrabber(
            bold_derivatives=bold_derivatives,
            anat_derivatives=anat_derivatives,
            subject_id=subject_id,
        ),
        name='bidssrc',
    )

    func_pre_desc = """
Functional data lesion specific data preprocessing

: For each of the {num_bold} BOLD runs found per subject (across all
tasks and sessions), the following lesion specific preprocessing was performed.
""".format(
        num_bold=len(bold_derivatives['bold_t1'])
    )

    func_preproc_wfs = []
    for bold_file, boldref, boldmask, confounds_file in zip(
        bold_derivatives['bold_t1'], 
        bold_derivatives['boldref_t1'],
        bold_derivatives['boldmask_t1'],
        bold_derivatives['confounds_file']):
        lesion_preproc_wf = init_lesion_preproc_wf(bold_file, boldref, boldmask, confounds_file)
        if lesion_preproc_wf is None:
            continue
        lesion_preproc_wf.inputs.inputnode.roi = roi["roi"][0]
        lesion_preproc_wf.__desc__ = func_pre_desc + (lesion_preproc_wf.__desc__ or "")
        
        # fmt:off
        workflow.connect([
            (bidssrc, lesion_preproc_wf, [
                ('t1w_preproc', 'inputnode.t1w_preproc'),
                ('t1w_mask', 'inputnode.t1w_mask'),
                ('t1w_dseg', 'inputnode.t1w_dseg'),
                ('t1w_tpms', 'inputnode.t1w_tpms'),
                ('template', 'inputnode.template'),
                ('anat2std_xfm', 'inputnode.anat2std_xfm'),
                ('std2anat_xfm', 'inputnode.std2anat_xfm'),
                # Undefined if freesurfer was not run
                ('t1w_aseg', 'inputnode.t1w_aseg')

            ]),
        ])
        # fmt:on
        func_preproc_wfs.append(lesion_preproc_wf)

    return workflow
