"""
Orchestrating the ROI-resampling workflow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: init_roi_preproc_wf
.. autofunction:: init_anat_derivatives_wf

"""
import os

import nibabel as nb
import numpy as np
from nipype.interfaces import utility as niu
from nipype.interfaces.fsl import Split as FSLSplit
from nipype.pipeline import engine as pe
from niworkflows.utils.connections import listify, pop_file

from ... import config
from ...interfaces import DerivativesDataSink

# BOLD workflows
from .outputs import init_anat_lesion_derivatives_wf
from .resample import (
    init_roi_std_trans_wf,
)

def init_roi_preproc_wf():
    """
    This workflow controls the lesion specific stages *fMRIStroke*.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from fmriprep.workflows.tests import mock_config
            from fmriprep import config
            from fmriprep.workflows.bold.base import init_func_preproc_wf
            with mock_config():
                bold_file = config.execution.bids_dir / "sub-01" / "func" \
                    / "sub-01_task-mixedgamblestask_run-01_bold.nii.gz"
                wf = init_func_preproc_wf(str(bold_file))

    Parameters
    ----------
    


    Inputs
    ------
    template
        List of templates to target
    anat2std_xfm
        List of transform files, collated with templates
    std2anat_xfm
        List of inverse transform files, collated with templates
    roi
        ROI mask

    Outputs
    -------
    roi_mask_t1
        ROI mask in T1w space
    roi_mask_std
        roi mask, resampled to template space


    """
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow
    from niworkflows.interfaces.nibabel import ApplyMask
    from niworkflows.interfaces.reportlets.registration import (
        SimpleBeforeAfterRPT as SimpleBeforeAfter,
    )
    from niworkflows.interfaces.utility import DictMerge, KeySelect

    mem_gb = {"filesize": 1, "resampled": 1, "largemem": 1}

    # Have some options handy
    omp_nthreads = config.nipype.omp_nthreads
    spaces = config.workflow.spaces
    output_dir = str(config.execution.output_dir)
    freesurfer = config.workflow.freesurfer

    # Extract BIDS entities
    layout = config.execution.layout

    # Take first file as reference
    ref_file = pop_file(bold_file)
    metadata = all_metadata[0]

    if os.path.isfile(ref_file):
        bold_tlen, mem_gb = _create_mem_gb(ref_file)

    wf_name = _get_wf_name(ref_file)
    config.loggers.workflow.debug(
        "Creating lesion processing workflow for <%s> (%.2f GB / %d TRs). "
        "Memory resampled/largemem=%.2f/%.2f GB.",
        ref_file,
        mem_gb["filesize"],
        bold_tlen,
        mem_gb["resampled"],
        mem_gb["largemem"],
    )
    # Build workflow
    workflow = Workflow(name=wf_name)
    workflow.__postdesc__ = """\
All resamplings can be performed with *a single interpolation
step* by composing all the pertinent transformations (i.e. co-registrations to anatomical and output spaces).
Gridded (volumetric) resamplings were performed using `antsApplyTransforms` (ANTs),
configured with Lanczos interpolation to minimize the smoothing
effects of other kernels [@lanczos].
"""

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "bold_t1",
                "boldref_t1",
                "boldmask_t1",
                "t1w_preproc",
                "t1w_mask",
                "t1w_dseg",
                "t1w_tpms",
                "t1w_aseg",
                "anat2std_xfm",
                "std2anat_xfm",
                "template",
                "confounds_file",
                "roi",
            ]
        ),
        name="inputnode",
    )
    inputnode.inputs.bold_t1 = bold_file
    inputnode.inputs.boldref_t1 = boldref_file
    inputnode.inputs.boldmask_t1 = boldmask_file
    inputnode.inputs.confounds_file = confounds_file

    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "roi_mask_t1",
                "roi_mask_std",
                "confounds_file",
                "confounds_metadata",
                "lagmaps",
            ]
        ),
        name="outputnode",
    )
    func_lesion_derivatives_wf = init_func_lesion_derivatives_wf(
        bids_root=layout.root,
        output_dir=output_dir,
        spaces=spaces,
    )
    func_lesion_derivatives_wf.inputs.inputnode.all_source_files = bold_file
    func_lesion_derivatives_wf.inputs.inputnode.source_file = bold_file

    
    anat_lesion_derivatives_wf = init_anat_lesion_derivatives_wf(
        bids_root=layout.root,
        output_dir=output_dir,
        spaces=spaces,
    )
    workflow.connect([
        (outputnode, func_lesion_derivatives_wf, [
            ("confounds_file", "inputnode.confounds_file"),
            ("confounds_metadata", "inputnode.confounds_metadata"),
            ("lagmaps", "inputnode.lagmaps"),
        ]),
        (inputnode, anat_lesion_derivatives_wf, [
            ("t1w_preproc", "inputnode.all_source_files"),
        ]),
    ])
    
    #LESION CONFOUNDS  #######################################################
    lesion_confounds_wf = init_confs_wf(
        mem_gb=mem_gb["largemem"],
        metadata=metadata,
        freesurfer=freesurfer,
        ncomp_method=config.workflow.ncomp_method,
        ica_method=config.workflow.ica_method,
        name="lesion_confounds_wf",
    )
    
    workflow.connect([
        # Connect bold_confounds_wf
        (inputnode, lesion_confounds_wf, [
            ("t1w_tpms", "inputnode.t1w_tpms"),
            ("t1w_mask", "inputnode.t1w_mask"),
            ("bold_t1", "inputnode.bold_t1"),
            ("boldref_t1", "inputnode.boldref_t1"),
            ("boldmask_t1", "inputnode.boldmask_t1"),
            ("confounds_file","inputnode.confounds_file"),
            ("roi", "inputnode.roi"),
        ]),
        
        (lesion_confounds_wf, outputnode, [
            ("outputnode.confounds_file", "confounds_file"),
            ("outputnode.confounds_metadata", "confounds_metadata"),
        ]),
    ])
    # HEMODYNAMIC WORKFLOW   #######################################################
    hemodynamics_wf = init_hemodynamic_wf(
        mem_gb=mem_gb["largemem"],
        metadata=metadata,
        maxlag=config.workflow.maxlag,
        name="hemodynamic_wf"
    )
    workflow.connect([
        # Connect bold_confounds_wf
        (inputnode, hemodynamics_wf, [
            ("bold_t1", "inputnode.bold_t1"),
            ("roi", "inputnode.roi"),
            ("t1w_preproc", "inputnode.t1w_preproc"),
            ("t1w_mask", "inputnode.t1w_mask"),
            ("t1w_tpms", "inputnode.t1w_tpms"),
            ("confounds_file","inputnode.confounds_file"),
            # undefined if freesurfer was not run
            ("t1w_aseg", "inputnode.t1w_aseg"),
        ]),
        (lesion_confounds_wf, hemodynamics_wf, [
            ("outputnode.boldmask", "inputnode.boldmask")
        ]),
        (hemodynamics_wf, outputnode, [
            ("outputnode.lagmaps", "lagmaps"),
        ]),
    ])
    
    # ROI RESAMPLING #######################################################

    if spaces.get_spaces(nonstandard=False, dim=(3,)):
        # Apply transforms in 1 shot
        roi_std_trans_wf = init_roi_std_trans_wf(
            mem_gb=mem_gb["resampled"],
            omp_nthreads=omp_nthreads,
            spaces=spaces,
            name="lesion_std_trans_wf",
            use_compression=not config.execution.low_mem,
        )
        # fmt:off
        workflow.connect([
            (inputnode, roi_std_trans_wf, [
                ("template", "inputnode.templates"),
                ("anat2std_xfm", "inputnode.anat2std_xfm"),
                ("roi", "inputnode.roi")
            ]),
            (roi_std_trans_wf, outputnode, [
                ("outputnode.roi_mask_std", "roi_mask_std"),
            ]),
        ])
        workflow.connect([
            (roi_std_trans_wf, anat_lesion_derivatives_wf, [
                ("outputnode.template", "inputnode.template"),
                ("outputnode.spatial_reference", "inputnode.spatial_reference"),
                ("outputnode.roi_mask_std", "inputnode.roi_mask_std"),
            ]),
        ])
        # fmt:on

    

    # REPORTING ############################################################

    # Fill-in datasinks of reportlets seen so far
    for node in workflow.list_node_names():
        if node.split(".")[-1].startswith("ds_report"):
            workflow.get_node(node).inputs.base_directory = output_dir
            workflow.get_node(node).inputs.source_file = ref_file

    return workflow


def _create_mem_gb(bold_fname):
    img = nb.load(bold_fname)
    nvox = int(np.prod(img.shape, dtype='u8'))
    # Assume tools will coerce to 8-byte floats to be safe
    bold_size_gb = 8 * nvox / (1024**3)
    bold_tlen = img.shape[-1]
    mem_gb = {
        "filesize": bold_size_gb,
        "resampled": bold_size_gb * 4,
        "largemem": bold_size_gb * (max(bold_tlen / 100, 1.0) + 4),
    }

    return bold_tlen, mem_gb


def _get_wf_name(bold_fname):
    """
    Derive the workflow name for supplied BOLD file.

    >>> _get_wf_name("/completely/made/up/path/sub-01_task-nback_bold.nii.gz")
    'func_preproc_task_nback_wf'
    >>> _get_wf_name("/completely/made/up/path/sub-01_task-nback_run-01_echo-1_bold.nii.gz")
    'func_preproc_task_nback_run_01_echo_1_wf'

    """
    from nipype.utils.filemanip import split_filename

    fname = split_filename(bold_fname)[1]
    fname_nosub = "_".join(fname.split("_")[1:])
    name = "lesion_preproc_" + fname_nosub.replace(".", "_").replace(" ", "").replace(
        "-", "_"
    ).replace("_bold", "_wf")

    return name

def extract_entities(file_list):
    """
    Return a dictionary of common entities given a list of files.

    Examples
    --------
    >>> extract_entities("sub-01/anat/sub-01_T1w.nii.gz")
    {'subject': '01', 'suffix': 'T1w', 'datatype': 'anat', 'extension': '.nii.gz'}
    >>> extract_entities(["sub-01/anat/sub-01_T1w.nii.gz"] * 2)
    {'subject': '01', 'suffix': 'T1w', 'datatype': 'anat', 'extension': '.nii.gz'}
    >>> extract_entities(["sub-01/anat/sub-01_run-1_T1w.nii.gz",
    ...                   "sub-01/anat/sub-01_run-2_T1w.nii.gz"])
    {'subject': '01', 'run': [1, 2], 'suffix': 'T1w', 'datatype': 'anat', 'extension': '.nii.gz'}

    """
    from collections import defaultdict

    from bids.layout import parse_file_entities

    entities = defaultdict(list)
    for e, v in [
        ev_pair for f in listify(file_list) for ev_pair in parse_file_entities(f).items()
    ]:
        entities[e].append(v)

    def _unique(inlist):
        inlist = sorted(set(inlist))
        if len(inlist) == 1:
            return inlist[0]
        return inlist

    return {k: _unique(v) for k, v in entities.items()}