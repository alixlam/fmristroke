"""
Orchestrating the BOLD-preprocessing workflow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: init_lesion_preproc_wf
.. autofunction:: init_func_lesion_derivatives_wf

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
from .confounds import init_confs_wf
from .outputs import init_func_lesion_derivatives_wf
from .registration import init_lesionplot_wf
from .denoise import init_denoise_wf

from .lagmaps import init_hemodynamic_wf

def init_lesion_preproc_wf(bold_file, boldref_file, boldmask_file, confounds_file, confounds_metadata):
    """
    This workflow controls the lesion specific stages *fMRIStroke*.

    Parameters
    ----------
    bold_file
        Path to NIfTI file 
    boldref_file
        Path to NIfTI file
    boldmask
        Path to NIfTI file
    confounds_file
        Path to tsv file
    confounds_metadata
        Path to json file describin confounds


    Inputs
    ------
    bold_t1
        BOLD series NIfTI file in anatomical space
    boldref_t1
        BOLD reference  NIfTI file in anatomical space
    boldmask_t1
        BOLD mask in anatomical space
    t1w_preproc
        Bias-corrected structural template image
    t1w_mask
        Mask of the skull-stripped template image
    t1w_dseg
        Segmentation of preprocessed structural image, including
        gray-matter (GM), white-matter (WM) and cerebrospinal fluid (CSF)
    t1w_aseg
        Segmentation of structural image, done with FreeSurfer.
    t1w_aparc
        Parcellation of structural image, done with FreeSurfer.
    t1w_tpms
        List of tissue probability maps in T1w space
    confounds_file
        confounds file
    roi
        ROI mask

    Outputs
    -------
    roi_mask_t1
        ROI mask in T1w space
    roi_mask_std
        roi mask, resampled to template space
    confounds
        TSV of confounds
    confounds_metadata
        Confounds metadata dictionary
    lagmaps
        lagmaps in T1w space


    """
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow
    from niworkflows.interfaces.nibabel import ApplyMask
    from niworkflows.interfaces.reportlets.registration import (
        SimpleBeforeAfterRPT as SimpleBeforeAfter,
    )
    from niworkflows.interfaces.utility import DictMerge, KeySelect

    img = nb.load(bold_file[0] if isinstance(bold_file, (list, tuple)) else bold_file)
    nvols = 1 if img.ndim < 4 else img.shape[3]

    mem_gb = {"filesize": 1, "resampled": 1, "largemem": 1}
    bold_tlen = 10

    # Have some options handy
    omp_nthreads = config.nipype.omp_nthreads
    spaces = config.workflow.spaces
    pipelines = config.workflow.pipelines
    output_dir = str(config.execution.output_dir)
    freesurfer = config.workflow.freesurfer

    # Extract BIDS entities and metadata from BOLD file(s)
    entities = extract_entities(bold_file)
    layout = config.execution.layout

    # Extract metadata
    all_metadata = [layout.get_metadata(fname) for fname in listify(bold_file)]

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
                "std2anat_xfm",
                "anat2std_xfm",
                "confounds_file",
                "confounds_metadata",
                "roi",
                "templates",
            ]
        ),
        name="inputnode",
    )
    inputnode.inputs.bold_t1 = bold_file
    inputnode.inputs.boldref_t1 = boldref_file
    inputnode.inputs.boldmask_t1 = boldmask_file
    inputnode.inputs.confounds_file = confounds_file
    inputnode.inputs.confounds_metadata = confounds_metadata

    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
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
        pipelines=pipelines,
    )
    
    func_lesion_derivatives_wf.inputs.inputnode.all_source_files = bold_file
    func_lesion_derivatives_wf.inputs.inputnode.source_file = bold_file

    
    workflow.connect([
        (outputnode, func_lesion_derivatives_wf, [
            ("confounds_file", "inputnode.confounds_file"),
            ("confounds_metadata", "inputnode.confounds_metadata"),
            ("lagmaps", "inputnode.lagmaps"),
            ("template", "inputnode.template"),
            ("spatial_reference", "inputnode.spatial_reference"),
            ("denoised_bold_t1", "inputnode.denoised_bold_t1"),
            ("denoised_bold_std", "inputnode.denoised_bold_std"),
            ("pipeline", "inputnode.pipeline")
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
            ("std2anat_xfm", "inputnode.std2anat_xfm"),
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
    
    # REGISTRATION PLOT WORKFLOW   ###################################################
    lesion_plot = init_lesionplot_wf(
        mem_gb=mem_gb["largemem"],
        name = "lesion_plot_wf"
    )
    workflow.connect([
        (inputnode, lesion_plot, [
            ("boldref_t1", "inputnode.boldref_t1"),
            ("t1w_preproc", "inputnode.t1w"),
            ("roi", "inputnode.t1w_roi"),
            ("t1w_mask", "inputnode.t1w_mask")
            ])
    ])

    # DENOISING WORKFLOW ########################################################
    denoising_wf = init_denoise_wf(
        mem_gb=mem_gb,
        omp_nthreads=omp_nthreads,
        metadata=metadata,
        pipelines=pipelines,
        spaces=spaces,
        name = "bold_denoise_wf"
    )
    workflow.connect([
        (inputnode, denoising_wf, [
            ("bold_t1", "inputnode.bold_t1"),
            ("anat2std_xfm", "inputnode.anat2std_xfm"),
            ("templates", "inputnode.templates"),
            ("confounds_metadata", "inputnode.confounds_metadata")
            ]),
        (lesion_confounds_wf, denoising_wf, [
            ("outputnode.confounds_file", "inputnode.confounds_file")
        ]),
        (denoising_wf, outputnode, [
            ("outputnode.template", "template"),
            ("outputnode.spatial_reference", "spatial_reference"),
            ("outputnode.denoised_bold_t1", "denoised_bold_t1"),
            ("outputnode.denoised_bold_std", "denoised_bold_std"),
            ("outputnode.pipeline", "pipeline")
        ])
    ])
    
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
