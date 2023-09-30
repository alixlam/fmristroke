"""
Orchestrating the BOLD-preprocessing workflow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: init_func_preproc_wf
.. autofunction:: init_func_derivatives_wf

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
from .outputs import init_func_lesion_derivatives_wf, init_anat_lesion_derivatives_wf
from .registration import 
from .resample import (
    init_roi_std_trans_wf,
)

def init_func_preproc_wf(bold_file, has_fieldmap=False):
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
    bold_file
        Path to NIfTI file (single echo) or list of paths to NIfTI files (multi-echo)
    has_fieldmap : :obj:`bool`
        Signals the workflow to use inputnode fieldmap files

    Inputs
    ------
    bold_file
        BOLD series NIfTI file
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
    template
        List of templates to target
    anat2std_xfm
        List of transform files, collated with templates
    std2anat_xfm
        List of inverse transform files, collated with templates

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
    if nvols <= 5 - config.execution.sloppy:
        config.loggers.workflow.warning(
            f"Too short BOLD series (<= 5 timepoints). Skipping processing of <{bold_file}>."
        )
        return

    mem_gb = {"filesize": 1, "resampled": 1, "largemem": 1}
    bold_tlen = 10

    # Have some options handy
    omp_nthreads = config.nipype.omp_nthreads
    spaces = config.workflow.spaces
    fmriprep_dir = str(config.execution.fmriprep_dir)

    # Extract BIDS entities and metadata from BOLD file(s)
    entities = extract_entities(bold_file)
    layout = config.execution.layout

    # Extract metadata
    all_metadata = [layout.get_metadata(fname) for fname in listify(bold_file)]

    # Take first file as reference
    ref_file = pop_file(bold_file)
    metadata = all_metadata[0]
    # get original image orientation

    echo_idxs = listify(entities.get("echo", []))
    
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

    # Find associated sbref, if possible
    overrides = {
        "suffix": "sbref",
        "extension": [".nii", ".nii.gz"],
    }

    if config.execution.bids_filters:
        overrides.update(config.execution.bids_filters.get('sbref', {}))
    sb_ents = {**entities, **overrides}
    sbref_files = layout.get(return_type="file", **sb_ents)

    sbref_msg = f"No single-band-reference found for {os.path.basename(ref_file)}."
    if sbref_files and "sbref" in config.workflow.ignore:
        sbref_msg = "Single-band reference file(s) found and ignored."
        sbref_files = []
    elif sbref_files:
        sbref_msg = "Using single-band reference file(s) {}.".format(
            ",".join([os.path.basename(sbf) for sbf in sbref_files])
        )
    config.loggers.workflow.info(sbref_msg)

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
                "bold_file",
                "t1w_preproc",
                "t1w_mask",
                "t1w_dseg",
                "t1w_tpms",
                "t1w_aseg",
                "t1w_aparc",
                "anat2std_xfm",
                "std2anat_xfm",
                "template",
                "anat_ribbon",
                "t1w2fsnative_xfm",
                "confounds_file",
                "t1_bold_xform",
                "bold_t1_xform",
                "t1w_roi",

            ]
        ),
        name="inputnode",
    )
    inputnode.inputs.bold_file = bold_file

    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "roi_mask_t1",
                "roi_mask_std",
                "confounds",
                "confounds_metadata",
                "lagmaps",
            ]
        ),
        name="outputnode",
    )

    # Generate a brain-masked conversion of the t1w
    t1w_brain = pe.Node(ApplyMask(), name="t1w_brain")

    func_lesion_derivatives_wf = init_func_lesion_derivatives_wf(
        bids_root=layout.root,
        output_dir=fmriprep_dir,
        spaces=spaces,
    )
    func_lesion_derivatives_wf.inputs.inputnode.all_source_files = bold_file
    
    anat_lesion_derivatives_wf = init_anat_lesion_derivatives_wf(
        bids_root=layout.root,
        output_dir=fmriprep_dir,
        spaces=spaces,
    )
    anat_lesion_derivatives_wf.inputs.inputnode.all_source_files = bold_file

    # fmt:off
    workflow.connect([
        (outputnode, func_lesion_derivatives_wf, [
            ("confounds", "inputnode.confounds"),
            ("confounds_metadata", "inputnode.confounds_metadata"),
            ("lagmaps", "inputnode.lagmaps"),
        ]),
        (outputnode, anat_lesion_derivatives_wf, [
            ("roi_mask_std", "inputnode.roi_mask_std"),
        ]),
    ])
    
    
    # fmt:on

    # get confounds
    lesion_confounds_wf = init_confs_wf(
        mem_gb=mem_gb["largemem"],
        metadata=metadata,
        freesurfer=freesurfer,
        ncomp_method=config.workflow.ncomp_method,
        ica_method=config.workflow.ica_method,
        name="lesion_confounds_wf",
    )
    
    # MAIN WORKFLOW STRUCTURE #######################################################
    # fmt:off
    workflow.connect([
        # Prepare masked T1w image
        (inputnode, t1w_brain, [("t1w_preproc", "in_file"),
                                ("t1w_mask", "in_mask")]),
        
        # Connect bold_confounds_wf
        (inputnode, lesion_confounds_wf, [
            ("t1w_tpms", "inputnode.t1w_tpms"),
            ("t1w_mask", "inputnode.t1w_mask"),
            ("bold", "inputnode.bold"),
            ("boldref", "inputnode.boldref"),
            ("confounds_file","inputnode.confounds_file"),
            ("t1w_roi", "inputnode.roi"),
            ("t1_bold_xform", "inputnode.t1_bold_xform")
        ]),
        
        (lesion_confounds_wf, outputnode, [
            ("outputnode.confounds_file", "confounds"),
            ("outputnode.confounds_metadata", "confounds_metadata"),
        ]),
    ])
    # 

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
                ("t1w_roi", "inputnode.roi")
            ]),
            (roi_std_trans_wf, outputnode, [
                ("outputnode.roi_mask_std", "roi_mask_std"),
            ]),
        ])
        workflow.connect([
            (roi_std_trans_wf, func_lesion_derivatives_wf, [
                ("outputnode.template", "inputnode.template"),
                ("outputnode.spatial_reference", "inputnode.spatial_reference"),
                ("outputnode.roi_mask_std", "inputnode.roi_mask_std"),
            ]),
        ])
        # fmt:on

    

    # REPORTING ############################################################
    ds_report_summary = pe.Node(
        DerivativesDataSink(desc="summary", datatype="figures", dismiss_entities=("echo",)),
        name="ds_report_summary",
        run_without_submitting=True,
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )

    ds_report_validation = pe.Node(
        DerivativesDataSink(desc="validation", datatype="figures", dismiss_entities=("echo",)),
        name="ds_report_validation",
        run_without_submitting=True,
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )

    # fmt:off
    workflow.connect([
        (summary, ds_report_summary, [("out_report", "in_file")]),
        (initial_boldref_wf, ds_report_validation, [("outputnode.validation_report", "in_file")]),
    ])
    # fmt:on

    # Fill-in datasinks of reportlets seen so far
    for node in workflow.list_node_names():
        if node.split(".")[-1].startswith("ds_report"):
            workflow.get_node(node).inputs.base_directory = fmriprep_dir
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
    name = "func_preproc_" + fname_nosub.replace(".", "_").replace(" ", "").replace(
        "-", "_"
    ).replace("_bold", "_wf")

    return name


def _to_join(in_file, join_file):
    """Join two tsv files if the join_file is not ``None``."""
    from niworkflows.interfaces.utility import JoinTSVColumns

    if join_file is None:
        return in_file
    res = JoinTSVColumns(in_file=in_file, join_file=join_file).run()
    return res.outputs.out_file


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


def get_img_orientation(imgf):
    """Return the image orientation as a string"""
    img = nb.load(imgf)
    return "".join(nb.aff2axcodes(img.affine))


def get_estimator(layout, fname):
    field_source = layout.get_metadata(fname).get("B0FieldSource")
    if isinstance(field_source, str):
        field_source = (field_source,)

    if field_source is None:
        import re
        from pathlib import Path

        from sdcflows.fieldmaps import get_identifier

        # Fallback to IntendedFor
        intended_rel = re.sub(r"^sub-[a-zA-Z0-9]*/", "", str(Path(fname).relative_to(layout.root)))
        field_source = get_identifier(intended_rel)

    return field_source