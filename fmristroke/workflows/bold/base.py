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
from .concatenate import init_concat_wf

# BOLD workflows
from .confounds import init_carpetplot_wf, init_confs_wf
from .connectivity import init_connectivity_wf, init_lesion_voxels_conn_wf
from .denoise import init_denoise_wf
from .lagmaps import init_hemodynamic_wf
from .outputs import init_func_lesion_derivatives_wf
from .registration import init_lesionplot_wf


def init_lesion_preproc_wf(
    bold_file, boldref_file, boldmask_file, confounds_file, confounds_metadata
):
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

    from ...utils.combine_runs import combine_run_source

    img = nb.load(
        bold_file[0] if isinstance(bold_file, (list, tuple)) else bold_file
    )
    nvols = 1 if img.ndim < 4 else img.shape[3]

    mem_gb = {"filesize": 1, "resampled": 1, "largemem": 1}
    bold_tlen = 10

    # Have some options handy
    omp_nthreads = config.nipype.omp_nthreads
    spaces = config.workflow.spaces
    pipelines = config.workflow.pipelines
    output_dir = str(config.execution.output_dir)
    freesurfer = config.workflow.freesurfer
    croprun = config.workflow.croprun
    session_level = config.execution.run_sessionlevel
    atlases = config.workflow.atlases
    conn_measure = config.workflow.conn_measure

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
                "roi_std",
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
                "template",
                "spatial_reference",
                "denoised_bold_t1",
                "denoised_bold_std",
                "pipeline",
                "atlases",
                "pipelines",
                "conn_measures",
                "conn_mat",
                "FCC",
                "lesion_conn",
            ]
        ),
        name="outputnode",
    )
    func_lesion_derivatives_wf = init_func_lesion_derivatives_wf(
        bids_root=layout.root,
        output_dir=output_dir,
        spaces=spaces,
        pipelines=pipelines,
        conn_measure=conn_measure,
        atlases=atlases,
    )

    func_lesion_derivatives_wf.inputs.inputnode.all_source_files = bold_file

    if not session_level:
        # fmt:off
        workflow.connect([
            (inputnode, func_lesion_derivatives_wf, [("bold_t1", "inputnode.source_file")],)
        ])
        # fmt:on

    else:
        # fmt:off
        workflow.connect([
            (inputnode, func_lesion_derivatives_wf, [
                (("bold_t1", combine_run_source), "inputnode.source_file",)],)
        ])
        # fmt:on

    # fmt:off
    workflow.connect([
        (outputnode, func_lesion_derivatives_wf, [
            ("confounds_file", "inputnode.confounds_file"),
            ("confounds_metadata", "inputnode.confounds_metadata"),
            ("lagmaps", "inputnode.lagmaps"),
            ("template", "inputnode.template"),
            ("spatial_reference", "inputnode.spatial_reference"),
            ("denoised_bold_t1", "inputnode.denoised_bold_t1"),
            ("denoised_bold_std", "inputnode.denoised_bold_std"),
            ("pipeline", "inputnode.pipeline"),
            ("pipelines", "inputnode.pipelines"),
            ("atlases", "inputnode.atlases"),
            ("conn_mat", "inputnode.conn_mat"),
            ("conn_measures", "inputnode.conn_measures"),
            ("lesion_conn", "inputnode.lesion_conn"),
            ("FCC", "inputnode.FCC")],),
    ])
    # fmt:on

    # COMBINE RUNS  ##########################################################
    concat_wf = init_concat_wf(
        mem_gb=mem_gb["largemem"],
        croprun=croprun,
        name="concat_runs_wf",
    )
    # fmt:off
    workflow.connect([
        (inputnode, concat_wf, [
            ("bold_t1", "inputnode.bold_t1"),
            ("boldref_t1", "inputnode.boldref_t1"),
            ("boldmask_t1", "inputnode.boldmask_t1"),
            ("confounds_file", "inputnode.confounds_file"),
            ("confounds_metadata", "inputnode.confounds_metadata"),],)
    ])
    # fmt:on

    # LESION CONFOUNDS  #######################################################
    lesion_confounds_wf = init_confs_wf(
        mem_gb=mem_gb["largemem"],
        metadata=metadata,
        freesurfer=freesurfer,
        ncomp_method=config.workflow.ncomp_method,
        ica_method=config.workflow.ica_method,
        session_level=session_level,
        name="lesion_confounds_wf",
    )

    # fmt:off
    workflow.connect([
        # Connect bold_confounds_wf
        (inputnode, lesion_confounds_wf, [
            ("t1w_tpms", "inputnode.t1w_tpms"),
            ("t1w_mask", "inputnode.t1w_mask"),
            ("roi", "inputnode.roi"),],),
        (concat_wf, lesion_confounds_wf, [
            ("outputnode.bold_t1", "inputnode.bold_t1"),
            ("outputnode.boldref_t1", "inputnode.boldref_t1"),
            ("outputnode.boldmask_t1", "inputnode.boldmask_t1"),
            ("outputnode.confounds_file", "inputnode.confounds_file"),],),
        (lesion_confounds_wf, outputnode, [
            ("outputnode.confounds_file", "confounds_file"),
            ("outputnode.confounds_metadata", "confounds_metadata"),],),
    ])
    # fmt:on

    # HEMODYNAMIC WORKFLOW   #################################################
    hemodynamics_wf = init_hemodynamic_wf(
        mem_gb=mem_gb["largemem"],
        metadata=metadata,
        maxlag=config.workflow.maxlag,
        name="hemodynamic_wf",
    )
    # fmt:off
    workflow.connect([
        # Connect bold_confounds_wf
        (inputnode, hemodynamics_wf, [
            ("roi", "inputnode.roi"),
            ("t1w_preproc", "inputnode.t1w_preproc"),
            ("t1w_mask", "inputnode.t1w_mask"),
            ("t1w_tpms", "inputnode.t1w_tpms"),
            ("std2anat_xfm", "inputnode.std2anat_xfm"),
            ("templates", "inputnode.template"),
            # undefined if freesurfer was not run
            ("t1w_aseg", "inputnode.t1w_aseg"),],),
        (concat_wf, hemodynamics_wf, [
            ("outputnode.bold_t1", "inputnode.bold_t1"),
            ("outputnode.confounds_file", "inputnode.confounds_file"),],),
        (lesion_confounds_wf, hemodynamics_wf, [
            ("outputnode.boldmask", "inputnode.boldmask")],),
        (hemodynamics_wf, outputnode, [
            ("outputnode.lagmaps", "lagmaps"),],),
    ])
    # fmt:on

    # DENOISING WORKFLOW #####################################################
    denoising_wf = init_denoise_wf(
        mem_gb=mem_gb["largemem"],
        omp_nthreads=omp_nthreads,
        metadata=metadata,
        pipelines=pipelines,
        spaces=spaces,
        name="bold_denoise_wf",
    )

    # fmt:off
    workflow.connect([
        (inputnode, denoising_wf, [
            ("anat2std_xfm", "inputnode.anat2std_xfm"),
            ("templates", "inputnode.templates"),],),
        (lesion_confounds_wf, denoising_wf, [
            ("outputnode.confounds_file", "inputnode.confounds_file")],),
        (concat_wf, denoising_wf, [
            ("outputnode.bold_t1", "inputnode.bold_t1"),],),
        (denoising_wf, outputnode, [
            ("outputnode.template", "template"),
            ("outputnode.spatial_reference", "spatial_reference"),
            ("outputnode.denoised_bold_t1", "denoised_bold_t1"),
            ("outputnode.denoised_bold_std", "denoised_bold_std"),
            ("outputnode.pipeline", "pipeline"),],),
    ])
    # fmt:on

    if session_level:
        # Get confounds descriptor from newly computed compcor confounds
        # fmt:off
        workflow.connect([
            (lesion_confounds_wf, denoising_wf, [
                ("outputnode.confounds_metadata_comp", "inputnode.confounds_metadata",)],),
        ])
        # fmt:on
    else:
        # fmt:off
        workflow.connect([
            (concat_wf, denoising_wf, [
                ("outputnode.confounds_metadata", "inputnode.confounds_metadata",)],),
        ])
        # fmt:on

    # CARPET PLOT IF SESSION LEVEL ##########################################
    if session_level:

        def _last(inlist):
            return inlist[-1]

        carpetplot_wf = init_carpetplot_wf(
            mem_gb=mem_gb["resampled"], metadata=metadata, name="carpetplot_wf"
        )
        carpetplot_select_std = pe.Node(
            KeySelect(fields=["std2anat_xfm"], key="MNI152NLin2009cAsym"),
            name="carpetplot_select_std",
            run_without_submitting=True,
        )
        # fmt:off
        workflow.connect([
            (inputnode, carpetplot_select_std, [
                ("std2anat_xfm", "std2anat_xfm"), ("templates", "keys")],),
            (carpetplot_select_std, carpetplot_wf, [
                ("std2anat_xfm", "inputnode.std2anat_xfm"),],),
            (concat_wf, carpetplot_wf, [
                ("outputnode.bold_t1", "inputnode.bold"),],),
            (lesion_confounds_wf, carpetplot_wf, [
                ("outputnode.boldmask", "inputnode.bold_mask"),],),
            (lesion_confounds_wf, carpetplot_wf, [
                ("outputnode.confounds_file", "inputnode.confounds_file",),
                ("outputnode.crown_mask", "inputnode.crown_mask"),
                (("outputnode.acompcor_masks", _last), "inputnode.acompcor_mask", ),],),
        ])
        # fmt:on

    # CONNECTIVITY WF ########################################################
    connectivity_wf = init_connectivity_wf(
        mem_gb=mem_gb["largemem"],
        atlases=atlases,
        pipelines=pipelines,
        conn_measure=conn_measure,
        name="connectivity_wf",
    )

    # fmt:off
    workflow.connect([
        (denoising_wf, connectivity_wf, [
            ("outputnode.template", "inputnode.templates"),
            ("outputnode.denoised_bold_std", "inputnode.bold_denoised",),
            ("outputnode.pipeline", "inputnode.pipelines"),],),
        (lesion_confounds_wf, connectivity_wf, [
            ("outputnode.boldmask", "inputnode.boldmask")],),
        (inputnode, connectivity_wf, [
            ("anat2std_xfm", "inputnode.anat2std_xfm"),
            ("roi_std", "inputnode.mask_lesion_std"),],),
        (connectivity_wf, outputnode, [
            ("outputnode.pipelines", "pipelines"),
            ("outputnode.atlases", "atlases"),
            ("outputnode.conn_measures", "conn_measures"),
            ("outputnode.conn_mat", "conn_mat"),
            ("outputnode.FCC", "FCC")])
    ])
    # fmt:on

    # LESION CONNECTIVITY ##################################################
    lesion_conn = init_lesion_voxels_conn_wf(
        mem_gb=mem_gb["largemem"],
        omp_nthreads=omp_nthreads,
        pipelines=pipelines,
        name="lesion_connectivity_wf",
    )

    # fmt:off
    workflow.connect([
        (denoising_wf, lesion_conn, [
            ("outputnode.denoised_bold_t1", "inputnode.denoised_bold_t1",),
            ("outputnode.pipeline", "inputnode.pipeline"),],),
        (inputnode, lesion_conn, [
            ("roi", "inputnode.t1w_mask_lesion"),
            ("t1w_preproc", "inputnode.t1w")],),
        (lesion_conn, outputnode, [
            ("outputnode.output_conn_roi", "lesion_conn"),])
    ])
    # fmt:on

    # REPORTING ############################################################
    if session_level:
        ref_file = combine_run_source(bold_file)

    # Fill-in datasinks of reportlets seen so far
    for node in workflow.list_node_names():
        if node.split(".")[-1].startswith("ds_report"):
            workflow.get_node(node).inputs.base_directory = output_dir
            workflow.get_node(node).inputs.source_file = ref_file

    # REGISTRATION PLOT WORKFLOW   ###########################################
    lesion_plot = init_lesionplot_wf(
        listify(boldref_file),
        mem_gb=mem_gb["largemem"],
        output_dir=output_dir,
        name="lesion_plot_wf",
    )
    # fmt:off
    workflow.connect([
        (inputnode, lesion_plot, [
            ("t1w_preproc", "inputnode.t1w"),
            ("roi", "inputnode.t1w_roi"),
            ("t1w_mask", "inputnode.t1w_mask"),],)
    ])
    # fmt:on

    return workflow


def _create_mem_gb(bold_fname):
    img = nb.load(bold_fname)
    nvox = int(np.prod(img.shape, dtype="u8"))
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

    """
    from nipype.utils.filemanip import split_filename

    fname = split_filename(bold_fname)[1]
    fname_nosub = "_".join(fname.split("_")[1:])
    name = "lesion_preproc_" + fname_nosub.replace(".", "_").replace(
        " ", ""
    ).replace("-", "_").replace("_bold", "_wf")

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
        ev_pair
        for f in listify(file_list)
        for ev_pair in parse_file_entities(f).items()
    ]:
        entities[e].append(v)

    def _unique(inlist):
        inlist = sorted(set(inlist))
        if len(inlist) == 1:
            return inlist[0]
        return inlist

    return {k: _unique(v) for k, v in entities.items()}
