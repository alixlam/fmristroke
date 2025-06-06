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
from .outputs import (
    init_connectivity_derivatives_wf,
    init_func_lesion_derivatives_wf,
)
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


    """
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow
    from niworkflows.interfaces.nibabel import ApplyMask
    from niworkflows.interfaces.reportlets.registration import (
        SimpleBeforeAfterRPT as SimpleBeforeAfter,
    )
    from niworkflows.interfaces.utility import DictMerge, KeySelect

    img = nb.load(
        bold_file[0] if isinstance(bold_file, (list, tuple)) else bold_file
    )
    nvols = 1 if img.ndim < 4 else img.shape[3]
    if (
        nvols * img.header["pixdim"][4] / 2 < config.workflow.maxlag
        and not config.execution.run_sessionlevel
    ):
        config.loggers.workflow.warning(
            f"Too short BOLD series with chosen maxlag (<= {config.workflow.maxlag} timepoints). Skipping processing of <{bold_file}>. To add this bold serie, consider decreasing the maxlag parameter or use the --session-level argument if you have multiple runs from the same session/task."
        )
        return

    mem_gb = {"filesize": 1, "resampled": 1, "largemem": 1}
    bold_tlen = 10

    # Have some options handy
    omp_nthreads = config.nipype.omp_nthreads
    spaces = config.workflow.spaces
    pipelines = config.workflow.pipelines
    output_dir = str(config.execution.output_dir)
    freesurfer = config.workflow.freesurfer
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
                "boldmask",
                "template",
                "spatial_reference",
                "bold_denoised_t1",
                "bold_denoised_std",
                "pipelines",
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

    # fmt:off
    workflow.connect([
        (inputnode, func_lesion_derivatives_wf, [("bold_t1", "inputnode.source_file")],)
    ])
    # fmt:on

    # fmt:off
    workflow.connect([
        (outputnode, func_lesion_derivatives_wf, [
            ("confounds_file", "inputnode.confounds_file"),
            ("confounds_metadata", "inputnode.confounds_metadata"),
            ("template", "inputnode.template"),
            ("spatial_reference", "inputnode.spatial_reference"),
            ("bold_denoised_t1", "inputnode.bold_denoised_t1"),
            ("bold_denoised_std", "inputnode.bold_denoised_std"),
            ("pipelines", "inputnode.pipelines"),],)
    ])
    # fmt:on

    # LESION CONFOUNDS  #######################################################
    lesion_confounds_wf = init_confs_wf(
        mem_gb=mem_gb["largemem"],
        metadata=metadata,
        freesurfer=freesurfer,
        ncomp_method=config.workflow.ncomp_method,
        ica_method=config.workflow.ica_method,
        name="lesion_confounds_wf",
    )

    # fmt:off
    workflow.connect([
        # Connect bold_confounds_wf
        (inputnode, lesion_confounds_wf, [
            ("t1w_tpms", "inputnode.t1w_tpms"),
            ("t1w_mask", "inputnode.t1w_mask"),
            ("roi", "inputnode.roi"),
            ("bold_t1", "inputnode.bold_t1"),
            ("boldref_t1", "inputnode.boldref_t1"),
            ("boldmask_t1", "inputnode.boldmask_t1"),
            ("confounds_file", "inputnode.confounds_file"),],),
        (lesion_confounds_wf, outputnode, [
            ("outputnode.boldmask", "boldmask"),
            ("outputnode.confounds_file", "confounds_file"),
            ("outputnode.confounds_metadata", "confounds_metadata"),],),
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
            ("templates", "inputnode.templates"),
            ("bold_t1", "inputnode.bold_t1"),
            ("confounds_metadata", "inputnode.confounds_metadata")],),
        (lesion_confounds_wf, denoising_wf, [
            ("outputnode.boldmask", "inputnode.boldmask"),
            ("outputnode.confounds_file", "inputnode.confounds_file")],),
        (denoising_wf, outputnode, [
            ("outputnode.template", "template"),
            ("outputnode.spatial_reference", "spatial_reference"),
            ("outputnode.bold_denoised_t1", "bold_denoised_t1"),
            ("outputnode.bold_denoised_std", "bold_denoised_std"),
            ("outputnode.pipeline", "pipelines"),],),
    ])
    # fmt:on

    """# CARPET PLOT DENOISED ##########################################
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
        (denoising_wf, carpetplot_wf, [
            ("outputnode.bold_denoised_t1", "inputnode.bold"),],),
        (lesion_confounds_wf, carpetplot_wf, [
            ("outputnode.boldmask", "inputnode.bold_mask"),],),
        (inputnode, carpetplot_wf, [
            ("confounds_file", "inputnode.confounds_file",),])
    ])
    # fmt:on"""

    # REPORTING ############################################################

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


def init_lesion_connectivity_wf(bold_file):
    """
    This workflow controls the connectivity specific stages *fMRIStroke*.

    Parameters
    ----------
    bold_file
        Path to NIfTI file


    Inputs
    ------
    bold_t1
        BOLD series NIFTI file in anatomical space
    bold_denoised_t1
        BOLD denoised series NIfTI file in anatomical space
    bold_denoised_std
        BOLD denoised series in standard space
    boldmask_t1
        BOLD mask in anatomical space
    templates
        List of templates that were applied as targets during
        spatial normalization.
    pipelines
        List of denoising pipelines
    anat2std_xfm
        transform from anatomical to standard space
    roi
        ROI mask
    gm_mask
        Gray-matter mask in T1w space

    Outputs
    -------
    conn_mat
        connectivity matrice
    pipelines
        Pipelines
    atlases
        Atlases
    conn_measure
        Conn_measure
    """
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow
    from niworkflows.interfaces.nibabel import ApplyMask
    from niworkflows.interfaces.reportlets.registration import (
        SimpleBeforeAfterRPT as SimpleBeforeAfter,
    )
    from niworkflows.interfaces.utility import DictMerge, KeySelect

    from ...utils.combine_runs import combine_run_source

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

    wf_name = _get_wf_name(ref_file, conn=True)

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

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "bold_t1",
                "bold_denoised_t1",
                "bold_denoised_std",
                "boldmask_t1",
                "t1w_preproc",
                "t1w_mask",
                "t1w_tpms",
                "t1w_aseg",
                "templates",
                "pipelines",
                "std2anat_xfm",
                "anat2std_xfm",
                "roi",
                "roi_std",
                "gm_mask",
            ]
        ),
        name="inputnode",
    )
    inputnode.inputs.bold_t1 = bold_file

    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "lagmaps",
                "atlases",
                "pipelines",
                "conn_measures",
                "conn_mat",
                "conn_mat_roi",
                "FCC",
                "FCC_roi",
                "lesion_conn",
            ]
        ),
        name="outputnode",
    )

    conn_derivatives_wf = init_connectivity_derivatives_wf(
        bids_root=layout.root,
        output_dir=output_dir,
        spaces=spaces,
        pipelines=pipelines,
        conn_measure=conn_measure,
        atlases=atlases,
    )

    conn_derivatives_wf.inputs.inputnode.all_source_files = bold_file

    if not session_level:
        # fmt:off
        workflow.connect([
            (inputnode, conn_derivatives_wf, [("bold_t1", "inputnode.source_file")],)
        ])
        # fmt:on

    else:
        # fmt:off
        workflow.connect([
            (inputnode, conn_derivatives_wf, [
                (("bold_t1", combine_run_source), "inputnode.source_file",)],)
        ])
        # fmt:on

    # fmt:off
    workflow.connect([
        (outputnode, conn_derivatives_wf, [
            ("lagmaps", "inputnode.lagmaps"),
            ("pipelines", "inputnode.pipelines"),
            ("atlases", "inputnode.atlases"),
            ("conn_mat", "inputnode.conn_mat"),
            ("conn_measures", "inputnode.conn_measures"),
            ("lesion_conn", "inputnode.lesion_conn"),
            ("FCC", "inputnode.FCC"),
            ("FCC_roi", "inputnode.FCC_roi"),
            ("conn_mat_roi", "inputnode.conn_mat_roi")],),
    ])
    # fmt:on

    # COMBINE RUNS  ##########################################################
    concat_wf = init_concat_wf(
        mem_gb=mem_gb["largemem"],
        croprun=croprun,
        session_level=session_level,
        name="concat_runs_wf",
    )

    # fmt:off
    workflow.connect([
        (inputnode, concat_wf, [
            ("boldmask_t1", "inputnode.boldmask_t1"),
            ("bold_denoised_t1", "inputnode.bold_denoised_t1"),
            ("bold_denoised_std", "inputnode.bold_denoised_std"),])
    ])
    # fmt:on

    # HEMODYNAMIC WORKFLOW   #################################################
    hemodynamics_wf = init_hemodynamic_wf(
        mem_gb=mem_gb["largemem"],
        metadata=metadata,
        maxlag=config.workflow.maxlag,
        name="hemodynamic_wf",
    )

    # select pipeline lag
    select_lag_pipeline = pe.Node(
        KeySelect(fields=["bold_denoised_t1"], key="Minimal"),
        name="select_pipeline_lag",
        run_without_submitting=True,
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
            ("outputnode.boldmask_t1", "inputnode.boldmask")]),
        (concat_wf, select_lag_pipeline, [
            ("outputnode.bold_denoised_t1", "bold_denoised_t1")]),
        (inputnode, select_lag_pipeline, [
            ("pipelines", "keys")]),
        (select_lag_pipeline, hemodynamics_wf, [
            ("bold_denoised_t1", "inputnode.bold_t1"),],),
        (hemodynamics_wf, outputnode, [
            ("outputnode.lagmaps", "lagmaps"),],),
    ])
    # fmt:on

    # CONNECTIVITY WF ########################################################
    connectivity_wf = init_connectivity_wf(
        mem_gb=mem_gb["largemem"],
        atlases=atlases,
        pipelines=pipelines,
        session_level=session_level,
        conn_measure=conn_measure,
        name="connectivity_wf",
    )

    # fmt:off
    workflow.connect([
        (inputnode, connectivity_wf, [
            ("templates", "inputnode.templates"),
            ("pipelines", "inputnode.pipelines"),
            ("anat2std_xfm", "inputnode.anat2std_xfm"),
            ("roi_std", "inputnode.mask_lesion_std"),],),
        (concat_wf, connectivity_wf, [
            ("outputnode.boldmask_t1", "inputnode.boldmask"),
            ("outputnode.bold_denoised_std", "inputnode.bold_denoised",),]),
        (connectivity_wf, outputnode, [
            ("outputnode.pipelines", "pipelines"),
            ("outputnode.atlases", "atlases"),
            ("outputnode.conn_measures", "conn_measures"),
            ("outputnode.conn_mat", "conn_mat"),
            ("outputnode.conn_mat_roi", "conn_mat_roi"),
            ("outputnode.FCC", "FCC"),
            ("outputnode.FCC_roi", "FCC_roi")])
    ])
    # fmt:on

    # LESION CONNECTIVITY ##################################################
    lesion_conn = init_lesion_voxels_conn_wf(
        mem_gb=mem_gb["largemem"],
        omp_nthreads=omp_nthreads,
        pipelines=pipelines,
        name="lesion_conn_wf",
    )

    # fmt:off
    workflow.connect([
        (inputnode, lesion_conn, [
            ("pipelines", "inputnode.pipelines"),
            ("roi", "inputnode.t1w_mask_lesion"),
            ("t1w_preproc", "inputnode.t1w")],),
        (concat_wf, lesion_conn, [
            ("outputnode.bold_denoised_t1", "inputnode.bold_denoised_t1"),]),
        (hemodynamics_wf, lesion_conn, [
            ("outputnode.gm_mask", "inputnode.gm_mask")]),
        (lesion_conn, outputnode, [
            ("outputnode.output_conn_roi", "lesion_conn"),]),
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


def _get_wf_name(bold_fname, conn=False):
    """
    Derive the workflow name for supplied BOLD file.

    """
    from nipype.utils.filemanip import split_filename

    fname = split_filename(bold_fname)[1]
    fname_nosub = "_".join(fname.split("_")[1:])

    if conn:
        name = "lesion_conn_" + fname_nosub.replace(".", "_").replace(
            " ", ""
        ).replace("-", "_").replace("_bold", "_wf")
    else:
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
