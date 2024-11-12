"""
Connectivity
^^^^^^^^^^^^

.. autofunction:: init_connectivity_wf

"""

import typing as ty

from nipype import Function
from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe

from ...interfaces import DerivativesDataSink
from ...utils.atlas import Atlases
from ...utils.pipelines import Pipelines


def init_connectivity_wf(
    mem_gb: float,
    atlases: Atlases,
    pipelines: Pipelines,
    session_level: bool,
    conn_measure: list = ["correlation"],
    name: str = "connectivity_wf",
):
    """
    Build a workflow to generate connectivity matrices.


    Parameters
    ----------
    mem_gb : :obj:`float`
        Size of BOLD file in GB - please note that this size
        should be calculated after resamplings that may extend
        the FoV
    conn_measure : :obj: `str`
        Measure to compute connectivity
    atlases : :py:class: `~fmristroke.utils.atlas.Atlases`
        A container for storing atlas for connectivity
    pipelines : :py:class:`~fmristroke.utils.pipelines.Pipelines`
        A container for storing denoising strategies.
    name : :obj:`str`
        Name of workflow (default: ``connectivity_wf``)

    Inputs
    ------
    bold_denoised
        BOLD image denoised in standard space
    boldmask
        fov mask
    mask_lesion_std
        Mask of the lesion in std space

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
    pipelines
        Pipelines
    atlases
        Atlases
    conn_measure
        Conn_measure
    """

    from niworkflows.engine.workflows import LiterateWorkflow as Workflow
    from niworkflows.interfaces.utility import KeySelect

    from ...interfaces.nilearn import Connectivity
    from ...interfaces.pyants import ApplyTransforms

    workflow = Workflow(name=name)

    workflow.__desc__ = f"""
    {'Then runs from the same session and task were concatenated and' if session_level else 'Then'} connectivity matrices were computed on the denoised BOLD signal using the following atlas : {atlases.get_atlases()} and the following
    connectivity measures : {conn_measure} using the Nilearn package  [@nilearn]."""

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "bold_denoised",
                "boldmask",
                "mask_lesion_std",
                "templates",
                "pipelines",
                "anat2std_xfm",
            ]
        ),
        name="inputnode",
    )

    iteratlas = pe.Node(
        niu.IdentityInterface(fields=["atlas"]),
        name="iteratlas",
    )
    iteratlas.iterables = [("atlas", atlases.atlases)]

    iterpipeline = pe.Node(
        niu.IdentityInterface(fields=["pipeline"]), name="iterpipe"
    )
    iterpipeline.iterables = [("pipeline", pipelines.get_pipelines())]

    iterconn = pe.Node(
        niu.IdentityInterface(fields=["conn_measure"]), name="iterconn"
    )
    iterconn.iterables = [("conn_measure", conn_measure)]

    # Get Atlas info
    atlas_info = pe.Node(
        niu.Function(
            function=_get_atlas_info,
            input_name=["in_atlas"],
            output_names=["atlas", "labels", "space", "mask_file"],
        ),
        name="atlas_info",
    )

    # Select correct space
    select_space = pe.Node(
        KeySelect(fields=["bold_denoised", "mask_lesion_std", "anat2std_xfm"]),
        name="select_space",
        run_without_submitting=True,
    )

    # Remove ROI voxels from Atlas ROI
    roi_resamp = pe.Node(
        niu.Function(function=_resample_to_img), name="roi_resamp"
    )
    roi_resamp.inputs.interpolation = "nearest"
    atlas_roi = pe.Node(niu.Function(function=_remove_roi), name="atlas_roi")

    # Bold mask in std space
    boldmask_tfm = pe.Node(
        ApplyTransforms(interpolation="multiLabel"), name="boldmask_std_tfm"
    )

    # atlas roi resamp
    atlas_roi_resamp = pe.Node(
        niu.Function(function=_resample_to_img), name="atlas_roi_resamp"
    )
    atlas_roi_resamp.inputs.interpolation = "nearest"

    # atlas resamp
    atlas_resamp = pe.Node(
        niu.Function(function=_resample_to_img), name="atlas_resamp"
    )
    atlas_resamp.inputs.interpolation = "nearest"

    # Compute connectivity
    select_pipeline = pe.Node(
        KeySelect(fields=["bold_denoised"]),
        name="select_pipeline",
        run_without_submitting=True,
    )
    connectivity_roi = pe.Node(Connectivity(), name="connectivity_roi")
    connectivity = pe.Node(Connectivity(), name="connectivity")

    # Compute FCC wo roi
    conn_FCC_roi = pe.Node(
        niu.Function(function=_get_FCC, output_names=["FCC"]), name="FCC_roi"
    )
    concat_FCC_roi = pe.Node(
        niu.Function(function=_concat_FCC, output_names=["FCC"]),
        name="concat_FCC_roi",
    )

    # Compute FCC
    conn_FCC = pe.Node(
        niu.Function(function=_get_FCC, output_names=["FCC"]), name="FCC"
    )
    concat_FCC = pe.Node(
        niu.Function(function=_concat_FCC, output_names=["FCC"]),
        name="concat_FCC",
    )

    # fmt:off
    workflow.connect(
        [
            (iteratlas, atlas_info, [("atlas", "in_atlas")]),
            (iterconn, connectivity_roi, [("conn_measure", "conn_measure")]),
            (iterconn, connectivity, [("conn_measure", "conn_measure")]),
            (inputnode, select_pipeline, [("bold_denoised", "bold_denoised"), ("pipelines", "keys")],),
            (iterpipeline, select_pipeline, [("pipeline", "key")]),
            (select_pipeline, select_space, [("bold_denoised", "bold_denoised")],),
            (inputnode, select_space, [("mask_lesion_std", "mask_lesion_std"),
                                    ("anat2std_xfm", "anat2std_xfm"),
                                    ("templates", "keys"),],),
            (atlas_info, select_space, [("space", "key")]),
            (select_space, roi_resamp, [("mask_lesion_std", "input_image")]),
            (atlas_info, roi_resamp, [("mask_file", "reference_image")]),
            (atlas_info, boldmask_tfm, [("mask_file", "reference_image")]),
            (select_space, boldmask_tfm, [("anat2std_xfm", "transforms")]),
            (inputnode, boldmask_tfm, [("boldmask", "input_image")]),
            (boldmask_tfm, connectivity_roi, [("output_image", "brain_mask")]),
            (boldmask_tfm, connectivity, [("output_image", "brain_mask")]),
            (roi_resamp, atlas_roi, [("out", "roi")]),
            (atlas_info, atlas_roi, [("mask_file", "in_atlas")]),
            (atlas_roi, atlas_roi_resamp, [("out", "input_image")]),
            (select_space, atlas_roi_resamp, [("bold_denoised", "reference_image")]),
            (atlas_roi_resamp, connectivity_roi, [("out", "atlas")]),
            (atlas_info, atlas_resamp, [("mask_file", "input_image")]),
            (select_space, atlas_resamp, [("bold_denoised", "reference_image")]),
            (atlas_resamp, connectivity, [("out", "atlas")]),
            (select_space, connectivity, [("bold_denoised", "input_image")]),
            (select_space, connectivity_roi, [("bold_denoised", "input_image")]),
            (connectivity, conn_FCC, [("output_conn", "conn_mat")]),
            (iterpipeline, conn_FCC, [("pipeline", "pipeline")]),
            (iteratlas, conn_FCC, [("atlas", "atlas")]),
            (iterconn, conn_FCC, [("conn_measure", "measure")]),
            (connectivity_roi, conn_FCC_roi, [("output_conn", "conn_mat")]),
            (iterpipeline, conn_FCC_roi, [("pipeline", "pipeline")]),
            (iteratlas, conn_FCC_roi, [("atlas", "atlas")]),
            (iterconn, conn_FCC_roi, [("conn_measure", "measure")]),
        ])
    # fmt:on

    # Handle outputs
    output_names = [
        "conn_mat",
        "FCC",
        "conn_mat_roi",
        "FCC_roi",
        "conn_measures",
        "pipelines",
        "atlases",
    ]
    output_connectivity = pe.Node(
        niu.IdentityInterface(
            fields=["conn_mat", "FCC", "conn_mat_roi", "FCC_roi"]
        ),
        name="output_connectivity",
    )
    output_measure = pe.JoinNode(
        niu.IdentityInterface(
            fields=[
                "conn_mat",
                "FCC",
                "conn_mat_roi",
                "FCC_roi",
                "conn_measures",
            ]
        ),
        joinfield=[
            "conn_mat",
            "FCC",
            "conn_mat_roi",
            "FCC_roi",
            "conn_measures",
        ],
        joinsource=iterconn,
        name="output_measure",
    )
    output_pipeline = pe.JoinNode(
        niu.IdentityInterface(
            fields=["conn_mat", "FCC", "conn_mat_roi", "FCC_roi", "pipelines"]
        ),
        joinfield=["conn_mat", "FCC", "conn_mat_roi", "FCC_roi", "pipelines"],
        joinsource=iterpipeline,
        name="output_pipeline",
    )
    output_atlas = pe.JoinNode(
        niu.IdentityInterface(
            fields=["conn_mat", "FCC", "conn_mat_roi", "FCC_roi", "atlases"]
        ),
        joinfield=["conn_mat", "FCC", "conn_mat_roi", "FCC_roi", "atlases"],
        joinsource=iteratlas,
        name="output_atlas",
    )
    outputnode = pe.Node(
        niu.IdentityInterface(fields=output_names), name="outputnode"
    )

    # fmt:off
    workflow.connect([
        (connectivity, output_connectivity, [("output_conn", "conn_mat")]),
        (conn_FCC, output_connectivity, [("FCC", "FCC")]),
        (connectivity_roi, output_connectivity, [("output_conn", "conn_mat_roi")]),
        (conn_FCC_roi, output_connectivity, [("FCC", "FCC_roi")]),
        (output_connectivity, output_measure, [("conn_mat", "conn_mat"),
                                                ("FCC", "FCC"),
                                                ("conn_mat_roi", "conn_mat_roi"),
                                                ("FCC_roi", "FCC_roi")]),
        (iterconn, output_measure, [("conn_measure", "conn_measures")]),
        (output_measure, output_pipeline, [("conn_mat", "conn_mat"),
                                            ("FCC", "FCC"),
                                            ("conn_mat_roi", "conn_mat_roi"),
                                            ("FCC_roi", "FCC_roi")]),
        (iterpipeline, output_pipeline, [("pipeline", "pipelines")]),
        (output_pipeline, output_atlas, [("conn_mat", "conn_mat"),
                                        ("FCC", "FCC"),
                                        ("conn_mat_roi", "conn_mat_roi"),
                                        ("FCC_roi", "FCC_roi")]),
        (iteratlas, output_atlas, [(("atlas", _get_atlas_name), "atlases")]),
        (output_atlas, concat_FCC, [("FCC", "in_FCC")]),
        (output_atlas, concat_FCC_roi, [("FCC_roi", "in_FCC")]),
        (concat_FCC, outputnode, [("FCC", "FCC")]),
        (concat_FCC_roi, outputnode, [("FCC", "FCC_roi")]),
        (output_atlas, outputnode, [("conn_mat", "conn_mat"),
                                    ("conn_mat_roi", "conn_mat_roi"),
                                    ("atlases", "atlases")]),
        (output_measure, outputnode, [("conn_measures", "conn_measures")]),
        (output_pipeline, outputnode, [("pipelines", "pipelines")]),
    ])

    return workflow


def init_lesion_voxels_conn_wf(
    mem_gb: float,
    omp_nthreads: int,
    pipelines: Pipelines,
    name: str = "bold_lesion_conn_wf",
):
    """
    Build a workflow to compute lesion to voxels connectivity.

    This workflow computes lesion to voxels connectivity according to the specified denoising strategies.

    Parameters
    ----------
    mem_gb : :obj:`float`
        Size of BOLD file in GB - please note that this size
        should be calculated after resamplings that may extend
        the FoV
    pipelines : :py:class:`~fmristroke.utils.pipelines.Pipelines`
        A container for storing denoising strategies.
    name : :obj:`str`
        Name of workflow (default: ``bold_lesion_conn_wf``)

    Inputs
    ------
    bold_denoised_t1
        BOLD image in T1 space, after the prescribed corrections (STC, HMC and SDC)
        and denoising.
    boldmask
        Bold fov mask in t1 space
    t1w_mask_lesion
        Mask of the lesion in T1w space
    gm_mask
        Gray matter mask
    t1w
        preproc T1w image


    Outputs
    -------
    output_conn_roi
        Lesion to voxels connectivity

    """
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow
    from niworkflows.interfaces.images import SignalExtraction
    from niworkflows.interfaces.utility import KeySelect

    from ...interfaces.nilearn import ROI2VoxelConnectivity
    from ...interfaces.reports import ROIConnPlot

    workflow = Workflow(name=name)
    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "bold_denoised_t1",
                "t1w_mask_lesion",
                "pipelines",
                "t1w",
                "gm_mask",
            ]
        ),
        name="inputnode",
    )

    iterpipeline = pe.Node(
        niu.IdentityInterface(fields=["pipeline"]), name="iterpipe"
    )
    iterpipeline.iterables = [("pipeline", pipelines.get_pipelines())]

    # Select pipeline
    select_pipeline = pe.Node(
        KeySelect(fields=["bold_denoised_t1"]),
        name="select_pipeline",
        run_without_submitting=True,
    )

    # ROI resamp
    select_tmp = pe.Node(niu.Select(index=0), name="select_resamp_temp")
    roi_resamp = pe.Node(
        niu.Function(function=_resample_to_img), name="roi_resamp"
    )
    roi_resamp.inputs.interpolation = "nearest"

    # GM Mask resamp
    gm_resamp = pe.Node(
        niu.Function(function=_resample_to_img), name="gm_resamp"
    )
    gm_resamp.inputs.interpolation = "nearest"

    # Extract ROI signal
    signal_roi = pe.Node(
        SignalExtraction(class_labels=["lesion_roi"]),
        name="signal_roi",
        mem_gb=mem_gb,
    )

    # ROI to voxels connectivity
    roi2voxelconn = pe.Node(ROI2VoxelConnectivity(), name="roi2voxelconn")

    # fmt:off
    workflow.connect([
        (inputnode, select_pipeline, [("bold_denoised_t1", "bold_denoised_t1"),
                                    ("pipelines", "keys")]),
        (iterpipeline, select_pipeline, [("pipeline", "key")]),
        (inputnode, roi_resamp, [("t1w_mask_lesion", "input_image")]),
        (inputnode, select_tmp, [("bold_denoised_t1", "inlist")]),
        (select_tmp, roi_resamp, [("out", "reference_image")]),
        (roi_resamp, signal_roi, [("out", "label_files")]),
        (select_pipeline, signal_roi, [("bold_denoised_t1", "in_file")]),
        (select_pipeline, roi2voxelconn, [("bold_denoised_t1", "input_image")]),
        (signal_roi, roi2voxelconn, [("out_file", "roi_ts")]),
        (inputnode, gm_resamp, [("gm_mask", "input_image")]),
        (select_pipeline, gm_resamp, [("bold_denoised_t1", "reference_image")]),
        (gm_resamp, roi2voxelconn, [("out", "brain_mask")]),
        (iterpipeline, roi2voxelconn, [("pipeline", "pipeline_name")]),
    ])
    # fmt:on
    output_names = ["pipeline", "output_conn", "output_img"]
    output_corr = pe.Node(
        niu.IdentityInterface(fields=output_names), name="output_corr"
    )
    join_out = pe.JoinNode(
        niu.IdentityInterface(
            fields=output_names,
        ),
        joinfield=output_names,
        joinsource=iterpipeline,
        name="join_conn",
    )

    concat_lesion_corr = pe.Node(
        niu.Function(function=_concat_tsv), name="concat_pipelines"
    )

    # Reporting
    lesion_conn_plot = pe.Node(
        ROIConnPlot(generate_report=True), name="lesion_conn_plot"
    )
    ds_report_lesioncorr = pe.Node(
        DerivativesDataSink(
            desc="lesioncorr",
            datatype="figures",
        ),
        run_without_submitting=True,
        name="ds_report_lesioncorr",
    )

    outputnode = pe.Node(
        niu.IdentityInterface(fields=["output_conn_roi"]), name="outputnode"
    )

    # fmt:off
    workflow.connect([
        (iterpipeline, output_corr, [("pipeline", "pipeline")]),
        (roi2voxelconn, output_corr, [("output_conn", "output_conn"),
                                    ("output_img", "output_img")]),
        (output_corr, join_out, [(f, f) for f in output_names]),
        (join_out, concat_lesion_corr, [("output_conn", "in_tsvs")]),
        (join_out, lesion_conn_plot, [("output_img", "input_images"),
                                      ("pipeline", "pipeline")]),
        (inputnode, lesion_conn_plot, [("t1w", "anat_img"),
                                        ("t1w_mask_lesion", "roi")]),
        (lesion_conn_plot, ds_report_lesioncorr, [("out_report", "in_file")]),
        (concat_lesion_corr, outputnode, [("out", "output_conn_roi")])
    ])
    # fmt:on

    return workflow


def _resample_to_img(input_image, reference_image, interpolation):
    from pathlib import Path

    from nilearn.image import resample_to_img

    new_img = resample_to_img(
        input_image, reference_image, interpolation=interpolation
    )
    out_name = Path("mask_resamp.nii.gz").absolute()
    new_img.to_filename(out_name)
    return str(out_name)


def _remove_roi(in_atlas, roi):
    from pathlib import Path

    import nibabel as nib

    roi = nib.load(roi).get_fdata()
    in_atlas = nib.load(in_atlas)
    new_at = in_atlas.get_fdata()
    new_at[roi == 1] = 0
    out_atlas = nib.Nifti1Image(new_at, in_atlas.affine)
    out_name = Path("atlas_roi.nii.gz").absolute()
    out_atlas.to_filename(out_name)
    return str(out_name)


def _get_atlas_info(in_atlas):
    atlas = in_atlas.atlas
    labels = in_atlas.labels
    space = in_atlas.space
    mask_file = in_atlas.mask_file
    return atlas, labels, space, mask_file


def _get_atlas_name(in_atlas):
    return in_atlas.atlas


def _get_FCC(conn_mat, atlas, measure, pipeline):
    import json
    from pathlib import Path

    import numpy as np

    from fmristroke.utils.metrics import compute_FCC

    conn_mat = np.load(conn_mat)
    FCC, _ = compute_FCC(conn_mat, atlas.labels)

    dic = {
        "pipeline": pipeline,
        "atlas": atlas.atlas,
        "measure": measure,
        "FCC": FCC,
    }
    json_object = json.dumps(dic, indent=4)
    out_name = Path("conn_FCC.json").absolute()

    with open(out_name, "w") as outfile:
        outfile.write(json_object)
    return str(out_name)
    # return dic


def _concat_FCC(in_FCC):
    import json
    from pathlib import Path

    import pandas as pd

    out = Path("FCC.tsv").absolute()

    flatten_FCC = [
        item for atlas in in_FCC for pipeline in atlas for item in pipeline
    ]
    jsons = []
    for FCC in flatten_FCC:
        with open(FCC, "r") as in_file:
            jsons.append(json.load(in_file))
    df = pd.DataFrame(jsons)

    df.to_csv(out, sep="\t", na_rep="n/a")
    return str(out)


def _concat_tsv(in_tsvs):
    from pathlib import Path

    import pandas as pd

    out = Path("lesion_conn.tsv").absolute()

    final_df = pd.DataFrame()
    for tsv in in_tsvs:
        df = pd.read_table(tsv, sep="\t")
        final_df = pd.concat((final_df, df), axis=1)
    final_df.to_csv(out, sep="\t", na_rep="n/a")
    return str(out)


def _get_pipeline_name(in_pipeline):
    return in_pipeline.pipeline
