"""
Denoising
^^^^^^^^^^

.. autofunction:: init_denoise_wf

"""

from __future__ import annotations

import typing as ty

from nipype import Function
from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from niworkflows.interfaces.fixes import (
    FixHeaderApplyTransforms as ApplyTransforms,
)

from ...config import DEFAULT_MEMORY_MIN_GB

if ty.TYPE_CHECKING:
    from niworkflows.utils.spaces import SpatialReferences

    from ..utils.pipelines import Pipelines


def init_denoise_wf(
    mem_gb: float,
    omp_nthreads: int,
    metadata: dict,
    pipelines: Pipelines,
    spaces: SpatialReferences,
    name: str = "bold_denoise_wf",
):
    """
    Build a workflow to denoise BOLD series.

    This workflow denoises the BOLD series according to the specified denoising strategies.

    Parameters
    ----------
    mem_gb : :obj:`float`
        Size of BOLD file in GB - please note that this size
        should be calculated after resamplings that may extend
        the FoV
    metadata : :obj:`dict`
        BIDS metadata for BOLD file
    pipelines : :py:class:`~fmristroke.utils.pipelines.Pipelines`
        A container for storing denoising strategies.
    spaces: py:class:`~niworkflows.utils.spaces.SpatialReferences`
        A container for storing, organizing, and parsing spatial normalizations. Composed of
        :py:class:`~niworkflows.utils.spaces.Reference` objects representing spatial references.
        Each ``Reference`` contains a space, which is a string of either TemplateFlow template IDs
        (e.g., ``MNI152Lin``, ``MNI152NLin6Asym``, ``MNIPediatricAsym``), nonstandard references
        (e.g., ``T1w`` or ``anat``, ``sbref``, ``run``, etc.), or a custom template located in
        the TemplateFlow root directory. Each ``Reference`` may also contain a spec, which is a
        dictionary with template specifications (e.g., a specification of ``{"resolution": 2}``
        would lead to resampling on a 2mm resolution of the space).
    name : :obj:`str`
        Name of workflow (default: ``bold_denoise_wf``)

    Inputs
    ------
    bold_t1
        BOLD image in T1 space, after the prescribed corrections (STC, HMC and SDC)
        when available.
    boldmask
        Bold fov mask in t1 space
    anat2std_xfm
        List of anatomical-to-standard space transforms generated during
        spatial normalization.
    templates
        List of templates that were applied as targets during
        spatial normalization.
    confounds_file
        TSV of all aggregated confounds.
    confounds_metadata
        JSON of confounds metadata file

    Outputs
    -------
    bold_denoised_t1
        BOLD image in T1 denoised
    bold_denoised_std
        BOLD image in std spaces denoised
    template
        Template identifiers synchronized correspondingly to previously
        described outputs.
    spatial reference
        Spatial reference
    pipeline
        Pipeline identifyers

    """

    from niworkflows.engine.workflows import LiterateWorkflow as Workflow

    from ...interfaces.confounds import SelectConfounds
    from ...interfaces.nilearn import Denoise
    from .resample import init_bold_std_trans_wf

    workflow = Workflow(name=name)

    workflow.__desc__ = f"""
    The BOLD time-series were denoised using the following strategies \n:
    """
    for pipeline in pipelines.pipelines:
        workflow.__desc__ += f"""*{pipeline.pipeline}: {pipeline._desc} \n"""

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "bold_t1",
                "boldmask",
                "anat2std_xfm",
                "templates",
                "confounds_file",
                "confounds_metadata",
            ]
        ),
        name="inputnode",
    )

    iterablesource = pe.Node(
        niu.IdentityInterface(fields=["pipeline"]), name="iterablesource2"
    )
    iterablesource.iterables = [("pipeline", pipelines.pipelines)]

    if spaces.get_spaces(nonstandard=False, dim=(3,)):
        bold_std_trans_wf = init_bold_std_trans_wf(
            mem_gb=mem_gb,
            omp_nthreads=omp_nthreads,
            spaces=spaces,
            name="bold_std_trans",
        )
        # fmt:off
        workflow.connect([
            (inputnode, bold_std_trans_wf, [
                ("templates", "inputnode.templates"),
                ("anat2std_xfm", "inputnode.anat2std_xfm"),
                ("bold_t1", "inputnode.bold_t1")
            ]),
        ])
        # fmt:on

    # Get pipeline info
    pipeline_info = pe.Node(
        niu.Function(
            function=_get_pipeline_info,
            input_name=["in_pipeline"],
            output_names=[
                "pipeline",
                "confounds_spec",
                "demean",
                "clean_specs",
            ],
        ),
        name="pipeline_info",
    )

    # Select and prepare confounds for given strategy
    select_confounds = pe.Node(SelectConfounds(), name="select_confounds")

    # Denoise image
    denoise_t1 = pe.Node(Denoise(), name="denoise_t1")
    if "RepetitionTime" in metadata:
        denoise_t1.inputs.tr = metadata["RepetitionTime"]

    denoise_std = pe.MapNode(
        Denoise(), iterfield=["input_image"], name="denoise_std"
    )
    if "RepetitionTime" in metadata:
        denoise_std.inputs.tr = metadata["RepetitionTime"]

    # fmt:off
    workflow.connect([
        (inputnode, select_confounds, [
            ("confounds_file", "confounds"),
            ("confounds_metadata", "confounds_metadata"),],),
        (inputnode, denoise_t1, [
            ("boldmask", "mask_img")]),
        (iterablesource, pipeline_info, [("pipeline", "in_pipeline")]),
        (pipeline_info, select_confounds, [
            ("pipeline", "pipeline"),
            ("confounds_spec", "confounds_spec"),
            ("demean", "demean"),],),
        (select_confounds, denoise_t1, [
            ("selected_confounds", "confounds_file")],),
        (select_confounds, denoise_std, [
            ("selected_confounds", "confounds_file")],),
        (pipeline_info, denoise_t1, [("clean_specs", "pipeline")]),
        (pipeline_info, denoise_std, [("clean_specs", "pipeline")]),
        (inputnode, denoise_t1, [("bold_t1", "input_image")]),
        (bold_std_trans_wf, denoise_std, [
            ("outputnode.bold_std", "input_image")],),
    ])
    # fmt:on

    output_names = ["bold_denoised_t1", "bold_denoised_std", "pipeline"]
    poutputnode = pe.Node(
        niu.IdentityInterface(fields=output_names), name="poutputnode2"
    )
    # fmt:off
    workflow.connect([
        # Connecting outputnode
        (iterablesource, poutputnode, [
            (("pipeline", _get_pipeline_name), "pipeline")]),
        (denoise_std, poutputnode, [("output_image", "bold_denoised_std")]),
        (denoise_t1, poutputnode, [("output_image", "bold_denoised_t1")]),
    ])
    # fmt:on

    # Connect parametric outputs to a Join outputnode
    outputnode_denoise = pe.JoinNode(
        niu.IdentityInterface(fields=output_names),
        name="outputnode_denoise",
        joinfield=output_names,
        joinsource="iterablesource2",
    )
    # fmt:off
    workflow.connect([
        (poutputnode, outputnode_denoise, [(f, f) for f in output_names]),
    ])
    # fmt:on

    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=output_names + ["template", "spatial_reference"]
        ),
        name="outputnode",
    )
    # fmt:off
    workflow.connect([
        (bold_std_trans_wf, outputnode, [("outputnode.template", "template"),
                                        ("outputnode.spatial_reference", "spatial_reference")]),
        (outputnode_denoise, outputnode, [(f, f) for f in output_names]),
    ])
    # fmt:on

    return workflow


def _get_pipeline_info(in_pipeline):
    pipeline = in_pipeline.pipeline
    confounds_spec = in_pipeline.confounds_spec
    demean = in_pipeline.demean
    clean_specs = in_pipeline.clean_specs
    return pipeline, confounds_spec, demean, clean_specs


def _get_pipeline_name(in_pipeline):
    return in_pipeline.pipeline
