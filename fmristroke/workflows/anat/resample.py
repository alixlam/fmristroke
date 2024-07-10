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


def init_roi_std_trans_wf(
    mem_gb: float,
    omp_nthreads: int,
    spaces: SpatialReferences,
    name: str = "lesion_std_trans_wf",
    use_compression: bool = True,
):
    """
    Sample lesion roi into standard space with a single-step resampling of the T1w lesion.

    .. important::
        This workflow provides two outputnodes.
        One output node (with name ``poutputnode``) will be parameterized in a Nipype sense
        (see `Nipype iterables
        <https://miykael.github.io/nipype_tutorial/notebooks/basic_iteration.html>`__), and a
        second node (``outputnode``) will collapse the parameterized outputs into synchronous
        lists of the output fields listed below.


    Parameters
    ----------
    mem_gb : :obj:`float`
        Size of BOLD file in GB
    omp_nthreads : :obj:`int`
        Maximum number of threads an individual process may use
    spaces : :py:class:`~niworkflows.utils.spaces.SpatialReferences`
        A container for storing, organizing, and parsing spatial normalizations. Composed of
        :py:class:`~niworkflows.utils.spaces.Reference` objects representing spatial references.
        Each ``Reference`` contains a space, which is a string of either TemplateFlow template IDs
        (e.g., ``MNI152Lin``, ``MNI152NLin6Asym``, ``MNIPediatricAsym``), nonstandard references
        (e.g., ``T1w`` or ``anat``, ``sbref``, ``run``, etc.), or a custom template located in
        the TemplateFlow root directory. Each ``Reference`` may also contain a spec, which is a
        dictionary with template specifications (e.g., a specification of ``{"resolution": 2}``
        would lead to resampling on a 2mm resolution of the space).
    name : :obj:`str`
        Name of workflow (default: `roi_std_trans_wf``)
    use_compression : :obj:`bool`
        Save registered roi masks as ``.nii.gz``

    Inputs
    ------
    anat2std_xfm
        List of anatomical-to-standard space transforms generated during
        spatial normalization.
    templates
        List of templates that were applied as targets during
        spatial normalization.

    Outputs
    -------
    roi_std
        roi, resampled to template space
    template
        Template identifiers synchronized correspondingly to previously
        described outputs.

    """
    from fmriprep.interfaces.maths import Clip
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow
    from niworkflows.interfaces.nibabel import GenerateSamplingReference
    from niworkflows.interfaces.utility import KeySelect
    from niworkflows.utils.spaces import format_reference

    from ...interfaces.pyants import ApplyTransforms

    workflow = Workflow(name=name)
    output_references = spaces.cached.get_spaces(nonstandard=False, dim=(3,))
    std_vol_references = [
        (s.fullname, s.spec)
        for s in spaces.references
        if s.standard and s.dim == 3
    ]

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "anat2std_xfm",
                "roi",
                "templates",
            ]
        ),
        name="inputnode",
    )

    iterablesource = pe.Node(
        niu.IdentityInterface(fields=["std_target"]), name="iterablesource"
    )
    # Generate conversions for every template+spec at the input
    iterablesource.iterables = [("std_target", std_vol_references)]

    split_target = pe.Node(
        niu.Function(
            function=_split_spec,
            input_names=["in_target"],
            output_names=["space", "template", "spec"],
        ),
        run_without_submitting=True,
        name="split_target",
    )

    select_std = pe.Node(
        KeySelect(fields=["anat2std_xfm"]),
        name="select_std",
        run_without_submitting=True,
    )

    select_tpl = pe.Node(
        niu.Function(function=_select_template),
        name="select_tpl",
        run_without_submitting=True,
    )

    gen_ref = pe.Node(GenerateSamplingReference(), name="gen_ref", mem_gb=0.01)

    mask_std_tfm = pe.Node(
        ApplyTransforms(interpolation="multiLabel"),
        name="mask_std_tfm",
        mem_gb=1,
    )

    # Write corrected file in the designated output dir
    mask_merge_tfms = pe.Node(
        niu.Merge(2),
        name="mask_merge_tfms",
        run_without_submitting=True,
        mem_gb=DEFAULT_MEMORY_MIN_GB,
    )

    merge_xforms = pe.Node(
        niu.Merge(4),
        name="merge_xforms",
        run_without_submitting=True,
        mem_gb=DEFAULT_MEMORY_MIN_GB,
    )
    # fmt:off
    workflow.connect([
        (inputnode, gen_ref, [("roi", "moving_image")]),
        (iterablesource, split_target, [("std_target", "in_target")]),
        (iterablesource, select_tpl, [("std_target", "template")]),
        (inputnode, select_std, [("anat2std_xfm", "anat2std_xfm"),
                                ("templates", "keys")]),
        (inputnode, mask_std_tfm, [("roi", "input_image")]),
        (split_target, select_std, [("space", "key")]),
        (select_std, merge_xforms, [("anat2std_xfm", "in1")]),
        (select_std, mask_merge_tfms, [("anat2std_xfm", "in1")]),
        (split_target, gen_ref, [(("spec", _is_native), "keep_native")]),
        (select_tpl, gen_ref, [("out", "fixed_image")]),
        (gen_ref, mask_std_tfm, [("out_file", "reference_image")]),
        (mask_merge_tfms, mask_std_tfm, [("out", "transforms")]),
    ])
    # fmt:on

    output_names = [
        "roi_mask_std",
        "spatial_reference",
        "template",
    ]

    poutputnode = pe.Node(
        niu.IdentityInterface(fields=output_names), name="poutputnode"
    )
    # fmt:off
    workflow.connect([
        # Connecting outputnode
        (iterablesource, poutputnode, [
            (("std_target", format_reference), "spatial_reference")]),
        (mask_std_tfm, poutputnode, [("output_image", "roi_mask_std")]),
        (select_std, poutputnode, [("key", "template")]),
    ])
    # fmt:on

    # Connect parametric outputs to a Join outputnode
    outputnode = pe.JoinNode(
        niu.IdentityInterface(fields=output_names),
        name="outputnode",
        joinsource="iterablesource",
    )
    # fmt:off
    workflow.connect([
        (poutputnode, outputnode, [(f, f) for f in output_names]),
    ])
    # fmt:on
    return workflow


def _split_spec(in_target):
    space, spec = in_target
    template = space.split(":")[0]
    return space, template, spec


def _select_template(template):
    from niworkflows.utils.misc import get_template_specs

    template, specs = template
    template = template.split(":")[0]  # Drop any cohort modifier if present
    specs = specs.copy()
    specs["suffix"] = specs.get("suffix", "T1w")

    # Sanitize resolution
    res = specs.pop("res", None) or specs.pop("resolution", None) or "native"
    if res != "native":
        specs["resolution"] = res
        return get_template_specs(template, template_spec=specs)[0]

    # Map nonstandard resolutions to existing resolutions
    specs["resolution"] = 2
    try:
        out = get_template_specs(template, template_spec=specs)
    except RuntimeError:
        specs["resolution"] = 1
        out = get_template_specs(template, template_spec=specs)

    return out[0]


def _is_native(in_value):
    return (
        in_value.get("resolution") == "native"
        or in_value.get("res") == "native"
    )
