"""
Genrerate Registration plot with roi
^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: init_lesionplot_wf

"""
from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from templateflow.api import get as get_template
from pathlib import Path


from fmriprep import config

from ...interfaces import DerivativesDataSink

def init_lesionplot_wf(
    mem_gb: float, name: str = "lesion_plot_wf"
):
    """
    Build a workflow to generate registration plots with lesions.

    Parameters
    ----------
    mem_gb : :obj:`float`
        Size of BOLD file in GB - please note that this size
        should be calculated after resamplings that may extend
        the FoV
    name : :obj:`str`
        Name of workflow (default: ``lesion_plot_wf``)

    Inputs
    ------
    bold_ref_t1
        Bold reference
    t1w
        t1w image
    t1w_roi
        Roi mask in T1w space 
    t1w_mask
        Mask of the skull-stripped template image

    Outputs
    -------
    out_registration_plot
        Path of the generated SVG file

    """
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow
    from ...interfaces.reports import RegisterLesionRPT
    from ...interfaces.pyants import ApplyTransforms


    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "boldref_t1",
                "t1w",
                "t1w_roi",
                "t1w_mask",
            ]
        ),
        name="inputnode",
    )

    outputnode = pe.Node(niu.IdentityInterface(fields=["out_regplot"]), name="outputnode")
    
    resample_bold = pe.Node(
        niu.Function(function=_resample_to_img), name="boldref_resamp")
    resample_bold.inputs.interpolation = "continuous"
    
    mask_img = pe.Node(
        niu.Function(function=_mask_img), name="mask_img")
    
    reg_plot = pe.Node(
        RegisterLesionRPT(),
        name="reg_plot"
    )
    ds_report_reg = pe.Node(
            DerivativesDataSink(datatype="figures", desc="reglesion", dismiss_entities=("echo",)),
            name='ds_report_reg',
            run_without_submitting=True,
            dismiss_entities=("space",)

        )

    workflow = Workflow(name=name)

    # fmt:off
    workflow.connect([
        (inputnode, resample_bold, [
                            ("boldref_t1", "in_img"),
                            ("t1w", "ref"),]),
        (inputnode, mask_img, [("t1w", "in_img"),
                                ("t1w_mask", "in_mask")]),
        (mask_img, reg_plot, [("out", "reference_file")]),
        (inputnode, reg_plot, [("t1w_roi", "roi")]),
        (resample_bold, reg_plot, [("out", "moving_image")]),
        (reg_plot, ds_report_reg, [("out_report", "in_file")]),
        (reg_plot, outputnode, [("out_report", "out_regplot")]),
    ])
    # fmt:on
    return workflow


def _resample_to_img(in_img, ref, interpolation):
    from pathlib import Path
    from nilearn.image import resample_to_img
    new_img = resample_to_img(in_img, ref, interpolation=interpolation)
    out_name = Path("mask_resamp.nii.gz").absolute()
    new_img.to_filename(out_name)
    return str(out_name)

def _mask_img(in_img, in_mask):
    from pathlib import Path
    from nilearn.image import threshold_img, load_img
    from nilearn.masking import apply_mask, unmask
    img= load_img(in_img)
    out_image = unmask(
                apply_mask(in_img, in_mask),
                in_mask,
            )
    out_name = Path("mask_resamp.nii.gz").absolute()
    out_image.to_filename(out_name)
    return str(out_name)
