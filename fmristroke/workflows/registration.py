"""
Genrerate Registration plot with roi
^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: init_lesionplot_wf

"""
from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from templateflow.api import get as get_template

from fmriprep import config

from ..interfaces import DerivativesDataSink

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
    from ..interfaces.reports import RegisterLesionRPT
    from ..interfaces.pyants import ApplyTransforms


    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "bold_ref",
                "t1w",
                "t1w_roi",
                "t1w_mask",
                "bold_t1_xform",
            ]
        ),
        name="inputnode",
    )

    outputnode = pe.Node(niu.IdentityInterface(fields=["out_regplot"]), name="outputnode")
    bold_t1 = pe.Node(ApplyTransforms(interpolation="lanczosWindowedSinc"), "bold_t1_trans")
    
    reg_plot = pe.Node(
        RegisterLesionRPT(),
        name="reg_plot"
    )
    ds_report_reg = pe.Node(
            DerivativesDataSink(datatype="figures", dismiss_entities=("echo",)),
            name='ds_report_reg',
            run_without_submitting=True,
        )

    workflow = Workflow(name=name)

    # fmt:off
    workflow.connect([
        (inputnode, bold_t1, [("bold_ref", "input_image"),
                                ("t1w", "reference_image"),
                                ("bold_t1_xform", "transforms")]),
        (inputnode, reg_plot, [("t1w_mask", "reference_image_mask"),
                            ("t1w", "reference_file"),
                            ("t1w_roi", "roi")]),
        (bold_t1, reg_plot, [("output_image", "moving_image")]),
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
