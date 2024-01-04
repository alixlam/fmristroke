"""
Connectivity
^^^^^^^^^^^^

.. autofunction:: init_connectivity_wf

"""

import typing as ty


from nipype import Function
from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe

from ...utils.atlas import Atlases
from ...utils.pipelines import Pipelines
from ...interfaces import DerivativesDataSink


def init_connectivity_wf(
    mem_gb: float,
    atlases: Atlases,
    pipelines: Pipelines,
    conn_measure: list = ['correlation'],
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
    """
    
    from ...interfaces.nilearn import Connectivity
    from ...interfaces.pyants import ApplyTransforms
    from niworkflows.interfaces.utility import KeySelect
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow



    workflow = Workflow(name=name)
    inputnode = pe.Node(niu.IdentityInterface(
        fields=[
            "bold_denoised",
            "boldmask",
            "mask_lesion_std",
            "templates",
            "pipelines",
            "anat2std_xfm"
        ]
    ), name="inputnode")
    
    iterablesource = pe.Node(niu.IdentityInterface(fields=["atlas", "conn_measure", "pipeline"]),name="iterablesource3")
    iterablesource.iterables = [("atlas", atlases.atlases),
                                ("conn_measure", conn_measure,),
                                ("pipeline", pipelines.get_pipelines())]

    # Get Atlas info
    atlas_info = pe.Node(
        niu.Function(function=_get_atlas_info,
        input_name=["in_atlas"],
        output_names=["atlas", "labels", "space", "mask_file"]),
        name="atlas_info"
    )

    # Select correct space
    select_space = pe.Node(
        KeySelect(fields=["bold_denoised", "mask_lesion_std", "anat2std_xfm"]),
        name="select_space",
        run_without_submitting=True,
    )

    # Remove ROI voxels from Atlas ROI
    roi_resamp = pe.Node(
        niu.Function(function=_resample_to_img), name="roi_resamp")
    roi_resamp.inputs.interpolation = 'nearest'   
    atlas_roi = pe.Node(
        niu.Function(function=_remove_roi), name="atlas_roi"
    )
    
    # Bold mask in std space
    boldmask_tfm = pe.Node(
        ApplyTransforms(interpolation="multiLabel"), name="boldmask_std_tfm"
    )
    
    # Compute connectivity
    select_pipeline = pe.Node(
        KeySelect(fields=["bold_denoised"]),
        name="select_pipeline",
        run_without_submitting=True
    )
    connectivity = pe.Node(Connectivity(), name="connectivity")
    
    # Save derivatives 
    ds_report_connectivity = pe.Node(
        DerivativesDataSink(
            desc='connectivity',
            suffix="mat",
            dismiss_entities=("echo"),
            space=None,
        ),
        name='ds_report_connectivity',
        run_without_submitting=True,
        mem_gb=0.1,
    )    
    workflow.connect([
        (iterablesource, atlas_info, [("atlas", "in_atlas")]),
        (iterablesource, connectivity, [("conn_measure", "conn_measure")]),
        (inputnode, select_pipeline, [("bold_denoised", "bold_denoised"),
                                    ("pipelines", "keys")]),
        (iterablesource, select_pipeline, [("pipeline", "key")]),
        (select_pipeline, select_space, [("bold_denoised", "bold_denoised")]),
        (inputnode, select_space, [("mask_lesion_std", "mask_lesion_std"),
                                    ("anat2std_xfm", "anat2std_xfm"),
                                    ("templates", "keys")]),
        (atlas_info, select_space, [("space", "key")]),
        (select_space, roi_resamp, [("mask_lesion_std", "input_image")]),
        (atlas_info, roi_resamp, [("mask_file", "reference_image")]),
        (atlas_info, boldmask_tfm, [("mask_file", "reference_image")]),
        (select_space, boldmask_tfm, [("anat2std_xfm", "transforms")]),
        (inputnode, boldmask_tfm, [("boldmask", "input_image")]),
        (boldmask_tfm, connectivity, [("output_image", "brain_mask")]),
        (roi_resamp, atlas_roi, [("out", "roi")]),
        (atlas_info, atlas_roi, [("mask_file", "in_atlas")]),
        (atlas_roi, connectivity, [("out", "atlas")]),
        (select_space, connectivity, [("bold_denoised", "input_image")]),
        (connectivity, ds_report_connectivity, [("output_conn", "in_file")]),
        (iterablesource, ds_report_connectivity, [("atlas", "atlas"),
                                            ("conn_measure",  "measure"),
                                            ("pipeline", "pipeline")])
    ])
    
    return workflow

def _resample_to_img(input_image, reference_image, interpolation):
    from pathlib import Path
    from nilearn.image import resample_to_img
    new_img = resample_to_img(input_image, reference_image, interpolation=interpolation)
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

def _get_conn_meta(in_atlas, measure):
    import json
    from pathlib import Path
    dic = {
        "conn_measure": measure,
        "atlas": in_atlas.atlas,
    }
    json_object = json.dumps(dic, indent=4)
    out_name = Path("conn_meta.json").absolute()

    with open(out_name, "w") as outfile:
        outfile.write(json_object)
    return str(out_name)

