"""
Combining runs
^^^^^^^^^^^^^^^

.. autofunction:: init_concat_wf

"""

from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe

from ...interfaces import DerivativesDataSink


def init_concat_wf(
    mem_gb: float,
    croprun: int,
    name: str = "concat_wf",
):
    """
    Build a workflow to combine runs from the same session and task if multiple runs are available.

    This workflow combines bold series, confounds and masks.
    It crops the beginning of the run if necesssary.

    Parameters
    ----------
    mem_gb : :obj:`float`
        Size of BOLD file in GB - please note that this size
        should be calculated after resamplings that may extend
        the FoV
    croprun : :obj: int`
        crop beginning of bold run by n volumes.
    name : :obj:`str`
        Name of workflow (default: ``concat_wf``)

    Inputs
    ------
    bold_denoised_t1
        List of BOLD images in T1w space from the same session and task
    bold_denoised_std
        List of BOLD images in T1w space from the same session and task
    boldmask_t1
        List of BOLD masks in T1w space from the same session and task
    
    Outputs
    -------
    bold_denoised_t1
        Concatenated BOLD denoised image in T1w
    bold_denoised_std
        Concatenated BOLD denoised image in std space
    boldmask_t1
        Bold mask for concatenated BOLD files
    """

    from niworkflows.engine.workflows import LiterateWorkflow as Workflow

    from ...interfaces.reports import SessionSummary

    workflow = Workflow(name=name)

    workflow.__desc__ = f"""Concatenate bold runs from same session and task
    """

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "bold_denoised_t1",
                "bold_denoised_std",
                "boldmask_t1",
            ]
        ),
        name="inputnode",
    )

    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "bold_denoised_t1",
                "bold_denoised_std",
                "boldmask_t1",
            ]
        ),
        name="outputnode",
    )

    # Concatenate BOLD series
    concat_bold_denoised_t1 = pe.Node(
        niu.Function(function=_concat_bold), name="concat_bold_denoised_t1"
    )
    concat_bold_denoised_t1.inputs.croprun = croprun
    
    concat_bold_denoised_std = pe.Node(
        niu.Function(function=_concat_bold), name="concat_bold_denoised_std"
    )
    concat_bold_denoised_std.inputs.croprun = croprun
    concat_bold_denoised_std.inputs.std = True
    

    # Intersect boldmasks
    concat_boldmask = pe.Node(
        niu.Function(function=_intersect_masks), name="inter_mask"
    )

    # Report number of runs included in new bold run
    report_concat = pe.Node(SessionSummary(), name="concat_summary")
    report_concat.inputs.croprun = croprun

    ds_report_concat = pe.Node(
        DerivativesDataSink(desc="concat", datatype="figures", space=None),
        name="ds_report_concat",
        run_without_submitting=True,
    )

    # fmt:off
    workflow.connect([
        # Concatenate all
        (inputnode, concat_bold_denoised_std, [("bold_denoised_std", "in_files")]),
        (inputnode, concat_bold_denoised_t1, [("bold_denoised_t1", "in_files")]),
        (inputnode, concat_boldmask, [("boldmask_t1", "in_files")]),

        #report
        (inputnode, report_concat, [("bold_t1", "bold_t1")]),
        (report_concat, ds_report_concat, [("out_report", "in_file")]),

        # outputs
        (concat_bold_denoised_std, outputnode, [("out", "bold_denoised_std")]),
        (concat_bold_denoised_t1, outputnode, [("out", "bold_denoised_t1")]),

        (concat_boldmask, outputnode, [("out", "boldmask_t1")]),
    ])

    return workflow


def _concat_bold(in_files, croprun, bold_std=False):
    from pathlib import Path

    import nibabel as nib
    from nilearn.image import concat_imgs

    if isinstance(in_files, str):
        in_files = [in_files]

    cropped_files = []
    for f in in_files:
        pipelines = []
        for i, pipe in enumerate(f):
            if std:
                space = []
                for j, s in enumerate(pipe):
                    nii = nib.load(s)
                    croped_file = nib.Nifti1Image(
                        nii.get_fdata()[:, :, :, croprun:], nii.affine
                    )
                    space.append(croped_file)
                pipelines.append(space)
            else:
                nii = nib.load(pipe)
                croped_file = nib.Nifti1Image(
                        nii.get_fdata()[:, :, :, croprun:], nii.affine
                    )
                pipelines.append([croped_file])
        cropped_files.append(pipelines)

    reshape_files = list(np.array(cropped_files).T)
    concatenated_img = [[concat_imgs(space) for space in pipe] for pipe in reshape_files]
    
    n_pipeline = len(concatenated_img)
    n_space = len(concatenated_img[0])
    
    out = []
    for i, pipe in concatenated_img:
        out_pipe = []
        for j, space in pipe:
            outname = Path(f"concat_bold_pipe_{i}_space_{j}.nii.gz").absolute()
            j.to_filename(outname)
            out_pipe.append(str(outname))
        out.append(out_pipe)
    return out


def _intersect_masks(in_files):
    from pathlib import Path

    from nilearn.masking import intersect_masks

    if isinstance(in_files, str):
        in_files = [in_files]
    out_name = Path("mask_inter.nii.gz").absolute()
    new_mask = intersect_masks(in_files, 1)
    new_mask.to_filename(out_name)
    return str(out_name)
