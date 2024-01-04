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
    bold_t1
        List of BOLD images in T1w space from the same session and task
    boldref_t1
        List of BOLD ref images in T1w space from the same session and task
    boldmask_t1
        List of BOLD masks in T1w space from the same session and task
    confounds_file
        List of confounds files
    confounds_metadata
        List of confounds metadata files

    Outputs
    -------
    bold_t1
        Concatenated BOLD image
    boldref_t1
        Bold ref for new concatenated files
    boldmask_t1
        Bold mask for concatenated BOLD files
    confounds_file
        Concatenated confounds
    confounds_metadata
        TO IMPLEMENT ??
    """

    from niworkflows.engine.workflows import LiterateWorkflow as Workflow

    from ...interfaces.reports import SessionSummary

    workflow = Workflow(name=name)

    workflow.__desc__ = f"""Concatenate bold runs from same session and task
    """

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "bold_t1",
                "boldref_t1",
                "boldmask_t1",
                "t1w_mask",
                "confounds_file",
                "confounds_metadata",
            ]
        ),
        name="inputnode",
    )

    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "bold_t1",
                "boldref_t1",
                "boldmask_t1",
                "t1w_mask",
                "confounds_file",
                "confounds_metadata",
            ]
        ),
        name="outputnode",
    )

    # Concatenate BOLD series
    concat_bold = pe.Node(
        niu.Function(function=_concat_bold), name="conca_bold"
    )
    concat_bold.inputs.croprun = croprun

    # Concatenate confounds
    concat_confounds = pe.Node(
        niu.Function(function=_concat_confounds), name="concat_confounds"
    )
    concat_confounds.inputs.croprun = croprun

    # Intersect boldmasks
    concat_boldmask = pe.Node(
        niu.Function(function=_intersect_masks), name="inter_mask"
    )

    # Confounds metadata and boldref
    # Boldref is used only in confounds_wf as reference for resampling
    # can select any boldref and it will be the same result
    # Confunds metadata is used when using compcor regresors in denoising
    # it will be computed again in session level pipelines
    select_boldref = pe.Node(niu.Select(index=0), name="select_boldref")

    select_confoundsmeta = pe.Node(
        niu.Select(index=0), name="select_confoundsmeta"
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
        (inputnode, concat_bold, [("bold_t1", "in_files")]),
        (inputnode, concat_confounds, [("confounds_file", "in_files")]),
        (inputnode, concat_boldmask, [("boldmask_t1", "in_files")]),
        (inputnode, select_boldref, [("boldref_t1", "inlist")]),
        (inputnode, select_confoundsmeta, [("confounds_metadata", "inlist")]),

        #report
        (inputnode, report_concat, [("bold_t1", "bold_t1")]),
        (report_concat, ds_report_concat, [("out_report", "in_file")]),

        # outputs
        (concat_bold, outputnode, [("out", "bold_t1")]),
        (concat_boldmask, outputnode, [("out", "boldmask_t1")]),
        (concat_confounds, outputnode, [("out", "confounds_file")]),
        (select_confoundsmeta, outputnode, [("out", "confounds_metadata")]),
        (select_boldref, outputnode, [("out", "boldref_t1")]),



    ])

    return workflow


def _concat_bold(in_files, croprun):
    from pathlib import Path

    import nibabel as nib
    from nilearn.image import concat_imgs

    # make sure files are sorted
    if isinstance(in_files, str):
        in_files = [in_files]
    in_files = sorted(in_files)

    if croprun != 0:
        list_runs = []
        for f in in_files:
            nii = nib.load(f)
            croped_file = nib.Nifti1Image(
                nii.get_fdata()[:, :, :, croprun:], nii.affine
            )
            list_runs.append(croped_file)
        new_img = concat_imgs(list_runs)
    else:
        new_img = concat_imgs(in_files)
    out_name = Path("concatbold.nii.gz").absolute()
    new_img.to_filename(out_name)
    return str(out_name)


def _concat_confounds(in_files, croprun):
    from pathlib import Path

    import pandas as pd

    if isinstance(in_files, str):
        in_files = [in_files]
    in_files = sorted(in_files)

    confs = []
    for f in in_files:
        conf_run = pd.read_csv(f, sep="\t")
        confs.append(conf_run.iloc[croprun:])
    concat = pd.concat(confs)
    # If concatenating runs remove compcor confounds
    if len(in_files) > 1:
        concat = concat.drop(
            [c for c in concat.columns if "comp" in c], axis=1
        )
    out_name = Path("concat_confounds.tsv").absolute()
    concat.to_csv(str(out_name), sep="\t")
    return str(out_name)


def _intersect_masks(in_files):
    from pathlib import Path

    from nilearn.masking import intersect_masks

    if isinstance(in_files, str):
        in_files = [in_files]
    out_name = Path("mask_inter.nii.gz").absolute()
    new_mask = intersect_masks(in_files, 1)
    new_mask.to_filename(out_name)
    return str(out_name)
