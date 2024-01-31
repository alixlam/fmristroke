"""
Hemodynamics
^^^^^^^^^^^^^

.. autofunction:: init_hemodynamic_wf

"""
from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe

from ...interfaces import DerivativesDataSink


def init_hemodynamic_wf(
    mem_gb: float,
    metadata: dict,
    maxlag: int = 10,
    name: str = "hemodynamic_wf",
):
    """
    Build a workflow to generate lagmaps.

    This workflow calculates the hemodynamic lag for stroke patients.
    as recommended by Siegel et al. (2017) ``https://pubmed.ncbi.nlm.nih.gov/28541130/``
    It uses the existing tool RapidTide that calculates a similarity function between a
    “probe” signal and every voxel of a BOLD fMRI dataset. It then determines the peak value,
    time delay, and width of the similarity function to determine when and how strongly that
    probe signal appears in each voxel.
    For details about the method visit : ``https://rapidtide.readthedocs.io/en/latest/``


    Parameters
    ----------
    mem_gb : :obj:`float`
        Size of BOLD file in GB - please note that this size
        should be calculated after resamplings that may extend
        the FoV
    metadata : :obj:`dict`
        BIDS metadata for BOLD file
    maxlag :obj: `int`
        Max lag to compute.
    name : :obj:`str`
        Name of workflow (default: ``hemodynamic_wf``)

    Inputs
    ------
    bold_t1
        BOLD image in T1w space, after the prescribed corrections (STC, HMC and SDC)
        when available.
    boldmask
        Bold fov mask
    roi
        Roi mask in T1w space
    t1w
        T1w image
    t1w_mask_lesion
        Mask of the lesion in T1w space
    t1w_mask
        Brain Mask of T1
    t1w_tpms
        List of tissue probability maps in T1w space
    t1w_aseg
        Segmentation of structural image, done with FreeSurfer.
    confounds_file
        TSV of all aggregated confounds.

    Outputs
    -------
    maps
        lagmap
    """

    from niworkflows.engine.workflows import LiterateWorkflow as Workflow
    from niworkflows.interfaces.morphology import BinaryDilation
    from niworkflows.interfaces.nibabel import Binarize
    from niworkflows.interfaces.utility import KeySelect
    from templateflow.api import get as get_template

    from ...interfaces.nilearn import Smooth
    from ...interfaces.pyants import ApplyTransforms
    from ...interfaces.rapidtide.rapidtide import RapidTide
    from ...interfaces.reports import HemodynamicsSummary, HemoPlot

    workflow = Workflow(name=name)

    workflow.__desc__ = f"""\\TO DO
        """

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "bold_t1",
                "boldmask",
                "roi",
                "t1w_preproc",
                "t1w_mask",
                "t1w_tpms",
                "t1w_aseg",
                "confounds_file",
                "std2anat_xfm",
                "template",
            ]
        ),
        name="inputnode",
    )
    outputnode = pe.Node(
        niu.IdentityInterface(fields=["lagmaps", "lagmap_metadata", "gm_mask"]),
        name="outputnode",
    )

    # Generate GM mask and interesct with bold fov
    # Select xform to "MNI152Lin2009cAsym" (always computed)
    select_mni_xfm = pe.Node(
        KeySelect(fields=["std2anat_xfm"], key="MNI152NLin2009cAsym"),
        name="select_mni_xfm",
        run_without_submitting=True,
    )

    # Warp segmentation into T1 space
    resample_parc_T1 = pe.Node(
        ApplyTransforms(
            input_image=str(
                get_template(
                    "MNI152NLin2009cAsym",
                    resolution=1,
                    desc="carpet",
                    suffix="dseg",
                    extension=[".nii", ".nii.gz"],
                )
            ),
            interpolation="multiLabel",
        ),
        name="resample_parc",
    )
    gm_mask = pe.Node(
        niu.Function(function=_gm_parcellation), name="GM_parcel"
    )

    intersect_mask = pe.Node(
        niu.Function(function=_intersect_masks), name="intersect_mask"
    )

    # Resample masks to bold resolution
    roi_resamp = pe.Node(
        niu.Function(function=_resample_to_img), name="roi_resamp"
    )
    roi_resamp.inputs.interpolation = "nearest"
    gm_mask_resamp = pe.Node(
        niu.Function(function=_resample_to_img), name="gm_resamp"
    )
    gm_mask_resamp.inputs.interpolation = "nearest"

    # generate lagmaps
    lagmaps = pe.Node(RapidTide(searchrange=maxlag), name="rapidtide")
    if "RepetitionTime" in metadata:
        lagmaps.inputs.repetitiontime = metadata["RepetitionTime"]

    # Mask for not significative correlations and Smooth lagmaps
    merge_maps = pe.Node(
        niu.Merge(2), name="merge_maps", run_without_submitting=True
    )
    smooth_lagmap = pe.MapNode(
        Smooth(),
        iterfield=["input_image"],
        name="smooth_maps",
    )
    mask_maps = pe.MapNode(
        niu.Function(function=_mask_img),
        iterfield=["input_image"],
        name="mask_maps",
    )

    # Extract metadata
    report_hemo = pe.Node(
        HemodynamicsSummary(),
        name="hemosummary",
    )

    select_lagmap = pe.Node(niu.Select(index=0), name="lagmap_select")
    select_corrmap = pe.Node(niu.Select(index=1), name="corrmap_select")

    # Mask T1w for report
    mask_img = pe.Node(niu.Function(function=_mask_img), name="mask_img")

    # Hemo report
    lagplot = pe.Node(
        HemoPlot(label="LagMap"),
        name="hemo_plot_lag",
    )
    corrplot = pe.Node(
        HemoPlot(label="CorrPlot", vmax=1), name="hemo_plot_corr"
    )

    ds_report_lagplot = pe.Node(
        DerivativesDataSink(desc="lagmap", datatype="figures"),
        name="ds_report_lagplot",
        run_without_submitting=True,
        dismiss_entities=("space",),
    )

    ds_report_corrplot = pe.Node(
        DerivativesDataSink(desc="corrmap", datatype="figures"),
        name="ds_report_corrplot",
        run_without_submitting=True,
        dismiss_entities=("space",),
    )

    ds_report_hemo = pe.Node(
        DerivativesDataSink(desc="hemo", datatype="figures"),
        name="ds_report_hemo",
        run_without_submitting=True,
        dismiss_entities=("space",),
    )
    # fmt:off
    workflow.connect([
        # Brain mask and GM mask
        (inputnode, select_mni_xfm, [("std2anat_xfm", "std2anat_xfm"),
                                     ("template", "keys")]),
        (inputnode, resample_parc_T1, [("t1w_preproc", "reference_image")]),
        (select_mni_xfm, resample_parc_T1, [("std2anat_xfm", "transforms")]),
        (resample_parc_T1, gm_mask, [("output_image", "segmentation")]),
        (gm_mask, gm_mask_resamp, [("out", "in_img")]),
        (inputnode, gm_mask_resamp, [("boldmask", "ref")]),
        (gm_mask_resamp, intersect_mask, [("out", "in1")]),
        (inputnode, intersect_mask, [("boldmask", "in2")]),
        (inputnode, roi_resamp, [("roi", "in_img"), 
                                ("bold_t1", "ref")]),

        # Hemodynamics
        (inputnode, lagmaps, [("bold_t1", "in_file"),
                                ("confounds_file", "confounds_file")]),
        (intersect_mask, lagmaps, [("out", "corrmask"),
                                    ("out", "globalmeaninclude")]),
        (roi_resamp, lagmaps, [("out", "globalmeanexclude")]),
        (lagmaps, merge_maps, [("output_lagmap", "in1"),
                                ("output_corrmap", "in2")]),
        (merge_maps, smooth_lagmap, [("out", "input_image")]),
        (lagmaps, mask_maps, [("output_corrfit_mask", "mask_file")]),
        (smooth_lagmap, mask_maps, [("output_image", "input_image")]),
        (mask_maps, report_hemo, [("out", "in_files")]),
        (lagmaps, report_hemo, [("output_corrfit_mask", "corrfit_mask")]),
        (gm_mask_resamp, report_hemo, [("out", "brain_mask")]),
        (inputnode, mask_img, [("t1w_preproc", "input_image"),
                                ("t1w_mask", "mask_file")]),
        (inputnode, report_hemo, [("t1w_aseg", "hemi_mask")]),
        (inputnode, lagplot, [("roi", "roi")]),
        (mask_img, lagplot, [("out", "anat_file")]),
        (inputnode, corrplot, [("roi", "roi")]),
        (mask_img, corrplot, [("out", "anat_file")]),
        (mask_maps, select_lagmap, [("out", "inlist")]),
        (mask_maps, select_corrmap, [("out", "inlist")]),
        (select_lagmap, lagplot, [("out", "lagmap_file")]),
        (select_corrmap, corrplot, [("out", "lagmap_file")]),

        # Set outputs
        (mask_maps, outputnode, [("out", "lagmaps")]),
        (gm_mask, outputnode, [("out", "gm_mask")]),
        (report_hemo, ds_report_hemo, [("out_report", "in_file")]),
        (lagplot, ds_report_lagplot, [("out_report", "in_file")]),
        (corrplot, ds_report_corrplot, [("out_report", "in_file")]),

    ])
    return workflow


def _mask_img(input_image, mask_file):
    from pathlib import Path

    from nilearn.masking import apply_mask, unmask

    new_file = unmask(apply_mask(input_image, mask_img=mask_file), mask_file)
    out_name = Path("masked_image.nii.gz").absolute()
    new_file.to_filename(out_name)
    return str(out_name)


def _intersect_masks(in1, in2):
    from pathlib import Path

    import nibabel as nib
    import numpy as np

    img = nib.load(in1)
    mskarr1 = np.asanyarray(img.dataobj, dtype=int) > 0
    mskarr2 = np.asanyarray(nib.load(in2).dataobj, dtype=int) > 0
    out = img.__class__(mskarr1 & mskarr2, img.affine, img.header)
    out.set_data_dtype("uint8")
    out_name = Path("mask_inter.nii.gz").absolute()
    out.to_filename(out_name)
    return str(out_name)


def _resample_to_img(in_img, ref, interpolation):
    from pathlib import Path

    from nilearn.image import resample_to_img

    new_img = resample_to_img(in_img, ref, interpolation=interpolation)
    out_name = Path("mask_resamp.nii.gz").absolute()
    new_img.to_filename(out_name)
    return str(out_name)


def _gm_parcellation(segmentation):
    """Generate GM mask"""
    from pathlib import Path

    import nibabel as nb
    import numpy as np

    img = nb.load(segmentation)

    lut = np.zeros((256,), dtype="uint8")
    lut[1:11] = 0  # WM+CSF
    lut[100:201] = 1  # Ctx GM
    lut[30:99] = 1  # dGM
    lut[255] = 0  # Cerebellum
    # Apply lookup table
    seg = lut[np.uint16(img.dataobj)]

    outimg = img.__class__(seg.astype("uint8"), img.affine, img.header)
    outimg.set_data_dtype("uint8")
    out_file = Path("segments.nii.gz").absolute()
    outimg.to_filename(out_file)
    return str(out_file)
