from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe

from ..interfaces import DerivativesDataSink

def init_hemodynamic_wf(
    mem_gb: float,
    metadata: dict,
    maxlag : int = 10, 
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
    


    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from fmristroke.workflows.lagmaps import init_hemodynamic_wf
            wf = init_hemodynamic_wf(
                mem_gb=1,
                metadata={},
            )

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
        Name of workflow (default: ``bold_confs_wf``)

    Inputs
    ------
    bold_t1
        BOLD image in T1w space, after the prescribed corrections (STC, HMC and SDC)
        when available.
    bold_mask
        Bold fov mask
    roi_t1w
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
    bold_t1_xform
        Affine matrix that maps the native BOLD space into alignement the T1w space
    confounds_file
        TSV of all aggregated confounds.

    Outputs
    -------
    maps
        lagmap and corrmap
    """
    
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow
    from ..interfaces.rapidtide.rapidtide import RapidTide
    from ..interfaces.pyants import ApplyTransforms
    from ..interfaces.nilearn import Smooth
    from niworkflows.interfaces.morphology import BinaryDilation
    from niworkflows.interfaces.nibabel import Binarize
    from ..interfaces.reports import HemodynamicsSummary, HemoPlot

    workflow = Workflow(name=name)
    
    workflow.__desc__ = f"""\TO DO
        """
    
    
    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "bold_t1",
                "bold_mask",
                "roi",
                "t1w",
                "t1w_mask",
                "t1w_tpms",
                "bold_t1_xform",
                "t1w_aseg",
                "confounds_file",
            ]
        ),
        name="inputnode",
    )
    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "lagmaps",
                "lagmap_metadata"]
        ),
        name="outputnode",
    )
    
    # Project bold fov mask into T1w space
    t1w_mask_fov_tfm = pe.Node(
        ApplyTransforms(interpolation="multiLabel"),
        name="t1w_mask_fov_tfm",
    )
    
    # Generate GM mask and interesct with bold fov
    gm_mask = pe.Node(niu.Select(index=0), name="gm_mask")
    gm_mask_bin = pe.Node(Binarize(thresh_low = 0.02), name="gm_mask_bin")
    gm_mask_dilate = pe.Node(BinaryDilation(), name="gm_mask_dilate")
    intersect_mask = pe.Node(niu.Function(function=_intersect_masks), name="intersect_mask")
    
    # Resample masks to bold resolution
    roi_resamp = pe.Node(
        niu.Function(function=_resample_to_img), name="roi_resamp")
    roi_resamp.inputs.interpolation = "nearest"
    gm_mask_resamp =  pe.Node(
        niu.Function(function=_resample_to_img), name="gm_resamp")
    gm_mask_resamp.inputs.interpolation = "nearest"

    
    # generate lagmaps
    lagmaps = pe.Node(
        RapidTide(searchrange=maxlag),
        name="rapidtide")
    if "RepetitionTime" in metadata:
        lagmaps.inputs.repetitiontime = metadata["RepetitionTime"]
    
    # Mask for not significative correlations and Smooth lagmaps
    merge_maps = pe.Node(
        niu.Merge(2), name="merge_maps", run_without_submitting=True
    )
    smooth_lagmap = pe.MapNode(
        Smooth(),
        iterfield=["input_image"],
        name="smooth_maps",)
    mask_maps = pe.MapNode(
        niu.Function(function=_mask_img),
        iterfield=['input_image'],
        name="mask_maps",
    )
    
    # Extract metadata 
    report_hemo = pe.Node(
        HemodynamicsSummary(),
        name="hemosummary",
    )
    
    select_lagmap = pe.Node(niu.Select(index=0), name="lagmap_select")
    select_corrmap =  pe.Node(niu.Select(index=1), name="corrmap_select")
    
    # Hemo report
    lagplot = pe.Node(
        HemoPlot(label="LagMap"),
        name="hemo_plot_lag",
    )
    corrplot = pe.Node(
        HemoPlot(label = "CorrPlot", vmax = 1),
        name="hemo_plot_corr"
    )
    
    ds_report_lagplot = pe.Node(
        DerivativesDataSink(desc="lagmap", datatype="figures"),
        name="ds_report_lagplot",
        run_without_submitting=True,
    )
    
    ds_report_corrplot = pe.Node(
        DerivativesDataSink(desc="corrmap", datatype="figures"),
        name="ds_report_corrplot",
        run_without_submitting=True,
    )
    
    ds_report_hemo = pe.Node(
        DerivativesDataSink(desc="hemo", datatype="figures"),
        name="ds_report_hemo",
        run_without_submitting=True,
    )
    # fmt:off
    workflow.connect([
        # Brain mask and GM mask
        (inputnode, t1w_mask_fov_tfm, [("bold_mask", "input_image"),
                                    ("t1w_mask", "reference_image"),
                                    ("bold_t1_xform", "transforms")]),
        (inputnode, gm_mask, [("t1w_tpms", "inlist")]),
        (gm_mask, gm_mask_bin, [("out", "in_file")]),
        (gm_mask_bin, gm_mask_dilate, [("out_file","in_mask")]),
        (gm_mask_dilate, intersect_mask, [("out_mask", "in1")]),
        (t1w_mask_fov_tfm, intersect_mask, [("output_image", "in2")]),
        (intersect_mask, gm_mask_resamp, [("out", "in_img")]),
        (inputnode, gm_mask_resamp, [("bold_t1", "ref")]),
        (inputnode, roi_resamp, [("roi", "in_img"), 
                                ("bold_t1", "ref")]),
        
        # Hemodynamics
        (inputnode, lagmaps, [("bold_t1", "in_file"),
                                ("confounds_file","confounds_file")]),
        (gm_mask_resamp, lagmaps, [("out", "corrmask"),
                                    ("out", "globalmeaninclude")]),
        (roi_resamp, lagmaps, [("out", "globalmeanexclude")]),
        (lagmaps, merge_maps, [("output_lagmap", "in1"),
                                ("output_corrmap", "in2")]),
        (merge_maps, smooth_lagmap, [("out", "input_image")]),
        (lagmaps, mask_maps, [("output_corrfit_mask", "mask_file")]),
        (smooth_lagmap, mask_maps, [("output_image", "input_image")]),
        (mask_maps, report_hemo, [("out", "in_files")]),
        (lagmaps, report_hemo, [("output_corrfit_mask", "corrfit_mask")]),
        (gm_mask_resamp, report_hemo, [("out","brain_mask")]),
        (inputnode, report_hemo, [("t1w_aseg", "hemi_mask")]),
        (inputnode, lagplot, [("roi", "roi"),
                            ("t1w", "anat_file")]),
        (inputnode, corrplot, [("roi", "roi"),
                            ("t1w", "anat_file")]),
        (mask_maps, select_lagmap, [("out", "inlist")]),
        (mask_maps, select_corrmap, [("out", "inlist")]),
        (select_lagmap, lagplot, [("out", "lagmap_file")]),
        (select_corrmap, corrplot, [("out", "lagmap_file")]),
        
        # Set outputs
        (mask_maps, outputnode, [("out", "lagmaps")]),
        (report_hemo, ds_report_hemo, [("out_report", "in_file")]),
        (lagplot, ds_report_lagplot, [("out_report", "in_file")]),
        (corrplot, ds_report_corrplot, [("out_report", "in_file")]),
    
    ])
    return workflow

def _mask_img(input_image, mask_file):
    from nilearn.masking import apply_mask, unmask
    from pathlib import Path
    new_file = unmask(
        apply_mask(
            input_image, mask_img = mask_file), mask_file)
    out_name = Path("masked_image.nii.gz").absolute()
    new_file.to_filename(out_name)
    return str(out_name)

def _intersect_masks(in1, in2):
    from pathlib import Path
    import numpy as np
    import nibabel as nib
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
