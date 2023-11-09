from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe

from ...interfaces import DerivativesDataSink

def init_confs_wf(
    mem_gb: float,
    metadata: dict,
    ica_method : str = 'canica',
    freesurfer: bool = False,
    ncomp_method: str or int = "varexp", 
    name: str = "confounds_wf",
):
    """
    Build a workflow to generate and write out confounding signals.

    This workflow calculates roi related confounds based on the ones obtained by fmriprep for a BOLD series from stroke patients, and aggregates them
    into a :abbr:`TSV (tab-separated value)` file, for use as nuisance
    regressors in a :abbr:`GLM (general linear model)`.
    The following confounds are added, with column headings in parentheses:

    #. Region-wise average signal excluding lesion signal (``csf_nolesion``, ``white_matter_nolesion``, ``combined_nolesion``)
    #. Region-wise average signal in roi (``lesion``)
    #. Region-wise average signal including roi (``csf_lesion``, ``combined_lesion``)
    #. ICA based lesion confounds


    Parameters
    ----------
    mem_gb : :obj:`float`
        Size of BOLD file in GB - please note that this size
        should be calculated after resamplings that may extend
        the FoV
    metadata : :obj:`dict`
        BIDS metadata for BOLD file
    name : :obj:`str`
        Name of workflow (default: ``bold_confs_wf``)

    Inputs
    ------
    bold
        BOLD image, after the prescribed corrections (STC, HMC and SDC)
        when available.
    t1w_mask
        Mask of the skull-stripped template image
    t1w_mask_lesion
        Mask of the lesion in T1w space
    t1w_tpms
        List of tissue probability maps in T1w space
    t1_bold_xform
        Affine matrix that maps the T1w space into alignment with
        the native BOLD space
    confounds_file
        TSV of all aggregated confounds computed by fmriprep

    Outputs
    -------
    confounds_file
        TSV of all aggregated confounds
    confounds_metadata
        Confounds metadata dictionary.
    ica_plot 
        Reportlet visualizing ICA components included in confounds files.
    """
    
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow
    from niworkflows.interfaces.images import SignalExtraction
    from niworkflows.interfaces.morphology import BinaryDilation, BinarySubtraction
    from niworkflows.interfaces.nibabel import ApplyMask, Binarize
    from niworkflows.interfaces.utility import TSV2JSON, DictMerge
    from niworkflows.interfaces.confounds import ExpandModel
    
    from fmriprep.interfaces.confounds import aCompCorMasks
    from ...interfaces.confounds import LesionMasks, ROIcomp, GatherConfounds
    from ...interfaces.nilearn import CanICAInterface
    from ...interfaces.reports import ICPlot
    
    workflow = Workflow(name=name)
    
    workflow.__desc__ = f"""\TO DO
        """
    
    
    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "bold_t1",
                "boldref_t1",
                "boldmask_t1",
                "t1w_mask",
                "roi",
                "t1w_tpms",
                "confounds_file",
            ]
        ),
        name="inputnode",
    )
    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "confounds_file",
                "confounds_metadata",
                "boldmask"]
        ),
        name="outputnode",
    )
    
    # Resample T1w mask into BOLD space and merge with BOLD brainmask
    t1w_mask_tfm = pe.Node(
        niu.Function(function=_resample_to_img), name="t1w_mask_tfm")
    t1w_mask_tfm.inputs.interpolation = 'nearest'
    t1w_mask_tfm_roi = pe.Node(
        niu.Function(function=_resample_to_img), name="t1w_mask_tfm_roi")
    t1w_mask_tfm_roi.inputs.interpolation = 'nearest'

    
    intersect_mask = pe.Node(niu.Function(function=_boldmsk_from_T1w), name="intersect_mask")
    
    # Generate csf, wm and combined masks
    acc_masks = pe.Node(aCompCorMasks(is_aseg=freesurfer), name="acc_masks")

    # Resample probseg maps in BOLD space via T1w-to-BOLD transform
    acc_msk_tfm = pe.MapNode(
        niu.Function(function=_resample_to_img), 
        iterfield=["input_image"],
        name="acc_msk_tfm",)
    acc_msk_tfm.inputs.interpolation = "continuous"
    
    acc_msk_brain = pe.MapNode(ApplyMask(), name="acc_msk_brain", iterfield=["in_file"])
    acc_msk_bin = pe.MapNode(Binarize(thresh_low=0.99), name="acc_msk_bin", iterfield=["in_file"])
    
    # Add lesion in masks
    acc_msk_bin_lesion = pe.MapNode(LesionMasks(), name="acc_msk_bin_lesion", iterfield=["in_file"])
    
    
    # Global and segment regressors
    signals_class_labels = [
        "lesion",
        "csf_lesion",
        "csf_nolesion",
        "white_matter_lesion",
        "white_matter_nolesion",
        "csf_wm_lesion",
        "csf_wm_nolesion",
    ]
        
    
    #Merge rois
    merge_rois = pe.Node(
        niu.Merge(2, ravel_inputs=True), name="merge_rois_csf", run_without_submitting=True
    )
    
    signals = pe.Node(
        SignalExtraction(class_labels=signals_class_labels), name="signals", mem_gb=mem_gb
    )
    
    ###### ICA ##### 
    # Automatically etimate number of components
    n_component = pe.Node(niu.Function(function = _nCompICA), name="nCompICA")
    n_component.inputs.method = ncomp_method
    
    # compute ICA
    canica = pe.Node(CanICAInterface(), name= "canICA")
    canica.inputs.algorithm = ica_method
    
    # Classify components
    lesion_components = pe.Node(ROIcomp(), name="lesion_components") 
    
    # Plot ROI ICs
    IC_plot = pe.Node(ICPlot(generate_report=True), name='ic_plots')
    
    ds_report_ica = pe.Node(DerivativesDataSink(
        desc="icaroi",
        datatype="figures",),
        run_without_submitting=True,
        name="ds_report_ica")
    
    ica_metadata_fmt = pe.Node(TSV2JSON(
        index_column="ica",
        drop_columns=None,
        output=None,
        enforce_case=False,),
        name = "ica_metadata_fmt",
        )
    
    lesion_conf_metadata_fmt = pe.Node(TSV2JSON(
        index_column="p",
        drop_columns=None,
        output=None,
        enforce_case=False,),
        name = "lesionconf_metadata_fmt",
        )
    
    mrg_conf_metadata = pe.Node(
        niu.Merge(2), name="merge_confound_metadata", run_without_submitting=True
    )
    
    mrg_conf_metadata2 = pe.Node(
        DictMerge(), name="merge_confound_metadata2", run_without_submitting=True
    )
    
    model_expand = pe.Node(
        ExpandModel(model_formula="(dd1(wm_nolesion + csf_lesion))^^2 + others"),
        name="model_expansion",
    )
    
    # Concatenate all confounds
    concat = pe.Node(GatherConfounds(), name= "concat", mem_gb=0.01, run_without_submitting=True)
    
    # fmt:off
    workflow.connect([
        # Brain masks
        (inputnode, t1w_mask_tfm, [("t1w_mask", "input_image"),
                                    ("boldmask_t1", "reference_image"),]),
        (inputnode, intersect_mask, [("bold_t1", "bold")]),
        (t1w_mask_tfm, intersect_mask, [("out", "T1wmask")]),
        
        # CSF, WM and combined masks with lesions
        (inputnode, acc_masks, [("t1w_tpms", "in_vfs"),
                                (("bold_t1", _get_zooms), "bold_zooms")]),
        (intersect_mask, acc_msk_tfm, [("out", "reference_image")]),
        (intersect_mask, acc_msk_brain, [("out", "in_mask")]),
        (acc_masks, acc_msk_tfm, [("out_masks", "input_image")]),
        (acc_msk_tfm, acc_msk_brain, [("out", "in_file")]),
        (acc_msk_brain, acc_msk_bin, [("out_file", "in_file")]),
        (acc_msk_bin, acc_msk_bin_lesion, [("out_file", "in_file")]),
        
        # Lesion masks
        (inputnode, t1w_mask_tfm_roi, [("roi", "input_image"),
                                        ("boldmask_t1", "reference_image"),]),
        (t1w_mask_tfm_roi, acc_msk_bin_lesion, [("out", "roi_file")]),

        # Mean region signals
        (t1w_mask_tfm_roi, merge_rois, [("out", "in1")]),
        (acc_msk_bin_lesion, merge_rois, [("out_masks", "in2")]),
        (inputnode, signals,[("bold_t1", "in_file")]),
        (merge_rois, signals, [("out", "label_files")]),


        # ICA
        (inputnode, n_component, [("bold_t1", "bold")]),
        (intersect_mask, n_component, [("out", "mask")]),
        (n_component, canica, [("out", "n_components")]),
        (inputnode, canica, [("bold_t1", "in_files")]),
        (intersect_mask, canica, [("out", "mask")]),
        (canica, lesion_components, [("components_img_zscored", "IC_components")]),
        (canica, lesion_components, [("ts", "ts")]),
        (t1w_mask_tfm_roi, lesion_components, [("out", "ROI_mask")]),
        

        # Collate computed confounds together
        (signals, concat, [("out_file", "region_signals")]),
        (lesion_components, concat, [("out_signal", "ic_roi_signals")]),
        (inputnode, concat, [("confounds_file", "confounds_file")]),
        
        # Expand model with derivatives
        (concat, model_expand, [("confounds_file", "confounds_file")]),
        
        # Confounds metadata
        (canica, ica_metadata_fmt, [("metadata_file", "in_file")]),
        (lesion_components, lesion_conf_metadata_fmt, [("metadata_file", "in_file")]),
        (ica_metadata_fmt, mrg_conf_metadata, [("output", "in1")]),
        (lesion_conf_metadata_fmt, mrg_conf_metadata, [("output", "in2")]),
        (mrg_conf_metadata, mrg_conf_metadata2, [("out", "in_dicts")]),
        
        # Set outputs
        (lesion_components, IC_plot, [("out_signal","in_ts")]),
        (inputnode, IC_plot, [("boldref_t1", "in_file")]),
        (canica, IC_plot, [("components_img_zscored", "in_ICs")]),
        (t1w_mask_tfm_roi, IC_plot, [("out", "in_mask")]),
        (model_expand, outputnode, [("confounds_file", "confounds_file")]),
        (mrg_conf_metadata2, outputnode, [("out_dict", "confounds_metadata")]),
        (intersect_mask, outputnode, [("out", "boldmask")]),
        (IC_plot, ds_report_ica, [("out_report", "in_file")])
    
    ])
    return workflow

def _boldmsk_from_T1w(T1wmask, bold):
    from nilearn.masking import compute_epi_mask, intersect_masks
    from pathlib import Path
    mask_fov = compute_epi_mask(bold)
    out = intersect_masks([T1wmask, mask_fov], 1)
    out_name = Path("mask_intersect.nii.gz").absolute()
    out.to_filename(out_name)
    return str(out_name)

def _nCompICA(bold, mask, method):
    assert method in ["varexp", "aic", "kic", "mdl"], "Method is not implemented"
    
    from nilearn.maskers import NiftiMasker
    from sklearn.utils.extmath import randomized_svd
    import numpy as np
    masker = NiftiMasker(mask, standardize=True, detrend=True, smoothing_fwhm=6)
    X  = masker.fit_transform(bold)
    
    # fit svd
    p, N = X.shape
    U, S, V = randomized_svd(X, n_components = p,)

    # Compute metrics
    
    # Explained variance
    if method == 'varexp':
        explained_variance_ = (S**2) / (p - 1)
        total_var = explained_variance_.sum()
        explained_variance_ratio_ = explained_variance_ / total_var
        cumsum_varexp = np.cumsum(explained_variance_ratio_)
        n_comp_varexp_90 = np.where(cumsum_varexp >= 0.90)[0][0] + 1
        
        return n_comp_varexp_90
    
    # other metrics
    else:
        metrics = {}
        aic = np.zeros(p - 1)
        kic = np.zeros(p - 1)
        mdl = np.zeros(p - 1)
        eigenvalues = S

        for k_idx, k in enumerate(np.arange(1, p)):
            LH = np.log(np.prod(np.power(eigenvalues[k:], 1 / (p - k))) / np.mean(eigenvalues[k:]))
            mlh = 0.5 * N * (p - k) * LH
            df = 1 + 0.5 * k * (2 * p - k + 1)
            aic[k_idx] = (-2 * mlh) + (2 * df)
            kic[k_idx] = (-2 * mlh) + (3 * df)
            mdl[k_idx] = -mlh + (0.5 * df * np.log(N))
            
        itc = np.row_stack([aic, kic, mdl])
        dlap = np.diff(itc, axis=1)
        
        for i, (criterion, id)  in enumerate(zip(itc, ['aic', 'kic', 'mdl'])) :
            a = np.where(dlap[i, :] > 0)[0] + 1
            if len(a) == 0 : 
                n = p
            else : 
                n = a[0]

        metrics[id] = n
        return metrics[method]

def _get_zooms(in_file):
    import nibabel as nb

    return tuple(nb.load(in_file).header.get_zooms()[:3])

def _resample_to_img(input_image, reference_image, interpolation):
    from pathlib import Path
    from nilearn.image import resample_to_img
    new_img = resample_to_img(input_image, reference_image, interpolation=interpolation)
    out_name = Path("mask_resamp.nii.gz").absolute()
    new_img.to_filename(out_name)
    return str(out_name)