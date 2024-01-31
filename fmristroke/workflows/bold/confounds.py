"""
Confounds generation
^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: init_confs_wf

"""
from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from templateflow.api import get as get_template

from ...interfaces import DerivativesDataSink


def init_confs_wf(
    mem_gb: float,
    metadata: dict,
    ica_method: str = "canica",
    freesurfer: bool = False,
    ncomp_method: str or int = "varexp",
    session_level: bool = True,
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
        TSV of all aggregated confouconf_corr_plotnds computed by fmriprep

    Outputs
    -------
    confounds_file
        TSV of all aggregated confounds
    confounds_metadata
        Confounds metadata dictionary.
    ica_plot
        Reportlet visualizing ICA components included in confounds files.
    """

    from fmriprep.interfaces.confounds import (
        FilterDropped,
        RenameACompCor,
        aCompCorMasks,
    )
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow
    from niworkflows.interfaces.confounds import ExpandModel
    from niworkflows.interfaces.images import SignalExtraction
    from niworkflows.interfaces.morphology import (
        BinaryDilation,
        BinarySubtraction,
    )
    from niworkflows.interfaces.nibabel import ApplyMask, Binarize
    from niworkflows.interfaces.patches import RobustACompCor as ACompCor
    from niworkflows.interfaces.patches import RobustTCompCor as TCompCor
    from niworkflows.interfaces.plotting import (
        CompCorVariancePlot,
        ConfoundsCorrelationPlot,
    )
    from niworkflows.interfaces.reportlets.masks import ROIsPlot
    from niworkflows.interfaces.utility import TSV2JSON, DictMerge

    from ...interfaces.confounds import GatherConfounds, LesionMasks, ROIcomp
    from ...interfaces.nilearn import CanICAInterface
    from ...interfaces.reports import ICPlot

    workflow = Workflow(name=name)

    workflow.__desc__ = f"""\\TO DO
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
                "boldmask",
                # if recomputed compcors (session level)
                "confounds_metadata_comp",
                "acompcor_masks",
                "crown_mask",
            ]
        ),
        name="outputnode",
    )

    # Resample T1w mask into BOLD space and merge with BOLD brainmask
    t1w_mask_tfm = pe.Node(
        niu.Function(function=_resample_to_img), name="t1w_mask_tfm"
    )
    t1w_mask_tfm.inputs.interpolation = "nearest"
    t1w_mask_tfm_roi = pe.Node(
        niu.Function(function=_resample_to_img), name="t1w_mask_tfm_roi"
    )
    t1w_mask_tfm_roi.inputs.interpolation = "nearest"

    intersect_mask = pe.Node(
        niu.Function(function=_boldmsk_from_T1w), name="intersect_mask"
    )

    # Generate csf, wm and combined masks
    acc_masks = pe.Node(aCompCorMasks(is_aseg=freesurfer), name="acc_masks")

    # Resample probseg maps in BOLD space via T1w-to-BOLD transform
    acc_msk_tfm = pe.MapNode(
        niu.Function(function=_resample_to_img),
        iterfield=["input_image"],
        name="acc_msk_tfm",
    )
    acc_msk_tfm.inputs.interpolation = "continuous"

    acc_msk_brain = pe.MapNode(
        ApplyMask(), name="acc_msk_brain", iterfield=["in_file"]
    )
    acc_msk_bin = pe.MapNode(
        Binarize(thresh_low=0.99), name="acc_msk_bin", iterfield=["in_file"]
    )

    # Add lesion in masks
    acc_msk_bin_lesion = pe.MapNode(
        LesionMasks(), name="acc_msk_bin_lesion", iterfield=["in_file"]
    )

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

    # Merge rois
    merge_rois = pe.Node(
        niu.Merge(2, ravel_inputs=True),
        name="merge_rois_csf",
        run_without_submitting=True,
    )

    signals = pe.Node(
        SignalExtraction(class_labels=signals_class_labels),
        name="signals",
        mem_gb=mem_gb,
    )

    ###### ICA #####
    # Automatically etimate number of components
    n_component = pe.Node(niu.Function(function=_nCompICA), name="nCompICA")
    n_component.inputs.method = ncomp_method

    # compute ICA
    canica = pe.Node(CanICAInterface(), name="canICA")
    canica.inputs.algorithm = ica_method

    # Classify components
    lesion_components = pe.Node(ROIcomp(), name="lesion_components")

    # Plot ROI ICs
    IC_plot = pe.Node(ICPlot(generate_report=True), name="ic_plots")

    ds_report_ica = pe.Node(
        DerivativesDataSink(
            desc="icaroi",
            datatype="figures",
        ),
        run_without_submitting=True,
        name="ds_report_ica",
    )

    ica_metadata_fmt = pe.Node(
        TSV2JSON(
            index_column="ica",
            drop_columns=None,
            output=None,
            enforce_case=False,
        ),
        name="ica_metadata_fmt",
    )

    lesion_conf_metadata_fmt = pe.Node(
        TSV2JSON(
            index_column="p",
            drop_columns=None,
            output=None,
            enforce_case=False,
        ),
        name="lesionconf_metadata_fmt",
    )

    mrg_confic_metadata = pe.Node(
        niu.Merge(2),
        name="merge_confound_metadata_ic",
        run_without_submitting=True,
    )

    mrg_confic_metadata2 = pe.Node(
        DictMerge(),
        name="merge_confound_metadata2_ic",
        run_without_submitting=True,
    )

    model_expand = pe.Node(
        ExpandModel(
            model_formula="(dd1(wm_nolesion + csf_lesion))^^2 + others"
        ),
        name="model_expansion",
    )

    # Redo CompCor if session level code adapted from fmriprep repo
    if session_level:
        acompcor = pe.Node(
            ACompCor(
                components_file="acompcor.tsv",
                header_prefix="a_comp_cor_",
                pre_filter="cosine",
                save_pre_filter=True,
                save_metadata=True,
                mask_names=["CSF", "WM", "combined"],
                merge_method="none",
                failure_mode="NaN",
            ),
            name="acompcor",
            mem_gb=mem_gb,
        )
        # Create the crown mask
        dilated_mask = pe.Node(BinaryDilation(), name="dilated_mask")
        subtract_mask = pe.Node(BinarySubtraction(), name="subtract_mask")

        crowncompcor = pe.Node(
            ACompCor(
                components_file="crown_compcor.tsv",
                header_prefix="edge_comp_",
                pre_filter="cosine",
                save_pre_filter=True,
                save_metadata=True,
                mask_names=["Edge"],
                merge_method="none",
                failure_mode="NaN",
                num_components=24,
            ),
            name="crowncompcor",
            mem_gb=mem_gb,
        )

        tcompcor = pe.Node(
            TCompCor(
                components_file="tcompcor.tsv",
                header_prefix="t_comp_cor_",
                pre_filter="cosine",
                save_pre_filter=True,
                save_metadata=True,
                percentile_threshold=0.02,
                failure_mode="NaN",
            ),
            name="tcompcor",
            mem_gb=mem_gb,
        )

        acompcor.inputs.variance_threshold = 0.5
        tcompcor.inputs.variance_threshold = 0.5

        # Set TR if present
        if "RepetitionTime" in metadata:
            tcompcor.inputs.repetition_time = metadata["RepetitionTime"]
            acompcor.inputs.repetition_time = metadata["RepetitionTime"]
            crowncompcor.inputs.repetition_time = metadata["RepetitionTime"]

        # Split aCompCor results into a_comp_cor, c_comp_cor, w_comp_cor
        rename_acompcor = pe.Node(RenameACompCor(), name="rename_acompcor")

        # CompCor metadata
        tcc_metadata_filter = pe.Node(
            FilterDropped(), name="tcc_metadata_filter"
        )
        acc_metadata_filter = pe.Node(
            FilterDropped(), name="acc_metadata_filter"
        )
        tcc_metadata_fmt = pe.Node(
            TSV2JSON(
                index_column="component",
                drop_columns=["mask"],
                output=None,
                additional_metadata={"Method": "tCompCor"},
                enforce_case=True,
            ),
            name="tcc_metadata_fmt",
        )
        acc_metadata_fmt = pe.Node(
            TSV2JSON(
                index_column="component",
                output=None,
                additional_metadata={"Method": "aCompCor"},
                enforce_case=True,
            ),
            name="acc_metadata_fmt",
        )
        crowncc_metadata_fmt = pe.Node(
            TSV2JSON(
                index_column="component",
                output=None,
                additional_metadata={"Method": "EdgeRegressor"},
                enforce_case=True,
            ),
            name="crowncc_metadata_fmt",
        )
        mrg_compcor = pe.Node(
            niu.Merge(3, ravel_inputs=True),
            name="mrg_compcor",
            run_without_submitting=True,
        )

        # extract tcompcor
        signals_tcompcor = pe.Node(
            SignalExtraction(class_labels=["tcompcor"]),
            name="signal_tcompcor",
            mem_gb=mem_gb,
        )

        # Generate metadata file
        mrg_conf_metadata = pe.Node(
            niu.Merge(3),
            name="merge_confound_metadata",
            run_without_submitting=True,
        )
        mrg_conf_metadata.inputs.in3 = {
            label: {"Method": "Mean"} for label in signals_class_labels
        }
        mrg_conf_metadata2 = pe.Node(
            DictMerge(),
            name="merge_confound_metadata2",
            run_without_submitting=True,
        )
        out_meta = pe.Node(niu.Function(function=_dict2json), name="out_meta")

        # Generate reportlet (ROIs)
        mrg_compcor = pe.Node(
            niu.Merge(3, ravel_inputs=True),
            name="mrg_compcor",
            run_without_submitting=True,
        )
        rois_plot = pe.Node(
            ROIsPlot(colors=["b", "magenta", "g"], generate_report=True),
            name="rois_plot",
            mem_gb=mem_gb,
        )

        ds_report_bold_rois = pe.Node(
            DerivativesDataSink(
                desc="rois", datatype="figures", dismiss_entities=("echo",)
            ),
            name="ds_report_bold_rois",
            run_without_submitting=True,
            mem_gb=0.01,
        )

        # Generate reportlet (CompCor)
        mrg_cc_metadata = pe.Node(
            niu.Merge(2),
            name="merge_compcor_metadata",
            run_without_submitting=True,
        )
        compcor_plot = pe.Node(
            CompCorVariancePlot(
                variance_thresholds=(0.5, 0.7, 0.9),
                metadata_sources=["tCompCor", "aCompCor", "crownCompCor"],
            ),
            name="compcor_plot",
        )

        ds_report_compcor = pe.Node(
            DerivativesDataSink(
                desc="compcorvar",
                datatype="figures",
                dismiss_entities=("echo",),
            ),
            name="ds_report_compcor",
            run_without_submitting=True,
            mem_gb=0.01,
        )

    # Concatenate all confounds
    concat = pe.Node(
        GatherConfounds(),
        name="concat",
        mem_gb=0.01,
        run_without_submitting=True,
    )

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
        (inputnode, signals, [("bold_t1", "in_file")]),
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
    ])

    if session_level:
        workflow.connect([
            # CompCor 
            (inputnode, acompcor, [("bold_t1", "realigned_file")]),
            (acc_msk_bin, acompcor, [("out_file", "mask_files")]),
            (acompcor, rename_acompcor, [("components_file", "components_file"),
                                        ("metadata_file", "metadata_file")]),
            # crownCompCor
            (intersect_mask, dilated_mask, [("out", "in_mask")]),
            (intersect_mask, subtract_mask, [("out", "in_subtract")]),
            (dilated_mask, subtract_mask, [("out_mask", "in_base")]),
            (inputnode, crowncompcor, [("bold_t1", "realigned_file"),]),
            (subtract_mask, crowncompcor, [("out_mask", "mask_files")]),
            (inputnode, tcompcor, [("bold_t1", "realigned_file"),
                                ("boldmask_t1", "mask_files")]),
            (inputnode, signals_tcompcor, [("bold_t1", "in_file")]),
            (tcompcor, signals_tcompcor, [("high_variance_masks", "label_files")]),
            (tcompcor, concat, [("pre_filter_file", "cosine")]),
            (signals_tcompcor, concat, [("out_file", "tcompcor")]),
            (rename_acompcor, concat, [("components_file", "acompcor")]),
            (crowncompcor, concat, [("components_file", "crowncompcor")]),

            # Reporting
            (inputnode, rois_plot, [("bold_t1", "in_file"),
                                ("boldmask_t1", "in_mask")]),
            (tcompcor, mrg_compcor, [("high_variance_masks", "in1")]),
            (acc_msk_bin, mrg_compcor, [(("out_file", _last), "in2")]),
            (subtract_mask, mrg_compcor, [("out_mask", "in3")]),
            (mrg_compcor, rois_plot, [("out", "in_rois")]),
            (rois_plot, ds_report_bold_rois, [("out_report", "in_file")]),
            (tcompcor, mrg_cc_metadata, [("metadata_file", "in1")]),
            (acompcor, mrg_cc_metadata, [("metadata_file", "in2")]),
            (crowncompcor, mrg_cc_metadata, [("metadata_file", "in3")]),
            (mrg_cc_metadata, compcor_plot, [("out", "metadata_files")]),
            (compcor_plot, ds_report_compcor, [("out_file", "in_file")]),

            # Outputs
            (tcompcor, tcc_metadata_filter, [("metadata_file", "in_file")]),
            (tcc_metadata_filter, tcc_metadata_fmt, [("out_file", "in_file")]),
            (rename_acompcor, acc_metadata_filter, [("metadata_file", "in_file")]),
            (acc_metadata_filter, acc_metadata_fmt, [("out_file", "in_file")]),
            (crowncompcor, crowncc_metadata_fmt, [("metadata_file", "in_file")]),
            (tcc_metadata_fmt, mrg_conf_metadata, [("output", "in1")]),
            (acc_metadata_fmt, mrg_conf_metadata, [("output", "in2")]),
            (crowncc_metadata_fmt, mrg_conf_metadata, [("output", "in3")]),
            (mrg_conf_metadata, mrg_conf_metadata2, [("out", "in_dicts")]),
            (mrg_conf_metadata2, out_meta, [("out_dict", "dicti")]),
            (out_meta, outputnode, [("out", "confounds_metadata_comp")]),
            (acc_msk_bin, outputnode, [("out_file", "acompcor_masks")]),
            (subtract_mask, outputnode, [("out_mask", "crown_mask")]),

        ])

    workflow.connect([
        # Expand model with derivatives
        (concat, model_expand, [("confounds_file", "confounds_file")]),

        # Confounds metadata
        (canica, ica_metadata_fmt, [("metadata_file", "in_file")]),
        (lesion_components, lesion_conf_metadata_fmt, [("metadata_file", "in_file")]),
        (ica_metadata_fmt, mrg_confic_metadata, [("output", "in1")]),
        (lesion_conf_metadata_fmt, mrg_confic_metadata, [("output", "in2")]),
        (mrg_confic_metadata, mrg_confic_metadata2, [("out", "in_dicts")]),

        # Set outputs
        (lesion_components, IC_plot, [("out_signal", "in_ts")]),
        (inputnode, IC_plot, [("boldref_t1", "in_file")]),
        (canica, IC_plot, [("components_img_zscored", "in_ICs")]),
        (t1w_mask_tfm_roi, IC_plot, [("out", "in_mask")]),
        (model_expand, outputnode, [("confounds_file", "confounds_file")]),
        (mrg_confic_metadata2, outputnode, [("out_dict", "confounds_metadata")]),
        (intersect_mask, outputnode, [("out", "boldmask")]),
        (IC_plot, ds_report_ica, [("out_report", "in_file")])

    ])

    return workflow


def init_carpetplot_wf(
    mem_gb: float, metadata: dict, name: str = "bold_carpet_wf"
):
    """
    Adapted from fmriprep `https://fmriprep.org/en/stable/`

    Build a workflow to generate *carpet* plots.

    Resamples the MNI parcellation (ad-hoc parcellation derived from the
    Harvard-Oxford template and others).

    Parameters
    ----------
    mem_gb : :obj:`float`
        Size of BOLD file in GB - please note that this size
        should be calculated after resamplings that may extend
        the FoV
    metadata : :obj:`dict`
        BIDS metadata for BOLD file
    name : :obj:`str`
        Name of workflow (default: ``bold_carpet_wf``)

    Inputs
    ------
    bold
        BOLD image, after the prescribed corrections (STC, HMC and SDC)
        when available.
    bold_mask
        BOLD series mask
    confounds_file
        TSV of all aggregated confounds
    std2anat_xfm
        ANTs-compatible affine-and-warp transform file
    crown_mask
        Mask of brain edge voxels
    acompcor_mask
        Mask of deep WM+CSF

    Outputs
    -------
    out_carpetplot
        Path of the generated SVG file

    """
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow

    from ...interfaces.confounds import FMRISummary
    from ...interfaces.pyants import ApplyTransforms

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "bold",
                "bold_mask",
                "confounds_file",
                "std2anat_xfm",
                "crown_mask",
                "acompcor_mask",
            ]
        ),
        name="inputnode",
    )

    outputnode = pe.Node(
        niu.IdentityInterface(fields=["out_carpetplot"]), name="outputnode"
    )

    # Carpetplot and confounds plot
    conf_plot = pe.Node(
        FMRISummary(
            tr=metadata["RepetitionTime"],
            confounds_list=[
                ("global_signal", None, "GS"),
                ("csf", None, "GSCSF"),
                ("white_matter", None, "GSWM"),
                ("std_dvars", None, "DVARS"),
                ("framewise_displacement", "mm", "FD"),
            ],
        ),
        name="conf_plot",
        mem_gb=mem_gb,
    )
    ds_report_bold_conf = pe.Node(
        DerivativesDataSink(
            desc="carpetplot",
            datatype="figures",
            extension="svg",
            dismiss_entities=("echo",),
        ),
        name="ds_report_bold_conf",
        run_without_submitting=True,
        mem_gb=0.01,
    )

    parcels = pe.Node(
        niu.Function(function=_carpet_parcellation), name="parcels"
    )

    # Warp segmentation into T1 space
    resample_parc = pe.Node(
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

    workflow = Workflow(name=name)

    # fmt:off
    workflow.connect([
        (inputnode, resample_parc, [("bold_mask", "reference_image")]),
        (inputnode, parcels, [("crown_mask", "crown_mask")]),
        (inputnode, parcels, [("acompcor_mask", "acompcor_mask")]),
        (inputnode, conf_plot, [("bold", "in_nifti"),
                                ("confounds_file", "confounds_file"),]),
        (inputnode, resample_parc, [("std2anat_xfm", "transforms")]),
        (resample_parc, parcels, [("output_image", "segmentation")]),
        (parcels, conf_plot, [("out", "in_segm")]),
        (conf_plot, ds_report_bold_conf, [("out_file", "in_file")]),
        (conf_plot, outputnode, [("out_file", "out_carpetplot")]),
    ])
    # fmt:on
    return workflow


def _boldmsk_from_T1w(T1wmask, bold):
    from pathlib import Path

    from nilearn.masking import compute_epi_mask, intersect_masks

    mask_fov = compute_epi_mask(bold)
    out = intersect_masks([T1wmask, mask_fov], 1)
    out_name = Path("mask_intersect.nii.gz").absolute()
    out.to_filename(out_name)
    return str(out_name)


def _nCompICA(bold, mask, method):
    assert method in [
        "varexp",
        "aic",
        "kic",
        "mdl",
    ], "Method is not implemented"

    import numpy as np
    from nilearn.maskers import NiftiMasker
    from sklearn.utils.extmath import randomized_svd

    masker = NiftiMasker(
        mask, standardize=True, detrend=True, smoothing_fwhm=6
    )
    X = masker.fit_transform(bold)

    # fit svd
    p, N = X.shape
    U, S, V = randomized_svd(
        X,
        n_components=p,
    )

    # Compute metrics

    # Explained variance
    if method == "varexp":
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
            LH = np.log(
                np.prod(np.power(eigenvalues[k:], 1 / (p - k)))
                / np.mean(eigenvalues[k:])
            )
            mlh = 0.5 * N * (p - k) * LH
            df = 1 + 0.5 * k * (2 * p - k + 1)
            aic[k_idx] = (-2 * mlh) + (2 * df)
            kic[k_idx] = (-2 * mlh) + (3 * df)
            mdl[k_idx] = -mlh + (0.5 * df * np.log(N))

        itc = np.row_stack([aic, kic, mdl])
        dlap = np.diff(itc, axis=1)

        for i, (criterion, id) in enumerate(zip(itc, ["aic", "kic", "mdl"])):
            a = np.where(dlap[i, :] > 0)[0] + 1
            if len(a) == 0:
                n = p
            else:
                n = a[0]

        metrics[id] = n
        return metrics[method]


def _get_zooms(in_file):
    import nibabel as nb

    return tuple(nb.load(in_file).header.get_zooms()[:3])


def _resample_to_img(input_image, reference_image, interpolation):
    from pathlib import Path

    from nilearn.image import resample_to_img

    new_img = resample_to_img(
        input_image, reference_image, interpolation=interpolation
    )
    out_name = Path("mask_resamp.nii.gz").absolute()
    new_img.to_filename(out_name)
    return str(out_name)


def _last(inlist):
    return inlist[-1]


def _dict2json(dicti):
    from json import dumps
    from pathlib import Path

    out_name = Path("confounds_metadata.json").absolute()
    out_name.write_text(dumps(dicti, sort_keys=True, indent=2))
    return str(out_name)


def _carpet_parcellation(segmentation, crown_mask, acompcor_mask):
    """Generate the union of two masks."""
    from pathlib import Path

    import nibabel as nb
    import numpy as np

    img = nb.load(segmentation)

    lut = np.zeros((256,), dtype="uint8")
    lut[100:201] = 1
    lut[30:99] = 2  # dGM
    lut[1:11] = 3  # WM+CSF
    lut[255] = 5  # Cerebellum
    # Apply lookup table
    seg = lut[np.uint16(img.dataobj)]
    seg[np.bool_(nb.load(crown_mask).dataobj)] = 6
    # Separate deep from shallow WM+CSF
    seg[np.bool_(nb.load(acompcor_mask).dataobj)] = 4

    outimg = img.__class__(seg.astype("uint8"), img.affine, img.header)
    outimg.set_data_dtype("uint8")
    out_file = Path("segments.nii.gz").absolute()
    outimg.to_filename(out_file)
    return str(out_file)
