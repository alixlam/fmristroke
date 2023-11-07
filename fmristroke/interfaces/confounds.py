
import os
import re
import shutil

import nibabel as nb
import numpy as np
import pandas as pd
from nipype import logging
from nipype.interfaces.base import (
    BaseInterfaceInputSpec,
    Directory,
    File,
    InputMultiObject,
    OutputMultiObject,
    SimpleInterface,
    TraitedSpec,
    isdefined,
    traits,
)
from nipype.utils.filemanip import fname_presuffix

from nilearn.interfaces.fmriprep.load_confounds_components import (
    _load_high_pass,
    _load_compcor,
    _load_global_signal,
    _load_wm_csf,
    _load_motion,
    _load_non_steady_state,
    _load_scrub,
)
from nilearn.interfaces.fmriprep.load_confounds_utils import _prepare_output, MissingConfound
from ..utils.confounds import _load_iclesion, _load_wm_csf_lesion

component_parameters = {
    "motion": ["motion"],
    "wm_csf": ["wm_csf"],
    "global_signal": ["global_signal"],
    "compcor": ["meta_json", "compcor", "n_compcor"],
    "scrub": ["scrub", "fd_threshold", "std_dvars_threshold"],
    "wm_csf_lesion": ["wm_csf"],
}

class _LesionMasksInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, desc="Input mask")
    roi_file = File(exists=True, desc="Input roi mask.")
    


class _LesionMasksOutputSpec(TraitedSpec):
    out_masks = OutputMultiObject(
        File(exists=True), desc="input mask intersected and combined with lesion masks, respectively"
    )


class LesionMasks(SimpleInterface):
    """Combine lesion mask with masks."""

    input_spec = _LesionMasksInputSpec
    output_spec = _LesionMasksOutputSpec

    def _run_interface(self, runtime):
        from ..utils.masks import subtract_roi, binary_union

        # mask + lesion
        mask_lesion = binary_union(self.inputs.in_file, self.inputs.roi_file)
        mask_lesionfile = self.inputs.in_file.replace(".nii.gz", "lesion.nii.gz")
        mask_lesion.to_filename(mask_lesionfile)

        
        # mask U lesion
        mask_nolesion = subtract_roi(in_mask= self.inputs.in_file, in_roi = self.inputs.roi_file)
        mask_nolesionfile = self.inputs.in_file.replace(".nii.gz", "nolesion.nii.gz")
        mask_nolesion.to_filename(mask_nolesionfile)
        
        
        self._results["out_masks"] = [mask_lesionfile, mask_nolesionfile]
        return runtime
    
class _ROIcompInputSpec(BaseInterfaceInputSpec):
    IC_components = File(exists=True, desc="Input ICA components")
    ts = File(exists=True, desc="Input timeseries of ICA components")
    ROI_mask = File(exists=True, desc="Mask of ROI")
    threshold = traits.Float(0.05, desc="Threshold to use for overlap", usedefault=True)


class _ROIcompOutputSpec(TraitedSpec):
    out_signal = OutputMultiObject(File(exists=True), desc="ts of ICs that overlap with ROI")
    metadata_file = File(exists=True, description= "Text file containing component metadata")

class ROIcomp(SimpleInterface):
    """Get ROI confounds signals from ICA"""

    input_spec = _ROIcompInputSpec
    output_spec = _ROIcompOutputSpec

    def _run_interface(self, runtime):
        from nilearn.image import resample_to_img, binarize_img, index_img, iter_img
        from sklearn.metrics import jaccard_score
        import numpy as np

        # Resample lesion
        lesion = resample_to_img(self.inputs.ROI_mask, index_img(self.inputs.IC_components, 1), interpolation='nearest').get_fdata()

        # Binarize spatial maps
        binary_spatial_maps = binarize_img(self.inputs.IC_components, 2)
        
        ICs = []
        J = []
        retained = []
        retained_ICS = []

        for IC, cur_img in enumerate(iter_img(binary_spatial_maps)):
            jaccard = jaccard_score(lesion.flatten(), cur_img.get_fdata().flatten())
            ICs.append(IC)
            J.append(jaccard)
            if jaccard > self.inputs.threshold:
                retained_ICS.append(IC)
                retained.append("True")
            else:
                retained.append("False")
            
        # save components time serie
        ts = np.load(self.inputs.ts).squeeze()
        roi_ts = ts[:, retained_ICS]
                
        output = np.vstack((np.array(['ica_lesion_'+str(i).zfill(3) for i in retained_ICS]).reshape(1,-1), roi_ts.astype(str)))
        self._results["out_signal"] = os.path.join(runtime.cwd, "components_ICs.txt")
        np.savetxt(self._results["out_signal"], output, fmt=b"%s", delimiter="\t")
        
        # Save metadata
        self._results["metadata_file"] = os.path.join(
            runtime.cwd, "components_metadata.tsv"
        )            
        with open(self._results["metadata_file"], "w") as f:
            f.write("\t".join(["Component", "Jaccard", "retained", "\n"]))
            for i in zip(np.array(ICs).astype(str), np.array(J).astype(str), retained):
                f.write(
                        "{0[0]}\t{0[1]}\t{0[2]}\n".format(i)
                    )
        
        return runtime
    
class _GatherConfoundsInputSpec(BaseInterfaceInputSpec):
    region_signals = File(exists=True, desc="Average time series from regions")
    ic_roi_signals = File(exists=True, desc="Signals from ICs overlapping with regions")
    confounds_file = File(exists=True, desc="Confounds file from fmriprep")

class _GatherConfoundsOutputSpec(TraitedSpec):
    confounds_file = File(exists=True, description="Confounds from fmriprep and confounds from roi")
    confounds_list = traits.List(traits.Str, desc='list of headers')

class GatherConfounds(SimpleInterface):
    """Gather confounds from fmriprep and confounds from roi"""

    input_spec = _GatherConfoundsInputSpec
    output_spec = _GatherConfoundsOutputSpec

    def _run_interface(self, runtime):
        combined_out, confounds_list = _gather_confounds(
            signals=self.inputs.region_signals,
            ic_roi=self.inputs.ic_roi_signals,
            confounds_fmriprep=self.inputs.confounds_file,
            newpath=runtime.cwd)
        self._results['confounds_file'] = combined_out
        self._results['confounds_list'] = confounds_list
        return runtime

class _SelectConfoundsInputSpec(BaseInterfaceInputSpec):
    pipeline = traits.Str(
        mandatory=True,
        desc="Denoising pipeline's name")
    confounds_spec = traits.Dict(
        mandatory=True,
        desc="Confounds to include"
    )
    demean =  traits.Bool(
        True,
        desc="Wether to demean confounds",
        usedefault=True,
    )
    confounds = File(
        exist=True,
        mandatory=True,
        desc="Confounds table")
    confounds_metadata = File(
        exist=True,
        mandatory=True,
        desc="Confounds description (aCompCor)")
    output_dir = Directory(
        exists=True,
        desc="Output path")

class _SelectConfoundsOutputSpec(TraitedSpec):
    selected_confounds = File(
        exists=True,
        desc="selected confounds table")
    sample_mask = traits.Either(
        None,
        File(),
        desc="shape: (number of scans - number of volumes removed,The index of the niimgs along time/fourth dimension for valid volumes for subsequent analysis."
    )


class SelectConfounds(SimpleInterface):
    """Filter confounds table according to denoising pipeline and prepare them for denoising.

    This interface reads raw confounds table and process it 
    retaining regressors of interest and creating additional regressors if 
    needed. This interface operates on single BIDS entity
    i.e. single subject, task and (optionally) session. 
    
    """
    input_spec = _SelectConfoundsInputSpec
    output_spec = _SelectConfoundsOutputSpec

    
    def _load_confounds(self):
        import pandas as pd
        import json
        with open(self.inputs.confounds_metadata, "rb") as f:   
            meta_json = json.load(f) 
        
        confounds_all = pd.read_csv(
            self.inputs.confounds, delimiter="\t", encoding="utf-8"
            )
        
        missing = {"confounds": [], "keywords": []}
        confounds_select, missing = _load_noise_component(confounds_all, "non_steady_state", missing, meta_json=meta_json)
        for component, params in self.inputs.confounds_spec.items():
            loaded_confounds, missing = _load_noise_component(
                confounds_all, component, missing, meta_json=meta_json, **params 
            )
            confounds_select = pd.concat([confounds_select, loaded_confounds], axis=1)
        
        if missing["confounds"] or missing["keywords"]:
            error_msg = (
                "The following keywords or parameters are missing: "
                + f" {missing['confounds']}"
                + f" {missing['keywords']}"
                + ". You may want to try a different denoising strategy."
            )
            raise ValueError(error_msg)

        return _prepare_output(confounds_select, self.inputs.demean)

    def _run_interface(self, runtime):
        sample_mask, selected_confounds = self._load_confounds()
        self._results["selected_confounds"] = os.path.join(runtime.cwd, f"selected_confounds_{self.inputs.pipeline}.tsv")
        selected_confounds.to_csv(self._results["selected_confounds"], sep='\t', index=False, na_rep='n/a')

        if sample_mask is not None:
            import numpy as np
            self._results["sample_mask"] = os.path.join(runtime.cwd, f"sample_mask_{self.inputs.pipeline}.txt")
            np.savetxt(self._results["sample_mask"], sample_mask)
        else:
            self._results["sample_mask"] = None
        return runtime



def _gather_confounds(
        signals=None, 
        ic_roi=None, 
        confounds_fmriprep=None,
        newpath=None):
    r"""
    Load confounds from multiple sources and concatenate them together horizontall
    Mostly based on gatherconfounds from fmriprep
    >>> pd.DataFrame({'Region signals': [0.1]}).to_csv('regions_signals.tsv', index=False, na_rep='n/a')
    >>> pd.DataFrame({'IC signals': [0.2]}).to_csv('ic_roi_signals.tsv', index=False, na_rep='n/a')
    >>> out_file, confound_list = _gather_confounds('regions_signals.tsv', 'ic_roi_signals.tsv')
    >>> confounds_list
    ['Region signals', 'IC signals']
    >>> pd.read_csv(out_file, sep='\s+', index_col=None,
    ...             engine='python')  # doctest: +NORMALIZE_WHITESPACE
        Region signals  IC signals
    0            0.1        0.2


    """
    def less_breakable(a_string):
        '''hardens the string to different envs (i.e., case insensitive, no whitespace, '#' '''
        return ''.join(a_string.split()).strip('#')

    # Taken from https://stackoverflow.com/questions/1175208/
    def camel_to_snake(name):
        s1 = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', name)
        return re.sub('([a-z0-9])([A-Z])', r'\1_\2', s1).lower()

    def _adjust_indices(left_df, right_df):
        # This forces missing values to appear at the beginning of the DataFrame
        # instead of the end
        index_diff = len(left_df.index) - len(right_df.index)
        if index_diff > 0:
            right_df.index = range(index_diff, len(right_df.index) + index_diff)
        elif index_diff < 0:
            left_df.index = range(-index_diff, len(left_df.index) - index_diff)

    all_files = []
    confounds_list = []
    for confound, name in (
        (signals, 'Global signals'),
        (ic_roi, 'IC ROI'),
        (confounds_fmriprep, 'fmriprep Confounds')
    ):
        if confound is not None and isdefined(confound):
            confounds_list.append(name)
            if os.path.exists(confound) and os.stat(confound).st_size > 0:
                all_files.append(confound)

    confounds_data = pd.DataFrame()
    for file_name in all_files:  # assumes they all have headings already
        try:
            new = pd.read_csv(file_name, sep="\t")
        except pd.errors.EmptyDataError:
            # No data, nothing to concat
            continue
        for column_name in new.columns:
            new.rename(
                columns={column_name: camel_to_snake(less_breakable(column_name))}, inplace=True
            )

        _adjust_indices(confounds_data, new)
        confounds_data = pd.concat((confounds_data, new), axis=1)

    if newpath is None:
        newpath = os.getcwd()

    combined_out = os.path.join(newpath, 'confounds.tsv')
    confounds_data.to_csv(combined_out, sep='\t', index=False, na_rep='n/a')

    return combined_out, confounds_list

def _load_noise_component(confounds_raw, component, missing, **kargs):
    """
    copied from Nilearn (https://nilearn.github.io)
    Load confound of a single noise component.

    Parameters
    ----------
    confounds_raw : :class:`pandas.DataFrame`
        The confounds loaded from the confounds file.

    component : :obj:`str`
        The noise component to be loaded. The item from the confounds list.

    missing : :obj:`dict`
        A dictionary of missing confounds and noise component keywords.

    kargs : :obj:`dict`
        Extra relevant parameters for the given `component`.

    Returns
    -------
    loaded_confounds : :class:`pandas.DataFrame`
        The confounds loaded from the confounds file for the given component.

    missing : :obj:`dict`
        A dictionary of missing confounds and noise component keywords.

    Raises
    ------
    MissingConfound
        If any of the confounds specified in the strategy are not found in the
        confounds file or confounds json file.
    """
    try:
        need_params = component_parameters.get(component)
        if need_params:
            params = {param: kargs.get(param) for param in need_params}
            loaded_confounds = eval(f"_load_{component}")(
                confounds_raw, **params
            )
        else:
            loaded_confounds = eval(f"_load_{component}")(
                confounds_raw
            )
    except MissingConfound as exception:
        missing["confounds"] += exception.params
        missing["keywords"] += exception.keywords
        loaded_confounds = pd.DataFrame()
    return loaded_confounds, missing

