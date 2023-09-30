
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
                
        output = np.vstack((np.array(retained_ICS).reshape(1,-1), roi_ts.astype(str)))
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


