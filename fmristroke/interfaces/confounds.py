
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
    pipeline = traits.Either(
        traits.Dict(), File,
        mandatory=True,
        desc="Denoising pipeline")
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
    selected_confounds_metadata = File(
        exists=True,
        desc="Summary of selected confounds"
    )


# class Confounds(SimpleInterface):
#     """Filter confounds table according to denoising pipeline.

#     This interface reads raw confounds table (fmriprep output) and process it 
#     retaining regressors of interest and creating additional regressors if 
#     needed. Additionally, it creates summary file containing all relevant 
#     information about confounds. This interface operates on single BIDS entity
#     i.e. single subject, task and (optionally) session. 
    
#     Summary contains fields:
#         'mean_fd': 
#             Mean framewise displacement.
#         'max_fd': 
#             Highest recorded framewise displacement.
#         'n_conf': 
#             Total number of confounds included.
#         'include':
#             Decision about subject inclusion in connectivity analysis based on
#             three criteria: (1) mean framewise displacement is lower than 
#             specified by the pipeline, (2) max framewise displacement did not 
#             exceed 5mm and (3) percentage of outlier scans did not exceed 20%. 
#             Note that if spikes strategy is not specified, include flag defaults
#             to True.
#         'n_spikes':
#             Number of outlier scans (only if spikes strategy is specified).
#         'perc_spikes':
#             Percentage of outlier scans (only if spikes strategy is specified).
#     """
#     input_spec = _SelectConfoundsInputSpec
#     output_spec = _SelectConfoundsOutputSpec

#     def _keep(self, regressor_names):
#         """
#         Copies selected regressors from confounds to selcted_confounds.
#         """
#         if regressor_names:
#             self.conf_prep = pd.concat((
#                 self.conf_prep,
#                 self.conf_raw[regressor_names]
#             ), axis=1)

#     def _load_tissue_signals(self):
#         tissue_regressors = []
#         for confound, setting in self.inputs.pipeline['confounds'].items():
            
#             if confound in ('white_matter', 'csf', 'global_signal'):
#                 for transform, include in setting.items():
#                     if transform == 'raw' and include:
#                         tissue_regressors.append(confound)
#                     elif include:
#                         tissue_regressors.append(f'{confound}_{transform}')
        
#         self._keep(tissue_regressors)

#     def _load_motion_parameters(self):
#         hmp_regressors = []
#         hmp_names = [f'{type_}_{axis}' 
#                     for type_ in ('trans', 'rot') 
#                     for axis in ('x', 'y', 'z')]

#         setting = self.inputs.pipeline['confounds']['motion']

#         for transform, include in setting.items():
#             if transform == 'raw' and include:
#                 hmp_regressors.extend(hmp_names)
#             elif include:
#                 hmp_regressors.extend(f'{hmp}_{transform}' for hmp in hmp_names)
        
#         self._keep(hmp_regressors)

#     def _load_compcors(self):
#         if not self.inputs.pipeline['confounds']['compcor']:
#             return

#         compcor_regressors = []
#         for mask in ('CSF', 'WM'):
#             acompcors = {
#                 (name, dict_['VarianceExplained']) 
#                 for name, dict_ in self.conf_json.items()
#                 if dict_.get('Retained') and dict_.get('Mask') == mask 
#                 }
#             acompcors = sorted(acompcors, key=lambda tpl: tpl[1], reverse=True)
#             acompcor_regressors.extend(acompcor[0] for acompcor in acompcors[:5])

#         self._keep(acompcor_regressors)
#     def _load_high_pass(self):
#         if not slef.inputs.pipeline['high_pass']:
#             return

#     def _load_non_steady_states(self):
         
        
#     def _create_spike_regressors(self):
#         if not self.inputs.pipeline['spikes']:
#             return

#         fd_th = self.inputs.pipeline['spikes']['fd_th']
#         dvars_th = self.inputs.pipeline['spikes']['dvars_th']

#         outliers = (self.conf_raw['framewise_displacement'] > fd_th) \
#                  | (self.conf_raw['std_dvars'] > dvars_th) 
#         outliers = list(outliers[outliers].index)

#         if outliers:
#             spikes = np.zeros((self.n_volumes, len(outliers)))
#             for i, outlier in enumerate(outliers):
#                 spikes[outlier, i] = 1.
                
#             conf_spikes = pd.DataFrame(
#                 data=spikes, 
#                 columns=[f'motion_outlier_{i:02}' for i in range(len(outliers))]
#                 )

#             self.conf_prep = pd.concat((
#                 self.conf_prep,
#                 conf_spikes,
#             ),
#             axis=1)
        
#         self.n_spikes = len(outliers)

#     def _create_summary_dict(self, subject: str, session: str, task: str, run: str):
#         self.conf_summary = {
#             'subject': subject,
#             'task': task,
#             'mean_fd': self.conf_raw["framewise_displacement"].mean(),
#             'max_fd': self.conf_raw["framewise_displacement"].max(),
#             'n_conf': len(self.conf_prep.columns),
#             'include': self._inclusion_check()
#         }

#         if self.inputs.pipeline['spikes']:
#             self.conf_summary['n_spikes'] = self.n_spikes 
#             self.conf_summary['perc_spikes'] = self.n_spikes / self.n_volumes * 100

#         if session:
#             self.conf_summary['session'] = session
#         if run:
#             self.conf_summary['run'] = run

#     def _inclusion_check(self):
#         '''Decide if subject should be included in connectivity analysis'''
#         if not self.inputs.pipeline['spikes']:
#             return True

#         mean_fd = self.conf_raw['framewise_displacement'].mean()
#         max_fd = self.conf_raw['framewise_displacement'].max()
#         fd_th = self.inputs.pipeline['spikes']['fd_th']

#         if mean_fd > fd_th or max_fd > 5 or self.n_spikes / self.n_volumes > 0.2:
#             return False
#         return True

#     def _run_interface(self, runtime):

#         # Setup useful properties
#         self.conf_raw = pd.read_csv(self.inputs.conf_raw, sep='\t')
#         with open(self.inputs.conf_json, 'r') as json_file:
#             self.conf_json = json.load(json_file)
#         self.n_volumes = len(self.conf_raw)
#         self.conf_prep = pd.DataFrame()

#         # entities
#         entities = parse_file_entities_with_pipelines(self.inputs.conf_raw)

#         # Create preprocessed confounds step-by-step
#         self._filter_motion_parameters()
#         self._filter_tissue_signals()
#         self._filter_acompcors()
#         self._create_spike_regressors()
#         self._create_summary_dict(
#             subject=entities.get('subject'), task=entities.get('task'),
#             session=entities.get('session'), run=entities.get('run'))

#         # Store output
#         entities['pipeline'] = self.inputs.pipeline['name']
#         conf_prep = join(self.inputs.output_dir, build_path(entities, self.conf_prep_pattern, False))
#         conf_summary = join(self.inputs.output_dir, build_path(entities, self.conf_summary_pattern, False))
#         self.conf_prep.to_csv(conf_prep, sep='\t', index=False, na_rep=0)
#         with open(conf_summary, 'w') as f:
#             json.dump(self.conf_summary, f)
#         self._results['conf_prep'] = conf_prep
#         self._results['conf_summary'] = conf_summary
#         return runtime



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


