# -*- coding: utf-8 -*-
"""
Nipype interfaces to canica and dictlearning in nilearn.decomposition
"""
import os

import numpy as np
from nilearn.decomposition import CanICA, DictLearning
from nipype.interfaces.base import (
    SimpleInterface,
    BaseInterfaceInputSpec,
    TraitedSpec,
    InputMultiPath,
    OutputMultiPath,
    traits,
    isdefined,
)

from nipype.utils.filemanip import fname_presuffix, split_filename




class _CanICAInputSpec(BaseInterfaceInputSpec):
    in_files = InputMultiPath(traits.File(desc="NifTI image file(s) from where to extract the data. \n"
                                        "If more than one, all should be spatially normalized.",
                                        exists=True, mandatory=True))
    mask = traits.File(desc="Mask to be used on data. If an instance of masker is passed, then its mask will be used.\n"
                            " If no mask is given, it will be computed automatically by a MultiNiftiMasker \n"
                            "with default parameters.",
                        exists=True)
    algorithm = traits.Enum(['canica', 'dictlearning'], desc="Desired nilearn ICA method.")
    confounds = traits.File(desc="CSV file path."
                                "This parameter is passed to nilearn.signal.clean. "
                                "Please see the related documentation for details.",
                            exists=True)
    n_components = traits.Int(desc="Number of components to extract.", default_value=20, usedefault=True)
    


class _CanICAOutputSpec(TraitedSpec):
    components_img = traits.File(desc="A nifti file with the reconstructed volume for each ICs.")
    components_img_zscored = traits.File(desc="A nifti file with the reconstructed volume for each ICs zscored.")
    ts = OutputMultiPath(traits.File(desc="Decomposition components shape: "
                                                "(number of scans, number of regions)"))
    metadata_file = traits.File(exists=True, desc= "text file containing ICA metadata")


class CanICAInterface(SimpleInterface):
    """ Nipype Interface to NiLearn methods to perform Canonical Independent Component Analysis.

    For more information look at: nilearn.decomposition.CanICA
    """
    input_spec = _CanICAInputSpec
    output_spec = _CanICAOutputSpec

    def _run_interface(self, runtime):
        from scipy.stats.mstats import zscore
        import numpy as np
        # init the estimator
        if self.inputs.algorithm == 'canica':
            self._estimator = CanICA(
                mask=self.inputs.mask,
                n_components=self.inputs.n_components,
                random_state=1,
                verbose=0,
            )

        elif self.inputs.algorithm == 'dictlearning':
            self._estimator = DictLearning(
                mask=self.inputs.mask,
                n_components=self.inputs.n_components,
                random_state=1,
                verbose=0,
            )

        # fit and transform
        ts = self._estimator.fit_transform(self.inputs.in_files)
        self._loadings = self._estimator.transform(self.inputs.in_files)
        
        self._results["components_img"] = fname_presuffix(
            self.inputs.in_files[0], suffix="_ICs", newpath=runtime.cwd,
        )
        self._results["components_img_zscored"] = fname_presuffix(
            self.inputs.in_files[0], suffix="_ICs_zscore", newpath=runtime.cwd,
        )
        self._results["ts"] = fname_presuffix(
            self.inputs.in_files[0], suffix="_ICsts" + ".npy", newpath=runtime.cwd, use_ext = False,
        )

        self._results["metadata_file"] = os.path.join(
            runtime.cwd, "ica_metadata.tsv"
        )
        
        z_scored_components = zscore(self._estimator.components_, axis = 1)
        z_scored_components_img = self._estimator.masker_.inverse_transform(z_scored_components)
        
        z_scored_components_img.to_filename(self._results["components_img_zscored"])
        self._estimator.components_img_.to_filename(self._results["components_img"])
        
        np.save(self._results["ts"], ts)

        with open(self._results["metadata_file"], "w") as f:
            f.write("\t".join(["algorithm", "n_Component", "variance_explained","\n"]))
            f.write(
                        "{}\t{}\t{}\n".format(self.inputs.algorithm, self.inputs.n_components, self._estimator.score(self.inputs.in_files, per_component=False))
                    )
        


        return runtime
    

class _SmoothInputSpec(BaseInterfaceInputSpec):
    input_image = traits.File(desc="NifTI image file(s) to smooth",
                                        exists=True, mandatory=True)
    fwhm = traits.Float(6, desc="Smoothing strength", usedefault= True)
    out_postfix = traits.Str("_smoothed", usedefault=True, desc="Postfix that is appended to all output files")
    output_image = traits.Str(desc="output file name", genfile=True, hash_files=False)

class _SmoothOutputSpec(TraitedSpec):
    output_image = traits.File(desc="Smoothed images")
    
class Smooth(SimpleInterface):
    """ Nipype Interface to NiLearn smooth image.

    For more information look at: nilearn.image.smooth_img
    """
    input_spec = _SmoothInputSpec
    output_spec = _SmoothOutputSpec
    
    def _gen_filename(self, name):
        if name == "output_image":
            output = self.inputs.output_image
            if not isdefined(output):
                _, name, ext = split_filename(self.inputs.input_image)
                output = name + self.inputs.out_postfix + ext
            return output
        return None
    
    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs["output_image"] = os.path.abspath(self._gen_filename("output_image"))
        return outputs
    
    def _run_interface(self, runtime):
        from nilearn.image import smooth_img
        self._results["output_image"] = os.path.abspath(self._gen_filename("output_image"))
        smoothed_image = smooth_img(self.inputs.input_image, self.inputs.fwhm)
        smoothed_image.to_filename(self._results["output_image"])
        return runtime
    