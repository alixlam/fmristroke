# -*- coding: utf-8 -*-
"""
Nipype interfaces to pyAnts
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
    File,
    InputMultiObject,
    isdefined,
)

from nipype.utils.filemanip import split_filename



class _ApplyTransformsInputSpec(BaseInterfaceInputSpec):
    input_image = File(
        mandatory=True,
        desc=("image to apply transformation to (generally a coregistered functional)"),
        exists=True,
    )
    output_image = traits.Str(
        "%s", desc="output file name", genfile=True, hash_files=False
    )
    out_postfix = traits.Str(
        "_trans",
        usedefault=True,
        desc=("Postfix that is appended to all output files (default = _trans)"),
    )
    reference_image = File(
        mandatory=True,
        desc="reference image space that you wish to warp INTO",
        exists=True,
    )
    interpolation = traits.Enum(
        "gaussian",
        "linear",
        "multiLabel", 
        "bSpline",
        "cosineWindowedSinc",
        "genericLabel",
        "lanczosWindowedSinc",
        "hammingWindowedSinc",
        "welchWindowedSinc",
        "nearestNeighbor",
        argstr="%s",
        usedefault=True,
    )
    transforms = InputMultiObject(
        traits.Either(File(exists=True), "identity"),
        argstr="%s",
        mandatory=True,
        desc="transform files: will be applied in reverse order. For "
        "example, the last specified transform will be applied first.",
    )    
    default_value = traits.Float(0.0, usedefault=True)
    


class _ApplyTransformsOutputSpec(TraitedSpec):
    output_image = File(exists=True, desc="Warped image")


class ApplyTransforms(SimpleInterface):
    """
    """
    input_spec = _ApplyTransformsInputSpec
    output_spec = _ApplyTransformsOutputSpec
    
    def _gen_filename(self, name):
        if name == "output_image":
            output = self.inputs.output_image
            if not isdefined(output):
                _, name, ext = split_filename(self.inputs.input_image)
                output = name + self.inputs.out_postfix + ext
            return output
        return None

    def _get_output_warped_filename(self):
        return (self._gen_filename("output_image"))

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs["output_image"] = os.path.abspath(self._gen_filename("output_image"))
        return outputs
    
    def _run_interface(self, runtime):
        from ..utils.ants import apply_transform_ants
        from niworkflows.utils.images import _copyxform
        from niworkflows import __version__
        self._results["output_image"] = os.path.abspath(self._gen_filename("output_image"))
        apply_transform_ants(fixed = self.inputs.reference_image, moving = self.inputs.input_image,
                                            transforms=self.inputs.transforms, interpolation=self.inputs.interpolation,
                                            outdir = self._results["output_image"],
                                            default_value = 0,)
        _copyxform(
            self.inputs.reference_image,
            os.path.abspath(self._gen_filename("output_image")),
            message="%s (niworkflows v%s)" % (self.__class__.__name__, __version__),
        )
        return runtime