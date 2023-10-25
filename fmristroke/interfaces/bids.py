"""Interfaces for handling BIDS-like neuroimaging structures."""
from collections import defaultdict
from contextlib import suppress
from json import dumps, loads
from pathlib import Path
import shutil
import os
from pkg_resources import resource_filename as _pkgres
import re

import nibabel as nb
import numpy as np

from nipype import logging
from nipype.interfaces.base import (
    traits,
    isdefined,
    Undefined,
    TraitedSpec,
    BaseInterfaceInputSpec,
    DynamicTraitedSpec,
    File,
    Directory,
    InputMultiObject,
    OutputMultiObject,
    Str,
    SimpleInterface,
)
from nipype.interfaces.io import add_traits
import templateflow as tf

class _BIDSDerivativeDataGrabberInputSpec(BaseInterfaceInputSpec):
    bold_derivatives = traits.Dict(Str, traits.Any)
    anat_derivatives = traits.Dict(Str, traits.Any)
    subject_id = Str()


class _BIDSDerivativeDataGrabberOutputSpec(TraitedSpec):
    out_dict = traits.Dict(desc="output data structure")
    bold_t1 = OutputMultiObject(desc="output functional images in anatomical space")
    boldref_t1 = OutputMultiObject(desc="output functional reference images in anatomical space")
    boldmask_t1 = OutputMultiObject(desc="output bold mask in t1 space")
    confounds_file = OutputMultiObject(desc="confounds timeseries")
    t1w_preproc = OutputMultiObject(desc="output T1w preprocessed images")
    t1w_mask = OutputMultiObject(desc="output T1w brain mask ")
    t1w_dseg = OutputMultiObject(desc="output T1w dseg")
    t1w_aseg = OutputMultiObject(desc="output T1w aseg")
    t1w_aparc = OutputMultiObject(desc="output T1w aparc")
    t1w_tpms = OutputMultiObject(desc="output T1w tissue probability masks")
    std_preproc = traits.List(desc="output T1w preprocessed images in std spaces")
    std_mask = traits.List(desc="output T1w brain mask in std spaces")
    std_dseg = traits.List(desc="output T1w dseg in std spaces")
    std_aseg = traits.List(desc="output T1w aseg in std spaces")
    std_aparc = traits.List(desc="output T1w aparc in std spaces")
    std_tpms = traits.List(desc="output T1w tissue probability masks in std spaces")
    anat2std_xfm = traits.List(desc=" Transform from antomical space to std spaces")
    std2anat_xfm = traits.List(desc=" transform from std spaces to anatomical space")
    template = traits.List(desc="Templates")
    t1w2fsnative_xfm = OutputMultiObject(desc="t1w to fs native transform")
    fsnative2t1w_xfm = OutputMultiObject(desc="fs native to t1w transform")
    surfaces = OutputMultiObject(desc="Surfaces from freesurfer")
    morphometrics = OutputMultiObject(desc="morphometrics (freesurfer output)")
    anat_ribbon = OutputMultiObject(desc="anotomical ribbon (freesurfer output)")



class BIDSDerivativeDataGrabber(SimpleInterface):
    """
    Collect derivative files from a BIDS directory structure.

    .. testsetup::

        >>> data_dir_canary()

    >>> bids_src = _BIDSDerivativeDataGrabber()
    >>> bids_src.inputs.bold_derivatives = bids_collect_bold_derivatives(
    ...     str(datadir / 'ds114'), '01')
    >>> bids_src.inputs.subject_id = '01'
    >>> res = bids_src.run()
    >>> res.outputs.t1w  # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    ['.../ds114/sub-01/ses-retest/anat/sub-01_ses-retest_T1w.nii.gz',
    '.../ds114/sub-01/ses-test/anat/sub-01_ses-test_T1w.nii.gz']

    """

    input_spec = _BIDSDerivativeDataGrabberInputSpec
    output_spec = _BIDSDerivativeDataGrabberOutputSpec

    def _run_interface(self, runtime):
        bids_dict = self.inputs.bold_derivatives
        bids_dict.update(self.inputs.anat_derivatives)
        self._results["out_dict"] = bids_dict
        self._results.update(bids_dict)
        return runtime
    
class _DerivativesDataSinkInputSpec(DynamicTraitedSpec, BaseInterfaceInputSpec):
    base_directory = traits.Directory(
        desc="Path to the base directory for storing data."
    )
    check_hdr = traits.Bool(True, usedefault=True, desc="fix headers of NIfTI outputs")
    compress = InputMultiObject(
        traits.Either(None, traits.Bool),
        usedefault=True,
        desc="whether ``in_file`` should be compressed (True), uncompressed (False) "
        "or left unmodified (None, default).",
    )
    data_dtype = Str(
        desc="NumPy datatype to coerce NIfTI data to, or `source` to"
        "match the input file dtype"
    )
    dismiss_entities = InputMultiObject(
        traits.Either(None, Str),
        usedefault=True,
        desc="a list entities that will not be propagated from the source file",
    )
    in_file = InputMultiObject(
        File(exists=True), mandatory=True, desc="the object to be saved"
    )
    meta_dict = traits.DictStrAny(desc="an input dictionary containing metadata")
    source_file = InputMultiObject(
        File(exists=False), mandatory=True, desc="the source file(s) to extract entities from")


class _DerivativesDataSinkOutputSpec(TraitedSpec):
    out_file = OutputMultiObject(File(exists=True, desc="written file path"))
    out_meta = OutputMultiObject(File(exists=True, desc="written JSON sidecar path"))
    compression = OutputMultiObject(
        traits.Either(None, traits.Bool),
        desc="whether ``in_file`` should be compressed (True), uncompressed (False) "
        "or left unmodified (None).",
    )
    fixed_hdr = traits.List(traits.Bool, desc="whether derivative header was fixed")
