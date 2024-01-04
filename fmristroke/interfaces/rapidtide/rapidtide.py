"""RapidTide interface
"""
import os

from nipype.interfaces.base import (
    CommandLineInputSpec,
    File,
    InputMultiObject,
    TraitedSpec,
    isdefined,
    traits,
)
from nipype.utils.filemanip import split_filename

from .base import RapidTideCommand


class RapidTideInputSpec(CommandLineInputSpec):
    in_file = File(
        argstr="%s",
        mandatory=True,
        desc=("image to apply rapidtide to"),
        exists=True,
        position=0,
    )

    output_path = traits.Str(
        "",
        argstr="%s",
        desc="output folder",
        position=1,
        hash_files=False,
        usedefault=True,
    )

    repetitiontime = traits.Float(
        desc="Repetition time (TR) of series - derived from image header if "
        "unspecified"
    )

    corrmask = File(
        argstr="--corrmask %s",
        desc="Only do correlations in nonzero voxels in mask",
        exists=True,
    )

    globalmeaninclude = File(
        argstr="--globalmeaninclude %s",
        desc=" Only use voxels mask for global regressor generation",
        exists=True,
    )

    globalmeanexclude = File(
        argstr="--globalmeanexclude %s",
        desc=" Voxels in mask are excluded for global regressor generation",
        exists=True,
    )

    spatialfilt = traits.Float(
        3,
        argstr="--spatialfilt %f",
        desc="Spatially filter fMRI data prior to analysis using GAUSSSIGMA in mm.",
        usedefault=True,
    )

    confounds_file = File(
        argstr="--motionfile %s",
        desc="Path to confounds file",
        exists=True,
    )

    searchrange = traits.Int(
        10,
        argstr="%d",
        desc="Search range for lags",
        usedefault=True,
    )

    filterband = traits.Str(
        "None",
        argstr="--filterband %s",
        desc="Filter data and regressors to specific band. Use “None” to disable filtering",
        usedefault=True,
    )
    filterfreq = traits.Str(
        "0.009 0.09",
        argstr="--filterfreqs %s",
        desc="Filter data and regressors to retain LOWERPASS to UPPERPASS. If –filterstopfreqs is not also specified, LOWERSTOP and UPPERSTOP will be calculated automatically.",
        usedefault=False,
    )

    noglm = traits.Bool(
        True,
        argstr="--noglm",
        desc="urn off GLM filtering to remove delayed regressor from each voxel",
        usedefault=True,
    )

    motdpos = traits.Bool(
        True,
        argstr="--motpos",
        desc="Toggle whether displacement regressors will be used in motion regression.",
        usedefault=True,
    )


class RapidTideOutputSpec(TraitedSpec):
    output_lagmap = File(exists=True, desc="output lagmap")
    output_corrmap = File(exists=True, desc="max correlation map")
    output_corrfit_mask = File(exists=True, desc="Mask where lag was computed")


class RapidTide(RapidTideCommand):
    """RapidTide,  This is the program that calculates a similarity function between a
    “probe” signal and every voxel of a BOLD fMRI dataset. It then determines the peak value,
    time delay, and wi dth of the similarity function to determine when and how strongly that
    probe signal appears in each voxel.


    """

    _cmd = "rapidtide"
    input_spec = RapidTideInputSpec
    output_spec = RapidTideOutputSpec

    def _format_output_path(self):
        _, name, _ = split_filename(self.inputs.in_file)
        name = name.split("_")[0]
        if self.inputs.output_path == "":
            return "{}/{}".format(os.getcwd(), name)
        else:
            return "{}/{}".format(self.inputs.output_path, name)

    def _gen_filenames(self):
        _, name, ext = split_filename(self.inputs.in_file)
        name = name.split("_")[0]
        output_lagmap = name + "_desc-maxtime_map" + ext
        output_corrmap = name + "_desc-maxcorr_map" + ext
        output_corrfit_mask = name + "_desc-corrfit_mask" + ext
        return output_lagmap, output_corrmap, output_corrfit_mask

    def _list_outputs(self):
        outputs = self._outputs().get()
        output_files = self._gen_filenames()
        outputs["output_lagmap"] = os.path.abspath(output_files[0])
        outputs["output_corrmap"] = os.path.abspath(output_files[1])
        outputs["output_corrfit_mask"] = os.path.abspath(output_files[2])
        return outputs

    def _format_arg(self, name, trait_spec, value):
        if name == "searchrange":
            return "--searchrange -{value} {value}".format(
                value=self.inputs.searchrange
            )
        if name == "output_path":
            return self._format_output_path()
        return super()._format_arg(name, trait_spec, value)
