"""
Nipype interfaces to metrics for evaluating denoising pipelines 
"""
import os

import numpy as np
from nipype.interfaces.base import (
    BaseInterfaceInputSpec,
    InputMultiPath,
    OutputMultiPath,
    SimpleInterface,
    TraitedSpec,
    isdefined,
    traits,
)
from nipype.utils.filemanip import fname_presuffix, split_filename

