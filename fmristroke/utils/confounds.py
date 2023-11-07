from nilearn.interfaces.fmriprep.load_confounds_utils import (
    MissingConfound,
    _add_suffix,
    _check_params,
    _find_confounds,
)
import pandas as pd

def _load_iclesion(confounds_raw):
    """Load the IC-Lesion regressors.

    Parameters
    ----------
    confounds_raw : pandas.DataFrame
        DataFrame of confounds.

    Returns
    -------
    pandas.DataFrame
        DataFrame of IC-Lesion regressors.

    Raises
    ------
    MissingConfound
        When IC-lesion regressors are not found, raise error as
    
    """
    ic_lesion_params = _find_confounds(confounds_raw, ["ica_lesion"])
    if not ic_lesion_params:
        raise MissingConfound(keywords=["ica_lesion"])
    return confounds_raw[ic_lesion_params]

def _load_wm_csf_lesion(confounds_raw, wm_csf):
    """Load the regressors derived from the white matter and CSF masks taking into account lesion mask.

    Parameters
    ----------
    confounds_raw : pandas.DataFrame
        DataFrame of confounds.

    wm_csf : str
        White matter and CSF strategy to use. Options are "basic",
        "derivatives", "power2", or "full".

    Returns
    -------
    pandas.DataFrame
        DataFrame of white matter and CSF regressors.

    Raises
    ------
    MissingConfound
        When white matter no lesion and CSF lesion regressors are not found, raise error as
        wm_csf is not a valid choice of strategy.
    """
    wm_csf_params = _add_suffix(["csf_lesion", "white_matter_nolesion"], wm_csf)
    if _check_params(confounds_raw, wm_csf_params):
        return confounds_raw[wm_csf_params]
    else:
        raise MissingConfound(keywords=["wm_csf"])
