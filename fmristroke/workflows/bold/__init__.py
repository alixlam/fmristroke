"""

Lesion preprocessing fMRI - BOLD signal workflows
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: fmristroke.workflows.bold.base
.. automodule:: fmristroke.workflows.bold.confounds
.. automodule:: fmristroke.workflows.bold.lagmaps
.. automodule:: fmristroke.workflows.bold.registration
.. automodule:: fmristroke.workflows.bold.denoise


"""

from .base import init_lesion_preproc_wf
from .confounds import init_confs_wf
from .registration import init_lesionplot_wf
from .lagmaps import init_hemodynamic_wf
from .denoise import init_denoise_wf

__all__ = [
    'init_lesion_preproc_wf',
    'init_confs_wf',
    'init_lesionplot_wf',
    'init_hemodynamic_wf',
    'init_denoise_wf',
]