*fMRIStroke*: A preprocessing pipeline for fMRI Data from Stroke patients 
=========================================================================
*fMRIStoke* is a BIDs application that runs on the outputs of *fmriprep*
(`www.fmriprep.org <https://www.fmriprep.org>`__) for the preprocessing of
task-based and resting-state functional MRI (fMRI) from stroke patients.

.. image:: https://readthedocs.org/projects/fmristroke/badge/?version=latest
    :target: https://fmristroke.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status
.. image:: https://github.com/alixlam/fmristroke/workflows/Pytest/badge.svg
    :target: https://github.com/alixlam/fmristroke/actions/Pytest
    :alt: github actions pytest


About
-----

.. figure:: https://github.com/alixlam/fmristroke/blob/main/docs/_static/pipeline.png?raw=true
    :target: https://github.com/alixlam/fmristroke/blob/main/docs/_static/pipeline.png

Motivation
~~~~~~~~~~
Stroke not only leads to structural damages of gray/white matter in affected patietns, 
but can also produce remote changes in structurally normal brain areas by a variety of different mechanisms [Siegel2017]_.
As a result, it is highly recommended, notably by [Siegel2017]_, to add specific quality checks and strategies to mitigate lesion specific confounds when dealing with
stroke data especially when doing :abbr:`FC analysis (functional connectivity)`.


To do this we propose *fMRIStroke* a functional magnetic resonance imaging (fMRI) data
quality checks and preprocessing pipeline tailored to stroke data. It is designed to provide an easily accessible,
interface that is robust to variations in scan acquisition
protocols and that requires minimal user input, while providing easily
interpretable and comprehensive reports.
It uses fmriprep (`www.fmriprep.org <https://www.fmriprep.org>`_) outputs derivatives to generate
new quality checks plots for stroke patients when lesion masks are available and
computes new confounds like signals in lesion masks, and ICA based confounds (as proposed in [Yourganov2017]_).


Added quality checks: 
~~~~~~~~~~~~~~~~~~~~~
-  **hemodynamics lagmap** using the Rapidtide_ python tool providing output reports on the hemodynamic lags present bold series.
-  **homotopic connectivity** if freesurfer reconstruction was run.
-  **Registration plots with lesion mask**

Added confounds:
~~~~~~~~~~~~~~~

- **lesion**: signal in lesion mask.
- **CSF lesion**: signal in CSF + lesion combined mask.
- **ICA_comp**: ICA based confounds [Yourganov2017]_.


Added outputs:
~~~~~~~~~~~~~~

- **ROI masks** in standardized space.
- **Denoised fMRI**: Denoised BOLD series using the provided pipelines.
- **Functional Connectivity**: Connectivity matrix using provided atlases and connectivity measures.


The *fMRIStroke* pipeline uses a combination of tools from well-known software
packages, including ANTs_,  FreeSurfer_, Rapidtide_ and Nilearn_

.. important::
  This pipeline was designed to run after fmriprep. Any other fMRI preprocessing tools might not provide the required derivatives for fMRIStroke to run properly. 


In summary this tool allows you to easily do the following:

- Generate preprocessing quality reports specific to stroke patients, with which the user can easily
  identify outliers.
- Receive verbose output concerning the stage of preprocessing for each
  subject, including meaningful errors.
- Automate and parallelize processing steps, which provides a significant
  speed-up from manual processing or shell-scripted pipelines.

More information and documentation can be found at
https://fmristroke.readthedocs.io/


Citation
--------
**Citation**.




Acknowledgements
----------------
This work makes great use of the work by the `NiPreps Community <https://www.nipreps.org>`__.
and the work done by `rapidtides authors <https://rapidtide.readthedocs.io/en/latest/>`__. 


References
----------

.. [Siegel2017] J. S. Siegel, G. L. Shulman, and M. Corbetta, Measuring functional connectivity in stroke: Approaches and considerations, J Cereb Blood Flow Metab, 2017.
     doi: `10.1177/0271678X17709198. <https://doi.org/10.1177/0271678X17709198>`_

.. [Yourganov2017] Yourganov, G., Fridriksson, J., Stark, B., Rorden, C., Removal of artifacts from resting-state fMRI data in stroke. Neuroimage Clin 2017.
     doi: `10.1016/j.nicl.2017.10.027 <https://doi.org/10.1016/j.nicl.2017.10.027>`_
