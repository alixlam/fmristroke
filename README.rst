*fMRIStroke*: A Preprocessing Pipeline for fMRI Data from Stroke patients 
=========================================================
*fMRIStoke* is an application that runs on the outputs of *fmriprep*
(`www.fmriprep.org <https://www.fmriprep.org>`__) for the preprocessing of
task-based and resting-state functional MRI (fMRI) from stroke patients.

.. image:: https://readthedocs.org/projects/fmriprep/badge/?version=latest
  :target: http://fmriprep.readthedocs.io/en/latest/?badge=latest
  :alt: Documentation Status


About
-----

.. image:: https://github.com/alixlam/fmristroke/blob/main/images/pipeline.png 

*fMRIStroke* is a functional magnetic resonance imaging (fMRI) data
quality checks pipeline that is designed to provide an easily accessible,
interface that is robust to variations in scan acquisition
protocols and that requires minimal user input, while providing easily
interpretable and comprehensive error and output reporting.
It uses fmriprep (`www.fmriprep.org <https://www.fmriprep.org>`__) outputs derivatives to generate
new quality checks plots for stroke patients when lesion masks are available (as recommended in [1]_) and
computes new confounds like signals in lesion masks, and ICA based confounds (as proposed in [2]_).

Added quality checks: 
~~~~~~~~~~~~~~~~~~~~~
- 1- *hemodynamics lagmap* using the *rapidtide* python tool (`https://rapidtide.readthedocs.io/en/latest/`__) providing
  output reports that are added to the fmriprep report.
- 2- *homotopic connectivity* if freesurfer reconstruction was run.
- 3- *parcellation homogeneity* if atlas is provided.

Added confounds:
~~~~~~~~~~~~~~~

- 1- *lesion*: signal in lesion mask.
- 2- *CSF lesion*: signal in CSF + lesion combined mask.
- 3- *ICA_comp*: ICA based confounds [2]_.

Added outputs:
~~~~~~~~~~~~~~

- 1- ROI masks in standardized space.

The *fMRIStroke* pipeline uses a combination of tools from well-known software
packages, including ANTs_ and FreeSurfer_.
This pipeline was designed to run after fmriprep.

This tool allows you to easily do the following:

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

.. [1] To  add 

.. [2] To add
