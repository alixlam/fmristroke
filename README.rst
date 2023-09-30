*fMRIStroke*: A Preprocessing Pipeline for fMRI Data from Stroke patients 
=========================================================
*fMRIStoke* is an application that runs on the outputs of *fmriprep*
(`www.fmriprep.org <https://www.fmriprep.org>`__) for the preprocessing of
task-based and resting-state functional MRI (fMRI) from stroke patients.

.. image:: https://img.shields.io/badge/docker-nipreps/fmriprep-brightgreen.svg?logo=docker&style=flat
  :target: https://hub.docker.com/r/nipreps/fmriprep/tags/
  :alt: Docker image available!

.. image:: https://readthedocs.org/projects/fmriprep/badge/?version=latest
  :target: http://fmriprep.readthedocs.io/en/latest/?badge=latest
  :alt: Documentation Status


About
-----
.. image:: https://github.com/oesteban/fmriprep/raw/f4c7a9804be26c912b24ef4dccba54bdd72fa1fd/docs/_static/fmriprep-21.0.0.svg


*fMRIPrep* is a functional magnetic resonance imaging (fMRI) data
preprocessing pipeline that is designed to provide an easily accessible,
state-of-the-art interface that is robust to variations in scan acquisition
protocols and that requires minimal user input, while providing easily
interpretable and comprehensive error and output reporting.
It performs basic processing steps (coregistration, normalization, unwarping,
noise component extraction, segmentation, skull-stripping, etc.) providing
outputs that can be easily submitted to a variety of group level analyses,
including task-based or resting-state fMRI, graph theory measures, and surface
or volume-based statistics.

.. note::

   *fMRIPrep* performs minimal preprocessing.
   Here we define 'minimal preprocessing'  as motion correction, field
   unwarping, normalization, bias field correction, and brain extraction.
   See the `workflows section of our documentation
   <https://fmriprep.readthedocs.io/en/latest/workflows.html>`__ for more details.

The *fMRIPrep* pipeline uses a combination of tools from well-known software
packages, including FSL_, ANTs_, FreeSurfer_ and AFNI_.
This pipeline was designed to provide the best software implementation for each
state of preprocessing, and will be updated as newer and better neuroimaging
software become available.

This tool allows you to easily do the following:

- Take fMRI data from raw to fully preprocessed form.
- Implement tools from different software packages.
- Achieve optimal data processing quality by using the best tools available.
- Generate preprocessing quality reports, with which the user can easily
  identify outliers.
- Receive verbose output concerning the stage of preprocessing for each
  subject, including meaningful errors.
- Automate and parallelize processing steps, which provides a significant
  speed-up from manual processing or shell-scripted pipelines.

More information and documentation can be found at
https://fmriprep.readthedocs.io/


Citation
--------
**Citation**.
Please acknowledge this work using the citation boilerplate that *fMRIPrep* includes
in the visual report generated for every subject processed.



Acknowledgements
----------------
This work is steered and maintained by the `NiPreps Community <https://www.nipreps.org>`__.
This work was supported by the Laura and John Arnold Foundation,
the NIH (grant NBIB R01EB020740, PI: Ghosh),
and NIMH (R24MH114705, R24MH117179, R01MH121867, PI: Poldrack)
