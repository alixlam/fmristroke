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

- 1- *hemodynamics lagmap* using the *rapidtide* python tool (`https://rapidtide.readthedocs.io/en/latest/`__) providing
  output reports that are added to the fmriprep report.
- 2- *homotopic connectivity* if freesurfer reconstruction was run.
- 3- *parcellation homogeneity* if atlas is provided.

Added confounds:

- 1- *lesion*: signal in lesion mask.
- 2- *CSF lesion*: signal in CSF + lesion combined mask.
- 3- *ICA_comp*: ICA based confounds [2]_.

Added outputs:

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

Usage Notes
===========
.. warning::
   This tool uses fMRIprep derivatives, it is made to be run after after fMRIprep and will most likely not work on derivatives from other preprocessing tools.


Execution and the BIDS format
-----------------------------
The *fMRIStroke* workflow takes as principal input the path of the dataset
that is to be processed and the path of the **fMRIprep** derivatives.
The input dataset is required to be in valid :abbr:`BIDS (Brain Imaging Data
Structure)` format, and it must include at least one T1w structural image and
a BOLD series.


The exact command to run *fMRIPRep* depends on the Installation_ method.
The common parts of the command follow the `BIDS-Apps
<https://github.com/BIDS-Apps>`_ definition.
Example: ::

    fmristroke data/bids_root/ out/ participant --fmriprep_dir derivatives/ -w work/


Command-Line Arguments
----------------------
    usage: fmristroke [-h]
              [--participant-label PARTICIPANT_LABEL [PARTICIPANT_LABEL ...]]
              [--bids-filter-file FILE] [--fmriprep_dir PATH]
              [--nprocs NPROCS] [--omp-nthreads OMP_NTHREADS]
              [--mem MEMORY_MB] [--low-mem] [--use-plugin FILE]
              [--reports-only] [--output-spaces [OUTPUT_SPACES ...]]
              [--ncomp_method {varexp,aic,kic,mdl}]
              [--ica_method {canica,dictlearning}] [--maxlag MAXLAG]
              [--fs-license-file FILE] [--freesurfer] [-w WORK_DIR]
              [--clean-workdir] [--resource-monitor] [--config-file FILE]
              [--write-graph] [-v] [--stop-on-first-crash]
              bids_dir output_dir {participant}

positional arguments:
~~~~~~~~~~~~~~~~~~~~~
  bids_dir              The root folder of a BIDS valid dataset (sub-XXXXX
                        folders should be found at the top level in this
                        folder).
  output_dir            The output path for the outcomes of preprocessing and
                        visual reports
  {participant}         Processing stage to be run, only "participant" in the
                        case of fMRIPrep (see BIDS-Apps specification).

options:
  -h, --help            show this help message and exit

Options for filtering BIDS queries:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  --participant-label PARTICIPANT_LABEL [PARTICIPANT_LABEL ...], --participant_label PARTICIPANT_LABEL [PARTICIPANT_LABEL ...]
                        A space delimited list of participant identifiers or a
                        single identifier (the sub- prefix can be removed)
                        (default: None)
  --bids-filter-file FILE
                        A JSON file describing custom BIDS input filters using
                        PyBIDS. (default: None)
  --fmriprep_dir PATH   Reuse the anatomical derivatives from another fMRIPrep
                        run or calculated with an alternative processing tool
                        (NOT RECOMMENDED). (default: None)

Options to handle performance:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  --nprocs NPROCS, --nthreads NPROCS, --n_cpus NPROCS, --n-cpus NPROCS
                        Maximum number of threads across all processes
                        (default: None)
  --omp-nthreads OMP_NTHREADS
                        Maximum number of threads per-process (default: None)
  --mem MEMORY_MB, --mem_mb MEMORY_MB, --mem-mb MEMORY_MB
                        Upper bound memory limit for fMRIPrep processes
                        (default: None)
  --low-mem             Attempt to reduce memory usage (will increase disk
                        usage in working directory) (default: False)
  --use-plugin FILE, --nipype-plugin-file FILE
                        Nipype plugin configuration file (default: None)

Options for performing only a subset of the workflow:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  --reports-only        Only generate reports, don't run workflows. This will
                        only rerun report aggregation, not reportlet
                        generation for specific nodes. (default: False)

Workflow configuration:
~~~~~~~~~~~~~~~~~~~~~~~
  --output-spaces [OUTPUT_SPACES ...]
                        Standard and non-standard spaces to resample
                        anatomical and functional images to. Standard spaces
                        may be specified by the form
                        ``<SPACE>[:cohort-<label>][:res-<resolution>][...]``,
                        where ``<SPACE>`` is a keyword designating a spatial
                        reference, and may be followed by optional, colon-
                        separated parameters. Non-standard spaces imply
                        specific orientations and sampling grids. Important to
                        note, the ``res-*`` modifier does not define the
                        resolution used for the spatial normalization. To
                        generate no BOLD outputs, use this option without
                        specifying any spatial references. For further
                        details, please check out
                        https://fmriprep.readthedocs.io/en/latest/spaces.html
                        (default: None)

Options relating to confounds:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  --ncomp_method {varexp,aic,kic,mdl}
                        method to estimate number of components for ICA lesion
                        confounds (default: varexp)
  --ica_method {canica,dictlearning}
                        Method to run ICA lesion (default: canica)

Specific options for hemodynmics analysis:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  --maxlag MAXLAG       Max lag for hemodynamic analysis (default: 10)

Specific options for FreeSurfer preprocessing:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  --fs-license-file FILE
                        Path to FreeSurfer license key file. Get it (for free)
                        by registering at
                        https://surfer.nmr.mgh.harvard.edu/registration.html
                        (default: None)
  --freesurfer          Was freesurfer run

Other options:
~~~~~~~~~~~~~~
  -w WORK_DIR, --work-dir WORK_DIR
                        Path where intermediate results should be stored
                        (default: /homes/a19lamou/fmristroke/work)
  --clean-workdir       Clears working directory of contents. Use of this flag
                        is not recommended when running concurrent processes
                        of fMRIPrep. (default: False)
  --resource-monitor    Enable Nipype's resource monitoring to keep track of
                        memory and CPU usage (default: False)
  --config-file FILE    Use pre-generated configuration file. Values in file
                        will be overridden by command-line arguments.
                        (default: None)
  --write-graph         Write workflow graph. (default: False)
  -v, --verbose         Increases log verbosity for each occurrence, debug
                        level is -vvv (default: 0)
  --stop-on-first-crash
                        Force stopping on first crash, even if a work
                        directory was specified. (default: False)


Troubleshooting
---------------
Logs and crashfiles are outputted into the
``<output dir>/fmriprep/sub-<participant_label>/log`` directory.
Information on how to customize and understand these files can be found on the
`nipype debugging <http://nipype.readthedocs.io/en/latest/users/debug.html>`_
page.



Usage
-----
.. image:: https://github.com/alixlam/fmristroke/blob/main/images/sub-02_ses-S0_task-MIpre_desc-flirtnobbrlesion_bold.svg



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
