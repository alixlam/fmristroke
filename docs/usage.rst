.. _Usage :

Usage Notes
===========
.. warning::
   This app runs on **fMRIprep derivatives**, fMRIprep needs to be run prior to running **fMRIStroke**. Derivatives obtained with other tools will most likely not work on this tool. 


Execution and the BIDS format
-----------------------------
The *fMRIStoke* workflow takes as principal input the path of the dataset
that is to be processed and the path to the fMRIprep derivatives. 
The input dataset must include at least one T1w structural image, an ROI mask of the stroke lesion and a BOLD series.

The exact command to run *fMRIStroke* depends on the Installation_ method.
The common parts of the command follow the `BIDS-Apps
<https://github.com/BIDS-Apps>`_ definition.
Example: ::

    fmristroke data/bids_root/ out/ participant -w work/ --fmriprep_dir derivatives/


ROI mask in BIDS format
~~~~~~~~~~~~~~~~~~~~~~~~
Lesion masks should be binary NIfTI images (damaged areas = 1, everywhere else = 0) in the same space and resolution as the T1w image, and use the ``_roi`` suffix, for example, ``sub-01_label-lesion_roi.nii.gz``. 
This file should be located in the ``sub-01/anat`` directory of the BIDS dataset. Since lesion masks are not yet part of the BIDS standard, it is also necessary to add a ``.bidsignore`` file in the root of the BIDS dataset to prevent the BIDS validator from raising errors.
Your ``.bidsignore`` file should contain the following line: ::

   *_roi.nii.gz

.. note::
   As mentioned by the `fMRIprep documentation <https://fmriprep.org/en/stable/workflows.html#cost-function-masking-during-spatial-normalization>`_ (v23.x.x):
   The lesion masking instructions in this section predate the release of BIDS Derivatives.
   As of BIDS 1.4.0, the recommended naming convention is::

       manual_masks/
       └─ sub-001/
          └─ anat/
             ├─ sub-001_desc-tumor_mask.nii.gz
             └─ sub-001_desc-tumor_mask.json

   In an upcoming version, fmriprep will search for lesion masks as pre-computed
   derivatives. Until this is supported, we will continue to look for the ``_roi`` suffix.



Command-Line Arguments
----------------------
.. argparse::
   :ref: fmristroke.cli.parser._build_parser
   :prog: fmristroke
   :nodefault:
   :nodefaultconst:


Troubleshooting
---------------
Logs and crashfiles are outputted into the
``<output dir>/sub-<participant_label>/log`` directory.
Information on how to customize and understand these files can be found on the
`nipype debugging <http://nipype.readthedocs.io/en/latest/users/debug.html>`_
page.

**Support and communication**.
The documentation of this project is found here: http://fmristroke.readthedocs.org/en/latest/.

All bugs, concerns and enhancement requests for this software can be submitted here:
https://github.com/alixlam/fmristroke/issues.

