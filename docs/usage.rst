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

    fmristroke data/bids_root/ out/ participant -w work/ --fmriprep-dir derivatives/


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

