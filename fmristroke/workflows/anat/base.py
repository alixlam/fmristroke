"""
Orchestrating the ROI-resampling workflow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: init_roi_preproc_wf
.. autofunction:: init_anat_derivatives_wf

"""
import os

import nibabel as nb
import numpy as np
from nipype.interfaces import utility as niu
from nipype.interfaces.fsl import Split as FSLSplit
from nipype.pipeline import engine as pe
from niworkflows.utils.connections import listify, pop_file

from ... import config
from ...interfaces import DerivativesDataSink

# BOLD workflows
from .outputs import init_anat_lesion_derivatives_wf
from .resample import (
    init_roi_std_trans_wf,
)

def init_roi_preproc_wf(name="roi_std_wf"):
    """
    This workflow controls the lesion specific stages *fMRIStroke*.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from fmriprep.workflows.tests import mock_config
            from fmriprep import config
            from fmriprep.workflows.bold.base import init_func_preproc_wf
            with mock_config():
                bold_file = config.execution.bids_dir / "sub-01" / "func" \
                    / "sub-01_task-mixedgamblestask_run-01_bold.nii.gz"
                wf = init_func_preproc_wf(str(bold_file))

    Inputs
    ------
    template
        List of templates to target
    anat2std_xfm
        List of transform files, collated with templates
    std2anat_xfm
        List of inverse transform files, collated with templates
    roi
        ROI mask

    Outputs
    -------
    roi_mask_t1
        ROI mask in T1w space
    roi_mask_std
        roi mask, resampled to template space


    """
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow
    
    mem_gb = {"filesize": 1, "resampled": 1, "largemem": 1}

    # Have some options handy
    omp_nthreads = config.nipype.omp_nthreads
    spaces = config.workflow.spaces
    output_dir = str(config.execution.output_dir)

    # Extract BIDS entities
    layout = config.execution.layout

    # Build workflow
    workflow = Workflow(name=name)
    workflow.__postdesc__ = """\
All resamplings can be performed with *a single interpolation
step* by composing all the pertinent transformations (i.e. co-registrations to anatomical and output spaces).
Gridded (volumetric) resamplings were performed using `antsApplyTransforms` (ANTs),
configured with Lanczos interpolation to minimize the smoothing
effects of other kernels [@lanczos].
"""

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "t1w_preproc",
                "anat2std_xfm",
                "std2anat_xfm",
                "template",
                "roi",
            ]
        ),
        name="inputnode",
    )

    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "roi_mask_t1",
                "roi_mask_std",
            ]
        ),
        name="outputnode",
    )
    
    anat_lesion_derivatives_wf = init_anat_lesion_derivatives_wf(
        bids_root=layout.root,
        output_dir=output_dir,
        spaces=spaces,
    )
    workflow.connect([
        (inputnode, anat_lesion_derivatives_wf, [
            ("t1w_preproc", "inputnode.all_source_files"),
        ]),
    ])
    
    # ROI RESAMPLING #######################################################

    if spaces.get_spaces(nonstandard=False, dim=(3,)):
        # Apply transforms in 1 shot
        roi_std_trans_wf = init_roi_std_trans_wf(
            mem_gb=mem_gb["resampled"],
            omp_nthreads=omp_nthreads,
            spaces=spaces,
            name="lesion_std_trans_wf",
            use_compression=not config.execution.low_mem,
        )
        # fmt:off
        workflow.connect([
            (inputnode, roi_std_trans_wf, [
                ("template", "inputnode.templates"),
                ("anat2std_xfm", "inputnode.anat2std_xfm"),
                ("roi", "inputnode.roi")
            ]),
            (roi_std_trans_wf, outputnode, [
                ("outputnode.roi_mask_std", "roi_mask_std"),
            ]),
        ])
        workflow.connect([
            (roi_std_trans_wf, anat_lesion_derivatives_wf, [
                ("outputnode.template", "inputnode.template"),
                ("outputnode.spatial_reference", "inputnode.spatial_reference"),
                ("outputnode.roi_mask_std", "inputnode.roi_mask_std"),
            ]),
        ])
        # fmt:on

    return workflow