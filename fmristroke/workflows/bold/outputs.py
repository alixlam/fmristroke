"""Writing out derivative files."""
from __future__ import annotations

import typing as ty

import numpy as np
from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe

from fmriprep import config
from fmriprep.config import DEFAULT_MEMORY_MIN_GB
from fmriprep.interfaces import DerivativesDataSink

from niworkflows.engine.workflows import LiterateWorkflow as Workflow
from smriprep.workflows.outputs import _bids_relative


if ty.TYPE_CHECKING:
    from niworkflows.utils.spaces import SpatialReferences

def init_func_lesion_derivatives_wf(
    bids_root: str,
    output_dir: str,
    spaces: SpatialReferences,
    name='func_lesion_derivatives_wf',
):
    """
    Set up a battery of datasinks to store derivatives in the right location.

    Parameters
    ----------
    bids_root : :obj:`str`
        Original BIDS dataset path.
    metadata : :obj:`dict`
        Metadata dictionary associated to the BOLD run.
    output_dir : :obj:`str`
        Where derivatives should be written out to.
    spaces : :py:class:`~niworkflows.utils.spaces.SpatialReferences`
        A container for storing, organizing, and parsing spatial normalizations. Composed of
        :py:class:`~niworkflows.utils.spaces.Reference` objects representing spatial references.
        Each ``Reference`` contains a space, which is a string of either TemplateFlow template IDs
        (e.g., ``MNI152Lin``, ``MNI152NLin6Asym``, ``MNIPediatricAsym``), nonstandard references
        (e.g., ``T1w`` or ``anat``, ``sbref``, ``run``, etc.), or a custom template located in
        the TemplateFlow root directory. Each ``Reference`` may also contain a spec, which is a
        dictionary with template specifications (e.g., a specification of ``{'resolution': 2}``
        would lead to resampling on a 2mm resolution of the space).
    name : :obj:`str`
        This workflow's identifier (default: ``func_lesion_derivatives_wf``).

    """

    workflow = Workflow(name=name)

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                'confounds_file',
                'confounds_metadata',
                'lagmaps',
                'source_file',
                'all_source_files',
            ]
        ),
        name='inputnode',
    )

    raw_sources = pe.Node(niu.Function(function=_bids_relative), name='raw_sources')
    raw_sources.inputs.bids_root = bids_root

    # Confounds
    ds_confounds = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
            desc='confounds2',
            suffix='timeseries',
            dismiss_entities=("echo",),
        ),
        name="ds_confounds",
        run_without_submitting=True,
        mem_gb=DEFAULT_MEMORY_MIN_GB,
    )

    # fmt:off
    workflow.connect([
        (inputnode, raw_sources, [('all_source_files', 'in_files')]),
        (inputnode, ds_confounds, [('source_file', 'source_file'),
                                    ('confounds_file', 'in_file'),
                                    ('confounds_metadata', 'meta_dict')]),
    ])
    # fmt:on

    if getattr(spaces, '_cached') is None:
        return workflow
    
    # Hemodynamics
    lagmap_select = pe.Node(niu.Select(index=0), name="lagmap_select")
    ds_lagmap = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
            desc='hemodyn',
            suffix='bold'
        ),
        name="ds_lagmap",
        run_without_submitting=True,
        mem_gb=DEFAULT_MEMORY_MIN_GB,
    )
    # fmt:off
    workflow.connect([
        (inputnode, lagmap_select, [('lagmaps', 'inlist')]),
        (inputnode, ds_lagmap, [('source_file', 'source_file')]),
        (lagmap_select, ds_lagmap, [('out', 'in_file')]),
    ])
    # fmt:on
    
    return workflow
