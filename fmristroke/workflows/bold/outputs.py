"""Writing out derivative files."""
from __future__ import annotations

import typing as ty

import numpy as np
from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe

from fmriprep import config
from fmriprep.config import DEFAULT_MEMORY_MIN_GB
from ...interfaces import DerivativesDataSink

from niworkflows.engine.workflows import LiterateWorkflow as Workflow
from smriprep.workflows.outputs import _bids_relative


if ty.TYPE_CHECKING:
    from niworkflows.utils.spaces import SpatialReferences
    from ...utils.pipelines import Pipelines

def init_func_lesion_derivatives_wf(
    bids_root: str,
    output_dir: str,
    spaces: SpatialReferences,
    pipelines: Pipelines,
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
    from niworkflows.interfaces.utility import KeySelect
    from niworkflows.interfaces.space import SpaceDataSource


    workflow = Workflow(name=name)

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                'confounds_file',
                'confounds_metadata',
                'lagmaps',
                'source_file',
                'all_source_files',
                'spatial_reference',
                'template',
                'denoised_bold_t1',
                'denoised_bold_std',
                'pipeline',
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
            desc='maxtime',
            suffix='map'
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
    
    # Denoising
    pipelinessource = pe.Node(
        niu.IdentityInterface(
            fields=["pipeline"]),
        name  = "pipelineinter"
        )
    pipelinessource.iterables = [("pipeline", pipelines.get_pipelines())]
    
    fields = ['denoised_bold_t1']
    select_pipeline = pe.Node(
        KeySelect(fields=fields),
        name='select_pipeline',
        run_without_submitting=True,
        mem_gb=DEFAULT_MEMORY_MIN_GB,
    )
    ds_denoised_t1 = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
            desc="denoised",
            suffix="bold",
            space="T1w",
            compress=True,
            dismiss_entities=("echo",),
        ),
        name="ds_denoised_t1",
        run_without_submitting=True,
        mem_gb=DEFAULT_MEMORY_MIN_GB
    )
    
    spacesourceiter = pe.Node(SpaceDataSource(), name='spacesourceiter', run_without_submitting=True)
    spacesourceiter.iterables= (
        'in_tuple',
        [(s.fullname, s.spec) for s in spaces.cached.get_standard(dim=(3,))],
    )
    
    fields_space = ['template', 'denoised_bold_std']
    fields = ['denoised_bold_std']
    select_pipeline_std = pe.Node(
        KeySelect(fields=fields),
        name='select_pipeline_std',
        run_without_submitting=True,
        mem_gb=DEFAULT_MEMORY_MIN_GB,
    )
    select_std = pe.Node(
            KeySelect(fields=fields_space),
            name='select_std',
            run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB,
    )
    
    ds_denoised_std = pe.Node(
        DerivativesDataSink(
                base_directory=output_dir,
                desc='denoised',
                suffix='bold',
                compress=True,
                dismiss_entities=("echo",),
            ),
            name='ds_denoised_std',
            run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB,
    )
    # fmt:off
    workflow.connect([
        # Denoised bold in T1w space
        (inputnode, ds_denoised_t1, [('source_file', 'source_file')]),
        (inputnode, select_pipeline, [('denoised_bold_t1', 'denoised_bold_t1'),
                                        ('pipeline', 'keys')]),
        (pipelinessource, select_pipeline, [('pipeline', 'key')]),
        (select_pipeline, ds_denoised_t1, [('denoised_bold_t1', 'in_file')]),
        (pipelinessource, ds_denoised_t1, [('pipeline', 'pipeline')]),
        
        # Denoised bold in std space
        (inputnode, ds_denoised_std, [('source_file', 'source_file')]),
        (inputnode, select_pipeline_std, [('denoised_bold_std', 'denoised_bold_std'),
                                            ('pipeline', 'keys')]),
        (pipelinessource, select_pipeline_std, [('pipeline', 'key')]),
        (inputnode, select_std, [('template', 'template'),
                                ('spatial_reference', 'keys')]),
        (spacesourceiter, select_std, [('uid', 'key')]),
        (select_pipeline_std, select_std, [('denoised_bold_std', 'denoised_bold_std')]),
        (select_std, ds_denoised_std, [('denoised_bold_std', 'in_file')]),
        (spacesourceiter, ds_denoised_std, [('space', 'space'),
                                        ('cohort', 'cohort'),
                                        ('resolution', 'resolution'),
                                        ('density', 'density')]),
        (pipelinessource, ds_denoised_std, [('pipeline', 'pipeline')])
    ])
    # fmt:on
    
    
    return workflow
