

#def init_denoise_wf(
#    mem_gb: float,
#    pipelines: Pipeline,
#    maxlag : int = 10, 
#    name: str = "hemodynamic_wf",
#):
#    """
#    Build a workflow to generate lagmaps.
#
#    This workflow calculates the hemodynamic lag for stroke patients.
#    as recommended by Siegel et al. (2017) ``https://pubmed.ncbi.nlm.nih.gov/28541130/``
#    It uses the existing tool RapidTide that calculates a similarity function between a 
#    “probe” signal and every voxel of a BOLD fMRI dataset. It then determines the peak value, 
#    time delay, and width of the similarity function to determine when and how strongly that 
#    probe signal appears in each voxel.
#    For details about the method visit : ``https://rapidtide.readthedocs.io/en/latest/``
#    
#
#
#    Workflow Graph
#        .. workflow::
#            :graph2use: orig
#            :simple_form: yes
#
#            from fmristroke.workflows.lagmaps import init_hemodynamic_wf
#            wf = init_hemodynamic_wf(
#                mem_gb=1,
#                metadata={},
#            )
#
#    Parameters
#    ----------
#    mem_gb : :obj:`float`
#        Size of BOLD file in GB - please note that this size
#        should be calculated after resamplings that may extend
#        the FoV
#    metadata : :obj:`dict`
#        BIDS metadata for BOLD file
#    maxlag :obj: `int`
#        Max lag to compute.
#    name : :obj:`str`
#        Name of workflow (default: ``bold_confs_wf``)
#
#    Inputs
#    ------
#    bold_t1
#        BOLD image in T1w space, after the prescribed corrections (STC, HMC and SDC)
#        when available.
#    boldmask
#        Bold fov mask
#    roi
#        Roi mask in T1w space 
#    t1w
#        T1w image
#    t1w_mask_lesion
#        Mask of the lesion in T1w space
#    t1w_mask
#        Brain Mask of T1
#    t1w_tpms
#        List of tissue probability maps in T1w space
#    t1w_aseg
#        Segmentation of structural image, done with FreeSurfer.
#    confounds_file
#        TSV of all aggregated confounds.
#
#    Outputs
#    -------
#    maps
#        lagmap and corrmap
#    """
#    return