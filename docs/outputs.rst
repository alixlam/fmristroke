.. _outputs:

------------------------
Outputs of *fMRIStroke*
------------------------
*fMRIStroke* outputs conform to the :abbr:`BIDS (brain imaging data structure)`
Derivatives specification (see `BIDS Derivatives`, along with the
upcoming `BEP 011` and `BEP 012`).
*fMRIStroke* generates three broad classes of outcomes:

1. **Visual QA (quality assessment) reports**:
   one :abbr:`HTML (hypertext markup language)` per subject, that allows the user a visual assessment of the quality of preprocessing.
   If the output dir is the fmriprep dir then the fmriprep report will be updated with fmriStroke figures.

2. **Derivatives (preprocessed data)** the input fMRI data ready for
   analysis, i.e., after the various preparation procedures
   have been applied.

3. **Confounds**: Confounds signals that can be utilized
   to run subsequent denoising steps.


Layout
------
Assuming fMRIStroke is invoked with::

    fmristroke <input_dir>/ <output_dir>/ participant --fmriprep-dir <derivatives_dir> [OPTIONS]

The outputs will be a `BIDS Derivatives`_ dataset of the form::

    <output_dir>/
      logs/
      sub-<label>/
      sub-<label>.html
      dataset_description.json
      .bidsignore

For each participant in the dataset,
a directory of derivatives (``sub-<label>/``)
and a visual report (``sub-<label>.html``) are generated.
``dataset_description.json`` is the metadata file in which fMRIPrep
records metadata recommended by the BIDS standard.


Visual Reports
--------------
*fMRIStroke* outputs summary reports, written to ``<output dir>/sub-<subject_label>.html``.
These reports make visual inspection of the results easy.
`View a sample report. <_static/SampleReport/sample_report.html>`_

   .. note::
       If the <output_dir> is the same as <fmirprep_dir> the fmriprep report will be updated with some visual report specific to lesions.
       Otherwise a report containing only lesion specific quality checks will be generated.

Coregistration
~~~~~~~~~~~~~~
The ROI mask is added to the coregistration plot.

.. figure:: _static/sub-027_ses-1week_task-rest_run-02_space-T1w_desc-reglesion_bold.svg


Hemodynamics
~~~~~~~~~~~~
Stroke disrupts the brain's vascular supply, not only within but also outside areas of infarction.
[Siegel2016]_ investigated temporal delays (lag) in resting state functional magnetic resonance imaging signals and found that significant hemodynamic lag was observed in 30% of stroke patients sub-acutely and
approximately 10% of patients showed lag at one-year post-stroke. Lags systematically alter measurements of fonctional connectivity from the affected nodes.
Thus [Siegel2017]_ recommends excluding patients with an hemodynamic lag > 1s.
Eiher way, lags should be kept in mind when doing :abbr:`FC analysis (functional connectivity)` of stroke patients.

Hemodynamic lag is determined here using cross-correlation with the global gray matter signal.
Visual report includes the lagmaps along with the max correlation map (Correlation corresponding to the lag in the lagmap).

.. figure:: _static/sub-024_task-rest_.svg


Derivatives of *fMRIStroke* (preprocessed data)
---------------------------------------------
Preprocessed, or derivative, data are written to
``<output dir>/sub-<subject_label>/``.
The `BIDS Derivatives`_ specification describes the naming and metadata conventions we follow.

Anatomical derivatives
~~~~~~~~~~~~~~~~~~~~~~
Anatomical derivatives are placed in each subject's ``anat`` subfolder::

  sub-<subject_label>/
    anat/
      sub-<subject_label>[_space-<space_label>]_label-lesion_roi-.nii.gz

Spatially-standardized derivatives are denoted with a space label,
such as ``MNI152NLin2009cAsym``, while derivatives in
the original ``T1w`` space omit the ``space-`` keyword.



Functional derivatives
~~~~~~~~~~~~~~~~~~~~~~
Functional derivatives are stored in the ``func/`` subfolder.
All derivatives contain ``task-<task_label>`` (mandatory) and ``run-<run_index>`` (optional), and
these will be indicated with ``[specifiers]``::

  sub-<subject_label>/
    func/
      sub-<subject_label>_[specifiers]_space-<space_label>_desc-lagmap.nii.gz


**Regularly gridded outputs (images)**.
Volumetric output spaces labels (``<space_label>`` above, and in the following) include
``T1w`` and ``MNI152NLin2009cAsym`` (default).


**Extracted confounding time series**.
For each :abbr:`BOLD (blood-oxygen level dependent)` run processed with *fMRIStroke*, an
accompanying *confounds* file will be generated. Thi confound file, contains both **fmriprep confounds** and additional lesion specific confounds.
Confounds_ are saved as a :abbr:`TSV (tab-separated value)` file::

  sub-<subject_label>/
    func/
      sub-<subject_label>_[specifiers]_desc-confounds_timeseries.tsv
      sub-<subject_label>_[specifiers]_desc-confounds_timeseries.json

These :abbr:`TSV (tab-separated values)` tables look like the example below,
where each row of the file corresponds to one time point found in the
corresponding :abbr:`BOLD (blood-oxygen level dependent)` time series::

  csf white_matter  global_signal std_dvars dvars framewise_displacement  t_comp_cor_00 t_comp_cor_01 t_comp_cor_02 t_comp_cor_03 t_comp_cor_04 t_comp_cor_05 a_comp_cor_00 a_comp_cor_01 a_comp_cor_02 a_comp_cor_03 a_comp_cor_04 a_comp_cor_05 non_steady_state_outlier00  trans_x trans_y trans_z rot_x rot_y rot_z
  682.75275 0.0 491.64752000000004  n/a n/a n/a 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 -0.00017029 -0.0  0.0
  669.14166 0.0 489.4421  1.168398  17.575331 0.07211929999999998 -0.4506846719 0.1191909139  -0.0945884724 0.1542023065  -0.2302324641 0.0838194238  -0.032426848599999995 0.4284323184  -0.5809158299 0.1382414008  -0.1203486637 0.3783661265  0.0 0.0 0.0207752 0.0463124 -0.000270924  -0.0  0.0
  665.3969  0.0 488.03  1.085204  16.323903999999995  0.0348966 0.010819676200000001  0.0651895837  -0.09556632150000001  -0.033148835  -0.4768871111 0.20641088559999998 0.2818768463  0.4303863764  0.41323714850000004 -0.2115232212 -0.0037154909000000004  0.10636180070000001 0.0 0.0 0.0 0.0457372 0.0 -0.0  0.0
  662.82715 0.0 487.37302 1.01591 15.281561 0.0333937 0.3328022893  -0.2220965269 -0.0912891436 0.2326688125  0.279138129 -0.111878887  0.16901660629999998 0.0550480212  0.1798747037  -0.25383302620000003  0.1646403629  0.3953613889  0.0 0.010164  -0.0103568  0.0424513 0.0 -0.0  0.00019174


Confounds
---------
The :abbr:`BOLD (blood-oxygen level dependent)` signal measured with fMRI is a mixture of fluctuations
of both neuronal and non-neuronal origin.
Neuronal signals are measured indirectly as changes in the local concentration of oxygenated hemoglobin.
Non-neuronal fluctuations in fMRI data may appear as a result of motion, scanner noise, physiological fluctuations (related to cardiac or respiratory effects) but also lesion specific artefacts. 

*Confounds* (or nuisance regressors) are variables representing these fluctuations with a potential
non-neuronal origin.
These non-neuronal fluctuations may drive spurious results in fMRI data analysis,
especially in functional connectivity analyses.
To minimize those confounding effects we can include
them as nuisance regressors and regress them out from
the fMRI data, also known as *denoising*.
There is currently no consensus on an optimal denoising strategy in the fMRI community.
Different strategies have been proposed, which achieve different compromises between
how much of the non-neuronal fluctuations are effectively removed, and how much of neuronal fluctuations
are damaged in the process.
The *fMRIPrep* pipeline generates a large array of possible confounds and the *fMRIStoke* pipeline adds to these confounds some lesion specific ones refer to [Yourganov2017]_ for more details.


.. danger::
   Similarly to fmriprep computed confounds: 
   Do not include all columns of ``~_desc-confounds_timeseries.tsv`` table
   into your design matrix or denoising procedure.
   Filter the table first, to include only the confounds (or components thereof)
   you want to remove from your fMRI signal.
   The choice of confounding variables depends on the analysis you want to perform and the data you have,
   it may be not straightforward as no gold standard procedure exists. **fMRIStroke** includes a denoising step in its pipeline
   using simple strategies proposed in the litterature. However, it is up to you to either use these denoised BOLD series, add a custom strategy or
   denoise using the ``~_desc-confounds_timeseries.tsv`` output independently of **fMRIStroke**. 


Confound regressors description
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Refer to `fmrirep doc <https://fmriprep.org/en/stable/outputs.html#confounds>`_ for more details about confounds and confounds regression.

**ICLesion confounds**.
:abbr:`ICLesion analysis (Independant Components Analysis based Lesion Noise Correction)` is a :abbr:`ICA (Independant component analysis)`,
based noise identification method.
In the method, independant components are calculated on the bold signal and components that overlap with an :abbr:`ROI (Region of Interest)`
that is unlikely to include signal related to neuronal activity, such as :abbr:`Lesion` masks are identified as potential noise component.
Signals extracted from ICA components can be further regressed out from the fMRI data with a
denoising procedure [Yourganov2017]_.

- ``IC_lesion_XX`` - additional noise components are calculated using :abbr:`ICA
  (ICA noise correction))`;

Each confounds data file will also have a corresponding metadata file
(``~desc-confounds.json``).
Metadata files contain additional information about columns in the confounds TSV file:

.. code-block:: json

    {
      "ica_lesion_06": {
        "Method": "canICA",
        "Retained": true,
        "jaccard": 0.06,
    }

For ICA decompositions, entries include:

  - ``Method``:  ICA method used.
  - ``Retained``: Indicates whether the component was saved in ``desc-confounds_timeseries.tsv``
    for use in denoising.
  - ``Jaccard``: Overlapping between spatial map of component and ROI mask. 


Confounds on the visual reports
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The visual reports provide several sections per task and run to aid designing
a denoising strategy for subsequent analysis.

Noise components computed during ICA decomposition are evaluated according
to the overlap of their spatial map with the ROI mask.
This is used by *fMRIStroke* to determine whether each component should be saved for
use in denoising operations: a component is saved if the jaccard index between ROI and binarized spatial map is > 5%.
*fMRIStroke* reports include a plot of the spatial map of each included component along with associated signal.

.. figure:: _static/sub-027_ses-1week_task-rest_run-02_space-T1w_desc-icaroi_bold.svg

Denoising
---------
As mentioned above there is no concensus on denoising strategy. However, **fmristroke** proposes some simple denoising pipelines [Yourganov2017]_ to preprocess your data.
By default **fmristroke** performs denoising using 2 common, and 2 lesion specific denoising pipelines. Each pipeline is defined as a single .json file.

Denoised BOLD series will be saved as::

  sub-<subject_label>/
    func/
      sub-<subject_label>_[specifiers]_space-<space_label>_pipeline-[Pipeline]_desc-denoised_bold.nii.gz


SimpleGS
~~~~~~~~~~
Denoising strategy based on regressing out: 24HMP - 24 head motion parameters including: 3 translations, 3 rotations, their temporal derivatives, and their quadratic terms,
8Phys - mean physiological signals from white matter (WM) and cerebrospinal fluid (CSF), their temporal derivatives, and quadratic terms,
and high pass filtering by adding discrete cosines transformation basis regressors to handle low-frequency signal drifts.

CompCorGS
~~~~~~~~~~
Denoising strategy based on regressing out: 24 head motion parameters including: 3 translations, 3 rotations, their temporal derivatives, and their quadratic terms, CompCor - Signals from CompCor,
and high pass filtering by adding discrete cosines transformation basis regressors to handle low-frequency signal drifts.",

SimpleLesionGS
~~~~~~~~~~~~~~~
Same as SimpleGS but with updating region signals (WM and CSF) with lesion mask.

ICLesionGS
~~~~~~~~~~
Same as SimpleLesionGS but adding IC_Lesion signals. [Yourganov2017]_

Adding custom strategy
~~~~~~~~~~~~~~~~~~~~~~
You can easily add a custom pipeline by creating a .json file. A file should follow the structure below.

.. code-block:: json

    {
    "name": "Name",
    "description": "Denoising strategy based on ...",
    "confounds": {
        "white_matter": {
            "raw": "False",
            "lesion": "False",
            "derivative1": "False",
            "power2": "False",
            "derivative1_power2": "False"
        },
        "csf": {
            "raw": "False",
            "lesion": "False",
            "derivative1": "False",
            "power2": "False",
            "derivative1_power2":  "False"
            },
        "global_signal": {
            "raw": "False",
            "derivative1": "False",
            "power2": "False",
            "derivative1_power2": "False"
            },
        "motion": {
            "raw": "False",
            "derivative1": "False",
            "power2": "False",
            "derivative1_power2": "False"
            },
        "high_pass": "False",
        "compcor": "False",
        "ic_lesion": "False"
        
    },
    "aroma": "False",
    "spikes": "False",
    "demean": "False"
    }


See implementation on :mod:`~fmristroke.workflows.bold.confounds.init_lesion_confs_wf`.



.. topic:: References

  .. [Yourganov2017] Yourganov, G., Fridriksson, J., Stark, B., Rorden, C., Removal of artifacts from resting-state fMRI data in stroke. Neuroimage Clin 2017.
     doi: `10.1016/j.nicl.2017.10.027 <https://doi.org/10.1016/j.nicl.2017.10.027>`_

  .. [Siegel2016] J. S. Siegel, A. Z. Snyder, L. Ramsey, G. L. Shulman, and M. Corbetta, The effects of hemodynamic lag on functional connectivity and behavior after stroke, J Cereb Blood Flow Metab 2016.
     doi: `10.1177/0271678X15614846. <http://journals.sagepub.com/doi/10.1177/0271678X15614846>`_

  .. [Siegel2017] J. S. Siegel, G. L. Shulman, and M. Corbetta, Measuring functional connectivity in stroke: Approaches and considerations, J Cereb Blood Flow Metab, 2017.
     doi: `10.1177/0271678X17709198. <https://doi.org/10.1177/0271678X17709198>`_

