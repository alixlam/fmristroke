---
title: 'fMRIStroke:  A preprocessing pipeline for fMRI Data from Stroke patients'
tags:
  - Python
  - fMRI
  - preprocessing
  - stroke
  - functional connectivity
authors:
  - name: Alix Lamouroux
    orcid: 
    equal-contrib: true
    affiliation: "1, 2" 

affiliations:
 - name: IMT Atlantique, Lab-STICC UMR CNRS 6285 F-29238, France
   index: 1
 - name: Inria
   index: 2
date: 20 December 2023
bibliography: paper.bib


---

# Summary

Functional Magnetic Resonance Imaging (fMRI) is a widely used neuroimaging technique for the analysis of neural activity and functional connectivity. Studying functional connectivity, which measures the temporal correlations between distinct brain regions, has yielded valuable insights into the intrinsic organization of the human brain and its alterations in response to diseases or injuries. 
Despite its utility, the Blood Oxygen Level Dependent (BOLD) signal obtained from fMRI data is inherently noisy and susceptible to various artifacts, posing challenges to the accuracy of connectivity analyses. Preprocessing and denoising play a crucial role in functional connectivity analysis by identifying and mitigating these nuisances, and reduce their impact on the data. This becomes particularly critical when dealing with stroke patients, given the added complexities associated with their neurological condition. Stroke induces disruptions in the brain's vascular supplies, resulting in structural damages to gray and white matter, along with remote changes in structurally normal brain areas through various mechanisms [@Siegel:2017]. Recognizing these challenges, it is highly recommended, especially by [@Siegel2017], to incorporate specific quality checks and strategies to address lesion-related confounds when working with stroke fMRI data. To address these considerations, we propose fMRIStroke, a tailored functional magnetic resonance imaging (fMRI) data quality checks and preprocessing pipeline designed specifically for stroke data.

# Statement of need
fMRI preprocessing is a crucial step of functional connectivity analysis based on fMRI. Several openly available tools already exist to automate most of the preprocessing steps, notably fMRIprep [@Esteban:2018]. However, these tools are not specific to fMRI data from stroke patients that require additional preprocessing steps and quality checks [@Siegel:2017].

To do this we propose *fMRIStroke* a functional magnetic resonance imaging (fMRI) data quality checks and preprocessing pipeline tailored to stroke data, designed to provide an interface requiring minimal user input, while providing easily interpretable and comprehensive reports with which the user can easily identify outliers.

It uses fmriprep (`www.fmriprep.org <https://www.fmriprep.org>`_) output derivatives to generate new quality checks plots for stroke patients (when lesion masks are available), computes new confounds, and does confounds regression (denoising) to provide a preprocessed fMRI ready for connectivty analysis.

Specifically the following quality checks are added: 
- hemodynamics lagmap
        Stroke can introduce altered blood flow patterns and a normal hemodynamic response cannot be assumed and as been found to be altered after stroke. [@Siegel:2016] investigated temporal delays (lag) in resting state functional magnetic resonance imaging signals and found thatsignificant hemodynamic lag was observed in 30% of stroke patients sub-acutely and approximately 10% of patients showed lag at one-year post-stroke. Lags systematically alter measurements of fonctional connectivity from the affected nodes, and need to be taken into consideration when doing FC analysis [@Siegel:2017]. FMRIStroke measures hemodynamic lag using cross-correlation with the global gray matter signal and includes a visual report including lagmaps using rapidTide an openly available python library [@ref].
- Registration plots with lesion mask.

the following confounds are added: 
- lesion: signal in lesion mask.

- CSF lesion: signal in CSF + lesion combined mask.

- ICA based confounds [@Yourganov2017]:
      Stroke can lead not only focal but also widespread changes in neural activity and it is possible that the presence of the lesion introduces a particular artefact into the fMRI data which is not removed by standard preprocessing techniques. [Yourganov:2018] proposes an ICA, based noise identification method. In the method, independant components are calculated on the bold signal and components that overlap with an ROI that is unlikely to include signal related to neuronal activity, such as Lesion masks are identified as potential noise component. Signals extracted from ICA components can be further regressed out from the fMRI data with a denoising procedure

and the additional outputs are prvided

- ROI masks in standardized space.

- Denoised fMRI: Denoised BOLD series using the provided pipelines.

![fMRIStroke pipeline.](https://github.com/alixlam/fmristroke/blob/main/docs/_static/fmristroke_pipeline.png)


# References