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
Functional connectivity analysis has yielded valuable insights into the intrinsic organization of the human brain and its alterations in response to diseases or injuries by examining the temporal correlations between different brain regions. One widely used technique for measuring this connectivity is Functional Magnetic Resonance Imaging (fMRI).
However, the Blood Oxygen Level Dependent (BOLD) signal obtained from fMRI data is inherently noisy and susceptible to various artifacts, compromising the accuracy and reliability of connectivity analyses. Careful preprocessing and denoising play important roles in functional connectivity analysis by identifying and mitigating these nuisances. This becomes particularly critical when dealing with stroke patients, given the added complexities associated with their neurological condition. Indeed, stroke induces disruptions in the brain's vascular supplies, resulting in structural damages to gray and white matter, along with remote changes in structurally normal brain areas through various mechanisms `[@Siegel:2017]`. Recognizing these challenges, it is highly recommended, especially by `[@Siegel2017]`, to incorporate specific quality checks and strategies to address lesion-related confounds when working with stroke fMRI data. To address these considerations, we propose fMRIStroke, a functional magnetic resonance imaging (fMRI) data quality checks and preprocessing pipeline designed specifically for stroke data.

# Statement of need
fMRI preprocessing is an integral step in functional connectivity analysis based on fMRI. Existing tools like fMRIprep `[@Esteban:2018]` have streamlined most of the preprocessing steps; however, they lack specificity for fMRI data from stroke patients that require additional steps and quality checks `[@Siegel:2017]`.

To this end we introduce fMRIStroke,a functional magnetic resonance imaging (fMRI) data quality checks and preprocessing pipeline tailored specifically for stroke data. It provides an interface requiring minimal user input while delivering easily interpretable and comprehensive reports.

Leveraging on fmriprep `[@Esteban:2018]` output derivatives, fMRIStroke generates new quality check plots for stroke patients, computes additional confounds, and executes confound regression (denoising), providing a preprocessed fMRI ready for connectivity analysis.

Specifically, novel quality checks include a hemodynamics lagmap. fMRIStroke assesses hemodynamic lag using cross-correlation with the global gray matter signal `[@ref]` using rapidTide an openly available python library `[@ref]`. As stroke can introduce altered blood flow patterns a normal hemodynamic response cannot be assumed. `[@Siegel:2016]` investigated temporal delays (lag) in resting state functional magnetic resonance imaging signals and found that significant hemodynamic lag was observed in 30% of stroke patients sub-acutely and approximately 10% of patients showed lag at one-year post-stroke. Lags systematically alter measurements of fonctional connectivity from the affected nodes, and need to be taken into consideration when doing FC analysis `[@Siegel:2017]`.

Additionalyy, signal in the lesion mask and signals in CSF and WM masks incorporating lesion masks are added to the confounds already computed by fMRIprep along with new ICA-based cofounds inspired by `@Yourganov:2018` proposed method. In their paper `@Yourganov:2018` point out that the lesion introduces a particular artefact into the fMRI data which is not removed by standard preprocessing techniques. To mitigate this effect, they propose an ICA-based noise identification method. Independant components are calculated on the bold signal and components that overlap with the lesion (that is unlikely to include signal related to neuronal activity), are identified as potential noise component. Signals extracted from ICA components can be further regressed out from the fMRI data with a denoising procedure.

Finally supplementary outputs include the lesion mask in standardized space and denoised fMRI. 

fMRIStroke pipeline is presented in more details in the following figure. 


![fMRIStroke pipeline.](https://github.com/alixlam/fmristroke/blob/main/docs/_static/fmristroke_pipeline.png)

# References
