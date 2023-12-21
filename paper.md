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
Functional connectivity analysis has yielded valuable insights into the intrinsic organization of the human brain and its alterations in response to diseases or injuries. One widely used technique for measuring this connectivity is Functional Magnetic Resonance Imaging (fMRI). However, fMRI is inherently noisy and susceptible to various artifacts, compromising the accuracy and reliability of connectivity analyses. This becomes particularly critical when dealing with stroke patients, given the added complexities associated with their neurological condition. Specific preprocessing and denoising are integral steps to identify the nuisance sources and mitigate their effect on functional connectivity analysis. 
To adress these consideration, we introduce fMRIStroke, a functional magnetic resonance imaging (fMRI) data quality checks and preprocessing pipeline tailored specifically for stroke data. fMRIStroke works as an add-on to typical preprocessing pipeline. Leveraging on a commonly used fMRI preprocessing tools outputs (fMRIprep), fMRIStroke generates new quality check plots, computes additional confounds, and executes confound regression (denoising), providing preprocessed fMRI ready for connectivity analysis. 



# Statement of need
Stroke, with its complex and varied impact on the brain's vascular supply, often results in damages that extend beyond the immediate lesion site, making functional connectivity analysis a promising tool to understand its immediate impact as well as network reorganizations during recovery. 

However, while the Blood Oxygen Level Dependent (BOLD) signal obtained from fMRI data is typically mixed with non-neuronal sources of variability [@ref] in healthy subjects it is even more true for stroke patients `[@Siegel:2017]`, compromising the reliability of connectivity analyses. 
Recognizing these challenges, it is highly recommended, especially by `@Siegel2017`, to incorporate specific quality checks and strategies to address lesion-related confounds when working with stroke fMRI data. 

Existing tools like fMRIprep `[@Esteban:2018]` have streamlined most of the preprocessing steps; however, they lack specificity for fMRI data from stroke patients. 

To this end, we propose fMRIStroke, a fMRI preprocessing pipeline designed specifically for stroke data that builds on the output derivatives of fMRIprep `[@Esteban:2018]`. It provides an interface requiring minimal user input while delivering easily interpretable and comprehensive reports. 

Specifically, novel quality checks include a hemodynamic lagmap. fMRIStroke assesses hemodynamic lag using cross-correlation with the global gray matter signal `[@ref]` using rapidTide an openly available python library `[@ref]`. As stroke can introduce altered blood flow patterns a normal hemodynamic response cannot be assumed. `[@Siegel:2016]` investigated temporal delays (lag) in resting state functional magnetic resonance imaging signals and found that significant hemodynamic lag was observed in 30% of stroke patients sub-acutely and approximately 10% of patients showed lag at one-year post-stroke. Lags systematically alter measurements of fonctional connectivity from the affected nodes, and need to be taken into consideration when doing FC analysis `[@Siegel:2017]`.

Additionalyy, signal in the lesion mask and signals in CSF and WM masks incorporating lesion masks are added to the confounds already computed by fMRIprep along with new ICA-based cofounds inspired by `@Yourganov:2018` proposed method. In their paper `@Yourganov:2018` point out that the lesion introduces a particular artefact into the fMRI data which is not removed by standard preprocessing techniques. To mitigate this effect, they propose an ICA-based noise identification method. Independant components are calculated on the bold signal and components that overlap with the lesion (that is unlikely to include signal related to neuronal activity), are identified as potential noise component. Signals extracted from ICA components can be further regressed out from the fMRI data with a denoising procedure.

Finally supplementary outputs include the lesion mask in standardized space and multiple denoised fMRI as there is currently no consensus in the fMRI community on an optimal denoising strategy that perform best.

fMRIStroke pipeline is presented in more details in the following figure. 


![fMRIStroke pipeline.](https://github.com/alixlam/fmristroke/blob/main/docs/_static/fmristroke_pipeline.png)

# References
