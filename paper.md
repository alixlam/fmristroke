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
Functional Magnetic Resonance Imaging (fMRI) is a commonly used neuroimaging technique to investigate neural activity. In particular, functional connectivity analysis, which examines the temporal correlation between different brain regions, has yielded important insights into the intrinsic organization of the human brain and its alterations due to diseases or injuries.
However, the Blood Oxygen Signal (BOLD) measured with fMRI data is inherently noisy and contains various artifacts that can compromise the accuracy and reliability of connectivity analyses.
Preprocessing and denoising fMRI data is particularly crucial for stroke patients due to several unique challenges associated with their neurological condition. Stroke, characterized by a sudden disruption of blood supply to the brain, not only leads to structural damages of gray/white matter in affected patients, but can also produce remote changes in structurally normal brain areas by a variety of different mechanisms [Siegel:2017]. Analyzing fMRI data from stroke patients poses additional complexities [Siegel:2017], making careful preprocessing even more essential. 


# Statement of need
TO DO : 
1 - Many tools exists but not for stroke
Many tools to preprocess fmri data exist notably fMRIprep, but not specific to stroke patients  

2 - Present fMRI stroke functionality 
To do this we propose *fMRIStroke* a functional magnetic resonance imaging (fMRI) data
quality checks and preprocessing pipeline tailored to stroke data. It is designed to provide an easily accessible, interface that requires minimal user input, while providing easily interpretable and comprehensive reports.
It generates preprocessing quality reports specific to stroke patients, with which the user can easily identify outliers.

Specifically, it uses fmriprep (`www.fmriprep.org <https://www.fmriprep.org>`_) output derivatives to generate new quality checks plots for stroke patients (when lesion masks are available), computes new confounds like signals in lesion masks, and ICA based confounds (as proposed in [Yourganov:2017]), and does confounds regression (denoising) to provide a preprocessed fMRI ready for connectiivty analysis.



Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements


# References