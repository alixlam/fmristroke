
Functional data lesion specific data preprocessing

: For each of the 2 BOLD runs found per subject (across all
tasks and sessions), the following lesion specific preprocessing was performed.

    Several additional 'stroke specific' confounding time series were calculated based on the
    *preprocessed BOLD*: Region-wise average signal excluding lesion signal, region-wise average signal in roi and
    region-wise average signal including roi.
    Additionally, a set of lesion related regressors are computed following the methods proposed by [@yourganov]. 
    Independant Components are calculated on the bold signal and components that overlap with an ROI that is unlikely to 
    include signal related to neuronal activity, such as lesion masks are identified as potential noise component. The remaining components are dropped from consideration.
        
    The BOLD time-series were denoised using the following strategies 
:
    *SimpleGS: 24 head motion parameters including: 3 translations, 3 rotations, their temporal derivatives, and their quadratic terms (Satterthwaite et al., 2013), 2Phys - mean physiological signals from white matter (WM) and cerebrospinal fluid (CSF), and high pass filtering by adding discrete cosines transformation basis regressors to handle low-frequency signal drifts. 
*Minimal: 24 head motion parameters including: 3 translations, 3 rotations, their temporal derivatives, and their quadratic terms (Satterthwaite et al., 2013) 
All resamplings can be performed with *a single interpolation
step* by composing all the pertinent transformations (i.e. co-registrations to anatomical and output spaces).
Gridded (volumetric) resamplings were performed using `antsApplyTransforms` (ANTs),
configured with Lanczos interpolation to minimize the smoothing
effects of other kernels [@lanczos].
Concatenate bold runs from same session and task
    
    Connectivity matrices are finally computed on the denoised BOLD signal in standard space using the following atlas : ['Scheaffer'] and the following
    connectivity measures : ['correlation'] using the Nilearn package  [@nilearn].
    Lagged cross-correlation analysis with reference to the global gray matter reference signal was performed for each voxel over 
    the range ±10 TR to generate hemodynamic lagmaps using the rapidtide library [@rapidtide].
        All resamplings can be performed with *a single interpolation
step* by composing all the pertinent transformations (i.e. co-registrations to anatomical and output spaces).
Gridded (volumetric) resamplings were performed using `antsApplyTransforms` (ANTs),
configured with Lanczos interpolation to minimize the smoothing
effects of other kernels [@lanczos].
