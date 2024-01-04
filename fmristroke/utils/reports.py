import nibabel as nib
import numpy as np


def get_hemoreport(lagmap, corrmap, corrfit_mask, brain_mask, hemimask=None):
    results = {}
    mask = nib.load(corrfit_mask).get_fdata()

    corrmap = nib.load(corrmap).get_fdata()
    lagmap = nib.load(lagmap).get_fdata()

    results["mean_lag"] = round(np.mean(np.abs(lagmap[mask > 0])), 2)
    results["mean_maxcorr"] = round(np.mean(corrmap[mask > 0]), 2)

    nvoxels = np.sum(nib.load(brain_mask).get_fdata())
    results["prop_voxel"] = round(1 - (nvoxels / np.sum(mask)))
    results["num_voxels"] = nvoxels

    if not hemimask:
        results["lh_mean_lag"] = None
        results["rh_mean_lag"] = None
    else:
        from copy import copy

        from nilearn.image import resample_to_img

        hemimask = resample_to_img(
            source_img=hemimask, target_img=brain_mask, interpolation="nearest"
        )
        # Labels from Freesurfer : 2 : Left cerebral wm, 3: Left cereabelar
        # cortex
        lh_mask = copy(hemimask).get_fdata()
        lh_mask[np.logical_and(lh_mask != 2, lh_mask != 3)] = 0
        lh_mask[np.logical_or(lh_mask == 2, lh_mask == 3)] = 1
        rh_mask = copy(hemimask).get_fdata()
        rh_mask[np.logical_and(rh_mask != 42, rh_mask != 43)] = 0
        rh_mask[np.logical_or(rh_mask == 42, rh_mask == 43)] = 1
        results["lh_mean_lag"] = np.round(
            np.mean(np.abs(lagmap[mask + lh_mask == 2])), 2
        )
        results["rh_mean_lag"] = np.round(
            np.mean(np.abs(lagmap[mask + rh_mask == 2])), 2
        )
    return results
