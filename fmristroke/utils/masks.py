import os

import nibabel as nib
import numpy as np


def subtract_roi(in_mask, in_roi):
    msk = nib.load(in_mask)
    mskdata = msk.get_fdata().astype(np.uint8)
    roi = np.invert(nib.load(in_roi).get_fdata().astype(np.uint8))
    mskdata = mskdata * roi
    msk.set_data_dtype(np.uint8)
    return msk


def binary_union(mask1, mask2):
    """Generate the union of two masks."""
    img = nib.load(mask1)
    mskarr1 = np.asanyarray(img.dataobj, dtype=int) > 0
    mskarr2 = np.asanyarray(nib.load(mask2).dataobj, dtype=int) > 0
    out = img.__class__(mskarr1 | mskarr2, img.affine, img.header)
    out.set_data_dtype("uint8")
    return out
