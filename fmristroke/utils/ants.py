from pathlib import Path

import ants
import nibabel as nib


def apply_transform_ants(
    fixed, moving, transforms, interpolation, outdir, **args
):
    fixed = ants.image_read(fixed)
    moving = ants.image_read(moving)
    transformed = ants.apply_transforms(
        fixed=fixed,
        moving=moving,
        transformlist=transforms,
        interpolator=interpolation,
        float=0,
        defaultvalue=0,
        **args
    )

    ants.image_write(transformed, outdir)

    return outdir
