# AUTO-GENERATED by tools/checkspecs.py - DO NOT EDIT
from ..confounds import FMRISummary


def test_FMRISummary_inputs():
    input_map = dict(
        confounds_file=dict(
            extensions=None,
        ),
        confounds_list=dict(),
        drop_trs=dict(
            usedefault=True,
        ),
        in_nifti=dict(
            extensions=None,
            mandatory=True,
        ),
        in_segm=dict(
            extensions=None,
        ),
        str_or_tuple=dict(),
        tr=dict(
            usedefault=True,
        ),
    )
    inputs = FMRISummary.input_spec()

    for key, metadata in list(input_map.items()):
        for metakey, value in list(metadata.items()):
            assert getattr(inputs.traits()[key], metakey) == value


def test_FMRISummary_outputs():
    output_map = dict(
        out_file=dict(
            extensions=None,
        ),
    )
    outputs = FMRISummary.output_spec()

    for key, metadata in list(output_map.items()):
        for metakey, value in list(metadata.items()):
            assert getattr(outputs.traits()[key], metakey) == value
