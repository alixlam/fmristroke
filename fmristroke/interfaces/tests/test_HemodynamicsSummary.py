# AUTO-GENERATED by tools/checkspecs.py - DO NOT EDIT
from ..reports import HemodynamicsSummary


def test_HemodynamicsSummary_inputs():
    input_map = dict(
        brain_mask=dict(
            extensions=None,
        ),
        corrfit_mask=dict(
            extensions=None,
        ),
        hemi_mask=dict(
            extensions=None,
        ),
        in_files=dict(),
    )
    inputs = HemodynamicsSummary.input_spec()

    for key, metadata in list(input_map.items()):
        for metakey, value in list(metadata.items()):
            assert getattr(inputs.traits()[key], metakey) == value


def test_HemodynamicsSummary_outputs():
    output_map = dict(
        out_report=dict(
            extensions=None,
        ),
    )
    outputs = HemodynamicsSummary.output_spec()

    for key, metadata in list(output_map.items()):
        for metakey, value in list(metadata.items()):
            assert getattr(outputs.traits()[key], metakey) == value