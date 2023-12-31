# AUTO-GENERATED by tools/checkspecs.py - DO NOT EDIT
from ..reports import HemoPlot


def test_HemoPlot_inputs():
    input_map = dict(
        anat_file=dict(
            extensions=None,
            mandatory=True,
        ),
        label=dict(),
        lagmap_file=dict(
            extensions=None,
            mandatory=True,
        ),
        out_report=dict(
            extensions=None,
            hash_files=False,
            usedefault=True,
        ),
        roi=dict(
            extensions=None,
        ),
        vmax=dict(),
    )
    inputs = HemoPlot.input_spec()

    for key, metadata in list(input_map.items()):
        for metakey, value in list(metadata.items()):
            assert getattr(inputs.traits()[key], metakey) == value


def test_HemoPlot_outputs():
    output_map = dict(
        out_report=dict(
            extensions=None,
        ),
    )
    outputs = HemoPlot.output_spec()

    for key, metadata in list(output_map.items()):
        for metakey, value in list(metadata.items()):
            assert getattr(outputs.traits()[key], metakey) == value
