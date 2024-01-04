"""Test parser."""
from argparse import ArgumentError
from contextlib import nullcontext

import pytest

from ... import config
from ...tests.test_config import _reset_config
from ..parser import _build_parser, parse_args

MIN_ARGS = ["data/", "out/", "participant"]


@pytest.mark.parametrize(
    "args,code",
    [
        ([], 2),
        (MIN_ARGS, 2),  # bids_dir does not exist
    ],
)
def test_parser_errors(args, code):
    """Check behavior of the parser."""
    with pytest.raises(SystemExit) as error:
        _build_parser().parse_args(args)

    assert error.value.code == code


@pytest.mark.parametrize("args", [MIN_ARGS, MIN_ARGS])
def test_parser_valid(tmp_path, args):
    """Check valid arguments."""
    datapath = tmp_path / "data"
    datapath.mkdir(exist_ok=True)
    args[0] = str(datapath)

    opts = _build_parser().parse_args(args)

    assert opts.bids_dir == datapath


def test_bids_filter_file(tmp_path, capsys):
    bids_path = tmp_path / "data"
    out_path = tmp_path / "out"
    bff = tmp_path / "filter.json"
    args = [str(bids_path), str(out_path), "participant", "--bids-filter-file", str(bff)]
    bids_path.mkdir()

    parser = _build_parser()

    with pytest.raises(SystemExit):
        parser.parse_args(args)

    err = capsys.readouterr().err
    assert "Path does not exist:" in err

    bff.write_text('{"invalid json": }')

    with pytest.raises(SystemExit):
        parser.parse_args(args)

    err = capsys.readouterr().err
    assert "JSON syntax error in:" in err
    _reset_config()

