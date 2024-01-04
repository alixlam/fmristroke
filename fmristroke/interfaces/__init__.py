# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
from json import loads
from pathlib import Path

from niworkflows.interfaces.bids import DerivativesDataSink as _DDSink
from pkg_resources import resource_filename as pkgrf

_pybids_spec = loads(Path(pkgrf("fmristroke", "data/bids.json")).read_text())
BIDS_DERIV_ENTITIES = _pybids_spec["entities"]
BIDS_DERIV_PATTERNS = tuple(_pybids_spec["default_path_patterns"])


class DerivativesDataSink(_DDSink):
    out_path_base = ""
    _config_entities = frozenset({e["name"] for e in BIDS_DERIV_ENTITIES})
    _config_entities_dict = BIDS_DERIV_ENTITIES
    _file_patterns = BIDS_DERIV_PATTERNS


__all__ = ("DerivativesDataSink",)
