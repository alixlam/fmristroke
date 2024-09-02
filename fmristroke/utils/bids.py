"""Utilities to handle BIDS derivatives inputs."""
from collections import defaultdict
from json import loads
from pathlib import Path

from bids import BIDSLayout
from bids.layout import Query
from bids.layout.writing import build_path
from pkg_resources import resource_filename as pkgrf

def collect_anat_derivatives(
    derivatives_dir, subject_id, std_spaces, freesurfer, spec=None, patterns=None
):
    """Gather existing derivatives and compose a cache."""
    if spec is None or patterns is None:
        _spec, _patterns = tuple(
            loads(Path(pkgrf("fmristroke", "data/io_spec_anat.json")).read_text()).values()
        )

        if spec is None:
            spec = _spec
        if patterns is None:
            patterns = _patterns

    derivs_cache = defaultdict(list, {})
    derivatives_dir = Path(derivatives_dir)
    layout = BIDSLayout(derivatives_dir, validate=False)
    
    for space in [None] + std_spaces:
        for k, q in spec["baseline"].items():
            q["subject"] = subject_id
            if space is not None:
                q["space"] = space
            item = layout.get(return_type="filename", **q)
            if space:
                derivs_cache["std_%s" % k] += item if len(item) == 1 else [item]
            else:
                derivs_cache["t1w_%s" % k] = item[0] if len(item) == 1 else item
    for space in std_spaces:
        for k, q in spec["std_xfms"].items():
            q["subject"] = subject_id
            q["from"] = q["from"] or space
            q["to"] = q["to"] or space
            item = layout.get(return_type="filename", **q)
            if not item:
                continue
            derivs_cache[k] += item
    derivs_cache = dict(derivs_cache)  # Back to a standard dictionary

    if freesurfer:
        for k, q in spec["surfaces"].items():
            q["subject"] = subject_id
            item = layout.get(return_type="filename", **q)
            if len(item) == 1:
                item = item[0]
            derivs_cache[k] = item
    derivs_cache["template"] = std_spaces
    return derivs_cache


def collect_bold_derivatives(
    derivatives_dir, subject_id, spec=None, bids_filters=None
):
    """Gather existing derivatives and compose a cache."""
    if isinstance(derivatives_dir, BIDSLayout):
        layout = derivatives_dir
    else:
        layout = BIDSLayout(str(derivatives_dir), validate=False)

    if spec is None:
        spec = loads(
            Path(pkgrf("fmristroke", "data/io_spec.json")).read_text()
        )

    bids_filters = bids_filters["bold"] if bids_filters else {}

    layout_get_kwargs = {"subject": subject_id, "session": Query.OPTIONAL}

    # add filters
    layout_get_kwargs.update(bids_filters)

    derivs_cache = {
        dtype: sorted(
            layout.get(**layout_get_kwargs, **query, invalid_filters="allow")
        )
        for dtype, query in spec.items()
    }

    return derivs_cache


def collect_roi_mask(bids_dir, subject_id):
    """Find roi mask and compose a cache."""
    if isinstance(bids_dir, BIDSLayout):
        layout = bids_dir
    else:
        layout = BIDSLayout(str(bids_dir), validate=False)

    query = {
        "datatype": "anat",
        "suffix": "roi",
        "subject": subject_id,
        "extension": [".nii", ".nii.gz"],
        "return_type": "file",
    }

    roi = {"roi": layout.get(**query)}
    return roi


def group_runs(bold):
    """Group runs from same session and same task"""
    import re
    from itertools import groupby

    def _grp_runs(x):
        if "_run-" not in x:
            return x
        run = re.search("_run-\\d*", x).group(0)
        return x.replace(run, "_run-?")

    grouped_runs = []
    for _, bold in groupby(bold, key=_grp_runs):
        bold = list(bold)
        grouped_runs.append(bold)
    return grouped_runs
