"""Utilities to handle BIDS derivatives inputs."""
from collections import defaultdict
from pathlib import Path
from json import loads
from pkg_resources import resource_filename as pkgrf
from bids.layout import Query
from bids import BIDSLayout

def collect_bold_derivatives(
    derivatives_dir, subject_id, spec=None, bids_filters=None 
):
    """Gather existing derivatives and compose a cache."""
    if isinstance(derivatives_dir, BIDSLayout):
        layout = derivatives_dir
    else:
        layout = BIDSLayout(str(derivatives_dir), validate=False)
    
    if spec is None:
        spec = loads(Path(pkgrf('fmristroke','data/io_spec.json')).read_text())
    
    bids_filters = bids_filters["bold"] if bids_filters else {}
    
    layout_get_kwargs = {
        'subject': subject_id,
        'session': Query.OPTIONAL
        }
    
    # add filters 
    layout_get_kwargs.update(bids_filters)
    
    derivs_cache = {
        dtype: sorted(layout.get(**layout_get_kwargs, **query))
        for dtype, query in spec.items()
    }

    return derivs_cache

def collect_roi_mask(
    bids_dir, subject_id
):
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
        "return_type": "file"
        }
        
    roi = {
        "roi": layout.get(**query)
    }
    return roi

def group_runs(bold):
    """Group runs from same session and same task"""
    from itertools import groupby
    import re
    
    def _grp_runs(x):
        if "_run-" not in x:
            return x
        run = re.search("_run-\\d*", x).group(0)
        return  x.replace(run, "_run-?")
    
    grouped_runs = []
    for _, bold in groupby(bold, key=_grp_runs):
        bold=list(bold)
        grouped_runs.append(bold)
    return grouped_runs