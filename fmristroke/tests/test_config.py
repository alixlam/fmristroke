"""Check the configuration module and file."""
import os
from pathlib import Path
from unittest.mock import patch

import pytest
from niworkflows.utils.spaces import format_reference
from pkg_resources import resource_filename as pkgrf
from toml import loads

from .. import config


def _reset_config():
    """
    Forcibly reload the configuration module to restore defaults.

    .. caution::
      `importlib.reload` creates new sets of objects, but will not remove
      previous references to those objects."""
    import importlib

    importlib.reload(config)


def test_reset_config():
    execution = config.execution
    setattr(execution, "bids_dir", "TESTING")
    assert config.execution.bids_dir == "TESTING"
    _reset_config()
    assert config.execution.bids_dir is None
    # Even though the config module was reset,
    # previous references to config classes
    # have not been touched.
    assert execution.bids_dir == "TESTING"
