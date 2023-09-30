"""The rapidtide module provides basic functions for interfacing with Rapidtide tools."""
import os 

from packaging.version import Version, parse

from nipype import logging
from nipype.interfaces.base import CommandLine, CommandLineInputSpec, traits, isdefined, PackageInfo

LOCAL_DEFAULT_NUMBER_OF_THREADS = 1

class Info(PackageInfo):
    version_cmd = (
        os.path.join(os.getenv("RAPIDTIDEPATH", ""), "rapidtide") + " --version"
    )

    @staticmethod
    def parse_version(raw_info):
        for line in raw_info.splitlines():
            if line.startswith("RapidTide Version: "):
                v_string = line.split()[2]
                break
        else:
            return None

        # -githash may or may not be appended
        v_string = v_string.split("-")[0]

        version = parse(v_string)

        return version.base_version

class RapidTideCommand(CommandLine):
    """Base class for ANTS interfaces"""

    _num_threads = LOCAL_DEFAULT_NUMBER_OF_THREADS

    @staticmethod
    def _format_xarray(val):
        """Convenience method for converting input arrays [1,2,3] to
        commandline format '1x2x3'"""
        return "x".join([str(x) for x in val])

    @property
    def version(self):
        return Info.version()