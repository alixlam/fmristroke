"""Combine runs utilities"""


def combine_run_source(in_files):
    """
    Create a new source name when
    combining multiple runs
    """
    import os

    from nipype.utils.filemanip import filename_to_list

    base, in_file = os.path.split(filename_to_list(in_files)[0])
    entities = [
        ent for ent in in_file.split("_") if not ent.startswith("run-")
    ]
    basename = "_".join(entities)
    return os.path.join(base, basename)
