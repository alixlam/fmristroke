import io
import os

from setuptools import find_packages, setup


def read(*paths, **kwargs):
    """
    Read the contents of a text file safely.

    >>> read("VERSION")
    '0.1.0'
    >>> read("README.md")
    ...
    """
    content = ""
    with io.open(
        os.path.join(os.path.dirname(__file__), *paths),
        encoding=kwargs.get("encoding", "utf8"),
    ) as open_file:
        content = open_file.read().strip()
    return content


def read_requirements(path):
    return [
        line.strip()
        for line in read(path).split("\n")
        if not line.startswith(('"', "#", "-", "git+"))
    ]


setup(
    name="fmristroke",
    description="Awesome preproceesing tool for stroke functional MRI using Nipype created by Alix L.",
    url="https://github.com/alixlam/fmristroke",
    version="0.1.0",
    author="Alix L.",
    author_email="alix.lamouroux@imt-atlantique.fr",
    packages=find_packages(exclude=["tests"]),
    data_files=[('data', ['data/*'])],
    install_requires=read_requirements("requirements.txt"),
    extras_require={"test": read_requirements("requirements-test.txt")},
    python_requires=">=3.8",
    entry_points = {'console_scripts': [
        'fmristroke=fmristroke.cli.run:main'
    ]}
)