# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# Kinetic Diagram Analysis Examples
#

import sys
from setuptools import setup, find_packages
import versioneer

needs_pytest = {"pytest", "test", "ptr"}.intersection(sys.argv)
pytest_runner = ["pytest-runner"] if needs_pytest else []

with open("README.md", "r") as handle:
    long_description = handle.read()

setup(
    # Self-descriptive entries which should always be present
    name="kda-examples",
    author="Nikolaus Awtrey",
    author_email="nawtrey@asu.edu",
    description="Kinetic Diagram Analysis Examples",
    long_description=long_description,
    long_description_content_type="text/markdown",
    keywords="science kinetics chemistry physics",
    packages=find_packages(),
    install_requires=["numpy", "networkx", "kda",],
    tests_require=["pytest"],
    include_package_data=True,
    setup_requires=[] + pytest_runner,
)
