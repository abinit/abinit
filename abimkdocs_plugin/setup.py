#!/usr/bin/env python
"""Setup script for abinit mkdocs plugins."""

from setuptools import setup


setup(
    name="abimkdocs-plugin",
    version="0.1",
    entry_points={
        'mkdocs.plugins': [
            'abimkdocs = abimkdocs_plugin:AbiMkdocsPlugin',
        ]
    }
)
