#!/usr/bin/python

"""
Setup script for pymethyl
"""

from distutils.core import setup, Extension

module1 = Extension('cFetch', sources=["pymethyl/cFetch.c"])

setup(name = "pymethyl",
	version = "0.2",
	author = "Greg Zynda",
	author_email="gjzynda@indiana.edu",
	license="GNU",
	description = "Indexing of MethylCoder output",
	packages = ["pymethyl"],
	ext_modules=[module1])
