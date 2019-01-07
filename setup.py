#!/usr/bin/python

"""
Setup script for Meth5py
"""

try:
	from setuptools import setup, Extension
except:
	from distutils.core import setup, Extension

setup(name = "Meth5py",
	version = "0.3",
	author = "Greg Zynda",
	author_email="zyndagj@gmail.com",
	license="BSD-3",
	description = "A class for converting BSMAPz methratio output into a fast and compressed hd5f format",
	packages = ["Meth5py"],
	requires = ["h5py"],
	tests_require = ['pydoc-markdown'],
	test_suite="tests")
