#!/usr/bin/python

"""
Setup script for Meth5py
"""

try:
	from setuptools import setup, Extension
except:
	from distutils.core import setup, Extension

setup(name = "Meth5py",
	version = "0.4.2",
	author = "Greg Zynda",
	author_email="zyndagj@gmail.com",
	license="BSD",
	url="https://github.com/zyndagj/Meth5py",
	description = "A class for converting BSMAPz methratio output into a fast and compressed hdf5 format",
	long_description_content_type='text/markdown',
	long_description=open('README.md','r').read(),
	packages = ["Meth5py"],
	requires = ["h5py","numpy"],
	test_suite="tests",
	classifiers = ["License :: OSI Approved :: BSD License"])
#tests_require = ['pydoc-markdown'],
