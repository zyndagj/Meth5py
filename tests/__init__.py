#!/usr/bin/python

import unittest, sys, os
try:
	from StringIO import StringIO
except:
	from io import StringIO
# Import path to test the CLI
#try: from unittest.mock import patch
#except: from mock import patch
# buffer for capturing log info
logStream = StringIO()
# Need to start logger BEFORE importing any pyPlateCalibrate code
import logging
FORMAT = "[%(levelname)s - %(filename)s:%(lineno)s - %(funcName)15s] %(message)s"
logging.basicConfig(stream=logStream, level=logging.DEBUG, format=FORMAT)
# Now import module
from Meth5py import Meth5py

class TestMeth5py(unittest.TestCase):
	def setUp(self):
		self.mr = os.path.join(os.path.dirname(__file__), 'test_meth.txt')
		self.h5 = os.path.join(os.path.dirname(__file__), 'test_meth.txt.h5')
		self.fa = os.path.join(os.path.dirname(__file__), 'test.fa')
	def tearDown(self):
		## Runs after every test function ##
		# Wipe log
		logStream.truncate(0)
		# Remove h5 file
		if os.path.exists(self.h5): os.remove(self.h5)
		## Runs after every test function ##
	def testIndexCreation(self):
		m5 = Meth5py(self.mr, self.fa, verbose=True)
	def testFaiReader(self):
		m5 = Meth5py(self.mr, self.fa, verbose=True)
		self.assertEqual(m5.sorted_chroms, ['Chr1','Chr2'])
		self.assertEqual(m5.chrom_dict, {'Chr1':20, 'Chr2':20})
	def testFetch(self):
		m5 = Meth5py(self.mr, self.fa, verbose=True)
		self.assertTrue(all(m5.fetch('Chr1', 10, 10)[0] == [2,0,10,20,1,1]))
		self.assertTrue(all(m5.fetch('Chr1', end=1)[0] == [-1]*6))
		self.assertTrue(all(map(all, m5.fetch('Chr1', 10, 11) == [[2,0,10,20,1,1],[2,0,12,20,1,1]])))
