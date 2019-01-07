#!/usr/bin/python

import subprocess as sp
import os.path
import h5py, logging, sys

class Meth5py:
	'''
	Class for converting BSMAPz methratio files to hdf5 for fast random access

	# Parameters
	methFile (str): Methylation input file
	faFile (str): Reference file
	h5File (str): previously created h5 file
	force (bool): Force the re-creation of h5 file
	verbose (bool): Enable verbose logging

	# Attributes
	logger (log): Logger instance
	contexts (tuple): Tuple of methylation contexts
	strands (tuple): Tuple of strands
	sorted_chroms (list): Sorted chromosomes in reference
	chom_dict (dict): Dictionary of chromosomes and their lengths

	# Example
	```
	from Meth5py import Meth5py
	m5 = Meth5py('tests/test_meth.txt', 'tests/test.fa')
	for record in m5.fetch('Chr1',10,11):
		print(record)
	m5.close()
	```
	'''
	def __init__(self, methFile='', faFile='', h5File='', force=False, verbose=False):
		'''
		Meth5py initializer
		'''
		# Init logger                                                                                       
		FORMAT = '[%(levelname)s - %(name)s.%(funcName)s] %(message)s'
		if verbose:
			self.log_level = 'DEBUG'
			logging.basicConfig(level=logging.DEBUG, format=FORMAT)
		else:
			self.log_level = 'INFO'
			logging.basicConfig(level=logging.INFO, format=FORMAT)
		self.logger = logging.getLogger('Meth5py')
		# Config static variables
		self.contexts = ('CG','CHG','CHH')
		self.strands = ('+', '-')
		if h5File and self._check_file(h5File):
			self.logger.debug("Loading previously created Meth5py file: %s"%(h5File))
			self.methBin = h5File
			self.H5 = h5py.File(self.methBin, 'r')
			self.sorted_chroms = sorted(self.H5.keys())
			self.logger.debug("Detected the following chromosomes:\n[%s]"%(', '.join(self.sorted_chroms)))
			self.chrom_dict = {chrom:self.H5[chrom].size[0] for chrom in self.sorted_chroms}
		else:
			# Check input files
			self.methFile = self._check_file(methFile)
			self.methBin = methFile + '.h5'
			# Check fasta and and fai
			self.fasta = self._check_file(faFile)
			try:
				self.fasta_fai = self._check_file(faFile+'.fai')
			except FileNotFoundError:
				try:
					self.logger.debug("Attempting to index with samtools")
					sp.check_call(['samtools','faidx',self.fasta])
					self.fasta_fai = self._check_file(faFile+'.fai')
				except sp.CalledProcessError:
					self.logger.error("Samtools was not present to index %s"%(self.fasta))
					raise FileNotFoundError
			# Get chromosomes and lengths
			self._readFAI(self.fasta_fai)
			# Create index if it doesn't exist
			if not os.path.exists(self.methBin) or force: self._makeIndex()
			# Open the index file
			self.H5 = h5py.File(self.methBin, 'r')
	def _check_file(self, file_name):
		'''
		Checks to see if a file exists

		# Parameters
		file_name (str): file path to check

		# Returns
		str: file name if it exists

		# Raises
		FileNotFoundError: If file doesn't exist
		'''
		if os.path.exists(file_name):
			return file_name
		self.logger.error("%s does not exist")
		raise FileNotFoundError
	def _readFAI(self, fai):
		"""
		Reads the sizes of each chromosome from the FASTA index

		# Parameters
		fai (str): path to reference fai

		# Attribues
		self.sorted_chroms (list): Sorted list of chromosome names
		self.chrom_dict (dict): Dictionary of chromosome lengths
		"""
		#FAI Format  http://www.biostars.org/p/1495/
                #chrName chrLen chrSeek lineBases lineLen
                #Chr1    30427671        6       79      80
                #Line len is bases+\n
		with open(fai, 'r') as FAI:
			lines = [line.rstrip('\n').split() for line in FAI]
		self.sorted_chroms = sorted([line[0] for line in lines])
		self.chrom_dict = {line[0]:int(line[1]) for line in lines}
	#def view(self, chrom, start = 1, end = -1, minCov = 5, maxCov = 200):
	#	pass
	def fetch(self, chrom, start = 1, end = -1, minCov = 0, maxCov = -1, index=True):
		'''
		Fetches a region of methylation ratios.
		
		If no start is given, 1-end is returned.
		If no end is given, start- is returned.

		-1 values are returned IF:

		 * depth < minCov
		 * depth > maxCov
		 * no reads for that position
		 * no methylation for that position
		
		# Parameters
		chrom (str): Chromosome name
		start (int): Start of region (1-indexed)
		end (int): End of region (1-indexed)
		minCov (int): Minimum coverage needed (Default: 0)
		maxCov (int): Maximum coverage allowed (Default: -1)
		index (bool): Retrun context and strand indices instead of strings

		# Returns
		list: [[context, strand, c, ct, g, ga], ...]
		'''
		if chrom not in self.sorted_chroms:
			self.logger.error("Chromosome %s not in %s"%(chrom, self.methBin))
			return []
		if end != -1 and end < start:
			self.logger.error("End coord %i is smaller than start %i"%(end, start))
			return []
		# samtools faidx test.fa Chr1:0-5 and Chr1:1-5 return the same thing
		startIndex = max(0, start-1)
		endIndex = self.chrom_dict[chrom] if end == -1 else end
		# get region
		ret = self.H5[chrom][startIndex:endIndex]
		if index:
			return ret
		else:
			return [(self.contexts[cI], self.strands[sI], c, ct, g, ga) for cI, sI, c, ct, g, ga in ret]
	def close(self):
		'''
		Closes the H5 file for cleaning instances of Meth5py

		# Attributes
		self.H5 (file): The file that is closed
		'''
		if self.H5:
			self.H5.close()
			del self.H5
	def _makeIndex(self):
		'''
		Reads the methratio file and creates an hd5f file of the following structure

		file[chrom_dataset][position] = [context_index, strand_index, c, ct, g, ga]

		# Attributes
		self.methFile (str): reads Methylation input file
		self.methBin (str): writes the structured binary output
		self.sorted_chroms (dict): reads the chomosomes to initializing h5 file
		self.chrom_dict (dict): reads chromosome lengths for initializing h5 file
		'''
		# Open file handles
		H5 = h5py.File(self.methBin, 'w')
		IM = open(self.methFile, 'r')
		# Create datasets
		for chrom in self.sorted_chroms:
			chrom_len = self.chrom_dict[chrom]
			# [context, strand, c, ct, g, ga]
			H5.create_dataset(chrom, (chrom_len, 6), compression='gzip', compression_opts=6, chunks=True, fillvalue=-1, dtype='i')
		# Skip the first line
		firstLine = IM.readline()
		#chr     pos     strand  context ratio   eff_CT_count    C_count CT_count        rev_G_count     rev_GA_count    CI_lower        CI_upper
		#test1   14      -       CHG     0.000   1.00    0       1       0       0       0.000   0.793
		for line in filter(lambda x: x[0] != "#", IM):
			tmp = line.split('\t')	
			chrom, pos, strandIndex = (tmp[0], int(tmp[1])-1, self.strands.index(tmp[2]))
			cIndex = self.contexts.index(tmp[3])
			c, ct, g, ga = map(int, tmp[6:10])
			# Set record
			H5[chrom][pos] = (cIndex, strandIndex, c, ct, g, ga)
		IM.close()
		H5.close()
