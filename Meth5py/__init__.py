#!/usr/bin/python

import subprocess as sp
import h5py, logging, sys, time, os
from array import array
import multiprocessing as mp
import numpy as np
try:
	from itertools import ifilter
except:
	ifilter = filter
global logger
logger = logging.getLogger('Meth5py')

class Meth5py:
	'''
	Class for converting BSMAPz methratio files to hdf5 for fast random access

	# Parameters
	methFile (str): Methylation input file
	faFile (str): Reference file
	h5File (str): previously created h5 file
	n_cores (int): The number of cores to use for indexing [0 all]
	force (bool): Force the re-creation of h5 file
	verbose (bool): Enable verbose logging

	# Attributes
	logger (log): Logger instance
	contexts (tuple): Tuple of methylation contexts
	strands (tuple): Tuple of strands
	sorted_chroms (list): Sorted chromosomes in reference
	chom_dict (dict): Dictionary of chromosomes and their lengths
	n_cores (int): The number of cores to use for indexing

	# Example
	```
	from Meth5py import Meth5py
	m5 = Meth5py('tests/test_meth.txt', 'tests/test.fa')
	for record in m5.fetch('Chr1',10,11):
		print(record)
	m5.close()
	```
	'''
	def __init__(self, methFile='', faFile='', h5File='', n_cores=0, force=False, verbose=False):
		'''
		Meth5py initializer
		'''
		# Init logger
		FORMAT = '[%(levelname)s - PID-%(process)d - %(name)s.%(funcName)s - %(msecs)d] %(message)s'
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
		self.n_cores = n_cores
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
		Reads the methratio file and creates an hdf5 file of the following structure

		file[chrom_dataset][position] = [context_index, strand_index, c, ct, g, ga]

		# Attributes
		self.methFile (str): reads Methylation input file
		self.methBin (str): writes the structured binary output
		self.sorted_chroms (dict): reads the chomosomes to initializing h5 file
		self.chrom_dict (dict): reads chromosome lengths for initializing h5 file
		'''
		if not self.n_cores or self.n_cores > 1:
			# Open file handles
			H5 = h5py.File(self.methBin, 'w')
			#IM = open(self.methFile, 'r')
			# Create datasets
			for chrom in self.sorted_chroms:
				chrom_len = self.chrom_dict[chrom]
				# [context, strand, c, ct, g, ga]
				H5.create_dataset(chrom, (chrom_len, 6), compression='lzf', chunks=True, fillvalue=-1, dtype='i')
			H5.close()
			NP = min(self.n_cores, mp.cpu_count()) if self.n_cores else mp.cpu_count()
			logger.debug("Creating index using %i cores"%(NP))
			# Create global variables
			global syncArray, chromArray, currentChrom
			syncArray = mp.RawArray('B', [0]*NP)
			chromArray = mp.RawArray('i', max(self.chrom_dict.values())*6)
			np_chromArray = np.frombuffer(chromArray, dtype='i')
			np_chromArray[:] = -1
			currentChrom = mp.RawArray('c', 30)
			currentChrom.value = get_first_chrom(self.methFile).encode('ascii')
			# Launch workers
			pool = mp.Pool(NP)
			pool.map(index_worker, [(self.methFile, self.methBin, self.chrom_dict, i, NP) for i in range(NP)])
			pool.close()
			pool.join()
			# Get an error if the shared mem is not deleted
			del syncArray, chromArray, currentChrom
		else:
			logger.debug("Creating index using 1 core")
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
			index_buffer = []
			value_buffer = []
			start_time, end_time = time.time(), time.time()
			current_chrom = False
			buffer_size = 1000
			for i, line in enumerate(ifilter(lambda x: x[0] != "#", IM)):
				tmp = line.split('\t')	
				chrom, pos, strandIndex = (tmp[0], int(tmp[1])-1, self.strands.index(tmp[2]))
				if not current_chrom: current_chrom = chrom
				cIndex = self.contexts.index(tmp[3])
				c, ct, g, ga = map(int, tmp[6:10])
				# Write if chrom changes
				if current_chrom != chrom:
					H5[current_chrom][index_buffer] = value_buffer
					end_time = time.time()
					elapsed_time = end_time - start_time
					#print("Processed %i records per second"%(len(index_buffer)/elapsed_time))
					#print("Finished %s"%(current_chrom))
					index_buffer, value_buffer = ([], [])
					current_chrom = chrom
					index_buffer.append(pos)
					value_buffer.append((cIndex, strandIndex, c, ct, g, ga))
					start_time = time.time()
				elif (i+1)%buffer_size == 0:
					index_buffer.append(pos)
					value_buffer.append((cIndex, strandIndex, c, ct, g, ga))
					H5[current_chrom][index_buffer] = np.array(value_buffer)
					#end_time = time.time()
					#elapsed_time = end_time - start_time
					#print("Processed %i records per second"%(len(index_buffer)/elapsed_time))
					index_buffer, value_buffer = ([], [])
					start_time = time.time()
				else:
					index_buffer.append(pos)
					value_buffer.append((cIndex, strandIndex, c, ct, g, ga))
			if index_buffer:
				H5[current_chrom][index_buffer] = value_buffer
				index_buffer, value_buffer = ([], [])
			IM.close()
			H5.close()

def get_first_chrom(meth_file):
	with open(meth_file,'r') as IM:
		line = IM.readline()
		line = IM.readline()
	chrom = line.split('\t')[0]
	return chrom

def index_worker(args):
	global currentChrom
	contexts = ('CG','CHG','CHH')
	strands = ('+', '-')
	meth_file, h5_file, chrom_dict, pid, NP = args
	if pid == 0:
		H5 = h5py.File(h5_file,'a')
		np_chromArray = np.frombuffer(chromArray, dtype='i')
	IM = open(meth_file,'r')
	short_sleep = 0.5
	long_sleep = 1
	assert(IM.readline().split('\t')[0] == 'chr')
	for i,line in ifilter(lambda x: x[0]%NP==pid, enumerate(IM)):
		tmp = line.rstrip('\n').split('\t')	
		chrom, pos, strandIndex = (tmp[0], int(tmp[1])-1, strands.index(tmp[2]))
		cIndex = contexts.index(tmp[3])
		c, ct, g, ga = map(int, tmp[6:10])
		# Handle write outs
		if currentChrom.value.decode('ascii') != chrom:
			syncArray[pid] = 0
			if pid != 0:
				#logger.debug("Waiting for chromosome %s to be written"%(currentChrom.value))
				while syncArray[pid] == 0:
					time.sleep(short_sleep)
					if syncArray[pid] == 1 and currentChrom.value.decode('ascii') != chrom:
						#logger.debug("Waiting for next chromosome")
						syncArray[pid] = 0
			else: # pid == 0
				while sum(syncArray) > 0:
					#logger.debug("MASTER waiting for others to finish writing "+str(list(syncArray)))
					time.sleep(short_sleep)
				clen = chrom_dict[currentChrom.value.decode('ascii')]
				assert(H5[currentChrom.value.decode('ascii')].shape == (clen,6))
				#H5[currentChrom.value][:,:] = np.reshape(chromArray[:clen*6], (clen,6))
				H5[currentChrom.value.decode('ascii')][:,:] = np.reshape(np_chromArray[:clen*6], (clen,6))
				np_chromArray[:] = -1
				assert(np.sum(chromArray) == -len(chromArray))
				#chromArray[:] = [-1]*len(chromArray)
				#logger.debug("Wrote %s. Updating global chrom to %s and removing barrier"%(currentChrom.value, chrom))
				currentChrom.value = chrom.encode('ascii')
				syncArray[:] = [1]*len(syncArray)
		chromArray[pos*6:(pos+1)*6] = (cIndex, strandIndex, c, ct, g, ga)
		#print("PID-%i %.2f wrote %s:%i-%i cn=%s"%(pid,time.time.time(), chrom, pos+1, pos+6, currentChrom.value))
	syncArray[pid] = 0
	if pid == 0:
		while sum(syncArray) > 0:
			#logger.debug("MASTER waiting for others to finish writing "+str(list(syncArray)))
			time.sleep(short_sleep)
		clen = chrom_dict[currentChrom.value.decode('ascii')]
		assert(H5[currentChrom.value.decode('ascii')].shape == (clen,6))
		#logger.debug("Wrote %s. The last chromosome"%(currentChrom.value))
		H5[currentChrom.value.decode('ascii')][:,:] = np.reshape(chromArray[:clen*6], (clen,6))
		H5.close()
	IM.close()
