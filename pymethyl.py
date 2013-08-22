#!/usr/bin/python

import sys
import os.path
import struct

class MethIndex:
	def __init__(self, methFile, faFile):
		self.methFile = methFile
		self.methIndex = methFile+'.idx'
		self.FA = faFile
		self.chromDict = {}
		self.readFAI()
		self.makeIndex()
	def readFAI(self):
		#FAI Format  http://www.biostars.org/p/1495/
		#chrName chrLen chrSeek lineBases lineLen
		#Chr1    30427671        6       79      80
		#Chr2    19698289        30812844        79      80
		#Line len is bases+\n
		fai = self.FA+'.fai'
		if not os.path.exists(fai): sys.exit("Please index the fasta file")
		IFAI = open(fai,'r')
		for line in IFAI:
			tmp = line.rstrip('\n').split()
			chrom = tmp[0]
			chromLen = int(tmp[1])
			self.chromDict[chrom] = chromLen
		IFAI.close()
	def makeIndex(self):
		#chr     pos     strand  context ratio   eff_CT_count    C_count CT_count        rev_G_count     rev_GA_count    CI_lower   CI_upper
		#Chr1    15      +       AACCC   1.000   1.00    1       1       8       8       0.207   1.000
		header = self.makeHeader()
		IM = open(self.methFile, 'r')
		OI = open(self.methIndex, 'wb')
		OI.write(self.makeHeader())
		curPos = 1
		line = IM.readline()
		order = []
		if line[:3] != 'chr': #check for header
			IM.seek(0)
		for line in IM:
			tmp = line.split('\t')
			chrom = tmp[0]
			pos = int(tmp[1])
			if pos < curPos: # add chrom order and make sure it prints for rest of chromosome
				curPos = 1
			C = int(tmp[6])
			CorT = int(tmp[7])
			while curPos < pos:
				writeBlank(OI)
				curPos += 1
			writeData(OI, C, CorT)
			curPos += 1
		IM.close()
		OI.close()
	def writeBlank(self,F):
		F.write('NNNN')
	def writeData(self,F,C,CorT):
		F.write(struct.pack('HH',C,CorT))
	def makeHeader(self):
		print sorted(self.chromDict)
		sys.exit()
		header = ""
		for k,v in self.chromDict.iteritems():
			header += k+" "+str(v)+" "
		header = header.rstrip(' ')+'\n'
		return header

def main():
	if len(sys.argv) != 3:
		sys.exit(' '.join((sys.argv[0],'meth.txt','fasta.fa')))
	inFile = sys.argv[1]
	inFA = sys.argv[2]
	MI = MethIndex(inFile, inFA)
	if os.path.exists(inFile+'.idx'):
		print "Index exists"
	else:
		print "Creating Meth Index"
	chromDict = readFAI(inFA)

if __name__ == "__main__":
	main()
