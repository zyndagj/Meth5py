#!/usr/bin/python

import sys
import os.path
import struct

class MethIndex:
	def __init__(self, methFile, faFile):
		self.methFile = methFile
		self.methBin = methFile+'.bin'
		self.methBinIndex = self.methBin+'.idx'
		self.FA = faFile
		self.chromDict = {}
		self.seekDict = {}
		self.readFAI()
		self.makeIndex()
	def fetch(self, chrom, start=-1, end=-1):
		#1 indexed
		if chrom not in self.seekDict:
			print "Not a real chromosome"
			return []
		if end < start:
			print "Bad coordinates"
			return []
		if start == -1:
			seekStart = self.seekDict[chrom]
		else:
			seekStart = self.seekDict[chrom]+(start-1)*4
		if end == -1:
			countEnd = self.chromDict[chrom]
		else:
			if start == -1:
				countEnd = end
			else:
				countEnd = end-start+1
		IB = open(self.methBin,'rb')
		IB.seek(seekStart)
		out = []
		for i in xrange(countEnd):
			vals = IB.read(4)
			tmp = struct.unpack('HH',vals)
			if tmp == (65535, 65535):
				out.append(-1)
			else:
				out.append(float(tmp[0])/float(tmp[1]))
		return out
	def readFAI(self):
		#FAI Format  http://www.biostars.org/p/1495/
		#chrName chrLen chrSeek lineBases lineLen
		#Chr1    30427671        6       79      80
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
		IM = open(self.methFile, 'r')
		OI = open(self.methBin, 'wb')
		curPos = 1
		line = IM.readline()
		order = []
		if line[:3] != 'chr': #check for header and get first chrom
			IM.seek(0)
			curChrom = IM.readline().split()[0]
			IM.seek(0)
		else:
			curChrom = IM.readline().split()[0]
			IM.seek(0)
			IM.readline()
		order.append(curChrom)
		for line in IM:
			tmp = line.split('\t')
			chrom = tmp[0]
			pos = int(tmp[1])
			C = int(tmp[6])
			CorT = int(tmp[7])
			if pos < curPos: # add chrom order and make sure it prints for rest of chromosome
				self.fillChrom(curChrom, curPos, OI)
				curPos = 1
				curChrom = chrom
				order.append(curChrom)
			while curPos < pos:
				self.writeBlank(OI)
				curPos += 1
			self.writeData(OI, C, CorT)
			curPos += 1
		self.fillChrom(curChrom, curPos, OI)
		IM.close()
		OI.close()
		self.makeBinIndex(order)
	def fillChrom(self, curChrom, curPos, F):
		curLim = self.chromDict[curChrom]
		for i in xrange(curLim-curPos+1):
			self.writeBlank(F)
	def writeBlank(self,F):
		# 65535 is the max vlue for H
		F.write('\xff\xff\xff\xff')
	def writeData(self,F,C,CorT):
		F.write(struct.pack('HH',C,CorT))
	def makeBinIndex(self, chromOrder):
		location = 0
		OBI = open(self.methBinIndex,'w')
		for chrom in chromOrder:
			size = self.chromDict[chrom]
			OBI.write(chrom+'\t')
			OBI.write(str(location)+'\n')
			self.seekDict[chrom] = location
			location = location+size*4
		OBI.close()

def main():
	if len(sys.argv) != 3:
		sys.exit(' '.join((sys.argv[0],'meth.txt','fasta.fa')))
	inFile = sys.argv[1]
	inFA = sys.argv[2]
	MI = MethIndex(inFile, inFA)
	sampleFetch(MI,'Chr1',s=1,e=20)
	sampleFetch(MI,'Chr1')
	sampleFetch(MI,'Chr1',e=20)
	sampleFetch(MI,'Chr1',s=18,e=20)

def sampleFetch(mi, c, s=-1, e=-1):
	ret = mi.fetch(c,start=s,end=e)
	print ret
	print "length:",len(ret)

if __name__ == "__main__":
	main()
