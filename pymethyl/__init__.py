#!/usr/bin/python

import sys
import os.path
import struct

class MethIndex:
    """
    Creates and queries methylation index for quick random access.
    """

    def __init__(self, methFile, faFile):
        """
        Index constructor
        
        Arguments
        =================================
        methFile        - Methylation input file from MethylCoder
        faFile          - Fasta file of reference (index required)
        """
        self.methFile = methFile
        self.methBin = methFile + '.bin'
        self.methBinIndex = self.methBin + '.idx'
        self.FA = faFile
        self.FAI = self.FA + '.fai'
        self.chromDict = {}
        self.seekDict = {}
        self.readFAI()
        if not os.path.exists(self.methBin) or not os.path.exists(self.methBinIndex):
            self.makeIndex()
        else:
            self.readBinIndex()

    def fetch(self, chrom, start = -1, end = -1):
        """
        Fetches a region of methylation ratios.
        
        If no start is given, 1-end is returned.
        If no end is given, start- is returned
        
        Agurments
        =================================
        chrom   - Chromosome
        start   - Start of region (1-indexed)
        end     - End of region (1-indexed)
        """
        if chrom not in self.seekDict:
            print 'Not a real chromosome'
            return []
        if end < start:
            print 'Bad coordinates'
            return []
        if start == -1:
            seekStart = self.seekDict[chrom]
        else:
            seekStart = self.seekDict[chrom] + (start - 1) * 4
        if end == -1:
            countEnd = self.chromDict[chrom]
        elif start == -1:
            countEnd = end
        else:
            countEnd = end - start + 1
        IB = open(self.methBin, 'rb')
        IB.seek(seekStart)
        out = []
        for i in xrange(countEnd):
            vals = IB.read(4)
            tmp = struct.unpack('HH', vals)
            if tmp == (65535, 65535):
                out.append(-1)
            else:
                out.append(float(tmp[0]) / float(tmp[1]))

        return out

    def readFAI(self):
        """
        Reads the sizes of each chromosome from the FASTA index
        """
        if not os.path.exists(self.FAI):
            sys.exit('Please index the fasta file')
        IFAI = open(self.FAI, 'r')
        for line in IFAI:
            tmp = line.rstrip('\n').split()
            chrom = tmp[0]
            chromLen = int(tmp[1])
            self.chromDict[chrom] = chromLen

        IFAI.close()

    def makeIndex(self):
        """
        Makes the .bin and .bin.idx methylation index files
        """
        print 'Making Index'
        IM = open(self.methFile, 'r')
        OI = open(self.methBin, 'wb')
        curPos = 1
        line = IM.readline()
        order = []
        if line[:3] != 'chr':
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
            if pos < curPos:
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
        """
        Writes data for rest of chromosome
        """
        curLim = self.chromDict[curChrom]
        for i in xrange(curLim - curPos + 1):
            self.writeBlank(F)

    def writeBlank(self, F):
        """
        Writes two blank values (65535) to .bin file.
        65535 is the largest USHORT.
        """
        F.write('\xff\xff\xff\xff')

    def writeData(self, F, C, CorT):
        """
        Writes two unsigned shorts to .bin file.
        """
        F.write(struct.pack('HH', C, CorT))

    def makeBinIndex(self, chromOrder):
        """
        Makes the bin index based on the order the chromosomes
        were written.
        """
        location = 0
        OBI = open(self.methBinIndex, 'w')
        for chrom in chromOrder:
            size = self.chromDict[chrom]
            OBI.write(chrom + '\t')
            OBI.write(str(location) + '\n')
            self.seekDict[chrom] = location
            location = location + size * 4

        OBI.close()

    def readBinIndex(self):
        """
        Loads the bin index into a seek dictionary.
        """
        IF = open(self.methBinIndex, 'r')
        for line in IF:
            tmp = line.rstrip('\n').split('\t')
            self.seekDict[tmp[0]] = int(tmp[1])

        IF.close()


def main():
    pass


if __name__ == '__main__':
    main()
