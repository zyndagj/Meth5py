PyMethyl
==========

This will parse [MethylCoder](https://github.com/brentp/methylcode) output and make a binary index for quick random access to methylation ratios based on position.

Installation
------------

```python
python setup.py install
```

or

```Shell
python setup.py install --user
```

Documentation
-------------

### Class - MethIndex
```
     |  Creates and queries methylation index for quick random access.
     |  
     |  Methods defined here:
     |  
     |  __init__(self, methFile, faFile)
     |      Index constructor
     |      
     |      Arguments
     |      =================================
     |      methFile                - Methylation input file from MethylCoder
     |      faFile            - Fasta file of reference (index required)
     |  
     |  fetch(self, chrom, start=-1, end=-1)
     |      Fetches a region of methylation ratios.
     |      
     |      If no start is given, 1-end is returned.
     |      If no end is given, start- is returned
     |      
     |      Agurments
     |      =================================
     |      chrom   - Chromosome
     |      start   - Start of region (1-indexed)
     |      end      - End of region (1-indexed)
     |  
     |  fillChrom(self, curChrom, curPos, F)
     |      Writes data for rest of chromosome
     |  
     |  makeBinIndex(self, chromOrder)
     |      Makes the bin index based on the order the chromosomes
     |      were written.
     |  
     |  makeIndex(self)
     |      Makes the .bin and .bin.idx methylation index files
     |  
     |  readBinIndex(self)
     |      Loads the bin index into a seek dictionary.
     |  
     |  readFAI(self)
     |      Reads the sizes of each chromosome from the FASTA index
     |  
     |  writeBlank(self, F)
     |      Writes two blank values (65535) to .bin file.
     |      65535 is the largest USHORT.
     |  
     |  writeData(self, F, C, CorT)
     |      Writes two unsigned shorts to .bin file.
```

## Usage

example/example.py
```python
from pymethyl import MethIndex

MI = MethIndex("test_meth.txt","test.fa")

print MI.fetch("Chr1",start=1,end=10)
print MI.fetch("Chr1",end=10)

print MI.fetch("Chr2",start=11,end=20) #20 bases long
print MI.fetch("Chr2",start=11)

print MI.fetch("Chr1",start=1,end=20) #Whole Chromosome
print MI.fetch("Chr1")
```
