#!/usr/bin/python

from pymethyl import MethIndex

MI = MethIndex("test_meth.txt","test.fa")

print MI.fetch("Chr1",start=1,end=10)
print MI.fetch("Chr1",end=10)

print MI.fetch("Chr2",start=11,end=20) #20 bases long
print MI.fetch("Chr2",start=11)

print MI.fetch("Chr1",start=1,end=20) #Whole Chromosome
print MI.fetch("Chr1")
