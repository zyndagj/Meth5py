# Meth5py
```python
Meth5py(self, methFile='', faFile='', h5File='', n_cores=0, force=False, verbose=False)
```

Class for converting BSMAPz methratio files to hdf5 for fast random access

__Parameters__

- __methFile (str)__: Methylation input file
- __faFile (str)__: Reference file
- __h5File (str)__: previously created h5 file
- __n_cores (int)__: The number of cores to use for indexing [0 all]
- __force (bool)__: Force the re-creation of h5 file
- __verbose (bool)__: Enable verbose logging

__Attributes__

- `logger (log)`: Logger instance
- `contexts (tuple)`: Tuple of methylation contexts
- `strands (tuple)`: Tuple of strands
- `sorted_chroms (list)`: Sorted chromosomes in reference
- `chom_dict (dict)`: Dictionary of chromosomes and their lengths
- `n_cores (int)`: The number of cores to use for indexing

__Example__

```
from Meth5py import Meth5py
m5 = Meth5py('tests/test_meth.txt', 'tests/test.fa')
for record in m5.fetch('Chr1',10,11):
 print(record)
m5.close()
```

## fetch
```python
Meth5py.fetch(self, chrom, start=1, end=-1, minCov=0, maxCov=-1, index=True)
```

Fetches a region of methylation ratios.

If no start is given, 1-end is returned.
If no end is given, start- is returned.

-1 values are returned IF:

 * depth < minCov
 * depth > maxCov
 * no reads for that position
 * no methylation for that position

__Parameters__

- __chrom (str)__: Chromosome name
- __start (int)__: Start of region (1-indexed)
- __end (int)__: End of region (1-indexed)
- __minCov (int)__: Minimum coverage needed (Default: 0)
- __maxCov (int)__: Maximum coverage allowed (Default: -1)
- __index (bool)__: Retrun context and strand indices instead of strings

__Returns__

`list`: [[context, strand, c, ct, g, ga], ...]

## close
```python
Meth5py.close(self)
```

Closes the H5 file for cleaning instances of Meth5py

__Attributes__

- `self.H5 (file)`: The file that is closed

