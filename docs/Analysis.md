# SortByBarcodes

This class is designed to create new Fastq files containing specific barcode sequences. The barcode list should be provided as input files in DNA sequence format. The barcode pattern is specified as a string based on regular expressions. Currently, only the Fastq format is supported for input data files.

## Example

```python
from genet import analysis as ans
ans.SortByBarcodes('./MyNGS.fastq', './Barcode_pattern.csv', 'TCGTATGCCGTCTTCTGCTTG[ATGC]{14}', n_cores=10)
```

The output file will be generated in the current working directory by default. If you want to save your output in another location, you can specify the 'output_path' option.

Example:

```python
ans.SortByBarcodes(seq_file='./MyNGS.fastq',
                   barcode_file='./Barcode_pattern.csv',
                   barcode_pattern='TCGTATGCCGTCTTCTGCTTG[ATGC]{14}',
                   output_name='My_sorted_data',
                   output_path='/extdata/mydir/results',
                   n_cores=20
                   )
```

### Sorting the Fastq File by Barcode List

Is it possible to simultaneously sort and write directly to the final file using the `with open` method? This approach would eliminate the need to create a separate temporary file, potentially resulting in significant time savings and reduced I/O overhead. It would be advisable to implement and validate this approach on a small scale first.