# Multialignment processing tool for pangenomes

## Features
* Conversion from MAF (Multiple Alignment Format) to PO (POAGraph representation file - check [See Lee, Grasso & Sharlow article](https://academic.oup.com/bioinformatics/article/18/3/452/236691/Multiple-sequence-alignment-using-partial-order) for details).
* Conversion from MAF to FASTA.
* Consensus generation from aligned sequences (MAF or PO input) - in a single iteration
* Consensus generation from aligned sequences (MAF or PO input) - iteratively
* POA Graph visualization (MAF or PO input)
* Handling Ebola sequences and group names as provided by [Ebola Portal](https://genome.ucsc.edu/ebolaPortal/)

## Work in progress
* Advanced iteartive consensus generation (tree-based)
* Handling Mycoplasma sequences names

## Dependencies
* [BioPython](http://biopython.org/wiki/Download)

## Running

### Example
python3 src/pangenome.py -f examples/entire_ebola/ebola_ncbi.maf -format maf -c -iter -hbmin 0.9 -data ebola

### Arguments description
Currently, all features are provided by module *mln*. There are a few options for different features available

usage: pangenome.py mln [-h] -f FILE -format FILE FORMAT [-m MERGE_BLOCKS]
                        [-fasta] [-c] [-iter] [-hbmin HBMIN]
                        [-min_comp MINCOMP] [-v] -data DATATYPE

optional arguments:

  -h, --help         show this help message and exit
  
  -f FILE            path to the MAF (Multiple Alignment Format) file or PO (POAGraph) file
  
  -format FILE FORMAT  maf or po
  
  -m MERGE_BLOCKS    default behaviour is to merge all blocks, where possible; provide MERGE_BLOCKS if special way of merging the blocks is required; pass [idx1:idx2, idx3:idx4] to choose range of blocks to be merged; IDs begin from 1.
  
  -fasta             generate FASTA files
  
  -c                 generate consensus
  
  -iter              if c: generate consensus iteratively
  
  -hbmin HBMIN       if c: HBMIN value for POA heaviest bundling alogrithm, min 0,
                     max 1
                     
  -min_comp MINCOMP  if c and iter: minimum compatibility between source and consensus to
                     match them
                     
  -v                 generate visualization
  
  -data DATATYPE     ebola or mycoplasma
 
## Development
TBA

## Tests

### Dependencies
* [Unittest](https://docs.python.org/3/library/unittest.html)
* [DDT](https://github.com/txels/ddt)

### Set up
To setup unitests in PyCharm:
Run -> Edit Configruations -> Add New Configruation ('+') -> 
Target: Path
Write path of the test file to run
Setup Working Directory to .../pangenome/tests

### Running
To run unittests (e.g. end_to_end_tests.py file) from command line:
python3 -m unittest end_to_end_tests.py

