[![Build Status](https://travis-ci.org/meoke/pangtree.svg?branch=master)](https://travis-ci.org/meoke/pangtree)

# PangtreeeBuild

This repository contains tool for multiple sequence alignment analysis. It implements the idea of pan-genome ([Ref. 1](https://doi.org/10.1093/bib/bbw089)) by representing the multialignment as a PO-MSA structure (Partial Order Alignment Graph - [Ref. 2](https://doi.org/10.1093/bioinformatics/btg109)). The main purpose of this software is to construct an *Affinity Tree* - a phylogenetic-like tree, with an agreed sequence (*consensus sequence*) assigned for each node. The result is saved in JSON file (see its schema in pangtree/pangtreebuild/serialization/affinity_tree_schema.json). Its content can be visualised using [PangtreeVis](https://github.com/meoke/pangtreevis).

This software is a part of the article: P.Dziadkiewicz, N.Dojer 'Getting insight into the pan-genome structure with Pangtree' that will be published soon in BMC Genomics.


## Getting Started

### Prerequisites

Running:
* [Python 3.8](http://python.org)
* [BioPython](https://biopython.org/)
* [numpy](http://www.numpy.org/)
* [jsonpickle](http://jsonpickle.github.io/)
* [networkx](https://networkx.github.io/)

Testing:
* [DDT](https://github.com/txels/ddt)


### Installing

```
pip install pangtreebuild
```


### Quick installation check

This line builds a pan-genome model for an example alignment of 160 Ebola virus sequences and saves it to a JSON file.

```python3 -m pangtreebuild --multialignment example_data/Ebola/multialignment.maf```



## Usage

1. Import package **pangtreebuild** to your Python program and use it according to the documentation.

or

2. Use **pangtreebuild** via command line with following arguments:

python3 -m pangtreebuild [args]

| Name  | CLI | Required | Description
| ------------- | ------------- | ------- | ----------
| Arguments affecting PO-MSA construction: |
| MULTIALIGNMENT  | --multialignment  | Yes | Path to the mulitalignment file (.maf or .po)
| METADATA | --metadata | No | Optional information about sequences in csv format. The only required column: \'seqid\' and its value must match multialignment files identifiers as described in *Sequence Naming Convention* (below). Example: example_data/Ebola/metadata.csv
| RAW_MAF | --raw_maf | No, default=False | Build PO-MSA without transforming multialignment (MAF file) to DAG. PO-MSA built in this way does not reflect real life sequences.
| FASTA_PROVIDER | --fasta_provider | No | Nucleotides source if any residues are missed in the multialignment file. Possible values: 'ncbi', 'file'. If not specified: MISSING_NUCLEOTIDE is used.
| MISSING_SYMBOL | --missing_symbol | No, default='?' | Symbol for missing nucleotides used if no FASTA_PROVIDER is given.
| CACHE | --cache | No | If set, sequences downloaded from NCBI are stored on local disc and reused between program calls, used if FASTA_PROVIDER is 'ncbi'
| FASTA_PATH | -fasta_path | Yes if FASTA_PROVIDER='FILE' | Path to fasta file or zipped fasta files with whole sequences present in multialignment, used if FASTA_PROVIDER is 'FILE'.
| Arguments affecting Affinity Tree construction: |
| AFFINITY | -affinity | No | Possible values: 'TREE' (default algorithm, descibed in Documentation.md), 'POA' (simplified version, based solely on Ref. 2)
| BLOSUM | --blosum | No, default=bin\blosum80.mat |  Path to the blosum filem. Blosum file must include MISSING_NUCLEOTIDE. |
| HBMIN | --hbmin | No, default=0.9 | 'POA' parameter. The minimum value of sequence compatibility to generated consensus.
| STOP | --stop | No, default=0.99 | 'TREE' parameter. Minimum value of compatibility in tree leaves.
| P | -p | No, default=1 | 'TREE' parameter. It changes the linear meaning of compatiblities during cutoff finding because the compatibilities are raised to the power o P. For P from range [0,1] it decreases distances between small compatibilities and increases distances between the bigger ones. For p > 1 it increases distances between small compatibilities and decreases distances between the bigger ones.
| Arguments affecting output generation: |
| OUTPUT_DIR | --output_dir, -o | No, default=timestamped folder in current working directory | Output directory path.
| OUTPUT_FULL | --output_full | No, default=False | Set, if list of pangenome nodes for sequences and consensuses should be included in pangenome.json.
| VERBOSE | --verbose, -v | No, default=False | Set if detailed log files must be produced.
| QUIET | --quiet, -q | No, default=False | Set to turn off console logging.
| FASTA | --output_fasta | No, default=False | Set to create fasta files with consensuses.
| PO | -output_po | No, default=False | Set to create po file with multialignment (without consensuses).

#### Sequence Naming Convention

[seqid].[anything after first dot is ignored]

### Example use cases
1. Build PO-MSA using default settings (transform to DAG, download missing nucleotides from NCBI) and save to .po file :
```
python -m pangtreebuild --multialignment example_data/Ebola/multialignment.maf -po

```
will produce:

- pangenome.json
- poagraph.po

2. Generate Affinity Tree, use metadata, detailed logging and default algorithm settings.
```
python3 -m pangtreebuild --multialignemnt example_data/Ebola/multialignment.maf -metadata example_data/Ebola/metadata.csv -affinity tree -v
```
will produce:

- pangenome.json
- details.log
- affinitytree/
    - tresholds.csv
    - *.po files from internal calls to poa software*


## Tests
```
python3 -m unittest discover -s pangtreebuild -p tests_*
```
or

```
nosetests pangtreebuild
```
## Authors
This software is developed with support of [OPUS 11 scientific project of National Science Centre:  Incorporating genomic variation information 
into DNA sequencing data analysis](https://www.mimuw.edu.pl/~dojer/rmg/)


## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Bibliography
1. [**Computational pan-genomics: status, promises and challenges**](https://doi.org/10.1093/bib/bbw089) 
The Computational Pan-Genomics Consortium. Briefings in Bioinformatics, Volume 19, Issue 1, January 2018, Pages 118–135.

2. [**Generating consensus sequences from partial order multiple sequence alignment graphs**](https://doi.org/10.1093/bioinformatics/btg109) C. Lee, Bioinformatics, Volume 19, Issue 8, 22 May 2003, Pages 999–1008
                        
