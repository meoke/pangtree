## Idea and algorithm description
[Pan-genome](https://en.wikipedia.org/wiki/Pan-genome) is a gene data structure being able to store multiple genomes & related data and be efficiently processed. This is a challenging bioinformatics task not only to design a pan-genome itself but also the algorithms it can be put into.

The idea of using partial order graphs as multiple sequence alignment representation (hereafter: *PO-MSA*) is introduced in [Ref. 2](https://doi.org/10.1093/bioinformatics/18.3.452). Such a graph can be constructed from an ordinary alignment in the following way:

![konstrukcja](docs/images/Poagraph.png "PO-MSA construction")

Such a representation stores information about multiple genome sequences but is concise and intuitive at the same time. C. Lee in the above mentioned article asks a question - what is the minimum set of paths that can be identified with homogeneous subgroups of the input sequences. It provides algorithm called Heaviest Bundle (implemented in software *POA*) to find the paths in a multialignment given as PO-MSA. PangtreeBuild extends the idea -by iterative application of the Heaviest Bundle algorithm, a hierarchic division of the input sequences is performed.

The division is represented as a tree called Consensus Tree. Each node of this tree has the following attributes assigned:
- list of sequences this node represents,
- an agreed sequences - consensus - of the assigned sequences,
- MinComp (minimum compatibility) value to represent how diverse this node is.


PangtreeBuild has the following features:
- PO-MSA construction from MAF ([Multiple Sequence Alignment](https://genome.ucsc.edu/FAQ/FAQformat.html#format5)) file or PO (example: data/simulated_small/multialingnment.po) with possible metadata adnotations provided in CSV file
- Consensus Tree construction using *poa* (single application of Heaviest Bundle) or *tree* (iterative application of Heaviest Bundle) algorithm
- saving PO-MSA to PO file
- saving PO-MSA and Consensus Tree to JSON file (this file can be visualised in [PangtreeVis](https://github.com/meoke/pangtreevis)


### PangtreeBuild flow diagram:
![pangenome_flow](docs/block_diagrams/Program.png "Pangenome flow")

### PO-MSA construction
PO-MSA construction from PO file is straightforward, as this format directly describes PO-MSA:
![poagraph_construction](docs/block_diagrams/Build_poagraph_from_po.png "Poagraph construction from po")

PO-MSA construction from MAF is trickier as this format does not ensure DAG and may not include all nucleotides/proteins from aligned sequences. Solution to the first problem is transforming the multialignment to DAG using [Mafgraph](https://github.com/anialisiecka/Mafgraph) and to the second one - complementing missing symbols from NCBI or local fasta files. PO-MSA consruction from MAF:

![poagraph_construction](docs/block_diagrams/Build_poagraph_from_maf.png "Poagraph construction from maf")

### Consensus Tree construction

Consensus tree algorithm diagram:

![consensus_tree_alg](docs/block_diagrams/Generate_consensus_tree.png "Consensus tree algorithm diagram")

Get children nodes details diagram:

![consensus_tree_alg_get_children](docs/block_diagrams/Generate_consensus_get_children_nodes.png "Get children nodes details diagram")

#### PO file format specification
This section originates from [POA software](https://sourceforge.net/projects/poamsa/) README.

_File intro:_
```
VERSION= Current version of POA,e.g. LPO.1.0

NAME=  Name of PO-MSA.  Defaults to name of 1st sequence in PO-MSA

TITLE=  Title of PO-MSA.  Defaults to title of 1st sequence in PO-MSA

LENGTH= Number of nodes in PO-MSA

SOURCECOUNT= Number of sequences in PO-MSA
```
_For each sequence in the PoaGraph:_
```
SOURCENAME= Name of sequence taken from FASTA sequence header

SOURCEINFO= Number of nodes in sequence 
            [Index of first node containing sequence] [Sequence weight] [Index of bundle containing sequence] [Title of sequence taken from FASTA sequence header]
```
_Example:_
```
SOURCENAME=GRB2_HUMAN

SOURCEINFO=217 10 0 3 GROWTH FACTOR RECEPTOR-BOUND PROTEIN 2 (GRB2 ADAPTOR PROTEIN)(SH2)
```
_For each node in the PO-MSA:_
```
Residue label:'L' delimited index list of other nodes with edges into node
              'S' delimited index list of sequences stored in each node
              'A' index of next node in same align ring
                     NB: align ring indices must form a cycle.
                     e.g. if two nodes 121 and 122 are aligned, then 
                     the line for node 121 indicates "A122", and
                     the line for node 122 indicates "A121".
``` 
_Example:_
```
F:L156L155L22S2S3S7A158
```

## Bibliography
1. [**Computational pan-genomics: status, promises and challenges**](https://doi.org/10.1093/bib/bbw089) 
The Computational Pan-Genomics Consortium. Briefings in Bioinformatics, Volume 19, Issue 1, January 2018, Pages 118–135.

2. [**Multiple sequence alignment using partial order graphs**](https://doi.org/10.1093/bioinformatics/18.3.452) Christopher Lee,  Catherine Grasso,  Mark F. Sharlow.
Bioinformatics, Volume 18, Issue 3, March 2002, Pages 452–464.

3. [**The computation of consensus patterns in DNA sequences**](https://doi.org/10.1016/0895-7177(93)90117-H) William H.E.DayF.R.McMorris. Mathematical and Computer Modelling
Volume 17, Issue 10, May 1993, Pages 49-52

                        
