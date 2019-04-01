# Pang

Tool for analysis and visualisation of multiple sequence alignment. It implements the idea of pan-genome ([Ref. 1](https://doi.org/10.1093/bib/bbw089)) by representing the multialginment as a graph and construction of a phylogenetic tree joined with an agreed sequence for every node.

[PL]
Narzędzie służące do analizy i wizualizacji uliniowienia wielu sekwencji genetycznych. Implementuje ideę pangenomeu ([Ref. 1](https://doi.org/10.1093/bib/bbw089)) poprzez grafową reprezentację multiuliniowienia oraz konstrukcję drzewa filogenetycznego z kompromisową sekwencją dla każdego węzła. 


## Getting Started

### Prerequisites

Running:
* [BioPython](https://biopython.org/)
* [Mafgraph](https://github.com/anialisiecka/Mafgraph)
* [numpy](http://www.numpy.org/)
* [jsonpickle](http://jsonpickle.github.io/)

Testing:
* [DDT](https://github.com/txels/ddt)


### Installing

```
TBA
```

### Quick installation check

```python3 -m pangenome --multialignment data/Fabricated/f.maf --metadata data/Fabricated/f_metadata.csv```

## Idea and algorithm description 

SECTION UNDER DEVELOPMENT
```
Konstrukcja grafu:
![konstrukcja](docs/images/pangraph_construcion.png "Pangraph construction")


##### Algorytm
Drzewo jest budowane w porządku Breadth First. Pseudo-Python-Code:

```Python
sequences = pangraph.all_sequences
nodes_to_process = [TreeNode(sequences)]
while nodes_to_process:
    subtree_root = nodes_to_process.pop()
    children_nodes = get_children(subtree_root)
    if len(children_nodes.sequences) == 1:
        break
    for child in children_nodes:
        subtree_root.children.append(child)
        
        if node_ready(child):
            nodes_to_process.append(child)
            
def node_ready(node):
    if node.min_compatibility <= stop:
        return True
    return False
    
def get_children(node):
    sequences = node.sequences
    nodes = []
    while sequences:
        consensus = poa(Pangraph(sequences))
        compatibilities = get_compatibilities(sequences, consensus)
        
        max_cutoff = find_max_cutoff(compatibilities)
        the_most_compatible_sequences = sequences > max_cutoff
        
        max_consensus = poa(Pangraph(the_most_compatible_sequences))
        compatibilities = get_compatibilities(node.sequences, max_consensus)**p
        
        node_cutoff = find_node_cutoff(compatibilities)
        compatible_sequences = sequences > node_cutoff
        
        nodes.append(TreeNode(compatible_sequences))
        sequences = sequences - compatible_sequences
    
    if re_consensus:
        nodes.move_sequences_if_needed()
    return nodes
    
```
```
Słowny opis podziału węzła **N** (odpowiada get_children, uruchamiane tylko gdy w węźle istnieje sekwencja 
o compatibility do consensusu w tym węźle o wartości niższej niż **STOP** ):

1. **level_guards** - pusta lista
2. Uruchom *poa* na wszystkich sekwencjach w węźle, weź consensus z o największej liczbie 
przypisanych sekwencji (wg *poa*) jako **C**.
3. Policz compatibility **C** z wszystkimi sekwencjami w tym węźle i podnieś wartości do potęgi **p**. 
4. Znajdź próg odcięcia **P1** wśród compatibilities policzonych w 3.:
    - Strategia MAX1 *(oryginalna)* [parametry: cutoff_search_range]
        - uporządkuj rosnąco compatibilities
        - znajdź największą różnicę występującą pomiędzy dwoma kolejnymi compatibility **Ci**, **Cj** na przedziale
         *cutoff_search_range*, gdzie przedział określa indeksy na uporządkowanej liście compatibilities
        - **P1** = **Cj**
    - Strategia MAX2 *(usunięcie *cutoff_search_range*, bo jego działanie jest tożsame z parametrem *stop*)* [brak parametrów]
        - uporządkuj rosnąco compatibilities
        - znajdź największą różnicę występującą pomiędzy dwoma kolejnymi compatibility **Ci**, **Cj**
        - **P1** = **Cj**
5. **max_sequences** - sekwencje, których compatibility do **C** przekracza **P1**
6. Uruchom *poa* na **max_sequences**, weź consensus o największej liczbie przypisanych sekwencji jako **C_MAX**.
7. Policz compatibility **C_MAX** z wszystkimi sekwencjami w tym węźle.
8. Znajdź próg odcięcia **P2** wśród posortowanych compatibilities policzonych w 7.:
    - Strategia NODE1 (oryginalna) [parametry: multiplier]
        - policz średnią odległość między compatibilities
        - uporządkuj compatibilities rosnąco
        - znajdź pierwsze takie **Ci**, **Cj**, pomiędzy którymi odległość jest większa niż
         średnia odległość * *multiplier*. Jeśli nie istnieją, ponów wyszukiwanie dla multiplier = 1. 
        - **P2** = **Cj**
    - Strategia NODE2 (z level guardem) [parametry: level_guards, multiplier]
        - IF lista **level_guards** jest pusta:
            - użyj NODE1
        - ELSE
            - guard = min(**level_guards**)
            - IF guard <= wszystkie compatibilities:
                - **P2** = min(compatibilities)
            - ELIF guard > wszystkie compatibilities:
                - użyj NODE1
            - ELSE
                - dodaj guard do compatibilities 
                - uporządkuj compatibilities rosnąco
                - policz średnią odległość między compatibilities
                - usuń guard z compatibilities (nie ma sensu zwracać go jako wynik, 
                jeśli nie było go oryginalnie wśród compatibilities)
                - wśród compatibilities nie większych niż guard znajdź pierwsze takie **Ci**, **Cj**, pomiędzy którymi odległość 
                jest większa niż średnia odległość * *multiplier*
                - jeśli nie ma takich **Ci**, **Cj**:
                    - **P2** = pierwsze compatibility większe niż guard
    - Strategia NODE3 (z level guardem, uproszczona) [parametry: level_guards]
        - IF lista **level_guards** jest pusta:
            - użyj strategii MAX2
        - ELSE
            - guard = min(**level_guards**)
            - IF guard <= wszystkie compatibilities:
                - **P2** = min(compatibilities)
            - ELSE
                - search_boundary = indeks pierwszego compatibility większego niż guard 
                albo max(compatibilities) (gdy guard > wszystkie compatibilities)
                - użyj strategii MAX2 na przedziale [0, search_boundary]
    - Strategia NODE4 (bardzo uproszczona) [brak parametrów]
        - użyj strategii MAX2
9. **node_sequences** - sekwencje, których compatibility do **C_MAX** (policzone w 7.) przekracza **P2**
10. Parametry uzyskanego węzła consensusowego (dziecka **N**):
    - ID consensusu **C_MAX**
    - sekwencje **node_sequences**
    - minimalne compatibility wśród compatibilities **node_sequences** do **C_MAX**
11. Z węzła **N** usuń **node_sequences**
12. Zapisz **P2** do **level_guards**
13. Jeśli pozostały jakieś sekwencje w węźle **N** - idź do 2.
14. Jeśli **re_consensus**:
    - Dla każdej sekwencji, która należała początkowo do **N**, sprawdź, czy spośród wartości compatibilities 
    do consensusów utworzonych przy podziale **N**, najwyższa jest ta, która jest związana z consensusem, 
    do którego ta sekwencja została przyporządkowana. 
    Jeśli nie, przyporządkuj tę sekwencję do consensusu, do którego comptibility jest najwyższe.
```

## Usage

1. Import package **pangenome** to your Python program. Check [API documentation]().

or

2. Use **pangenome** from command line with following arguments:

python3 -m pangenome [args]

| Name  | CLI | Required | Description
| ------------- | ------------- | ------- | ----------
| Arguments affecting poagraph build process: |
| MULTIALIGNMENT  | --multialignment, -m  | Yes | Path to the mulitalignment file (.maf or .po)
| DATATYPE  | --datatype  | No, default = 'n' | Possible values: 'n' (nucleotides), 'p' (proteins).
| METADATA | --metadata | No | Optional information about sequences in csv format. The only required column: \'seqid\' and its value must match multialignment files identifiers as described in *Sequence Naming Convention* (below). Example: data\Ebola\ebola_metadata.csv
| RAW_MAF | -raw_maf | No, default=False | Build poagraph without transforming multialignment (maf) to DAG. Poagraph build in this way does not reflect real life sequences.
| FASTA_COMPLEMENTATION | -fasta_complementation | No, default=NCBI | Nucleotides/proteins source if any are missed in the multialignment. Possible values: 'NCBI', 'FILE', 'NO'
| MISSING_NUCLEOTIDE | -missing_n | No, default='?' | Symbol for missing nucleotides, used if FASTA_COMPLEMENTATION is 'NO'.
| EMAIL | -email | Yes if FASTA_COMPLEMENTATION='NCBI' | E-mail address for NCBI API, used if FASTA_COMPLEMENTATION is 'NCBI'.
| CACHE | -cache | No, default='Yes' | If True, sequences downloaded from NCBI are stored on local disc and reused between program calls, used if Fasta Complementation Option is 'NCBI'
| FASTA_FILE | -fasta_source_file | Yes if FASTA_COMPLEMENTATION='FILE' | Path to fasta file or zipped fasta files with whole sequences present in multialignment, used if FASTA_COMPLEMENTATION is 'FILE'.
| Arguments affecting consensuses tree algorithm: |
| CONSENSUS | -consensus | No, default='TREE' | Possible values: 'TREE' (tree algorithm), 'POA' (poa algorithm)
| BLOSUM | --blosum | No, default=bin\blosum80.mat |  Path to the blosum file which is used in consensus algorithm. Blosum file must include MISSING_NUCLEOTIDE. |
| HBMIN | -hbmin | No, defaUlt=0.9 | 'POA' parameter. The minimum value of sequence compatibility to generated consensus.
| STOP | -stop | No, default=0.99 | 'TREE' parameter. Minimum value of compatibility in tree leaves.
| MAX | -max | No, default=MAX2 | 'TREE' parameter. Max cutoff finding strategy. Available values: 'MAX1', 'MAX2'.
| NODE | -node | No, default=NODE3 | 'NODE' parameter. Node cutoff finding strategy. Available values: 'NODE1', 'NODE2', 'NODE3', 'NODE4'
| R | -r | No, default=[0,1] | 'MAX1' parameter. Specifies what part of sorted capabilities should be searched for node cutoff. Format: '[a, b]' where a, b in [0, 1] and a < b. 
| MULTIPLIER | -multiplier | No, default=1 | 'NODE1' and 'NODE2' parameter. It controls the size of gaps for node cutoff. The greater it is, the more granular the tree is.
| P | -p | No, default=1 | 'TREE' parameter. It changes the linear meaning of compatiblities during cutoff finding because the compatibilities are raised to the power o P. For p from range [0,1] it decreases distances between small compatibilities and increases distances between the bigger ones.For p > 1 it increases distances between small compatibilities and decreases distances between the bigger ones.
| Arguments affecting output format: |
| OUTPUT_DIR | --output, -o | No, default=timestamped folder in current working directory | Output directory path.
| VERBOSE | --verbose, -v | No, default=False | Set if detailed log files must be produced.
| QUIET | --quiet, -q | No, default=False | Set to turn off console logging.
| FASTA | --fasta | No, default=False | Set to create fasta files with consensuses.
| PO | -po | No, default=False | Set to create po file with multialignment (without consensuses).
| INCLUDE_PATHS | -output_with_paths | No, default=False | Set if output json should include paths (it significantly increases file size).

#### Sequence Naming Convention

[anything before first dot is ignored].[everything after first dot (also other dots) is interpreted as seqid]

### Example use cases
1. Build poagraph using default settings (transform to DAG, download missing nucleotides from NCBI) and save to .po file :
```
python3 -m data\Ebola\Ebola.maf -po
```
will produce:

- pangenome.json
- poagraph.po

2. Generate consensuses tree, use metadata, detailed logging and default algorithm settings.
```
python3 -m data\Ebola\Ebola.maf -metadata data\Ebola\Ebola.maf' -consensus tree -v
```
will produce:

- pangenome.json
- details.log
- consensus/
    - tresholds.csv
    - *.po files from internal calls to poa software*


## Development

### Documentation
```
TBA
```
### Running the tests
```
python -m unittest discover -s tests -p '*_test.py'
```

### Contributing
```
TBA
```


## Authors
This software is developed with support of [OPUS 11 scientific project of National Science Centre:  Incorporating genomic variation information 
into DNA sequencing data analysis](https://www.mimuw.edu.pl/~dojer/rmg/)


## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Bibliography
1. [**Computational pan-genomics: status, promises and challenges**](https://doi.org/10.1093/bib/bbw089) 
The Computational Pan-Genomics Consortium. Briefings in Bioinformatics, Volume 19, Issue 1, January 2018, Pages 118–135.

2. [**Multiple sequence alignment using partial order graphs**](https://doi.org/10.1093/bioinformatics/18.3.452) Christopher Lee,  Catherine Grasso,  Mark F. Sharlow.
Bioinformatics, Volume 18, Issue 3, March 2002, Pages 452–464.


                        
