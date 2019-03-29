# Pang


Narzędzie służące do analizy i wizualizacji uliniowienia wielu sekwencji genetycznych. Implementuje ideę pangenomeu ([Ref. 1](https://doi.org/10.1093/bib/bbw089)) poprzez grafową reprezentację multiuliniowienia oraz konstrukcję drzewa filogenetycznego wraz z kompromisową sekwencją dla każdego węzła. 

Tool for analysis and visualisation of multiple sequence alignment. It implements the idea of pan-genome ([Ref. 1](https://doi.org/10.1093/bib/bbw089)) by representing the multialginment as a graph and construction of a phylogenetic tree joined with an agreed sequence for every node.


## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

What things you need to install the software and how to install them

To use Pang via command line:
* [BioPython](https://biopython.org/)
* [numpy](http://www.numpy.org/)
* [jsonpickle](http://jsonpickle.github.io/)
* [Mafgraph](https://github.com/anialisiecka/Mafgraph)
* [DDT](https://github.com/txels/ddt)
* [pandas](https://pandas.pydata.org/)
* [networkx](https://networkx.github.io/)

To use also visualisations (via web browser):
* [plotly](https://plot.ly)
* [flask](http://flask.pocoo.org/)
* [dash, dash-core-components, dash-html-components, dash-table](https://dash.plot.ly/)


### Installing

```
TBA setup
```

#### Quick installation check - using terminal
```
python3 pang/main.py --multialignment ../examples/Fabricated/f.maf -- metadata ../examples/Fabricated/f_metadata.csv -consensus tree -v
```
#### Quick installation check - using web browser
Run:
```
python3 run_dash_app.py
```
Open web browser (Google Chrome is recommended) and go to http://127.0.0.1:8056/

## Running the tests

python -m unittest discover -s tests -p '*_test.py'

### Break down into end to end tests


### And coding style tests


## Deployment



## Built With


## Contributing


## Versioning


## Authors


## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

### Bibliography
1. [**Computational pan-genomics: status, promises and challenges**](https://doi.org/10.1093/bib/bbw089) 
The Computational Pan-Genomics Consortium. Briefings in Bioinformatics, Volume 19, Issue 1, January 2018, Pages 118–135 

* Hat tip to anyone whose code was used
* Inspiration
* etc


Pang - narzędzie służące do analizy i wizualizacji uliniowienia wielu sekwencji genetycznej. Implementuje ideę pangenomeu [

Główne funkcjonalności:
* 

## Wymagania
* [BioPython](https://biopython.org/)
* [numpy](http://www.numpy.org/)
* [jsonpickle](http://jsonpickle.github.io/)
* [Mafgraph-todo]()
* [Unittest](https://docs.python.org/3/library/unittest.html)
* [DDT](https://github.com/txels/ddt)

graphviz
graphviz-dev
## Użycie programu
``
  **-h, --help**            show this help message and exit

  --multialignment MULTIALIGNMENT, -m MULTIALIGNMENT
                        Path to the mulitalignment file. Accepted formats:
                        .maf, .po.

  --datatype DATATYPE   Input type: 'n' for nucleotides, 'p' for protieins.
  
  --metadata METADATA   Path to the csv file with genomes specification.
                        See... examples\Ebola\ebola_metadata.csv
  
  --blosum BLOSUM       Path to the BLOSUM matrix used in consensus generation
                        algorithm.If fasta_complementation option is NO and a
                        custom symbol is provided, the matrix specified here
                        must include this symbol.If fasta_complementation
                        option is NO and a custom symbol is not provided, the
                        matrix specified here must include symbol '?' as this
                        is the default symbol for missing nucleotide.
  
  --output OUTPUT, -o OUTPUT
                        Output directory path.
  
  -fasta                Set if fasta files for consensuses must be produced.
  
  -output_po            Set if po file with entire pangraph (without any
                        consensuses) must be produced.
  
  -consensus CONSENSUS  Set if consensus must be generated. Values to choose:
                        'simple' or 'tree'.
  
  -hbmin HBMIN          Simple POA algorithm parameter. The minimum value of
                        sequence compatibility to generated consensus.
  
  -r R R                Tree POA algorithm parameter.Specify what part of
                        sorted capabilities should be searched for node
                        cutoff. E.g. [0.2,0.8]
  
  -multiplier MULTIPLIER
                        Tree POA algorithm parameter.Cutoff value for node
                        parameter. The greater it is, the more granular the
                        tree is.
  
  -stop STOP            Tree POA algorithm parameter.Value of node
                        compatibility above which the node is no more split.
  
  -re_consensus         Tree POA algorithm parameter.Set if after producing
                        children nodes, sequences should be moved to siblings
                        nodes if compatibility to its consensus is higher.
  
  -not_dag              Pangraph building from maf file parameter.Set if the
                        maf content must not be transformed to DAG when
                        building pangraph. Pangraph that was build in this way
                        provides consensuses tree the consensuses do not
                        reflect the real life sequences.
  
  -fasta_complementation FASTA_COMPLEMENTATION
                        Pangraph building from maf file parameter. Ignored
                        when -not_dag parameter is set.Maf file usually
                        contains not full sequences but only parts of them,
                        aligned to each other. To build an exact pangraph the
                        full sequences must be retrieved from: ncbi or local
                        file system. Don't use this argument if you want the
                        pangraph to be build without full sequences.Pass
                        "ncbi" if you want to download the lacking fragments
                        from ncbiPass "local" if you want to use fasta from
                        local file system.
  
  -email EMAIL          E-mail address requiered when Fasta Complementation
                        Option is "NCBI" as using Entrez API obligates the
                        user to pass e-mail address.
  
  -cache                Used if Fasta Complementation Option is "NCBI" Stores
                        sequences downloaded from NCBI on local disc.They are
                        reused between uses of this program.
  
  -missing_n MISSING_N  If fasta_complementation is NO, a custom symbol for
                        missing nucleotides can be specified.Make sure it is
                        included in BLOSUM matrix you use.
  
  --fasta_source_file FASTA_SOURCE_FILE, -f FASTA_SOURCE_FILE
                        ZIP archive with fasta files used to complement
                        missing parts of sequences in maf file.
  
  -p P                  Tree consensus algorithm parameter.When deciding about
                        consensus node split, the compatibilities are raised
                        to the power o p.It enables to change the linear
                        meaing of compatibility values.For p from range [0,1]
                        it decreases distances between small compatibilities
                        and increases distances between the bigger ones.For p
                        > 1 it increases distances between small
                        compatibilities and decreases distances between the
                        bigger ones.
  
  -max MAX              Specify which strategy - MAX1 or MAX2 use for finding
                        max cutoff (see details in README.md)
  
  -node NODE            Specify which strategy - NODE1 (1), NODE2 (2), NODE3
                        (3) or NODE4 (4) use for finding max cutoff (see
                        details in README.md)
  
  -v, --verbose         Set if detailed log files must be produced.
  
  -q, --quiet           Set to turn off console logging .
  
  -output_with_nodes    Set if output json should include nodes (it
                        significantly increases file size).

``
## Przykłady
python3 pang examples/Fabricated/f.maf -d examples/Fabricated/f_metadata.json -consensus tree
                        
## Opisy funkcjonalności
### Konstrukcja poagraphu
Input: plik .maf (Multialignment Alignment Format) lub .po

Output: plik .po 

Konstrukcja grafu:
![konstrukcja](docs/images/pangraph_construcion.png "Pangraph construction")

### Wizualizacja
TBA
### Generowanie consensusów
#### Simple
TBA
#### Tree

### Założenia
##### Input
TBA

##### Output:
TBA

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
