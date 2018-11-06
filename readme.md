#Pang

## Wymagania
* [BioPython](https://biopython.org/)
* [numpy](http://www.numpy.org/)
* [jsonpickle](http://jsonpickle.github.io/)
* [Mafgraph-todo]()
* [Unittest](https://docs.python.org/3/library/unittest.html)
* [DDT](https://github.com/txels/ddt)
## Użycie programu
usage: pang [-h] --multialignment MULTIALIGNMENT --data DATA [--output OUTPUT]
            [-fasta] [-vis] [-consensus {simple,tree}] [-hbmin HBMIN] [-r R R]
            [-multiplier MULTIPLIER] [-stop STOP]

Consensus generation and visulization of Pangenome

optional arguments:

  **-h, --help**            show this help message and exit
  
  **--multialignment MULTIALIGNMENT, -m MULTIALIGNMENT**
                        Path to the mulitalignment file. Accepted formats:
                        .maf, .po.
                        
  **--data DATA, -d DATA**  Path to the json file with genomes specification.
                        See... examples\Ebola\ebola_metadata.json
                        
  **--output OUTPUT, -o OUTPUT**
                        Output directory path.
                        
  **-fasta**                Set if fasta files must be produced.
  
  **-vis**                  Set if visualization must be produced.
  
  **-consensus {simple,tree}**
                        Set if consensus must be generated. Algorithms to
                        choose: 'simple' or 'tree'.
                        
  **-hbmin HBMIN**          Simple POA algorithm parameter. The minimum value of
                        sequence compatibility to generated consensus
                        
  **-r R R**                Tree POA algorithm parameter.Specify what part of
                        sorted capabilities should be searched for node
                        cutoff. E.g. [0.2,0.8]
                        
  **-multiplier MULTIPLIER**
                        Tree POA algorithm parameter. The greater it is, the more granular the tree is.
                        
  **-stop STOP**            Tree POA algorithm parameter.Value of node
                        compatibility above which the node is no more split.
                        
                        

## Przykłady
python3 pang -m examples/Fabricated/f.maf -d examples/Fabricated/f_metadata.json -consensus tree
                        
## Opisy funkcjonalności
### Konstrukcja poagraphu
Input: plik .maf (Multialignment Alignment Format)

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
- *poagraph* - kosntruowany przez program z pliku .maf
- *cutoff_search_range* - Wśród posortowanych wartości compatibilities poszukiwany jest próg odcięcia, czyli taka wartość compatibility, poniżej której sekwencje nie są przydzielane do danego węzła. *cutoff_search_range* określa, w której części listy wartość będzie poszukiwana. Np. dla [0.2, 0.5, 0.8, 0.9, 0.95] i *cutoff_search_range* = [0.5, 1], próg odcięcia będzie szukany wśród [0.8, 0.9, 0.95]
- *multiplier* - Próg odcięcia jest liczony wg procedury: oblicz średnią odległość między wartościami compatibilities i znajdź pierwszą odległość, która tę średnią*multiplier przekracza.
- *stop* - compatibility węzła, przy osiągnięciu którego węzeł nie jest już dzielony

##### Output:
- lista consensusów
- drzewo, w którym każdy węzeł zawiera:
    - ID consensusu, któremu odpowiada ten węzeł
    - listę sekwencji, które zostały zakwalifikowane do tego węzła (consensusu)
    - listę compatibilities sekwencji pangraphu do consensusu tego węzła
    - minimalne compatibility wśród sekwencji zakwalifikowanych do tego węzła

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
        
        if not_ready(child):
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
        compatibilities = get_compatibilities(node.sequences, max_consensus)
        
        node_cutoff = find_node_cutoff(compatibilities)
        compatible_sequences = sequences > node_cutoff
        
        nodes.append(TreeNode(compatible_sequences))
        sequences = sequences - compatible_sequences
    return nodes
    
def find_max_cutoff(compatibilities):
    sc = sorted(compatibilities)
    search_range_values = get_search_range(sc, cutoff_search_range)
    max_d = get_max_distance(sorted(search_range_values))
    (c1, c2) = get_first_pair_with_distance_greater_than_max_d(search_range_values, max_d)
    return c2
    
def find_node_cutoff(compatibilites):
    sc = sorted(compatibilities)
    avg_d = get_average_distance(compatibilites)
    d = avg_d * multiplier
    (c1, c2) = get_first_pair_with_distance_greater_than_d(sc, d)
    return c2
```

Słowny opis podziału węzła (odpowiada get_children):


1. Uruchom *poa* na wszystkich sekwencjach w tym węźle, weź consensus o największej liczbie przypisanych sekwencji jako **C**.
2. Policz compatibility **C** z wszystkimi sekwencjami w tym węźle.
3. Znajdź próg odcięcia **P1** wśród compatibilities policzonych w 2.:
    - weź compatibilities z przedziału *cutoff_search_range*
    - znajdź największą różnicę występującą pomiędzy dwoma kolejnymi compatibility **Ci**, **Cj**
    - **P1** = **Cj**
4. **max_sequences** - sekwencje, których compatibility do **C** przekracza **P1**
5. Uruchom *poa* na **max_sequences**, weź consensus o największej liczbie przypisanych sekwencji jako **C_MAX**.
6. Policz compatibility **C_MAX** z wszystkimi sekwencjami w tym węźle.
7. Znajdź próg odcięcia **P2** wśród posortowanych compatibilities policzonych w 6.:
    - policz średnią odległość między compatibilities
    - znajdź takie **Ci**, **Cj** pomiędzy odległość jest większa niż średnia odległość * *multiplier*
    - **P2** = **Cj**
8. **node_sequences** - sekwencje, których compatibility do **C_MAX** (policzone w 6.) przekracza **P2**
9. Parametry uzyskanego węzła:
    - ID consensusu **C_MAX**
    - sekwencje **node_sequences**
    - minimalne compatibility wśród compatibilities **node_sequences** do **C_MAX**
