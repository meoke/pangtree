#Pang

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
usage: pang [-h] --multialignment MULTIALIGNMENT --data DATA [--output OUTPUT]
            [-fasta] [-vis] [-consensus {simple,tree}] [-hbmin HBMIN] [-r R R]
            [-multiplier MULTIPLIER] [-stop STOP] [-re_consensus]

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
                        

  **-re_consensus**         Tree POA algorithm parameter.Set if after producing
                        children nodes, sequences should be moved to siblings
                        nodes if compatibility to its consensus is higher.
                        

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
- *re_consensus* - przełącznik. Jeśli ma wartość True, to po wygenerowaniu dzieci danego węzła dokonywane jest sprawdzenie, czy sekwencje nie powinny przynależeć do innego węzła spośród rodzeństwa, niż zostały pierwotnie przypisane.

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
        compatibilities = get_compatibilities(node.sequences, max_consensus)
        
        node_cutoff = find_node_cutoff(compatibilities)
        compatible_sequences = sequences > node_cutoff
        
        nodes.append(TreeNode(compatible_sequences))
        sequences = sequences - compatible_sequences
    
    if re_consensus:
        nodes.move_sequences_if_needed()
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