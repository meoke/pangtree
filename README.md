# pangenome

## POA
~/Repos/pangenome/bin/poa -read_msa $input -hb -po $output -hb ~/Repos/pangenome/bin/blosum80.mat -hbmin 0.06

## Development

### Tests
The framework Unittest is used here.

To setup unitests in PyCharm:
Run -> Edit Configruations -> Add New Configruation ('+') -> 
Target: Path
Write path of the test file to run
Setup Working Directory to .../pangenome/tests

To run unittests from command line:
python3 -m unittest end_to_end_tests.py

## Instalation

#### BioPython
[BioPython Download and Instalation](http://biopython.org/wiki/Download)

#### DDT
[Download and Instalation](https://github.com/txels/ddt)

## Moduły programu:
### Konwersja maf (Multiple Alignment Format) do pox

#### Input
maf file

#### Output
- master file - blocks desciption
- pox files without consensus (if -c not specified) or pox files with consensus (if -c specified)
- html file (if -v specified)

The amount of pox files generated depends on -m option. It equals number of blocks produced by the converter which can be equal or less than the original. Some blocks can be merged by the converter. One can specify which blocks should be merged in the following way:
Pass 'idx1:idx2,idx3,idx4:idx5' as -m argument to indicate range of blocks to merge or \'all\' to merge all blocks where possible. IDs begin from 1.

#### Usage
./pg.py -f [maf_file] -v
python3 ./src/pg.py convert -f ./files/ebola.maf


## poa software

### Multialignment of fasta files 
./poa -read_fasta [...]/POA_article_single_block.fna  -clustal [..]/POA_article_single_block2.aln blosum80.mat











- wygenerować consensusy dla danego bloku przy użyciu programu poa, zaznaczyć je odpowiednim kolorem i powinno się pokreywac 
z grupami sekwencji w tym uliniowieniu

Aligning the unalignable: bacteriophage whole genome sequencing
 
- na ile dany consensus różni się od tej grupy, np. na ilu % pozycje się różnią


Running Tests
Required  packages: ddt:
sudo apt-get install python3-ddt