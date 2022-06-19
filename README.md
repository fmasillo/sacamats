# sacamats
Suffix sorting algorithm for set of similar strings

## Installation
```sh
git clone https://github.com/fmasillo/sacamats.git
cd sacamats
make
```

## Usage

This tool expects only one argument from command line, namely a txt file consisting of two rows. The first row should be the path to the reference sequence file and the sequence must be terminated with a dollar sign. The second row should be the path to the collection of sequences, either in FASTA format or all in one line separated by %.

Command example:
'''sh
./sacamats list_of_files.txt
'''

List of files example:
/home/Desktop/reference.fa
/home/Desktop/collection.fa

## Citation

Please, if you use this tool in an academic setting cite the following paper (to appear in WABI2022):

Zsuzsanna Lipták, Francesco Masillo, Simon J. Puglisi. Suffix sorting via matching statistics. In Proc. Workshop in Algorithm for Bioinformatics (WABI), 2022



