# sacamats
Algorithm for constructing the generalized suffix array of a collection of highly similar strings.

## Installation

```sh
git clone https://github.com/fmasillo/sacamats.git
git submodule update --init --recursive
cd sacamats
make
```

## Usage

This tool expects only one argument from command line, namely a txt file consisting of two rows. The first row should be the path to the reference sequence file and the sequence must be terminated with a dollar sign. The second row should be the path to the collection of sequences, either in FASTA format or all in one line separated by %.

Command example:
```sh
./sacamats list_of_files.txt
```

List of files example:
```
/data/reference.fa
/data/collection.fa
```

## Citation

Please, if you use this tool in an academic setting cite the following paper (to appear in WABI2022):

```
  @inproceedings{LiptakMP2022,
    author    = {Zsuzsanna Lipták and Francesco Masillo and Simon J. Puglisi},
    title     = {Suffix sorting via matching statistics},
    booktitle = {22nd International Workshop on Algorithms in Bioinformatics, {WABI}
                 2022, September 5-7, 2022, Potsdam, Germany},
    year      = {2022},
  }
```


