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
If the user wants to save the Generalized Suffix Array to (binary) file, he can specify the output file via adding the ```-o``` flag to the usual command, e.g.:

```sh
./sacamats list_of_files.txt -o outputGSA
```

## Citation

Journal paper:

Zsuzsanna Lipták, Francesco Masillo, Simon J. Puglisi. Suffix Sorting via Matching Statistics. Algorithms Mol. Biol. 2024: https://almob.biomedcentral.com/articles/10.1186/s13015-023-00245-z

```
@article{LiptakMP24,
  author       = {{\relax Zs}uzsanna Lipt{\'{a}}k and
                  Francesco Masillo and
                  Simon J. Puglisi},
  title        = {Suffix sorting via matching statistics},
  journal      = {Algorithms Mol. Biol.},
  volume       = {19},
  number       = {1},
  pages        = {11},
  year         = {2024},
}
```

Conference paper:

Zsuzsanna Lipták, Francesco Masillo, Simon J. Puglisi. Suffix Sorting via Matching Statistics. WABI 2022: https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.WABI.2022.20

```
  @inproceedings{LiptakMP22,
  author       = {{\relax Zs}uzsanna Lipt{\'{a}}k and
                  Francesco Masillo and
                  Simon J. Puglisi},
  title        = {Suffix Sorting via Matching Statistics},
  booktitle    = {22nd International Workshop on Algorithms in Bioinformatics, {WABI}
                  2022, September 5-7, 2022, Potsdam, Germany},
  series       = {LIPIcs},
  volume       = {242},
  pages        = {20:1--20:15},
  publisher    = {Schloss Dagstuhl - Leibniz-Zentrum f{\"{u}}r Informatik},
  year         = {2022},
}
```


