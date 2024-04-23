# sacamats
Algorithm for constructing the generalized suffix array of a collection of highly similar strings.

## Installation

```sh
git clone https://github.com/fmasillo/sacamats.git
cd sacamats
make
```

## Usage

This tool expects only one argument from command line, namely a txt file consisting of two rows. The first row should be the path to the reference sequence file. The second row should be the path to the collection of sequences in FASTA format.

By calling the executable with ```-h``` a help message is displayed.

```sh
./sacamats -h

Usage: ./sacamats [options] <input filename>
<input filename> is the name of the file containing paths to the reference sequence (in the first line) and to the collection file (in the second line).
  Options: 
        -p      read only a prefix of the file expressed in number of characters, def. whole file
        -t      number of threads to use, def. max available
        -o      basename for the output files, def. <input filename>
        -h      prints this help
```

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

Journal paper:

Zsuzsanna Lipták, Francesco Masillo, Simon J. Puglisi. Suffix Sorting via Matching Statistics. Algorithms Mol. Biol. 19(1):11, 2024: https://almob.biomedcentral.com/articles/10.1186/s13015-023-00245-z

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

Zsuzsanna Lipt´ak, Francesco Masillo, and Simon J. Puglisi. Suffix sorting via matching statistics. In Proc. of the 22nd International Workshop on Algorithms in Bioinformatics (WABI 2022), volume 242 of LIPIcs, pages 20:1–20:15. Schloss Dagstuhl - Leibniz-Zentrum f¨ur Informatik, 2022: https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.WABI.2022.20

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


