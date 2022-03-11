# header_renamer
A Python script for renaming fasta file headers using regex pattern matching.

## Authors
[CSynodinos](https://github.com/CSynodinos)

## Installation

```bash
  curl https://raw.githubusercontent.com/CSynodinos/rename-fasta-headers/master/header_renamer.py -o header_renamer.py
```

## Dependencies
    1) Biopython
    2) Pandas 

This fasta header renaming script works by looking for a pattern in every header of the .fasta file.
Once found, it switches the id and description of that header according to a specified id and description
respectively. The pattern to find, new id and new description are specified with a .csv file that has the following
structure:
```bash
    pattern,new_id,new_description
    header,foo,bar
    header2,foo2,bar2
                                                    .   .   .
```

## Example
```bash
    >>> python3 header_renamer.py -i yourfasta.fasta -cv patterns.csv  
```

For more information regarding all the options available:
```bash
    >>> python3 header_renamer.py -h
```

## Features

- [x] Input file format to rename:
    - [x] .fasta

- [x] Input pattern file format:
    - [x] .csv
