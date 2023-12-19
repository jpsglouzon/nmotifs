# n-motifs model #

**Represent RNA secondary structures** of arbitrary size uncovers **structural patterns** that can provide a better understanding of **RNA functions**. The n-motifs representation of RNA secondary structure computes a vector representation of secondary structures.

* Download **[here](https://github.com/jpsglouzon/nmotifs/releases)** the ubuntu 64bit executable and follow the 'how to use it' instructions to compare secondary structures. 

* Or follow the instructions below to compile and use the n-motifs program on almost every computer platform.

## How to compile the n-motifs program ? ##

* Download the source code **[here](https://github.com/jpsglouzon/nmotifs/zipball/master)** and unzip.

* Open the terminal, `cd path_to_supernmotifs_program` to access the super-n-motifs program folder then compile it by running the command `make`.

* The executable file named 'supernmotifs' can be found in `path_to_supernmotifs_program`.

## How to use it? ##

* **Computer n-motif representation** by calling: 
```
/path_to_supernmotifs_program/nmotifs -i fileInDb -o folderOfResults
```

* The Super-n-motifs program takes as input a file of **rna secondary structures** in **dot-bracket** format (-i fileInDb):
```
>RNA1
GCCCCGCUGAUGAGGUCAGGGAAAACCGAAAGUGUCGACUCUACGGGGC
((((((.......((((......))))...((((....)))).))))))
```
It ouputs a **n-motifs matrix** and **n-motifs position matrix** by default (-o folderOfResults) .
For further options: `./path_to_supernmotifs_program/nmotifs -h`

* The n-motifs model supports circular RNA, pseudoknots and g-quadruplexes. **Circular RNAs** require to add `c_` to the header of each RNA :  `>c_RNA1`. Base pairs involved in **pseudoknots** are typically represented by special characters `{}`, `<>`, `[]`, and alphabets such as `Aa` or `Bb` : `..AA..aa`. Interacting guanines in **g-quadruplexes** are represented by `+` : `..((..++.++.++..++.)).`.

## How to cite the n-motifs model ? ##

This n-motifs is part of the super-n-motif model.


Glouzon JS, Perreault JP, Wang S. The super-n-motifs model: a novel
alignment-free approach for representing and comparing RNA secondary structures. Bioinformatics. 2017 Jan 14. pii: btw773. doi: 10.1093/bioinformatics/btw773

## Licence ##

The n-motifs model is released under the terms of the GNU GPL licence. For further informations, please see the LICENCE file of the repository.


