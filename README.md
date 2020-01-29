# SquiggleKit

#### A toolkit for accessing and manipulating nanopore signal data

**Full documentation:** https://psy-fer.github.io/SquiggleKitDocs/

**publication:** [SquiggleKit: A toolkit for manipulating nanopore signal data](https://doi.org/10.1093/bioinformatics/btz586)

**Pre-print:** [SquiggleKit: A toolkit for manipulating nanopore signal data](https://doi.org/10.1101/549741)


## Coming changes

**Fast5_fetcher:** merge single files into multi-fast5 files

**SquigglePull:** python3, read from multi-fast5

**SquigglePlot:** python3, read from multi-fast5, image size args, arg clean-up

**Segmenter:** dynamic file formats and more stability

**MotifSeq:** Improved background modelling, custom modelling, RNA specific tools, custom alignment methods

## Overview

| **Tool**      | **Category**  | **Description** |
| ------------- | ------------- | --------------- |
| [Fast5_fetcher](#fast5_fetcher) | `File management`   | Fetches fast5 files given a filtered input list |
| [SquigglePull](#squigglepull)  | `Signal extraction` | Extracts event or raw signal from data files    |
| [SquigglePlot](#squiggleplot)  | `Signal visualisation` | Visualisation tool for signal data       |
| [Segmenter](#segmenter)     | `Signal analysis`   | Finds adapter stall, and homopolymer regions    |
| [MotifSeq](#motifseq)      | `Signal analysis`   | Finds nucleotide sequence motifs in signal, i.e.“Ctrl+F” |


<img src="img/kitpic.png" width="700">


## Requirements

Following a self imposed guideline, most things written to handle nanopore data or bioinformatics in general, will use as little 3rd party libraries as possible, aiming for only core libraries, or have all included files in the package.

In the case of `fast5_fetcher.py` and `batch_tater.py`, only core python libraries are used. So as long as **Python 2.7+** is present, everything should work with no extra steps.

There is one catch. Everything is written primarily for use with **Linux**. Due to **MacOS** running on Unix, so long as the GNU tools are installed (see below), there should be minimal issues running it. **Windows** however may require more massaging. The Windows-Subsystem-Linux must be installed. Follow the instructions [here](https://itsfoss.com/install-bash-on-windows/) to do this.

SquiggleKit tools were not made to be executable to allow for use with varying python environments on various operating systems. To make them executable, add `#!` paths, such as `#!/usr/bin/env python2.7` as the first line of each of the files, then add the SquiggleKit directory to the PATH variable in `~/.bashrc`,  `export PATH="$HOME/path/to/SquiggleKit:$PATH"`

### Install

    git clone https://github.com/Psy-Fer/SquiggleKit.git

Use pip for python 2 and pip3 for python 3. User environments may vary.

for
`fast5_fetcher.py`, `SquigglePull.py`, `segmenter.py`:

- numpy
- matplotlib
- h5py
- sklearn
- ont_fast5_api

```
pip install numpy h5py sklearn matplotlib
```

for `MotifSeq.py`:

- all of the above
- scipy
- scrappie
- mlpy 3.5.0 (only use pip3 in python 3 - see below)

```
pip install scipy scrappie
```


#### Installing mlpy:


##### Python2

- Download the [Files](https://sourceforge.net/projects/mlpy/files/)
- Install [Instructions](http://mlpy.sourceforge.net/docs/3.5/install.html)

##### Python3.6 (Python3.7 not working because of a cython/gsl issue)

```
pip3 install machine-learning-py
```



## Quick start

#### fast5_fetcher

If using MacOS, and NOT using homebrew, install it here:

[**homebrew installation instructions**](https://brew.sh/)

then install gnu-tar with:

    brew install gnu-tar

#### Building the index (not required for multi_fast5)

How the index is built depends on which file structure you are using. It will work with both tarred and un-tarred file structures. Tarred is preferred. (zip and other archive methods are being investigated)

##### - Raw structure (not preferred)

```bash
for file in $(pwd)/reads/*/*;do echo $file; done >> name.index

gzip name.index
```

##### - Local basecalled structure

```bash
for file in $(pwd)/reads.tar; do echo $file; tar -tf $file; done >> name.index

gzip name.index
```

##### - Parallel basecalled structure

```bash
for file in $(pwd)/fast5/*fast5.tar; do echo $file; tar -tf $file; done >> name.index
```

If you have multiple experiments, then cat them all together and gzip.

```bash
for file in ./*.index; do cat $file; done >> ../all.name.index

gzip all.name.index
```

##### Basic use on a local computer

**using a filtered paf file as input:**

```bash
python fast5_fetcher.py -p my.paf -s sequencing_summary.txt.gz -i name.index.gz -o ./fast5
```


### SquigglePull


**All raw data:**

    python SquigglePull.py -rv -p ~/data/test/reads/1/ -f all > data.tsv

**Positional event data:**

    python SquigglePull.py -ev -p ./test/ -t 50,150 -f pos1 > data.tsv

### SquigglePlot

**Plot individual fast5 file:**

    python SquigglePlot.py -i ~/data/test.fast5

**Plot files in path**

    python SquigglePlot.py -p ~/data/ --plot_colour -g

**Plot first 2000 data points of each read from signal file and save at 300dpi pdf:**

    python SquigglePlot.py -s signals.tsv.gz --plot_colour teal -n 2000 --dpi 300 --no_show o--save test.pdf --save_path ./test/plots/

### Segmenter

**Identify any segments in folder and visualise each one**

Use `f` to full screen a plot, and `ctrl+w` to close a plot and move to the next one.

    python segmenter.py -p ./test/ -v

**Stall identification**

    python segmenter.py -s signals.tsv.gz -ku -j 100 > signals_stall_segments.tsv


### MotifSeq

**Find kmer motif:**



fasta format for model:

>my_kmer.fa

    >my_kmer_name
    ATCGATCGCTATGCTAGCATTACG

find the best match to that k-mer in the signal:

    python MotifSeq.py -s signals.tsv -i my_kmer.fa > signals_kmer.tsv

<p align="left"><img src="img/MotifSeq_fig.jpg" alt="MotifSeq" width="100%" height="100%"></p>

## Limitations

k-mer length should not really be below 12nt, below this things get hairy based on modelling

The p-values and hit probabilities provided are based on loose modelling of negative background scores for a number of k-mers. It is currently only modelled on R9.4 model, not R10 or RNA.


## Acknowledgements

I would like to thank the members of my lab, Shaun Carswell, Kirston Barton, Hasindu Gamaarachchi, Kai Martin, Tansel Ersavas, Brent Bevear, Jillian Hammond, and Martin Smith, from the Genomic Technologies team from the [Garvan Institute](https://www.garvan.org.au/) for their feedback on the development of these tools.

## License

[The MIT License](https://opensource.org/licenses/MIT)
