# SquiggleKit
####  A toolkit for accessing and manipulating nanopore signal data

## File management
### fast5_fetcher

https://github.com/Psy-Fer/fast5_fetcher

## Signal extraction
### SquigglePull

#### Coming changes:

- different file inputs
- different outputs based on file inputs
- multiprocessing
- integrated targeting
- file format formalisation (happy to have feedback)

#### Usage

    usage: SquigglePull.py [-h] [-p PATH] [-t TARGET] [-f {pos1,all}] [-r | -e]
                           [-v] [-s]

    SquigglePull - extraction of raw/event signal from Oxford Nanopore fast5 files

    optional arguments:
      -h, --help            show this help message and exit
      -p PATH, --path PATH  Top directory path of fast5 files
      -t TARGET, --target TARGET
                            Target information as comma delimited string
                            structured by format type
      -f {pos1,all}, --form {pos1,all}
                            Format of target information
      -r, --raw             Target raw signal
      -e, --event           Target event signal
      -v, --verbose         Engage higher output verbosity
      -s, --scale           Scale signal output for comparison




**All raw data:**

    python SquigglePull.py -rv -p ~/data/test/reads/1/ -f all > data.tsv

**Positional event data:**

    python SquigglePull.py -ev -p ./test/ -t 50,150 -f pos1 > data.tsv

## Signal analysis
### Tool 1
coming soon...

### MotifSeq
coming soon...
