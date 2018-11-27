# SquiggleKit

#### A toolkit for accessing and manipulating nanopore signal data

## File management

### fast5_fetcher

A tool for fetching nanopore fast5 files to save time and simplify downstream analysis

<https://github.com/Psy-Fer/fast5_fetcher>

## Signal extraction

### SquigglePull

#### Coming changes:

-   different file inputs
-   different outputs based on file inputs
-   multiprocessing
-   integrated targeting
-   file format formalisation (happy to have feedback)

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

#### Use on HPC

When working on an HPC, here is the Sun Grid Engine script for extracting the files, assuming using array jobs for each tarball folder of fast5s. ${SGE_TASK_ID} is an integer > 1 so files should be named something like name.1.tar, name.2.tar, etc.

**spul.sge**

```bash
    source ~/venv2714/bin/activate
    BLAH=/dir/to/fast5s/name.${SGE_TASK_ID}.tar
    mkdir ${TMPDIR}/fast5

    echo "extracting..." >&2
    # extracts files wit no folder structure
    tar -xf ${BLAH} --transform='s/.*\///' -C ${TMPDIR}/fast5/
    echo "extraction complete!" >&2
    echo "Number of files:" >&2
    ls ${TMPDIR}/fast5/ | wc -l >&2

    time python /path/to/SCuigglePull.py -rv -p ${TMPDIR}/fast5/ -f all > ${TMPDIR}/signal.${SGE_TASK_ID}.tsv

    echo "copying data..." >&2

    cp ${TMPDIR}/SP_signal.${SGE_TASK_ID}.tsv /my/output/dir/

    echo "done!" >&2
```

## Signal analysis

### Segmenter

coming soon...

Stall and homopolymer identification in signal space.

### MotifSeq

coming soon...

**The Ctrl+F for signal**.
MotifSeq is used to identify any stretch of nucleotides, or other pattern, converted to signal. This is best coupled with Segmenter to achieve high accuracy targeting of regions in squiggle space.

## Acknowledgements

I would like to thank my lab (Shaun Carswell, Kirston Barton, Hasindu Gamaarachchi, Kai Martin) in Genomic Technologies team from the [Garvan Institute](https://www.garvan.org.au/) for their feedback on the development of these tools.

## License

[The MIT License](https://opensource.org/licenses/MIT)
