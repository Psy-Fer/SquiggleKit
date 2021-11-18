import os
import sys
import gzip
import io
import subprocess
import traceback
import argparse
import platform
'''

    James M. Ferguson (j.ferguson@garvan.org.au)
    Genomic Technologies
    Garvan Institute
    Copyright 2018

    fast5_fetcher is designed to help manage fast5 file data storage and organisation.
    It takes 3 files as input: fastq/paf/flat, sequencing_summary, index

    --------------------------------------------------------------------------------------
    version 0.0   - initial
    version 0.2   - added argparser and buffered gz streams
    version 0.3   - added paf input
    version 0.4   - added read id flat file input
    version 0.5   - pppp print output instead of extracting
    version 0.6   - did a dumb. changed x in s to set/dic entries O(n) vs O(1)
    version 0.7   - cleaned up a bit to share and removed some hot and steamy features
    version 0.8   - Added functionality for un-tarred file structures and seq_sum only
    version 1.0   - First release
    version 1.1   - refactor with dicswitch and batch_tater updates
    version 1.1.1 - Bug fix on --transform method, added OS detection
    version 1.2.0 - Added file trimming to fully segment selection

    TODO:
        - Python 3 compatibility
        - autodetect file structures
        - autobuild index file - make it a sub script as well
        - Consider using csv.DictReader() instead of wheel building
        - options to build new index of fetched fast5s
        - Add in verbosity for stderr outputs
        - Add sam file compatibility for filtering
        - Do extra checks for \r with mac/windows files
        - Memory reduction using hashed binary indexes of larger files
        - individual file flag
        - Multi-fast5 support
        - multiprocessing

    -----------------------------------------------------------------------------
    MIT License

    Copyright (c) 2017 James Ferguson

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:
    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
'''


class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)


def main():
    '''
    Main function to control each major control function.
    '''
    parser = MyParser(
        description="fast_fetcher - extraction of specific nanopore fast5 files")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-q", "--fastq",
                       help="fastq.gz for read ids")
    group.add_argument("-p", "--paf",
                       help="paf alignment file for read ids")
    group.add_argument("-f", "--flat",
                       help="flat file of read ids")
    parser.add_argument("--OSystem", default=platform.system(),
                        help="running operating system - leave default unless doing odd stuff")
    parser.add_argument("-s", "--seq_sum",
                        help="sequencing_summary.txt.gz file")
    parser.add_argument("-i", "--index",
                        help="index.gz file mapping fast5 files in tar archives")
    parser.add_argument("-o", "--output",
                        help="output directory for extracted fast5s")
    parser.add_argument("-t", "--trim", action="store_true",
                        help="trim files as if standalone experiment, (fq, SS)")
    parser.add_argument("-l", "--trim_list",
                        help="list of file names to trim, comma separated. fastq only needed for -p and -f modes")
    parser.add_argument("-x", "--prefix", default="default",
                        help="trim file prefix, eg: barcode_01, output: barcode_01.fastq, barcode_01_seq_sum.txt")
    # parser.add_argument("-t", "--procs", type=int,
    #                    help="Number of CPUs to use - TODO: NOT YET IMPLEMENTED")
    parser.add_argument("-z", "--pppp", action="store_true",
                        help="Print out tar commands in batches for further processing")
    parser.add_argument("--move", action="store_true",
                        help="Move (mv) instead of copy (cp) files. ***DANGER: USE AT OWN RISK***")
    args = parser.parse_args()

    # print help if no arguments given
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    # TODO: hide with verbosity
    print >> sys.stderr, "Starting things up!"

    p_dic = {}
    if args.pppp:
        # TODO: hide with verbosity
        print >> sys.stderr, "PPPP state! Not extracting, exporting tar commands"


    # Handles the trimmming of fastq or sequencing summary files depending on input
    trim_pass = False
    if args.trim:
        SS = False
        FQ = False
        if args.trim_list:
            A = args.trim_list.split(',')
            for a in A:
                if "fastq" in a:
                    FQ = a
                elif "txt" in a:
                    SS = a
                else:
                    print >> sys.stderr, "Unknown trim input. detects 'fastq' or 'txt' for files. Input:", a
        else:
            print >> sys.stderr, "No extra files given. Compatible with -q fastq input only"

        if args.fastq:
            FQ = args.fastq
        if args.seq_sum:
            SS = args.seq_sum

        # final check for trimming
        if FQ and SS:
            trim_pass = True
            print >> sys.stderr, "Trim setting detected. Writing to working direcory"
        else:
            print >> sys.stderr, "Unable to verify both fastq and sequencing_summary files. Please check filenames and try again. Exiting..."
            sys.exit()

    # Do the actual trimming depending on input type - non destructive, so can do it first
    ids = []
    if args.fastq:
        ids = get_fq_reads(args.fastq)
        if trim_pass:
            trim_SS(args, ids, SS)
    elif args.paf:
        ids = get_paf_reads(args.paf)
        if trim_pass:
            trim_both(args, ids, FQ, SS)
    elif args.flat:
        ids = get_flat_reads(args.flat)
        if trim_pass:
            trim_both(args, ids, FQ, SS)
    if not ids and trim_pass:
        filenames, ids = get_filenames(args.seq_sum, ids)
        trim_both(args, ids, FQ, SS)
    else:
        filenames, ids = get_filenames(args.seq_sum, ids)

    paths = get_paths(args.index, filenames)
    # TODO: hide with verbosity
    print >> sys.stderr, "extracting..."
    # TODO: place multiprocessing pool here
    for p, f in paths:
        # if -z option, get file paths for command lists
        if args.pppp:
            if p in p_dic:
                p_dic[p].append(f)
            else:
                p_dic[p] = [f]
            continue
        else:
            try:
                extract_file(args, p, f)
            except:
                traceback.print_exc()
                print >> sys.stderr, "Failed to extract:", p, f
    # For each .tar file, write a file with the tarball name as filename.tar.txt
    # and contains a list of files to extract - input for batch_tater.py
    if args.pppp:
        # TODO: check for naming and dynamically make tater_master filename
        with open("tater_master.txt", 'w') as m:
            for i in p_dic:
                fname = "tater_" + i.split('/')[-1] + ".txt"
                m_entry = "{}\t{}".format(fname, i)
                fname = args.output + "/tater_" + i.split('/')[-1] + ".txt"
                m.write(m_entry)
                m.write('\n')
                with open(fname, 'w') as f:
                    for j in p_dic[i]:
                        f.write(j)
                        f.write('\n')
    # TODO: hide with verbosity
    print >> sys.stderr, "done!"


def dicSwitch(i):
    '''
    A switch to handle file opening and reduce duplicated code
    '''
    open_method = {
        "gz": gzip.open,
        "norm": open
    }
    return open_method[i]


def get_fq_reads(fastq):
    '''
    read fastq file and extract read ids
    quick and dirty to limit library requirements - still bullet fast
    reads 4 lines at a time and on 1st line, split and get id, dropping the @ symbol
    '''
    c = 0
    read_ids = set()
    if fastq.endswith('.gz'):
        f_read = dicSwitch('gz')
    else:
        f_read = dicSwitch('norm')
    with f_read(fastq, 'rb') as fq:
        if fastq.endswith('.gz'):
            fq = io.BufferedReader(fq)
        for line in fq:
            c += 1
            line = line.strip('\n')
            if c == 1:
                idx = line.split()[0][1:]
                read_ids.add(idx)
            elif c >= 4:
                c = 0
    return read_ids


def get_paf_reads(reads):
    '''
    Parse paf file to pull read ids (from minimap2 alignment)
    First field is the readID
    '''
    read_ids = set()
    if reads.endswith('.gz'):
        f_read = dicSwitch('gz')
    else:
        f_read = dicSwitch('norm')
    with f_read(reads, 'rb') as fq:
        if reads.endswith('.gz'):
            fq = io.BufferedReader(fq)
        for line in fq:
            line = line.strip('\n')
            line = line.split()
            read_ids.add(line[0])
    return read_ids


def get_flat_reads(filename):
    '''
    Parse a flat file separated by line breaks \n
    ie, one readID per line
    TODO: make @ symbol check once, as they should all be the same
    '''
    read_ids = set()
    check = True
    if filename.endswith('.gz'):
        f_read = dicSwitch('gz')
    else:
        f_read = dicSwitch('norm')
    with f_read(filename, 'rb') as fq:
        if filename.endswith('.gz'):
            fq = io.BufferedReader(fq)
        for line in fq:
            line = line.strip('\n')
            if check:
                if line[0] == '@':
                    x = 1
                else:
                    x = 0
                check = False
            idx = line[x:]
            read_ids.add(idx)
    return read_ids


def trim_SS(args, ids, SS):
    '''
    Trims the sequencing_summary.txt file to only the input IDs
    '''
    if args.prefix:
        pre = args.prefix + "_seq_sum.txt"
    else:
        pre = "trimmed_seq_sum.txt"
    head = True
    if SS.endswith('.gz'):
        f_read = dicSwitch('gz')
    else:
        f_read = dicSwitch('norm')
    # make this compatible with dicSwitch
    with open(pre, "w") as w:
        with f_read(SS, 'rb') as sz:
            if SS.endswith('.gz'):
                sz = io.BufferedReader(sz)
            for line in sz:
                if head:
                    w.write(line)
                    head = False
                    continue
                l = line.split()
                if l[1] in ids:
                    w.write(line)


def trim_both(args, ids, FQ, SS):
    '''
    Trims the sequencing_summary.txt and fastq files to only the input IDs
    '''
    # trim the SS
    trim_SS(args, ids, SS)
    if args.prefix:
        pre = args.prefix + ".fastq"
    else:
        pre = "trimmed.fastq"

    # trim the fastq
    c = 0
    P = False
    if FQ.endswith('.gz'):
        f_read = dicSwitch('gz')
    else:
        f_read = dicSwitch('norm')
    with open(pre, "w") as w:
        with f_read(FQ, 'rb') as fq:
            if FQ.endswith('.gz'):
                fq = io.BufferedReader(fq)
            for line in fq:
                c += 1
                if c == 1:
                    if line.split()[0][1:] in ids:
                        P = True
                        w.write(line)
                elif P and c < 4:
                    w.write(line)
                elif c >= 4:
                    if P:
                        w.write(line)
                    c = 0
                    P = False


def get_filenames(seq_sum, ids):
    '''
    match read ids with seq_sum to pull filenames
    uses set to remove possible duplicates, as well as faster "if in" checks
    sets are hashed, lists are not, so O(N) vs O(1) complexity
    '''
    # for when using seq_sum for filtering, and not fq,paf,flat
    ss_only = False
    if not ids:
        ss_only = True
        ids = set()
    head = True
    files = set()
    if seq_sum.endswith('.gz'):
        f_read = dicSwitch('gz')
    else:
        f_read = dicSwitch('norm')
    with f_read(seq_sum, 'rb') as sz:
        if seq_sum.endswith('.gz'):
            sz = io.BufferedReader(sz)
        for line in sz:
            if head:
                head = False
                continue
            line = line.strip('\n')
            line = line.split()
            if ss_only:
                files.add(line[0])
                ids.add(line[1])
            else:
                if line[1] in ids:
                    files.add(line[0])
    return files, ids


def get_paths(index_file, filenames, f5=None):
    '''
    Read index and extract full paths for file extraction
    This could be done better with byte indexing
    '''
    tar = False
    paths = []
    c = 0
    if index_file.endswith('.gz'):
        f_read = dicSwitch('gz')
    else:
        f_read = dicSwitch('norm')
    # detect normal or tars
    with f_read(index_file, 'rb') as idz:
        if index_file.endswith('.gz'):
            idz = io.BufferedReader(idz)
        for line in idz:
            line = line.strip('\n')
            c += 1
            if c > 10:
                break
            if line.endswith('.tar'):
                tar = True
                break
    # extract paths
    with f_read(index_file, 'rb') as idz:
        if index_file.endswith('.gz'):
            idz = io.BufferedReader(idz)
        for line in idz:
            line = line.strip('\n')
            if tar:
                if line.endswith('.tar'):
                    path = line
                elif line.endswith('.fast5'):
                    f = line.split('/')[-1]
                    if f in filenames:
                        paths.append([path, line])
                else:
                    continue
            else:
                if line.endswith('.fast5'):
                    f = line.split('/')[-1]
                    if f in filenames:
                        paths.append(['', line])
                else:
                    continue

    return paths


def extract_file(args, path, filename):
    '''
    Do the extraction.
    I was using the tarfile python lib, but honestly, it sucks and was too volatile.
    if you have a better suggestion, let me know :)
    That --transform hack is awesome btw. Blows away all the leading folders. use
    cp for when using untarred structures. Not recommended, but here for completeness.

    --transform not working on MacOS. Need to use gtar
    Thanks to Kai Martin for picking that one up!

    TODO: move OS detection to top, and use to help compatibility for edge cases.
          works pretty well on all systems so far, but all about dat support!
    '''
    OSystem = ""
    OSystem = args.OSystem
    save_path = args.output
    if path.endswith('.tar'):
        if OSystem in ["Linux", "Windows"]:
            cmd = "tar -xf {} --transform='s/.*\///' -C {} {}".format(
                path, save_path, filename)
        elif OSystem == "Darwin":
            cmd = "gtar -xf {} --transform='s/.*\///' -C {} {}".format(
                path, save_path, filename)
        else:
            print >> sys.stderr, "Unsupported OSystem, trying Tar anyway, OS:", OSystem
            cmd = "tar -xf {} --transform='s/.*\///' -C {} {}".format(
                path, save_path, filename)
    else:
        if args.move:
            cmd = "mv {} {}".format(filename, save_path)
        else:
            cmd = "cp {} {}".format(filename, save_path)

        subprocess.call(cmd, shell=True, executable='/bin/bash')


if __name__ == '__main__':
    main()
