import os
import sys
import gzip
import subprocess
import traceback
import argparse
import platform
import shutil
from ont_fast5_api.conversion_tools import single_to_multi_fast5
from ont_fast5_api.conversion_tools.single_to_multi_fast5 import create_multi_read_file
# batch_conver_single_to_multi(input_path, output_folder, filename_base, batch_size, threads, recursive, follow_symlinks, target_compression)
from ont_fast5_api.multi_fast5 import MultiFast5File
from ont_fast5_api.multi_fast5 import Fast5File
for i in sys.argv:
    if i == '-m':
        from ont_fast5_api.conversion_tools import multi_to_single_fast5
        from ont_fast5_api.conversion_tools.multi_to_single_fast5 import create_single_f5
        from ont_fast5_api.multi_fast5 import MultiFast5File
'''

    James M. Ferguson (j.ferguson@garvan.org.au)
    Genomic Technologies
    Garvan Institute
    Copyright 2018

    fast5_fetcher is designed to help manage fast5 file data storage and organisation.
    It takes as input: fastq/paf/flat, sequencing_summary, (index)

    ----------------------------------------------------------------------------
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
    version 1.2.1 - Added to SquiggleKit - Adopting its versioning
    version 1.3.0 - Added Multifast5 and python3 support
    ----------------------------------------------------------------------------

    TODO:
        - autodetect file structures
        - autobuild index file - make it a sub script as well - multifast5 negates this
        - Consider using csv.DictReader() instead of wheel building
        - options to build new index of fetched fast5s
        - Add sam file compatibility for filtering - flat file minus headers
        - Do extra checks for \r with mac/windows files
        - Memory reduction using hashed binary indexes of larger files
        - individual file flag
        - multiprocessing
        - header using # with autodetect for metadata
        - reslease a 2.0.0 with major re-write and cleanup with multiprocessing

    ----------------------------------------------------------------------------
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
    VERSION = "1.3.0"

    parser = MyParser(
        description="fast_fetcher - extraction of specific nanopore fast5 files")
    inputs = parser.add_mutually_exclusive_group()
    index = parser.add_mutually_exclusive_group()
    trim = parser.add_mutually_exclusive_group()
    inputs.add_argument("-q", "--fastq",
                       help="fastq for read ids")
    inputs.add_argument("-p", "--paf",
                       help="paf alignment file for read ids")
    inputs.add_argument("-f", "--flat",
                       help="flat file of read ids")
    parser.add_argument("-s", "--seq_sum", required=True,
                        help="sequencing_summary.txt.gz file")
    index.add_argument("-i", "--index",
                        help="index.gz file mapping fast5 files in tar archives")
    index.add_argument("-m", "--multi_f5",
                        help="path to multi-fast5 files")
    parser.add_argument("-c", "--f5_format", default="multi", choices=["multi", "single"],
                        help="output fast5 file format")
    parser.add_argument("--threshold", default=4000, type=int,
                        help="threshold number for amount of reads in a single multifast5 output file")
    parser.add_argument("-o", "--output", required=True,
                        help="output directory for extracted fast5s")
    trim.add_argument("-t", "--trim", action="store_true",
                        help="trim files as if standalone experiment, (fq, SS)")
    parser.add_argument("-l", "--trim_list",
                        help="list of file names to trim, comma separated. fastq only needed for -p and -f modes")
    parser.add_argument("-x", "--prefix", default="trimmed",
                        help="trim file prefix, eg: barcode_01, output: barcode_01.fastq, barcode_01_seq_sum.txt")
    trim.add_argument("-D", "--seq_sum_1D2", choices=["first", "second", "both"],
                        help="Sequencing summary file is from 1D^2 basecalling")
    # parser.add_argument("-P", "--procs", type=int,
    #                     help="Number of CPUs to use - Only available for Multi-fast5 files")
    parser.add_argument("-z", "--pppp", action="store_true",
                        help="Print out tar commands in batches for further processing - advanced use/batch_tater.py")
    parser.add_argument("--OSystem", default=platform.system(),
                        help="running operating system - leave default unless doing odd stuff")
    parser.add_argument("-V", "--version", action="store_true",
                        help="Print version information")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="engage higher level of verbosity for troubleshooting")
    args = parser.parse_args()
    # --------------------------------------------------------------------------
    # print help if no arguments given
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    # print metadata
    if args.version:
        sys.stderr.write("SquiggleKit fast5_fetcher: {}\n".format(VERSION))
        sys.exit(1)

    if args.verbose:
        sys.stderr.write("Verbose mode active - dumping info to stderr\n")
        sys.stderr.write("SquiggleKit fast5_fetcher: {}\n".format(VERSION))
        sys.stderr.write("args: {}\n".format(args))

    if args.multi_f5:
        if args.verbose:
            sys.stderr.write("Multi-fast5 mode detected in mode: {}\n".format(args.f5_format))

    if args.seq_sum_1D2:
        if args.verbose:
            sys.stderr.write("1D^2 sequencing summary mode on. Getting {} read(s)\n".format(args.seq_sum_1D2))


    # --------------------------------------------------------------------------
    # do checks on arguments
    # check inputs, index, and output are present (no need for index if multi_f5)
    # check it can find all the files and paths before continuing
    if not args.fastq and not args.flat and not args.paf:
        if not args.seq_sum:
            sys.stderr.write("\nNo input detected.\nPlease indicate fastq/paf/flat, or sequencing_summary file fast5 selection.\n\n\n")
            parser.print_help(sys.stderr)
            sys.exit(1)
        else:
            sys.stderr.write("\nSequencing_summary file only mode detected. Extracting all files found in sequencing_summary file.\n\n\n")
    if not args.output:
        sys.stderr.write("\nNo output detected.\nPlease indicate a path to place the extracted fast5 files.\n\n\n")
        parser.print_help(sys.stderr)
        sys.exit(1)
    else:
        if not os.path.isdir(args.output):
            os.mkdir(args.output)
            if args.verbose:
                sys.stderr.write("Output folder '{}' created\n".format(args.output))
    # --------------------------------------------------------------------------

    if args.verbose:
        sys.stderr.write("Checks passed!\nStarting things up!\n")

    p_dic = {}
    if args.pppp and args.verbose:
        sys.stderr.write("PPPP state! Not extracting, exporting tar commands.\n")
    if args.pppp and args.multi_f5:
        sys.stderr.write("PPPP state and multi_f5 not yet supported.\n")
        sys.exit(1)

    # --------------------------------------------------------------------------
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
                    sys.stderr.write("Unknown trim input. detects 'fastq' or 'txt' for files. Input: {}\n".format(a))
        else:
            sys.stderr.write("No extra files given. Compatible with -q fastq input only.\n")

        if args.fastq:
            FQ = args.fastq
        if args.seq_sum:
            SS = args.seq_sum

        # final check for trimming
        if FQ and SS:
            trim_pass = True
            sys.stderr.write("Trim setting detected. Writing to working directory.\n")
        else:
            sys.stderr.write("Unable to verify both fastq and sequencing_summary files. Please check filenames and try again. Exiting...\n")
            sys.exit(1)

    # --------------------------------------------------------------------------
    # Do the actual trimming depending on input type - non destructive, so can do it first
    # get the ids
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
    # --------------------------------------------------------------------------
    # get fast5 file names using seq_sum
    if not ids and trim_pass:
        if args.multi_f5:
            filenames, ids = get_filenames_multi_f5(args.seq_sum, ids)
            trim_both(args, ids, FQ, SS)
        else:
            filenames, ids = get_filenames(args.seq_sum, ids)
            trim_both(args, ids, FQ, SS)
    else:
        if args.multi_f5:
            if args.verbose:
                sys.stderr.write("Getting multi-fast5 info...\n")
            if args.seq_sum_1D2:
                #hardcoded for both right now
                filenames, ids = get_1D2_filenames_multi_f5(args.seq_sum, ids)
            else:
                filenames, ids = get_filenames_multi_f5(args.seq_sum, ids)
        else:
            filenames, ids = get_filenames(args.seq_sum, ids)
    # --------------------------------------------------------------------------

    if not filenames:
        sys.stderr.write("no filenames list built, check inputs\n")

    if not ids:
        sys.stderr.write("No ids list built, check inputs\n")
    # --------------------------------------------------------------------------
    if args.multi_f5:
        if args.verbose:
            sys.stderr.write("Extracting reads from multi-fast5 files...\n")
        m_paths = get_paths_multi_f5(args.multi_f5, filenames)
        if not m_paths:
            sys.stderr.write("No file paths built\n")
        multi_f5_handler(args, m_paths, filenames)

    else:
        if args.verbose:
            sys.stderr.write("Doing legacy extraction...\n")
        paths = get_paths(args.index, filenames)

        # TODO: place multiprocessing pool here
        # if converting single to multi, initialise needed variables
        if args.f5_format == 'multi':
            count = 0
            i = 0
            outfile = "multi_"+str(i)
            os.mkdir("fast5_fetcher_temp")
            tmp_path = os.path.abspath("fast5_fetcher_temp")
            save_path = args.output
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
                    sys.stderr.write("Failed to extract: {} {}\n".format(p, f))
                else:
                    if args.f5_format == 'multi':
                        # increase count if file is sucessfully extracted and add details to mapping file
                        count += 1
                        with open(os.path.join(save_path, "filename_mapping.txt"), 'a') as out_sum:
                            out_sum.write("{}\t{}\n".format(outfile, f))
                        # if threshold is reached, all files in folder are converted into single multi fast5 file and folder and count are reset
                        if count == args.threshold:
                            s2m(tmp_path, save_path, outfile, None)
                            shutil.rmtree(tmp_path)
                            i += 1
                            outfile = "multi_"+str(i)
                            os.mkdir("fast5_fetcher_temp")
                            tmp_path = os.path.abspath("fast5_fetcher_temp")
                            count = 0
        # if there are any single fast5 left, convert them into a multi fast5 file
        if args.f5_format == 'multi':
            if count != 0:
                s2m(tmp_path, save_path, outfile, None)
            shutil.rmtree(tmp_path)

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
    if args.verbose:
        # do a check some files were actually extracted
        sys.stderr.write("done!\n\nScript end.\n")


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
    with f_read(fastq, 'rt') as fq:
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
    with f_read(reads, 'rt') as fq:
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
    with f_read(filename, 'rt') as fq:
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
        with f_read(SS, 'rt') as sz:
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
        with f_read(FQ, 'rt') as fq:
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
    with f_read(seq_sum, 'rt') as sz:
        for line in sz:
            if head:
                head = False
                continue
            line = line.strip('\n')
            line = line.split()
            # add 1D^2 logic here
            if ss_only:
                files.add(line[0])
                ids.add(line[1])
            else:
                if line[1] in ids:
                    files.add(line[0])
    return files, ids

def get_filenames_multi_f5(seq_sum, ids):
    '''
    Match read ids with seq_sum to pull file names when packed in multifast5 format
    '''
    # for when using seq_sum for filtering, and not fq,paf,flat
    ss_only = False
    if not ids:
        ss_only = True
        ids = set()
    head = True
    files = {}
    f5_idx=1
    r_idx=2
    if seq_sum.endswith('.gz'):
        f_read = dicSwitch('gz')
    else:
        f_read = dicSwitch('norm')
    with f_read(seq_sum, 'rt') as sz:
        for line in sz:
            line = line.strip('\n')
            line = line.split()
            if head:
                head = False
                for i in range(len(line)):
                    if line[i] == "filename_fast5":
                        f5_idx=i
                    elif line[i] == "read_id":
                        r_idx=i
                continue
            # add 1D^2 logic here
            if ss_only:
                if line[f5_idx] not in files:
                    files[line[f5_idx]] = []
                files[line[f5_idx]].append(line[r_idx])
                ids.add(line[r_idx])
            else:
                if line:
                    if line[2] in ids:
                        if line[f5_idx] not in files:
                            files[line[f5_idx]] = []
                        files[line[f5_idx]].append(line[r_idx])
    return files, ids


def get_1D2_filenames_multi_f5(seq_sum, ids):
    '''
    Match read ids with seq_sum to pull file names when packed in multifast5 format
    '''
    # for when using seq_sum for filtering, and not fq,paf,flat
    ss_only = False
    if not ids:
        ss_only = True
        ids = set()
    head = True
    files = {}
    if seq_sum.endswith('.gz'):
        f_read = dicSwitch('gz')
    else:
        f_read = dicSwitch('norm')
    with f_read(seq_sum, 'rt') as sz:
        for line in sz:
            if head:
                head = False
                continue
            line = line.strip('\n')
            line = line.split('\t')
            if ss_only:
                if line[0] not in files:
                    files[line[0]] = []
                files[line[0]].append(line[2])
                ids.add(line[2])
                if line[1] not in files:
                    files[line[1]] = []
                files[line[1]].append(line[3])
                ids.add(line[3])
            else:
                if line[2] in ids:
                    if line[0] not in files:
                        files[line[0]] = []
                    files[line[0]].append(line[2])
                    if line[1] not in files:
                        files[line[1]] = []
                    files[line[1]].append(line[3])
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
    with f_read(index_file, 'rt') as idz:
        for line in idz:
            line = line.strip('\n')
            c += 1
            if c > 10:
                break
            if line.endswith('.tar'):
                tar = True
                break
    # extract paths
    with f_read(index_file, 'rt') as idz:
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


def get_paths_multi_f5(path, filenames):
    '''
    Build the path list for the multi_fast5 files given a top dir, recursive
    '''
    filepaths = []
    for dirpath, dirnames, files in os.walk(path):
        for fast5 in files:
            if fast5.endswith('.fast5'):
                if fast5 in filenames:
                    filepaths.append([os.path.join(dirpath, fast5), fast5])
    return filepaths


def convert_multi_to_single(input_file, read_list, output_folder):
    '''
    Pull the exact read out of the file.
    '''
    results = [os.path.basename(input_file)]
    try:
        with MultiFast5File(input_file, 'r') as multi_f5:
            read_ids =  set(multi_f5.get_read_ids())
            for query_read in read_list:
                if query_read in read_ids:
                    try:
                        read = multi_f5.get_read(query_read)
                        output_file = os.path.join(output_folder, "{}.fast5".format(query_read))
                        create_single_f5(output_file, read)
                        results.append(os.path.basename(output_file))
                    except:
                        traceback.print_exc()
                        sys.stderr.write("{}\n\tFailed to copy read '{}' from {}\n".format("convert_multi_to_single", query_read, input_file))
                else:
                    sys.stderr.write("{}\n\tFailed to find read '{}' in {}\n".format("convert_multi_to_single", query_read, input_file))
    except:
        traceback.print_exc()
        sys.stderr.write("{}\n\tFailed to copy files from: {}\n".format("convert_multi_to_single", input_file))
    finally:
        return results


def m2s(f5_path, read_list, save_path):
    '''
    Open a single multi-fast5 file and extract 1 single read
    TODO: sort out the file mapping file writing
    '''
    result = convert_multi_to_single(f5_path, read_list, save_path)
    with open(os.path.join(save_path, "filename_mapping.txt"), 'a') as out_sum:
        M = result[0]
        for S in result[1:]:
           out_sum.write("{}\t{}\n".format(M,S))

def s2m(f5_path, save_path, output_file, target_compression):
    '''
    Combine single fast5 files in the f5_path dir into 1 multi fast5 file, saved as output_file in the save_path
    '''
    filenames = []
    results = []
    output = None
    output_file = output_file + ".fast5"
    if target_compression == "gzip":
        output_file = output_file + ".gz"

    for dirpath, dirnames, files in os.walk(f5_path):
        for fast5 in files:
            if fast5.endswith('.fast5'):
                filenames.append(os.path.join(dirpath, fast5))
    if filenames:
        results, output = create_multi_read_file(filenames, os.path.join(save_path, output_file), target_compression)

def multi_f5_handler(args, m_paths, filenames):
    '''
    handle multi-fast5 Files
    '''
    save_path = args.output

    if args.f5_format == "multi":
        threshold = args.threshold
        count = 0
        index = 0
        i = 0
        outfile = "multi_"+str(i)
        os.mkdir("fast5_fetcher_temp")
        tmp_path = os.path.abspath("fast5_fetcher_temp")
        for p, f in m_paths:
            readIDs = filenames[f]
            length = len(readIDs)
            count += length
            if count >= threshold:
                # add reads to reach the exact threshold
                index = length - (count - threshold)
                convert_multi_to_single(p, readIDs[:index], tmp_path)
                # convert the single files to a single multi
                s2m(tmp_path, save_path, outfile, None)
                # write to the mapping file
                with open(os.path.join(save_path, "filename_mapping.txt"), 'a') as out_sum:
                    for ID in readIDs[:index]:
                        out_sum.write("{}\t{}\n".format(outfile, ID))
                del readIDs[:index]
                shutil.rmtree(tmp_path)

                # if the number of reads left is enough to create another complete file, do so until you cannot
                while len(readIDs) >= threshold:
                    # create a new temp folder
                    os.mkdir("fast5_fetcher_temp")
                    tmp_path = os.path.abspath("fast5_fetcher_temp")
                    # add reads to new folder and create a new multi file from it
                    i += 1
                    outfile = "multi_"+str(i)
                    convert_multi_to_single(p, readIDs[:threshold], tmp_path)
                    s2m(tmp_path, save_path, outfile, None)
                    with open(os.path.join(save_path, "filename_mapping.txt"), 'a') as out_sum:
                        for ID in readIDs[:threshold]:
                            out_sum.write("{}\t{}\n".format(outfile, ID))
                    del readIDs[:threshold]
                    #delete temp folder
                    shutil.rmtree(tmp_path)

                # excess reads that are not enough to create a new file go in a new temp folder
                os.mkdir("fast5_fetcher_temp")
                tmp_path = os.path.abspath("fast5_fetcher_temp")
                i += 1
                outfile = "multi_"+str(i)
                # add remaining reads, if any, to a new folder
                if readIDs:
                    convert_multi_to_single(p, readIDs, tmp_path)
                    with open(os.path.join(save_path, "filename_mapping.txt"), 'a') as out_sum:
                        for ID in readIDs:
                            out_sum.write("{}\t{}\n".format(outfile, ID))
                count = len(readIDs)
            else:
                # add reads to folder
                convert_multi_to_single(p, readIDs, tmp_path)
                # add details to mapping file
                with open(os.path.join(save_path, "filename_mapping.txt"), 'a') as out_sum:
                    for ID in readIDs:
                        out_sum.write("{}\t{}\n".format(outfile, ID))

        # check if there are remaining files
        if os.listdir(tmp_path):
            # convert remaining files to multi
            s2m(tmp_path, save_path, outfile, None)
        # remove temp folder
        shutil.rmtree(tmp_path)

    elif args.f5_format == "single":
        with open(os.path.join(args.output, "filename_mapping.txt"), 'w') as out_sum:
            out_sum.write("multi_read_file\tsingle_read_file\n")
        for p, f in m_paths:
            try:
                read_list = filenames[f]
                m2s(p, read_list, save_path)
            except:
                traceback.print_exc()
                sys.stderr.write("Failed to extract from multi_fast5: {}\n".format(p))

    else:
        sys.stderr.write("Something went wrong, check --f5_format: {}\n".format(args.f5_format))



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
    if args.f5_format == "multi":
        save_path = os.path.abspath("fast5_fetcher_temp")

    if path.endswith('.tar'):
        if OSystem in ["Linux", "Windows"]:
            cmd = "tar -xf {} --transform='s/.*\///' -C {} {}".format(
                path, save_path, filename)
        elif OSystem == "Darwin":
            cmd = "gtar -xf {} --transform='s/.*\///' -C {} {}".format(
                path, save_path, filename)
        else:
            sys.stderr.write("Unsupported OSystem, trying Tar anyway, OS: {}\n".format(OSystem))
            cmd = "tar -xf {} --transform='s/.*\///' -C {} {}".format(
                path, save_path, filename)
    else:
        cmd = "cp {} {}".format(filename, save_path)
    subprocess.call(cmd, shell=True, executable='/bin/bash')


if __name__ == '__main__':
    main()
