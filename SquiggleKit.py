import os
import sys
import argparse
import traceback
import numpy as np
import slow5
#import pod5
import h5py
import time

"""
SquiggleKit python re-write 2023

re-writing SquiggleKit to be up to date with the latest file formats and methods

Supporting:
 - pod5
 - slow5 (prefered)

"""

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)


def main():
    parser = MyParser(
        description="SquiggleKit - set of tools for processing and manipulating Oxford Nanopore raw signal data",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subcommand = parser.add_subparsers(help='subcommand --help for help messages', dest="command")

    # main options
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Engage higher output verbosity")

    # sub-module for SquigglePlot
    plot = subcommand.add_parser('plot', help='plot the signal data after (optional) conversion to pA',
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # SquigglePlot sub-module options
    plot_group = plot.add_mutually_exclusive_group()
    #   need to make it so that raw_signal flag cannot be applied with the signal file
    #   currently cannot support it unless output file contains digitisation, range and offest values (extra info mode)

    # plot_group.add_argument("-p", "--f5_path",
    #                     help="Fast5 top dir")
    # plot_group.add_argument("-s", "--signal",
    #                     help="Extracted signal file from SquigglePull")
    # plot_group.add_argument("-i", "--ind", nargs="+",
    #                     help="Individual fast5 file/s.")
    plot.add_argument("-r", "--readID",
                        help="Individual readID to extract from a multifast5 file")
    #plot.add_argument("--single",action="store_true",
    #                    help="single fast5 files")
    # plot.add_argument("-t", "--type", action="store", default="auto", choices=["auto", "single", "multi"], help="Specify the type of files provided. Default is autodetection which enables a mix of single and multifast5 files.")
    # plot.add_argument("--head", action="store_true",
    #                     help="Header present in signal or flat file")
    plot.add_argument("--raw_signal",action="store_true",
                        help="Plot raw signal instead of converting to pA")
    plot.add_argument("-n", "--Num",
                        help="Section of signal to look at - -n 2000 or -n 100,1500")
    plot.add_argument("--lim_hi", type=int, default=1200,
                        help="Upper limit for signal outliers")
    plot.add_argument("--lim_low", type=int, default=0,
                        help="Lower limit for signal outliers")
    # Arguments for now, but best way forward will probably be a config file
    plot.add_argument("--plot_colour", default='grey',
                        help="Colour of signal plot, takes any pyplot entry: k,r,b,g,red,blue,etc...")
    plot.add_argument("--save",
                        help="Save file readname_saveArg.pdf --save saveArg.pdf, use png, etc for other file types")
    plot.add_argument("--save_path",
                        help="Save filepath")
    plot.add_argument("--no_show", action="store_true",
                        help="Do not show plot (used for saving many)")
    plot.add_argument("--dpi", type=int, default=100,
                        help="Change DPI for publication figs, eg: --dpi 300")

    # sub-module for segmenter
    segment = subcommand.add_parser('segmenter', help='find obvious structural regions in squiggle data',
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # segmenter sub-module options
    segmenter_group = segment.add_mutually_exclusive_group()
    # segmenter_group.add_argument("-f", "--f5f",
    #                    help="File list of fast5 paths")
    # segmenter_group.add_argument("-p", "--f5_path",
    #                    help="Fast5 top dir")
    # segmenter_group.add_argument("-s", "--signal",
    #                    help="Extracted signal file from squigglePull")
    # segment.add_argument("--single",action="store_true",
    #                     help="single fast5 files")
    segment.add_argument("-n", "--Num", type=int, default=0,
                        help="Section of signal to look at - default 0=all")
    segment.add_argument("-e", "--error", type=int, default=5,
                        help="Allowable error in segment algorithm")
    segment.add_argument("-c", "--corrector", type=int, default=50,
                        help="Window size for increasing total error correction - better long segment detection")
    segment.add_argument("-w", "--window", type=int, default=150,
                        help="Minimum segment window size to be detected")
    segment.add_argument("-d", "--seg_dist", type=int, default=50,
                        help="Maximum distance between 2 segments to be merged into 1")
    segment.add_argument("-t", "--std_scale", type=float, default=0.75,
                        help="Scale factor of STDev about median")
    segment.add_argument("-v", "--view", action="store_true",
                        help="view each output")
    segment.add_argument("-g", "--gap", action="store_true",
                        help="Turn on gap distance for stall to polyTAil")
    segment.add_argument("-b", "--gap_dist", type=int, default=3000,
                        help="Maximum distance between stall and polyTAil segment - for 10X/dRNA")
    segment.add_argument("-k", "--stall", action="store_true",
                        help="Turn on stall detection - must be present")
    segment.add_argument("-u", "--test", action="store_true",
                        help="Run Tests")
    segment.add_argument("-l", "--stall_len", type=float, default=0.25,
                        help="Minimum percentage of minimum window segment for initial stall segment")
    segment.add_argument("-j", "--stall_start", type=int, default=300,
                        help="Maximum distance for start of stall segment to be detected")
    segment.add_argument("-lim_hi", "--lim_hi", type=int, default=900,
                        help="Upper limit for signal outlier scaling")
    segment.add_argument("-lim_low", "--lim_low", type=int, default=0,
                        help="Lower limit for signal outlier scaling")
    segment.add_argument("--raw_signal",action="store_true",
                        help="Plot raw signal instead of converting to pA")

    # sub-module for MotifSeq
    motifseq = subcommand.add_parser('motifseq', help='MotifSeq - the Ctrl+f for signal. Signal-level local alignment of sequence motifs',
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    motifseq_group = motifseq.add_mutually_exclusive_group()
    motifseq_mods = motifseq.add_mutually_exclusive_group()
    # motifseq_group.add_argument("-f", "--f5f",
    #                    help="File list of fast5 paths")
    # motifseq_group.add_argument("-p", "--f5_path",
    #                    help="Fast5 top dir")
    # motifseq_group.add_argument("-s", "--signal",
    #                    help="Extracted signal file from SquigglePull")
    motifseq.add_argument("-l", "--scale", default="medmad", choices=["zscale", "medmad"],
                       help="scaling/normalisation factor to use")
    motifseq_mods.add_argument("-i", "--fasta_input",
                        help="fasta file to be converted to simulated signal by scrappy")
    motifseq.add_argument("--scrappie_model", default="squiggle_r94", choices=['squiggle_r94', 'squiggle_r94_rna', 'squiggle_r10'],
                        help="model to use with fasta_input for conversion")
    motifseq_mods.add_argument("-m", "--model",
                        help="custom multiline .tsv of signal to search for - see docs - name{tab}60{tab}435...")
    motifseq.add_argument("-x", "--sig_extract", action="store_true",
                        help="Extract signal of match")
    motifseq.add_argument("--slope", type=float, default=2.90,
                        help="[Experimental] slope")
    motifseq.add_argument("--intercept", type=float, default=-9.6,
                        help="[Experimental] intercept")
    motifseq.add_argument("--std_const", type=float, default=0.08468,
                        help="[Experimental] standard deviation constant")
    motifseq.add_argument("-v", "--view", action="store_true",
                        help="view each output")
    motifseq.add_argument("--save",
                        help="save path for images")
    motifseq.add_argument("--img", default="png",
                        help="Type of image to save. png, jpeg, pdf, svg, etc. (default: png)")
    motifseq.add_argument("-scale_hi", "--scale_hi", type=int, default=1200,
                        help="Upper limit for signal outlier scaling")
    motifseq.add_argument("-scale_low", "--scale_low", type=int, default=0,
                        help="Lower limit for signal outlier scaling")



    # collect args
    args = parser.parse_args()

    # print help if no arguments given
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    if args.verbose:
        sys.stderr.write("Verbose mode on. Starting timer.\n")
        start_time = time.time()

    elif args.command == "plot":
        if len(sys.argv) == 2:
            plot.print_help(sys.stderr)
            sys.exit(1)
        squigglePlot(args)

    elif args.command == "segmenter":
        if len(sys.argv) == 2:
            segment.print_help(sys.stderr)
            sys.exit(1)
        segmenter(args)


    elif args.command == "motifseq":
        if len(sys.argv) == 2:
            motifseq.print_help(sys.stderr)
            sys.exit(1)
        motifseq(args)

    else:
        parser.print_help(sys.stderr)
        sys.exit(1)
    
    if args.verbose:
        end_time = time.time() - start_time
        sys.stderr.write("Time taken: {}\n".format(end_time))

    sys.stderr.write("Done\n")




def squigglePlot():
    return

def segmenter():
    return

def motifseq():
    return