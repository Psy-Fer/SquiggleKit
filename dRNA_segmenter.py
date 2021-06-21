import numpy as np
import sys
import pandas as pd
import argparse


'''
    James M. Ferguson (j.ferguson@garvan.org.au)
    Genomic Technologies
    Garvan Institute
    Copyright 2019


    dRNA DNA barcode extraction.
    looks for drop in signal and get's it as a segment

    output structure:
    fast5, readID, start, stop

    MIT License

    Copyright (c) 2018 James Ferguson

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
    main function
    '''

    parser = MyParser(
        description="dRNA_segmenter - cut out adapter region of dRNA signal")
    #group = parser.add_mutually_exclusive_group()
    parser.add_argument("-s", "--signal",
                        help="Signal file")
    parser.add_argument("-c", "--start_col", type=int, default="4",
                        help="start column for signal")

    args = parser.parse_args()

    # print help if no arguments given
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    # arguments...put this into something better for Tansel
    sig = args.signal         # signal file
    w = 2000

    with open(sig, 'rt') as s:
        for read in s:
            read = read.strip('\n')
            read = read.split('\t')
            f5 = read[0]
            readID = read[1]
            sig = scale_outliers(np.array([int(i) for i in read[args.start_col:]], dtype=int))

            s = pd.Series(sig)
            t = s.rolling(window=w).mean()
            mn = t.mean()
            std = t.std()
            # might need to tighten the bounds a little more
            # top = mn + (std*0.5)
            bot = mn - (std*0.5)

            # main algo

            begin  = False
            seg_dist = 1500
            hi_thresh = 200000
            lo_thresh = 2000

            start = 0
            end = 0
            segs = []
            count = -1
            for i in t:
                count += 1
                if i < bot and not begin:
                    start = count
                    begin = True
                elif i < bot:
                    end = count
                elif i > bot and begin:
                    if segs and start - segs[-1][1] < seg_dist:
                        segs[-1][1] = end
                    else:
                        segs.append([start, end])
                    start = 0
                    end = 0
                    begin = False
                else:
                    continue

            x, y = 0, 0
            for a, b in segs:
                if b - a > hi_thresh:
                    continue
                if b - a < lo_thresh:
                    continue
                x, y = a - 1000, b - 1000
                print("{}\t{}\t{}\t{}".format(f5, readID, x, y))
                break


def scale_outliers(squig):
    ''' Scale outliers to within m stdevs of median '''
    k = (squig > 0) & (squig < 1200)
    return squig[k]


if __name__ == '__main__':
    main()
