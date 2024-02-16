import numpy as np
import sys
import pandas as pd
import argparse
import pyslow5 as slow5
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['figure.figsize'] = [20.0, 12.0]



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
    parser.add_argument("-f", "--slow5",
                        help="slow5 file")
    parser.add_argument("-c", "--start_col", type=int, default="4",
                        help="start column for signal")
    parser.add_argument("-p", "--plot", action="store_true",
                        help="Live plot each segment")

    args = parser.parse_args()

    # print help if no arguments given
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    # arguments...put this into something better for Tansel
    sig_file = args.signal         # signal file
    # w = 2000
    t_start = 1000
    t_end = 5000

    if args.slow5:
        s5 = slow5.Open(args.slow5, 'r')
        for read in s5.seq_reads():
            readID = read['read_id']
            sig = scale_outliers(read['signal'])

            prev = False  # previous string
            error = 5
            no_err_thresh = 2500
            err = 0       # total error
            prev_err = 0  # consecutive error
            c = 0         # counter
            # w = args.corrector        # window to increase total error thresh
            w = 1200
            window = 100
            # seg_dist = args.seg_dist  # distance between 2 segs to be merged as one
            seg_dist = 1200
            start = 0     # start pos
            end = 0       # end pos
            segs = []     # segments [(start, stop)]
            adapter_found = False

            median = np.median(sig[t_start:t_end])
            stdev = np.std(sig[t_start:t_end])
            top = median + (stdev * 0.8)
            count = -1
            for a in sig:
                count += 1
                if a < top: # If datapoint is within range
                    # print("in range")
                    if not prev:
                        start = count
                        prev = True
                        err = 0
                    c += 1 # increase counter
                    if prev_err:
                        prev_err = 0
                    # if current window longer than detect limit, and corrector, and is divisible by corrector
                    if c >= window and c >= w and not c % w:
                        err -= 1 # drop current error count by 1
                # not within range
                else:
                    # winthin segment and less than error
                    if prev and err < error:
                        c += 1
                        if count >= no_err_thresh:
                            err += 1
                            prev_err += 1
                        if c >= window and c >= w and not c % w:
                            err -= 1
                    # within segment, above error, and greater than window
                    elif prev:
                        if c >= window:
                            # go back to where error stretch began for accurate cutting
                            end = count - prev_err
                            prev = False
                            # if segs very close, merge them
                            if segs and start - segs[-1][1] < seg_dist:
                                segs[-1][1] = end
                            else:
                                # save segment
                                segs.append([start, end])
                            c = 0
                            err = 0
                            prev_err = 0
                        # within segment but not long enough
                        else:
                            prev = False
                            c = 0
                            err = 0
                            prev_err = 0
                    elif segs and count - segs[-1][1] > seg_dist:
                        break_point = count
                        prev = False
                        c = 0
                        err = 0
                        prev_err = 0
                        adapter_found = True
                        break
                    else:
                        continue
                if adapter_found:
                    break
                
                if not adapter_found:
                    break_point = count
                    # print("{}\t{}\t{}".format(readID, ".", "."))

            for a, b in segs:
                x, y = a, b
                print("{}\t{}\t{}".format(readID, x, y))
                break
            




        # for read in s5.seq_reads():
        #     readID = read['read_id']
        #     # if readID != "1e94e97d-8256-4a42-abe3-52ecefc2674e":
        #     #     continue
        #     sig = scale_outliers(read['signal'])

        #     # s = pd.Series(sig)
        #     # t = s.rolling(window=w).mean()
        #     # mn = s.mean()
        #     # std = t.std()
        #     mn = np.mean(sig)
        #     std = np.std(sig)
        #     # this controls the line the signal has to break to start/end a segment
        #     # bot = mn - (std*0.5)
        #     bot = mn - (std*0.1)
        #     # bot = mn

        #     # main algo

        #     begin  = False
        #     # Any segments less than this value are merged into the current segment
        #     seg_dist = 1500
        #     # no final segment can be longer than this
        #     hi_thresh = 200000
        #     # no final segment can be smaller than this
        #     lo_thresh = 2000

        #     start = 0
        #     end = 0
        #     segs = []
        #     count = -1
        #     for i in sig:
        #         count += 1
        #         # if the signal is below the line and a segment has not started yet, start a segment
        #         if i < bot and not begin:
        #             start = count
        #             begin = True
        #         # if we are in a segment, and it's still below the line, continue a segment
        #         elif i < bot:
        #             end = count
        #         # if the signal goes above the line while in a segment, we cut the end of the segment and do some tests.
        #         elif i > bot and begin:
        #             # Check to see if this segment is within seg_dist of last segment, if so, extend the last segment to the end value of this segment
        #             if segs and start - segs[-1][1] < seg_dist:
        #                 segs[-1][1] = end
        #             # otherwise, we have our segment
        #             else:
        #                 segs.append([start, end])
        #             # reset flags
        #             start = 0
        #             end = 0
        #             begin = False
        #         else:
        #             continue
            
        #     # this filters too long or short segments
        #     # and does some adjustments because of the window length from the rolling average calculations
        #     x, y = 0, 0
        #     for a, b in segs:
        #         if b - a > hi_thresh:
        #             continue
        #         if b - a < lo_thresh:
        #             continue
        #         # if a <= 2000:
        #         #     # x, y = a, b - 1000
        #         #     x, y = 0, b - 1000
        #         # else:
        #         #     # x, y = a - 1000, b - 1000
        #         #     x, y = a - 2000, b - 2000
        #         x, y = a, b
        #         print("{}\t{}\t{}".format(readID, x, y))
        #         break
        
            if args.plot:
                fig = plt.figure(1)
                ax = fig.add_subplot(111)
                fig.suptitle("readID: {}\nstart: {}, stop: {}\nbot: {}".format(readID, x, y, round(top, 2)), fontsize=16)
                ax.axvline(x=x, color='m')
                ax.axvline(x=y, color='m')
                ax.axvline(x=break_point, color='g')
                ax.axvline(x=t_start, color='b')
                ax.axvline(x=t_end, color='b')
                ax.axhline(y=top, color='b')
                ax.axvspan(x, y, alpha=0.5, color='orange')
                plt.plot(sig, color='k')
                plt.show()
                plt.clf()



    else:
        with open(sig_file, 'rt') as s:
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
