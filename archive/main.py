from flask import Flask, render_template, request, redirect, Markup
import os
import sys
import time
import h5py
import traceback
import numpy as np
from bokeh.plotting import figure
from bokeh.resources import CDN
from bokeh.embed import file_html
from bokeh.models import Title, HoverTool, ColumnDataSource, FreehandDrawTool, BoxEditTool, BoxAnnotation, CustomJS, Rect, Spacer
from bokeh.models.widgets.buttons import AbstractButton, Toggle
import json

from tornado.ioloop import IOLoop
from threading import Thread
from bokeh.embed import server_document
from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource, Slider
from bokeh.server.server import Server
from bokeh.themes import Theme


app = Flask(__name__)

signal = []
ut = None
lt = None
err = 5
err_win = 50
min_win = 150
max_merge = 50
stdev_scale = 0.75
stall_len = 0.25

@app.route('/test', methods=['GET'])
def bkapp_page():
    script = server_document('http://localhost:5006/bkapp')
    return render_template("embed.html", script=script, template="Flask")

@app.route("/",methods = ["POST", "GET"])
def home():
    if request.method == "POST":
        f5_path = request.form.get('f5_path')
        type = request.form.get('type')
        if not os.path.isdir(f5_path):
            return render_template("home.html", error=True, f5_path=f5_path)
        else:
            return render_template("loading.html", f5_path=f5_path, type=type)

    return render_template("home.html")

@app.route("/results")
def results(f5_path=None):

    f5_path = request.args['f5_path']
    type = request.args['type']
    if not os.path.isdir(f5_path):
        return render_template("error.html")
    if request.args['processing'] == '0':
        if os.path.isfile(f5_path+"/data_"+type+".tsv"):
            exception = "Data file for already exists"
            return render_template("exception.html", f5_path=f5_path, type=type, exception=exception)

    if request.args['processing'] == '1':
        os.remove(f5_path+"/data_"+type+".tsv")
        return render_template("loading.html", f5_path=f5_path, type=type)

    count = 0
    with open(os.path.join(f5_path, "data_"+type+".tsv"), 'a') as out_sum:
        for dirpath, dirnames, files in os.walk(f5_path):
            for fast5 in files:
                if fast5.endswith('.fast5'):
                    fast5_file = os.path.join(dirpath, fast5)
                    # extract data from file
                    data, multi = extract_f5_all(fast5_file, request.args['type'])
                    #print data to a single file
                    if not multi:
                        count += 1
                        ar = map(str, data['raw'])

                        out_sum.write('{}\t{}\t{}\n'.format(
                                fast5, data['readID'], '\t'.join(ar)))
                    else:
                        for read in data:
                            count += 1
                            ar = map(str, data[read]['raw'])

                            out_sum.write('{}\t{}\t{}\n'.format(
                                    fast5, data[read]['readID'], '\t'.join(ar)))
    return render_template("results.html", f5_path=f5_path, type=type, count=count)

@app.route("/view_graphs")
def view():

    f5_path = request.args['f5_path']
    type = request.args['type']
    read = request.args.get('read_id')
    id = ''

    script = []

    if read is None:
        read = ""

    reads = []
    sig = None
    segs = None

    if not os.path.isfile(f5_path+"/data_"+type+".tsv"):
        return render_template("error.html")

    with open(f5_path+"/data_"+type+".tsv", 'rt') as data:
        for num, l in enumerate(data):
            l = l.strip('\n')
            l = l.split('\t')
            readID = l[1]
            reads.append(l[1])
            if read == readID:
                fast5 = l[0]
                if "." in l[4]:
                    sig = np.array([float(i) for i in l[4:]], dtype=float)
                else:
                    sig = np.array([int(i) for i in l[4:]], dtype=int)

    graph = dict()
    if sig is not None:
        global signal
        signal = sig
        print(signal)
        Thread(target=bk_worker).start()
        id = str(read)

        script = server_document('http://localhost:5006/bkapp')
        return render_template("view_graphs.html", f5_path=f5_path, type=type, id=id, script=script, count=len(reads), reads=reads)
    else:
        error = "The signal was unable to be found for "+read+" :(."
        return render_template("error.html", error=error)

@app.route("/delete")
def delete():
    f5_path = request.args['f5_path']
    type = request.args['type']
    if os.path.isfile(f5_path+"/data_"+type+".tsv"):
        os.remove(f5_path+"/data_"+type+".tsv")
    return redirect("/")

def extract_f5_all(filename, type):
    '''
    inputs:
        filepath/name
        args from command line
    does:
        open fast5 files, extract whole signal and read data and converts to pA by default
    Returns:
        dic for further processing/printing
    '''
    f5_dic = {}
    multi = True
    raw = False
    if type == "raw":
        raw = True
    with h5py.File(filename, 'r') as hdf:

        reads = list(hdf.keys())
        if 'read' not in reads[1]:
            multi = False

        # single fast5 files
        if not multi:
            f5_dic = {'raw': [], 'seq': '', 'readID': '',
                    'digitisation': 0.0, 'offset': 0.0, 'range': 0.0,
                    'sampling_rate': 0.0}
            # extract the data
            try:
                c = list(hdf['Raw/Reads'].keys())
                for col in hdf['Raw/Reads/'][c[0]]['Signal'][()]:
                    f5_dic['raw'].append(int(col))
                f5_dic['readID'] = hdf['Raw/Reads/'][c[0]].attrs['read_id'].decode()
                digitisation = hdf['UniqueGlobalKey/channel_id'].attrs['digitisation']
                offset = hdf['UniqueGlobalKey/channel_id'].attrs['offset']
                range = float("{0:.2f}".format(hdf['UniqueGlobalKey/channel_id'].attrs['range']))

                # convert to pA
                if not raw:
                    f5_dic['raw'] = np.array(f5_dic['raw'], dtype=int)
                    f5_dic['raw'] = convert_to_pA_numpy(f5_dic['raw'], digitisation, range, offset)
                    f5_dic['raw'] = np.round(f5_dic['raw'], 2)

            except:
                traceback.print_exc()
                sys.stderr.write("extract_fast5_all():failed to extract raw signal or fastq from {}".format(filename))
                f5_dic = {}

        # multi fast5 files
        else:
            for read in reads:
                f5_dic[read] = {'raw': [], 'seq': '', 'readID': '',
                                'digitisation': 0.0, 'offset': 0.0, 'range': 0.0,
                                'sampling_rate': 0.0}

                # extract the data
                try:
                    for col in hdf[read]['Raw/Signal'][()]:
                        f5_dic[read]['raw'].append(int(col))
                    f5_dic[read]['readID'] = hdf[read]['Raw'].attrs['read_id'].decode()
                    digitisation = hdf[read]['channel_id'].attrs['digitisation']
                    offset = hdf[read]['channel_id'].attrs['offset']
                    range = float("{0:.2f}".format(hdf[read]['channel_id'].attrs['range']))

                    # convert to pA
                    if not raw:
                        f5_dic[read]['raw'] = np.array(f5_dic[read]['raw'], dtype=int)
                        f5_dic[read]['raw'] = convert_to_pA_numpy(f5_dic[read]['raw'], digitisation, range, offset)
                        f5_dic[read]['raw'] = np.round(f5_dic[read]['raw'], 2)

                except:
                    traceback.print_exc()
                    sys.stderr.write("extract_fast5_all():failed to read readID: {}".format(read))

    return f5_dic, multi

def convert_to_pA_numpy(d, digitisation, range, offset):
    raw_unit = range / digitisation
    return (d + offset) * raw_unit

def scale_outliers(sig, max, min):
    ''' Scale outliers to within m stdevs of median '''
    ''' Remove outliers that don't fit within the specified bounds '''
    k = (sig > min) & (sig < max)
    return sig[k]

def get_segs(sig, error, error_win, min_win, max_merge, std_scale, stall_len):
    '''
    Get segments from signal
    This works by running through the signal and finding regions that are above
    the bot and below the top parameters, with some error tollerance, for a
    minimum window of length.
    '''

    mn = sig.min()
    mx = sig.max()
    mean = np.mean(sig)
    median = np.median(sig)
    # use this with outlier rejection to fix stdev thresholds
    stdev = np.std(sig)
    top = median + (stdev * std_scale)
    bot = median - (stdev * std_scale)

    # parameter tuning visualisation
    # TODO: Put tuning plots here

    # this is the algo. Simple yet effective
    prev = False  # previous string
    err = 0       # total error
    prev_err = 0  # consecutive error
    c = 0         # counter
    w = error_win        # window to increase total error thresh
    seg_dist = max_merge # distance between 2 segs to be merged as one
    start = 0     # start pos
    end = 0       # end pos
    segs = []     # segments [(start, stop)]
    left = []
    right = []
    for i in range(len(sig)):
        a = sig[i]
        if a < top and a > bot: # If datapoint is within range
            if not prev:
                start = i
                prev = True
            c += 1 # increase counter
            w += 1 # increase window corrector count
            if prev_err:
                prev_err = 0
            if c >= min_win and c >= w and not c % w: # if current window longer than detect limit, and corrector, and is divisible by corrector
                err -= 1 # drop current error count by 1
        else:
            if prev and err < error:
                c += 1
                err += 1
                prev_err += 1
                if c >= min_win and c >= w and not c % w:
                    err -= 1
            elif prev and (c >= min_win or not segs and c >= min_win * stall_len):
                end = i - prev_err # go back to where error stretch began for accurate cutting
                prev = False
                if segs and start - segs[-1][1] < seg_dist: # if segs very close, merge them
                    segs[-1][1] = end
                else:
                    segs.append([start,end])
                    left.append(start)
                    right.append(end) # save segment
                c = 0
                err = 0
                prev_err = 0
            elif prev:
                prev = False
                c = 0
                err = 0
                prev_err = 0
            else:
                continue

    if segs:
        return left, right
    else:
        return False

def bkapp(doc):
    global signal
    global ut
    global lt
    global show_segs
    show_segs = False
    ut = 0
    lt = 0
    if signal.any():
        ut = max(signal)
        lt = min(signal)
        source = ColumnDataSource(data={
            'signal'    : signal,
            'position'  : list(range(0,len(signal)))
            })

        p = figure()

        p.line('position','signal', source=source)

        p.add_tools(HoverTool(
            tooltips=[
                ('signal', '@signal'),
                ('position', '@position'),
            ],
            formatters={
                'signal'    : 'printf',
                'position'    : 'printf'
            },
            mode='vline'
        ))

        renderer = p.multi_line([[1,1]], [[1,1]], line_width=4, alpha=0.5, color='green')
        draw_tool = FreehandDrawTool(renderers=[renderer])
        p.add_tools(draw_tool)

        src = ColumnDataSource({
            'x':[1,1,1], 'y':[1,1,1], 'width':[1,1,1], 'height':[1,1,1]
        })
        box_renderer = p.rect('x', 'y', 'width', 'height', fill_alpha=0.4, fill_color='orange', line_color='orange', source=src)
        box_draw_tool = BoxEditTool(renderers=[box_renderer], empty_value=1, num_objects = 5)
        p.add_tools(box_draw_tool)

        ut_slider = Slider(start=lt, end=max(signal), value=max(signal), name='upper_thresh', step=1, title="Upper Threshold")
        lt_slider = Slider(start=min(signal), end=ut, value=min(signal), name='lower_thresh', step=1, title="Lower Threshold")

        def ut_callback(attr, old, new):
            global signal
            global ut
            global lt
            ut = new
            new_signal = scale_outliers(signal, ut, lt)
            source.data = {
                        'signal'    : new_signal,
                        'position'  : list(range(0,len(new_signal)))
                        }
            update_segs()

        def lt_callback(attr, old, new):
            global signal
            global ut
            global lt
            lt = new
            new_signal = scale_outliers(signal, ut, lt)
            source.data = {
                        'signal'    : new_signal,
                        'position'  : list(range(0,len(new_signal)))
                        }
            update_segs()

        ut_slider.on_change('value', ut_callback)
        lt_slider.on_change('value', lt_callback)

        segments = ColumnDataSource(data={
            'top'     : [1,1],
            'bottom'  : [1,1],
            'left'    : [1,1],
            'right'   : [1,1]
            })

        button = Toggle(label="View Segments", sizing_mode="scale_width")

        def segment_handler(new):
            global show_segs
            show_segs = new
            if not new:
                segments.data = {
                    'top'     : [1,1],
                    'bottom'  : [1,1],
                    'left'    : [1,1],
                    'right'   : [1,1]
                            }
            update_segs()

        button.on_click(segment_handler)

        err_slider = Slider(start=0, end=20, value=5, name='error', step=1, title="Allowable Error")
        err_win_slider = Slider(start=0, end=100, value=50, name='err_win', step=1, title="Error Window Size")
        min_win_slider = Slider(start=0, end=500, value=150, name='min_win', step=1, title="Minimum Window Size")
        max_merge_slider = Slider(start=0, end=100, value=50, name='max_merge', step=1, title="Max Merge Distance")
        stdev_scale_slider = Slider(start=0, end=5, value=0.75, name='stdev_scale', step=0.01, title="Standard Deviation Scale Factor")
        stall_len_slider = Slider(start=0, end=5, value=0.25, name='stall_len', step=0.01, title="Stall Length")

        p.quad(top='top',bottom='bottom',left='left',right='right',source=segments,fill_alpha=0.5,fill_color='pink',line_color='pink')

        def err_callback(atrr, old, new):
            global err
            err = new
            update_segs()

        def err_win_callback(atrr, old, new):
            global err_win
            err_win = new
            update_segs()

        def min_win_callback(atrr, old, new):
            global min_win
            min_win = new
            update_segs()

        def max_merge_callback(atrr, old, new):
            global max_merge
            max_merge = new
            update_segs()

        def stdev_scale_callback(atrr, old, new):
            global stdev_scale
            stdev_scale = new
            update_segs()

        def stall_len_callback(atrr, old, new):
            global stall_len
            stall_len = new
            update_segs()

        def update_segs():
            #need to take into account the modified signal- somehow access it?
            global err
            global err_win
            global min_win
            global max_merge
            global stdev_scale
            global stall_len
            global ut
            global lt
            global show_segs
            left = None
            right = None
            if show_segs:
                sig = scale_outliers(signal, ut, lt)
                if sig.any():
                    left, right = get_segs(sig, err, err_win, min_win, max_merge, stdev_scale, stall_len)
                if left is not None and right is not None:
                    segments.data = {
                                'top'     :  np.full(len(left),1000),
                                'bottom'  :  np.full(len(left),0),
                                'left'    :  left,
                                'right'   :  right
                                }
                else:
                    segments.data = {
                        'top'     : [1,1],
                        'bottom'  : [1,1],
                        'left'    : [1,1],
                        'right'   : [1,1]
                                }


        err_slider.on_change('value', err_callback)
        err_win_slider.on_change('value', err_win_callback)
        min_win_slider.on_change('value', min_win_callback)
        max_merge_slider.on_change('value', max_merge_callback)
        stdev_scale_slider.on_change('value', stdev_scale_callback)
        stall_len_slider.on_change('value', stall_len_callback)

        doc.add_root(row(column(Spacer(height=10), ut_slider, lt_slider, Spacer(height=10), button, err_slider, err_win_slider, min_win_slider, max_merge_slider, stdev_scale_slider, stall_len_slider, Spacer(height=10), sizing_mode="stretch_height"), p, sizing_mode="stretch_both"))
        doc.theme = Theme(filename="theme.yaml")

def bk_worker():
    # Can't pass num_procs > 1 in this configuration. If you need to run multiple
    # processes, see e.g. flask_gunicorn_embed.py
    print("I, the bk_worker, am being run")
    server = Server({'/bkapp': bkapp}, io_loop=IOLoop(), allow_websocket_origin=["127.0.0.1:8080"])
    server.start()
    server.io_loop.start()


if __name__ == "__main__":
    print('Please open the page http://127.0.0.1:8080 to access the SquiggleKit Web Application')
    app.run(port="8080", debug=True)
