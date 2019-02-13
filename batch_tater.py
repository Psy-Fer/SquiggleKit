import sys
import subprocess
'''
    Potato scripting engaged.

    James M. Ferguson (j.ferguson@garvan.org.au)
    Genomic Technologies
    Garvan Institute
    Copyright 2018

    batch_tater.py takes list/s of files to extract, and speeds it up a bit, by only opening
    one tar file at a time and extracting what is needed.

    To run on sun grid engine using array jobs as a hacky way of doing multiprocessing.
    Also, helps check when things go wrong, and easy to relaunch failed jobs.
    Some things left in from running on some tasty nanopore single cell data.

    sge file:

    source ~/work/venv2714/bin/activate

    FILE=$(ls ./fast5/ | sed -n ${SGE_TASK_ID}p)
    BLAH=fast5/${FILE}

    mkdir ${TMPDIR}/fast5

    time python batch_tater.py tater_master.txt ${BLAH} ${TMPDIR}/fast5/

    echo "size of files:" >&2
    du -shc ${TMPDIR}/fast5/ >&2
    echo "extraction complete!" >&2
    echo "Number of files:" >&2
    ls ${TMPDIR}/fast5/ | wc -l >&2

    echo "copying data..." >&2

    tar -cf ${TMPDIR}/f5f.${SGE_TASK_ID}.tar --transform='s/.*\///' ${TMPDIR}/fast5/*.fast5
    cp ${TMPDIR}/f5f.${SGE_TASK_ID}.tar ./clean_f5s/

    CMD:

    CMD="qsub -cwd -V -pe smp 1 -N batchCln -S /bin/bash -t 1-10433 -tc 80 -l mem_requested=20G,h_vmem=20G,tmp_requested=20G ../batch.sge"

    Launch:

    echo $CMD && $CMD


    stats:

    fastq: 27491304
    mapped: 11740093
    z mode time: 10min
    batch_tater total time: 21min
    per job time: ~28s
    number of CPUs: 100
'''

# being lazy and using sys.argv...i mean, it is pretty lit
master = sys.argv[1]
tar_list = sys.argv[2]
save_path = sys.argv[3]

# this will probs need to be changed based on naming convention
# I think i was a little tired when I wrote this
list_name = tar_list.split('/')[-1]

PATH = 0

# not elegent, but gets it done
with open(master, 'r') as f:
    for l in f:
        l = l.strip('\n')
        l = l.split('\t')
        if l[0] == list_name:
            PATH = l[1]
            break

# for stats later and easy job relauncing
print >> sys.stderr, "extracting:", tar_list
# do the thing. That --transform hack is awesome. Blows away all the leading folders.
if PATH:
    cmd = "tar -xf {} --transform='s/.*\///' -C {} -T {}".format(
        PATH, save_path, tar_list)
    subprocess.call(cmd, shell=True, executable='/bin/bash')

else:
    print >> sys.stderr, "PATH not found! check index nooblet"
    print >> sys.stderr, "inputs:", master, tar_list, save_path
