#!/usr/bin/env python

# ============================================================================ #
# Copyright (c) Tanya Golubchik                                                #
# golubchi@well.ox.ac.uk                                                       #
# March 2017                                                                   #
# ============================================================================ #

# ============================================================================ #
# Import Modules                                                               #
# ============================================================================ #

from __future__ import division
import os, sys, argparse
import subprocess as sp
from Bio import SeqIO
from collections import deque

# ============================================================================ #
# GLOBAL VARIABLES                                                             #
# ============================================================================ #
_args = None
_progname=os.path.basename(sys.argv[0])

# ============================================================================ #
# LOGGING                                                                      #
# ============================================================================ #
def loginfo(s):
    sys.stderr.write('  Info: {0}\n'.format(s))
def logerr(s):
    sys.stderr.write('  Warning: {0}\n'.format(s))
def stoperr(s, errcode=1):
    errword = 'Finished' if not errcode else 'Error'
    sys.stderr.write('  {0}: {1}\n'.format(errword, s))
    sys.exit(errcode)

# ============================================================================ #
# PROGRAM USAGE                                                                #
# ============================================================================ #
def Initialise():
    '''
    Parse command-line arguments.
    '''
    global _args
    parser = argparse.ArgumentParser(description=__doc__, epilog='Example: {progname} -i /path/to/my/fastq.gz -x 9606'.format(progname=_progname))
    parser.add_argument( '-i', required=True, nargs='*', help='Fastq[.gz] files to filter. Files must NOT contain blank lines. If gzipped, filenames must end in .gz. Each file will produce output with the suffix _filt.fastq in the current working directory.' )
    parser.add_argument( '-k', default=None, required=True, help='Path to corresponding Kraken[.gz] output file, or other whitespace-separated text file where the second column is the read name and the third is the TaxID.' )
    parser.add_argument( '-r', default='', help='NCBI TaxID(s) to retain. Use comma-separated list without spaces for multiple TaxIDs). Default is none.' )
    parser.add_argument( '-x', default='', help='NCBI TaxID(s) to exclude. Use comma-separated list without spaces for multiple TaxIDs). Default is none. To remove all reads marked as human set to 9606.' )
    parser.add_argument( '--rT', default='', help='Text names of top-level taxa, to retain this and all taxa below it. Requires "--lineagefile". Eg. "--rT Bacteria,Viruses" to retain all bacterial and viral sequences.' )
    parser.add_argument( '--xT', default='', help='Text names of top-level taxa, to exclude this and all taxa below it. Requires "--lineagefile". Eg. "--xT \"Homo sapiens,Fungi\"" to exclude all human and fungal sequences.' )
    parser.add_argument( '--suffix', default='filt', help='Suffix to append to filtered files. Default is "filt".' )
    parser.add_argument( '--lineagefile', default='lineages-2018-03-12.csv.gz', help='Path to CSV file containing lineages of all NCBI taxa. Default is "lineages-2018-03-12.csv.gz".' )
    _args = parser.parse_args()
    return

def Clean_Commandline():
    '''
    Print errors, and exit if necessary, on bad input data.
    '''
    _args.o = []
    for inpath in _args.i:
        if not os.path.isfile(inpath):
            stoperr('Unable to open FastQ file {0}.'.format(inpath))
        outstem = os.path.basename(inpath.split('.gz')[0]) if inpath.endswith('.gz') else os.path.basename(inpath)
        outpath = '{0}_{1}.fastq'.format(os.path.splitext(outstem)[0], _args.suffix)
        _args.o.append(outpath)
    loginfo('Output files {}'.format(_args.o))
    if len(_args.i) != len(_args.o):
        stoperr('Could not create output paths for all given input files.')
    if not os.path.isfile(_args.k):
        stoperr('Unable to open Kraken file {0}.'.format(_args.i))
    retain_list, exclude_list = deque(), deque()
    if _args.rT or _args.xT:
        if not os.path.isfile(_args.lineagefile):
            stoperr('Unable to open lineage file {0}.'.format(_args.lineagefile))
        else:
            taxa = dict.fromkeys(_args.rT.split(','), retain_list)
            taxa.update(dict.fromkeys(_args.xT.split(','), exclude_list)) # if same taxid in both lists, exclusion takes precedence over retention
            if '' in taxa: 
                del taxa['']
            cmd = 'zcat' if _args.lineagefile.endswith('.gz') else 'cat'
            handle = sp.Popen((cmd, _args.lineagefile), bufsize=8192, stdout=sp.PIPE).stdout
            line = handle.readline()
            while line:
                for taxname, taxlist in taxa.iteritems():
                    if ',{},'.format(taxname) in line:
                        taxlist.append(line.split(',', 1)[0])
                line = handle.readline()
    try:
        _args.x = (frozenset(_args.x.split(',')) | frozenset(exclude_list)) - frozenset([''])
    except:
        stoperr('TaxID(s) {0} invalid.'.format(_args.x))
    try:
        _args.r = (frozenset(_args.r.split(',')) | frozenset(retain_list)) - frozenset([''])
    except:
        stoperr('TaxID(s) {0} invalid.'.format(_args.r))
    if not _args.r and not _args.x:
        stoperr('Nothing to do. Exiting.')
    return


# ============================================================================ #
# DATA PROCESSING FUNCTIONS                                                    #
# ============================================================================ #

def Filter_Reads(reads_to_exclude, reads_to_keep, inpath, out_h):
    '''
    Collect the four lines belonging to the next read and filter on readname.
    ASSUMES NO BLANK LINES IN INPUT.
    '''
    cmd = 'zcat' if inpath.endswith('.gz') else 'cat'
    handle = sp.Popen((cmd, inpath), bufsize=8192, stdout=sp.PIPE).stdout
    line = handle.readline()
    num_reads = 0
    while line:
        readheader, seq, qh, qual = line, handle.readline(), handle.readline(), handle.readline()
        line = handle.readline()
        readname = readheader[1:].split('/')[0].split()[0]
        to_keep = True if not reads_to_keep or (reads_to_keep and (readname in reads_to_keep)) else False
        to_exclude = True if (reads_to_exclude and (readname in reads_to_exclude)) else False
        if to_keep and not to_exclude:
            out_h.write(''.join((readheader, seq, qh, qual)))
            num_reads+=1
    return num_reads

# ============================================================================ #
# MAIN                                                                         #
# ============================================================================ #

if __name__ == '__main__':
    Initialise()
    Clean_Commandline()
    reads_to_exclude, reads_to_keep = deque(), deque()
    cmd = 'zcat' if _args.k.endswith('.gz') else 'cat'
    handle = sp.Popen((cmd, _args.k), bufsize=8192, stdout=sp.PIPE).stdout
    line = handle.readline()
    while line:
        readname, taxid = line.split()[1:3]
        if taxid in _args.x:
            reads_to_exclude.append(readname)
        elif taxid in _args.r:
            reads_to_keep.append(readname)
        line = handle.readline()
    reads_to_exclude, reads_to_keep = frozenset(reads_to_exclude), frozenset(reads_to_keep)
    loginfo('Excluding {} reads ({} taxa specified).'.format(len(reads_to_exclude), len(_args.x)) )
    loginfo('Retaining {} reads ({} taxa specified).'.format(len(reads_to_keep) if reads_to_keep else 'all other', len(_args.r)) )
    for (inpath, outpath) in zip(_args.i, _args.o):
        with open(outpath, 'w', buffering=8192) as out_h:
            num_reads = Filter_Reads(reads_to_exclude, reads_to_keep, inpath, out_h)
        loginfo('Wrote {0} reads to {1}.'.format(num_reads, outpath))

# ============================================================================ #
