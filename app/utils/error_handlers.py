import os
import subprocess as sp
from collections import deque

from app.utils.shell_cmds import loginfo, stoperr, read_line

def error_handler_filter_keep_reads(argies):
    '''
    Print errors, and exit if necessary, on bad input data.
    Extend incl/exlc ID lists from input names.
    '''
    retain_list, exclude_list, argies["o"] = deque(), deque(), []
    '''Check in path, set up output path'''
    for inpath in argies["input_file"]:
        if not os.path.isfile(inpath):
            stoperr(f'Unable to open FastQ file {inpath}.')
        outstem = os.path.basename(inpath.split('.gz')[0]) if inpath.endswith('.gz') else os.path.basename(inpath)
        '''Append suffix "filt" to output file'''
        outpath = f'{os.path.splitext(outstem)[0]}_filt.fastq'
        argies["o"].append(outpath)
    loginfo('Output files {}'.format(argies["o"]))

    '''Check input files'''
    if len(argies["input_file"]) != len(argies["o"]):
        stoperr('Could not create output paths for all given input files.')
    if not os.path.isfile(argies["kraken"]):
        stoperr(f'Unable to open Kraken file {argies["kraken"]} for input {argies["input_file"]}.')

    '''Check lineage file; iterate over linF, get ret/excl IDs from names'''
    if argies["RetainNames"] or argies["ExcludeNames"]:
        if not os.path.isfile(argies["LineageFile"]):
            stoperr(f'Unable to open lineage file {argies["LineageFile"]}.')
        else:
            taxa = dict.fromkeys(argies["RetainNames"].split(','), retain_list)
            '''If same taxid in both lists, exclusion takes precedence over retention'''
            taxa.update(dict.fromkeys(argies["ExcludeNames"].split(','), exclude_list)) 
            if '' in taxa: 
                del taxa['']
            cmd = 'zcat' if argies["LineageFile"].endswith('.gz') else 'cat'
            handle = sp.Popen((cmd, argies["LineageFile"]), bufsize=8192, stdout=sp.PIPE).stdout
            line = read_line(handle)
            while line:
                for taxname, taxlist in taxa.items():
                    if f',{taxname},' in line:
                        taxlist.append(line.split(',', 1)[0])
                line = read_line(handle)

    '''Check NCBI TaxIDs to retain/exclude'''
    try:
        argies["ExcludeIds"] = (frozenset(argies["ExcludeIds"].split(',')) | frozenset(exclude_list)) - frozenset([''])
    except:
        stoperr(f'TaxID(s) {argies["ExcludeIds"]} invalid.')
    try:
        argies["retain"] = (frozenset(argies["RetainIds"].split(',')) | frozenset(retain_list)) - frozenset([''])
    except:
        stoperr(f'TaxID(s) {argies["RetainIds"]} invalid.')
    if not argies["ExcludeIds"] and not argies["RetainIds"]:
        stoperr('Nothing to do. Exiting.')

    return argies["o"], argies["ExcludeIds"], argies["RetainIds"]

