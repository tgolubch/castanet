import subprocess as sp
import sys 

def shell(args, executable='/bin/bash'):
    '''Call Bash shell with input string as argument'''
    sp.run(args, text=True,shell=True, executable=executable)

def loginfo(s):
    '''Log Info statements to stderr'''
    sys.stderr.write(f'  Info: {s}\n')

def logerr(s):
    '''Log Warning/Error statements to stderr'''
    sys.stderr.write(f'  Warning: {s}\n')

def stoperr(s, errcode=1):
    '''Call sys exit on finish or err'''
    status = 'Finished' if not errcode else 'Error'
    sys.stderr.write(f'  {status}: {s}\n')
    raise SystemError(f'  {status}: {s}\nErrcode: {errcode}')

def read_line(h):
    '''Parse line from input file'''
    return h.readline().decode("utf-8")