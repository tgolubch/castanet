import subprocess as sp

def shell(args, executable='/bin/bash'):
    '''Call Bash shell with input string as argument'''
    sp.run(args, text=True,shell=True, executable=executable)
