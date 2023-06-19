import os


def make_exp_dir(ExpName):
    '''Checks if dir exists; creates if not'''
    if not os.path.exists(f"experiments/{ExpName}"):
        os.makedirs(f"experiments/{ExpName}")
