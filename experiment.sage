'''
Generic tools to run experiments in Sage.
'''

import glob
from multiprocessing.pool import Pool

def run_all_parallel(glob_pattern, experiment_function)
    instances = glob.glob(glob_pattern)
    pool = Pool()
    pool.map(experiment_function, instances)
