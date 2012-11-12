#! /usr/bin/env python
# -*- coding: utf-8 -*-


"""
=========================
Multiprocessing Solutions
=========================

:Authors:
    Moritz Emanuel Beber
:Date:
    2012-11-06
:Copyright:
    Copyright(c) 2012 Jacobs University of Bremen. All rights reserved.
:File:
    03_multi_tetra.py

"""


import sys
import os
import multiprocessing
import timeit

from glob import glob
from Bio import SeqIO

from collections import defaultdict
from copy_reg import pickle
from types import MethodType


def _pickle_method(method):
    func_name = method.im_func.__name__
    obj = method.im_self
    cls = method.im_class
    return _unpickle_method, (func_name, obj, cls)

def _unpickle_method(func_name, obj, cls):
    for cls in cls.mro():
        try:
            func = cls.__dict__[func_name]
        except KeyError:
            pass
        else:
            break
    return func.__get__(obj, cls)


def read_record(filename):
    """
    Read a sequence from a file containing a single record.
    """
    ext = os.path.splitext(filename)[1].lower()
    if ext[1] == "f":
        frmt = "fasta"
    elif ext[1] == "g":
        frmt = "genbank"
    else:
        frmt = "fasta"
    with open(filename, "r") as file_h:
        record = SeqIO.read(file_h, frmt)
    return record

def sub_freq(seq, sub_sz):
    freq = defaultdict(int)
    for i in range(len(seq) - sub_sz + 1):
        freq[seq[i: i + sub_sz]] += 1
    return dict(freq)

def parallel_freq((seq, sub_sz)):
    freq = defaultdict(int)
    for i in range(len(seq) - sub_sz + 1):
        freq[seq[i: i + sub_sz]] += 1
    return dict(freq)

if __name__ == '__main__':
    num_cpu = int(sys.argv[1])
    chunk_sz = int(sys.argv[2])
    files = sorted(glob("../Sessions/Data/03_Sequence_Classification/*"))
    sequences = [str(read_record(filename).seq) for filename in files\
            if not "random" in filename]
    print("loaded sequences")
    # linear approach
    setup = """
from __main__ import sub_freq
from __main__ import sequences
    """
#    print "linear:",\
#            min(timeit.repeat(stmt="[sub_freq(seq, 4) for seq in sequences]",
#            setup=setup, repeat=3, number=10)), "s"
    # parallel approach
    pickle(MethodType, _pickle_method, _unpickle_method)
    pool = multiprocessing.Pool(num_cpu)
    setup = """
from __main__ import parallel_freq
from __main__ import sequences
from __main__ import pool
    """
    print "parallel:",\
            min(timeit.repeat(stmt="pool.map(parallel_freq, [(seq, 4) for seq in sequences])",
            setup=setup, repeat=3, number=10)), "s"

