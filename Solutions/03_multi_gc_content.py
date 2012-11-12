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
    03_multi_gc_content.py

"""


import sys
import os
import multiprocessing
import timeit

from glob import glob
from Bio import SeqIO

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

def gc_content(seq, pool, chunk_sz):
    try:
        return sum(pool.map(seq.count, ["G", "g", "C", "c", "S", "s"],
            chunksize=chunk_sz)) * 100.0 / len(seq)
    except ZeroDivisionError:
        return 0.0

if __name__ == '__main__':
    num_cpu = int(sys.argv[1])
    chunk_sz = int(sys.argv[2])
    files = sorted(glob("Sessions/Data/03_Sequence_Classification/*"))
    record = read_record(files[0])
    print("loaded sequence")
    # linear approach
    setup = """
from Bio.SeqUtils import GC
from __main__ import record
    """
    print "linear:", min(timeit.repeat(stmt="GC(record.seq)", setup=setup,
        repeat=3, number=10)), "s"
    # parallel approach
    pickle(MethodType, _pickle_method, _unpickle_method)
    pool = multiprocessing.Pool(num_cpu)
    setup = """
from __main__ import gc_content
from __main__ import record
from __main__ import pool
from __main__ import chunk_sz
    """
    print "parallel:", min(timeit.repeat(stmt="gc_content(record.seq, pool, chunk_sz)",
        setup=setup, repeat=3, number=10)), "s"

