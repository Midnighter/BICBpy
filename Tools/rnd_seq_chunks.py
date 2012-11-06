#! /usr/bin/env python
# -*- coding: utf-8 -*-


"""
======================================
Configure IPython Notebook Environment
======================================

:Authors:
    Moritz Emanuel Beber
    Nikolaus Sonnenschein
:Date:
    2010-10-25
:Copyright:
    Copyright(c) 2010 Jacobs University of Bremen. All rights reserved.
:File:
    rnd_seq_chunks.py

Reads a number of separate sequence records in fasta or genbank format and
constructs a single fasta file with random chunks of those sequences.
Additionally, a dict with the correct mapping from random chunk to sequence
record is pickled.
"""


import sys
import os
import random
import uuid
import cPickle as pickle

from datetime import date
from Bio import SeqIO
# there was some reorganisation of the Bio package structure
try:
    from Bio.Seq import SeqRecord
except ImportError:
    from Bio.SeqRecord import SeqRecord


def random_chunks(records, num_chunks, chunk_sz):
    """Construct random sequence chunks."""
    rnd_chunks = list()
    mapping = dict()
    for rec in records:
        end = len(rec.seq) - chunk_sz
        for c in range(num_chunks):
            start = random.randint(0, end)
            rnd_seq = rec.seq[start: start + chunk_sz]
            rnd_id = str(uuid.uuid1())
            rnd_chunks.append(SeqRecord(seq=rnd_seq, id=rnd_id))
            mapping[rnd_id] = rec.id
    random.shuffle(rnd_chunks)
    return (rnd_chunks, mapping)

def read_record(filename):
    """
    Read a sequence from a file containing a single record.
    """
    ext = os.splitext(filename)[1].lower()
    if ext[1] == "f":
        frmt = "fasta"
    elif ext[1] == "g":
        frmt = "genbank"
    else:
        frmt = "fasta"
    with open(filename, "r") as file_h:
        record = SeqIO.read(file_h, frmt)
    return record

def main(seq_files, output, num_chunks=10, chunk_sz=1000):
    num_chunks = int(num_chunks)
    chunk_sz = int(chunk_sz)
    today = date.today()
    seq_output = "{0}_{1:%Y-%m-%d}_sz_{2:d}.fasta".format(output, today, chunk_sz)
    map_output = "{0}_{1:%Y-%m-%d}_sz_{2:d}.pkl".format(output, today, chunk_sz)
    records = [read_record(filename) for filename in seq_files]
    (rnd_chunks, mapping) = random_chunks(records, num_chunks, chunk_sz)
    with open(seq_output, "w") as file_h:
        SeqIO.write(rnd_chunks, file_h, "fasta")
    with open(map_output, "wb") as file_h:
        pickle.dump(mapping, file_h, protocol=pickle.HIGHEST_PROTOCOL)


if __name__ == '__main__':
    if len(sys.argv) < 5:
        print "Usage:\npython {0} <num chunks: int> <chunk size: int> "\
                "<sequence record: path> [sequence record 2] [...] "\
                "<output: path>".format(sys.argv[0])
        sys.exit(2)
    main(sys.argv[3:-2], sys.argv[-1], sys.argv[1], sys.argv[2])

