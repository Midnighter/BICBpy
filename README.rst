BICBpy - Bioinformatics and Computational Biology Tutorial with Python
======================================================================

:Authors:
    Moritz Emanuel Beber
    Marc-Thorsten Hütt
    Nikolaus Sonnenschein
:Institution:
    Jacobs University, Bremen, Germany


Purpose:
--------
- Undergrad teaching
- Python fundamentals and standard library (urllib2, pickle)
- Biopython: Genomic sequences
- NetworkX: Moritz?
- Numpy and Scipy: Solving ODEs

Tutorials:
----------

1. General introduction (using IPython + a Beamer)
    1.1 Python?
    1.2 IPython shell
    1.3 Neat example
    1.4 Syntax
        1.4.1 Indentation
        1.4.2 Assignments
        1.4.3 Functions
        1.4.4 Control structures
    1.5 Everything is an object ...
        1.5.1 Introspection

    1.6. Solve as many python riddles as you can (www.pythonchallenge.com)

2. Working with sequences
    2.1 Read in a genomic nucleotide sequence from a fasta/genbank file (see
    alo SeqIO.read and SeqIO.parse)
    2.2 Determine the GC content (G+CA+T+G+CÂ100 ) of the sequence (see
    Seq.count)
    2.3 Determine all tetranucleotide frequencies of the sequence (a python
    dictionary maybe useful) and plot them (see list_plot)
    2.4 Read in another genomic sequence and compare it to the previous
    sequence using both GC content and tetranucleotide frequencies (see
    scatter_plot and pearsonr)
    2.5 Read in the E. coli, Marinobacter aquaeoli and S. aureus sequences and
    determine their tetranucleotide frequency profiles
    2.6 Use these profiles to predict the origin of the anonymous sequence
    chunks contained in the random_chunks.fasta file. They all have been
    obtained from either of the previously mentioned genomes.


3. NetworkX + GraphTheory + Protein-Protein interaction networks

4. ODEs

Installation:

Resources:
----------

Biopython documentation
http://www.biopython.org/wiki/Documentation

NetworkX documentation
http://networkx.lanl.gov/contents.html

References:

[1] http://en.wikipedia.org/wiki/Pearson_product-moment_correlation_coefficient

Reading Material:

Bassi S, 2007 A Primer on Python for Life Science Researchers. PLoS Comput Biol
3(11): e199. doi:10.1371/journal.pcbi.0030199

Code Like a Pythonista: Idiomatic Python
http://python.net/%7Egoodger/projects/pycon/2007/idiomatic/handout.html
