#!env/bin/python

import datetime
from optparse import OptionParser
from cartdb.main import main

# Version
ver = '0.1.0'

# Command line argument parsing
descr = 'CARTdb v'+ver
parser = OptionParser(usage='RefSeqDB/env/bin/cartdb <options>', version=ver, description=descr)
parser.add_option('-i', default=None, dest='input', action='store', help="Input file name")
parser.add_option('-n', default=None, dest='ncbi', action='store', help="RefSeqDB output file with NCBI interim mapping data")
parser.add_option('-u', default=None, dest='ucsc', action='store', help="RefSeqDB output file with UCSC mapping data")
parser.add_option('-g', default=None, dest='hgnc', action='store', help="HGNC ID to Gene Symbol dictionary file")
parser.add_option('-o', default='output', dest='output', action='store', help="Output file name prefix [default value: %default]")
(options, args) = parser.parse_args()

# Welcome message
print '\n' + '='*100
now = str(datetime.datetime.now())
now = now[:now.find('.')]
print 'CARTdb v{} started: {}\n'.format(ver, now)

main(ver, options)

now = str(datetime.datetime.now())
now = now[:now.find('.')]
print '\nFinished: {}'.format(now)
print '='*100 + '\n'