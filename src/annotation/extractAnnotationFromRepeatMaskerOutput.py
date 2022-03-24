#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""extractAnnotationFromRepeatMaskerOutput.py
DESCRIPTION MISSING
"""

__author__ = "Marc W Schmid"
__version__ = "0"

import argparse
import sys
import logging
import textwrap
import operator
import gzip
import re
logging.basicConfig(format="=== %(levelname)s === %(asctime)s === %(message)s", level=logging.DEBUG, datefmt='%Y-%m-%d %H:%M:%S')

def myopen(fileName, mode="r"):
    """open either a regular or a compressed file"""
    if fileName.endswith(".gz"):
        return gzip.open(fileName, mode=mode)
    else:
        return open(fileName, mode=mode)

def extractAnnotationFromRepeatMaskerOutput(inputFileName):
    outstring = '\t'.join(["queryID", "queryStart", "queryEnd", "queryLeft", "repeatClassOrFamily", "swScore"]) + '\n'
    sys.stdout.write(outstring)
    with myopen(inputFileName) as infile:
        skip = infile.readline()
        skip = infile.readline()
        skip = infile.readline()
        for line in infile:
            fields = line.strip().split()
            (swScore, percDiv, percDel, percIns, queryID, queryStart, queryEnd, queryLeft, something, matchingRepeat, repeatClassOrFamily, repeatStart, repeatEnd, repeatLeft, ID) = fields[0:15]
            queryLeft = re.sub(r"\(|\)", "", queryLeft)
            outstring = '\t'.join([queryID, queryStart, queryEnd, queryLeft, repeatClassOrFamily, swScore]) + '\n'
            sys.stdout.write(outstring)
    pass

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        prog="extractAnnotationFromRepeatMaskerOutput.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""\
            DESCRIPTION MISSING.
            
            The GFF files may contain headers specified with #.
                        
            The resulting table is printed to std::out.
            """))

    parser.add_argument("-v", "--version", action="version",
                        version='%(prog)s {0}'.format(__version__))
    
    parser.add_argument("repeatMaskerOutput", type=str,
                        help="""
                        A '*.out' file from repeatMasker.
                        """)

    args = parser.parse_args()
    
    extractAnnotationFromRepeatMaskerOutput(args.repeatMaskerOutput)
