import sys
#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__ = "pmoreno"
__date__ = "$May 30, 2011 7:22:44 PM$"

if __name__ == "__main__":
    #print "Hello World"
    hmmscanRes = open(sys.argv[1], 'r')

    line = hmmscanRes.readline()

    while not line.startswith("Scores for complete sequence"):
        line = hmmscanRes.readline()

    line = hmmscanRes.readline()
    line = hmmscanRes.readline()
    line = hmmscanRes.readline()
    line = hmmscanRes.readline()
    import re

    # regRes = re.search('\s+(\S+)\s+(\S+)\s+(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)', line)
    # evalue = regRes.group
    tokens = re.split('\s+', line)
    evalue = tokens[1]
    score = tokens[2]
    model = tokens[9]

    print model+"\t"+evalue+"\t"+score

