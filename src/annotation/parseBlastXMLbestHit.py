#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
usage:
python pSofia_parseBlastXMLbestHit.py xmlFile

two files will be written:
.genes (single best hit)
.map (gi to description)
"""

import sys
from Bio.Blast import NCBIXML
import gzip

try:
    xmlFile = sys.argv[1]
except:
    print >> sys.stderr, __doc__
    sys.exit(1)

def myopen(fileName, mode="r"):
    """open either a regular or a compressed file"""
    if fileName.endswith(".gz"):
        return gzip.open(fileName, mode=mode)
    else:
        return open(fileName, mode=mode)

class blastHitFromXML:
    
    def __init__(self, queryName, subjectName, hsp):
        self.q_id = queryName
        self.q_s = hsp.query_start
        self.q_e = hsp.query_end
        self.s_id = subjectName
        self.s_s = hsp.sbjct_start
        self.s_e = hsp.sbjct_end
        self.perc_id = float(hsp.identities)/hsp.align_length
        self.al_len = hsp.align_length
        self.mm = hsp.align_length - hsp.identities
        self.go = hsp.gaps
        self.e_val = hsp.expect
        self.bit_score = hsp.bits
        self.identities = hsp.identities
    
    def __str__(self):
        return '\t'.join([str(self.q_s), str(self.q_e), self.s_id, str(self.s_s), str(self.s_e), str(self.perc_id), str(self.bit_score)])

class queryToHits:
    
    def __init__(self, query, printMode = "all"):
        self.printMode = printMode
        self.query = query
        self.hits = []
        self.hitNamesToIndex = {}
        self.hitNamesToBitScoreSum = {}
        self.hitNamesToIdentitySum = {}
        self.hitNamesToLengthSum = {}
    
    def __str__(self):
        if self.printMode == "oneBest":
            return self.printOneBest()
        elif self.printMode == "allBest":
            return self.printAllBest()
        elif self.printMode == "simpleOneBest":
            return self.printSimpleOneBest()
        elif self.printMode == "simpleAllBest":
            return self.printSimpleAllBest()
        else:
            return self.printAll()

    def printOneBest(self):
        self.sortHits()
        bestSubjects = self.getBestSubjects()
        toPrint = [str(self.hits[i]) for i in self.hitNamesToIndex[bestSubjects[0]]]
        return self.query + '\t' + str(self.hits[0])
        
    def printAllBest(self):
        self.sortHits()
        bestSubjects = self.getBestSubjects()
        toPrint = []
        for bestSubject in bestSubjects:
            toPrint.extend([str(self.hits[i]) for i in self.hitNamesToIndex[bestSubject]])
        return '\n'.join([self.query + '\t' + x for x in toPrint])
        
    def printAll(self):
        self.sortHits()
        toPrint = [str(x) for x in self.hits]
        return '\n'.join([self.query + '\t' + x for x in toPrint])

    def printSimpleOneBest(self):
        self.sortHits()
        bestSubjects = self.getBestSubjects()
        bestSubject = bestSubjects[0]
        return '\t'.join([self.query, bestSubject, str(self.hitNamesToBitScoreSum[bestSubject])])

    def printSimpleAllBest(self):
        self.sortHits()
        bestSubjects = self.getBestSubjects()
        toPrint = ['\t'.join([self.query, bestSubject, str(self.hitNamesToBitScoreSum[bestSubject])]) for bestSubject in bestSubjects]
        return '\n'.join(toPrint)

    def sortHits(self):
        self.hits.sort(lambda x,y: cmp(y.bit_score, x.bit_score))
        self.hitNamesToIndex = {}
        for i, x in enumerate(self.hits):
            try:
                self.hitNamesToIndex[x.s_id].append(i)
            except KeyError:
                self.hitNamesToIndex[x.s_id] = [i]
    
    def getBestSubjects(self):
        maxScore = max([val for (key, val) in self.hitNamesToBitScoreSum.items()])
        bestSubjects = [key for (key, val) in self.hitNamesToBitScoreSum.items() if val == maxScore]
        return bestSubjects
        
    def add(self, hit):
        try:
            self.hitNamesToIndex[hit.s_id].append(len(self.hits))
            self.hitNamesToBitScoreSum[hit.s_id] += hit.bit_score
            self.hitNamesToIdentitySum[hit.s_id] += hit.identities
            self.hitNamesToLengthSum[hit.s_id] += hit.al_len
        except KeyError:
            self.hitNamesToIndex[hit.s_id] = [len(self.hits)]
            self.hitNamesToBitScoreSum[hit.s_id] = hit.bit_score
            self.hitNamesToIdentitySum[hit.s_id] = hit.identities
            self.hitNamesToLengthSum[hit.s_id] = hit.al_len
        self.hits.append(hit)
    
class queryCollection:
    
    def __init__(self, printMode = "all"):
        self.data = {}
        self.printMode = printMode

    def __str__(self):
        out = [str(val) for (key, val) in self.data.items()]
        return '\n'.join(out)
    
    def add(self, hit):
        try:
            self.data[hit.q_id].add(hit)
        except KeyError:
            self.data[hit.q_id] = queryToHits(hit.q_id, self.printMode)
            self.data[hit.q_id].add(hit)

def parseXML(xmlFile):
    entryCounter = 0
    genesOut = queryCollection("simpleAllBest")
    mapOut = {}
    
    with myopen(xmlFile) as infile:
        entries = NCBIXML.parse(infile)
        for e in entries:
            entryCounter += 1
            if (entryCounter % 10000) == 0:
                print >> sys.stderr, "processed %d entries" % entryCounter
            queryName = e.query #split...
            for al in e.alignments:
                # this were thomas files
                #fields = al.hit_def.split('|')
                #subjectName = fields[1] # that's the gi number
                #subjectDescription = fields[4].strip()
                subjectName = al.hit_id
                fields = al.hit_def.split('|') # some have MANY species in the description, just take the first and note down how many others
                if len(fields) > 1:
                    subjectDescription = ''.join([fields[0], " and ", str(len(fields)-1), " others"])
                else:
                    subjectDescription = fields[0]
                mapOut[subjectName] = subjectDescription
                for hsp in al.hsps:
                    hit = blastHitFromXML(queryName, subjectName, hsp)
                    genesOut.add(hit)
    return (genesOut, mapOut)

(parsedData, mapData) = parseXML(xmlFile)
mapFile = xmlFile.replace(".gz", "") + ".map"
genesFile = xmlFile.replace(".gz", "") + ".genes"
with open(genesFile, 'wb') as outfile:
    print >> outfile, parsedData
with open(mapFile, 'wb') as outfile:
    for (key, val) in mapData.items():
        print >> outfile, '\t'.join([key, val])

