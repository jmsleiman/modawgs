#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2014 Joseph M. Sleiman

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.

from lxml import etree
import models

def importFromSeedXML(seedFileSource):
	'''
	importFromSeedXML(seedFileName):
	seedFileName: a reference to the opened file that will be analyzed (usually "seed.xml" after running Seed).

	The following is the XML structure of Seed, although there are some meta tags that are also used:

	<motifs> : List of all the motifs															0
		<motif id="0"> : List of a particular motif												0 0
		<seq>GGCTCTNNNNAGAGCC</seq> : basepair sequence found in RNA							0 1
		<sec>((((((....))))))</sec> : structural sequence found in RNA							0 2
		<match id="0"> : which string is this?													0 3
			<offset>30</offset> : location from position 0: ie this is the 30th nucleotide		0 3 0
			<energy>0.0</energy> : bch info														0 3 1
			<seq>GGCTCTTTTCAGAGCC</seq> : nucleotide sequence									0 3 2
			<sec>((((((....))))))</sec> : link/structure sequence								0 3 3
		</match>
	This function reads in the XML output from Seed and constructs objects to hold the data.
	'''
	tree = etree.ElementTree(file=seedFileSource)
	motifs = tree.getroot()
	
	# store all the motifs in a list for better access
	
	genericMotifList = list()
	specificMotifList = list()

	for motif in motifs[:-1]:
		currentMotif = models.RNA2DMotif(int(motif.attrib["id"]), motif[0].text, motif[1].text)
		genericMotifList.append(currentMotif)
		
		for x in xrange(2, len(motif)):
			specificMotifList.append(models.RNA2DMatch(currentMotif, int(motif[x].attrib["id"]), motif[x][0].text, motif[x][1].text, motif[x][2].text, motif[x][3].text))
	
	return (genericMotifList, specificMotifList)
	
def importFromDREMEXML(dremeFileName):
	'''
	importFromDREMEXML(dremeFileName):
	dremeFileName: a reference to the opened file that will be analyzed (refered to in this program as "dreme.xml").

	The following is the lxml dump of the structure:
		motif.attrib: {'seq': 'CTCTTTTC', 'evalue': '1.3e-005', 'pvalue': '6.4e-010', 'nsites': '23', 'n': '1', 'p': '23', 'length': '8', 'unerased_evalue': '1.3e-005', 'id': 'm01'}
			position.attrib:
			{'i': '1', 'A': '0.000000', 'C': '1.000000', 'T': '0.000000', 'G': '0.000000'}
			{'i': '2', 'A': '0.000000', 'C': '0.000000', 'T': '1.000000', 'G': '0.000000'}
			{'i': '3', 'A': '0.000000', 'C': '1.000000', 'T': '0.000000', 'G': '0.000000'}
			{'i': '4', 'A': '0.000000', 'C': '0.000000', 'T': '1.000000', 'G': '0.000000'}
			{'i': '5', 'A': '0.000000', 'C': '0.000000', 'T': '1.000000', 'G': '0.000000'}
			{'i': '6', 'A': '0.000000', 'C': '0.000000', 'T': '1.000000', 'G': '0.000000'}
			{'i': '7', 'A': '0.000000', 'C': '0.000000', 'T': '1.000000', 'G': '0.000000'}
			{'i': '8', 'A': '0.000000', 'C': '1.000000', 'T': '0.000000', 'G': '0.000000'}
			match info:  {'p': '23', 'evalue': '1.3e-005', 'pvalue': '6.4e-010', 'seq': 'CTCTTTTC', 'n': '1'}

		motif.attrib: {'seq': 'GAGCCAC', 'evalue': '1.3e-003', 'pvalue': '6.5e-008', 'nsites': '22', 'n': '1', 'p': '20', 'length': '7', 'unerased_evalue': '1.3e-003', 'id': 'm02'}
			position.attrib:
			{'i': '1', 'A': '0.000000', 'C': '0.000000', 'T': '0.000000', 'G': '1.000000'}
			{'i': '2', 'A': '1.000000', 'C': '0.000000', 'T': '0.000000', 'G': '0.000000'}
			{'i': '3', 'A': '0.000000', 'C': '0.000000', 'T': '0.000000', 'G': '1.000000'}
			{'i': '4', 'A': '0.000000', 'C': '1.000000', 'T': '0.000000', 'G': '0.000000'}
			{'i': '5', 'A': '0.000000', 'C': '1.000000', 'T': '0.000000', 'G': '0.000000'}
			{'i': '6', 'A': '1.000000', 'C': '0.000000', 'T': '0.000000', 'G': '0.000000'}
			{'i': '7', 'A': '0.000000', 'C': '1.000000', 'T': '0.000000', 'G': '0.000000'}
			match info:  {'p': '20', 'evalue': '1.3e-003', 'pvalue': '6.5e-008', 'seq': 'GAGCCAC', 'n': '1'}
		
		Which is essentially telling us that on the motif.attrib line, here is the motif's pattern, and on the position.attrib, here is a PSSM to describe the motif.
		
	'''
	tree = etree.ElementTree(file=dremeFileName)
	motifs = tree.getroot()
	
	motifList = list()
	
	for motif in motifs[1]:
		currentMotif = models.RNA1DMotif(motif.attrib["seq"], motif.attrib["evalue"], motif.attrib["pvalue"], motif.attrib["nsites"], motif.attrib["n"], motif.attrib["p"], motif.attrib["length"], motif.attrib["unerased_evalue"], motif.attrib["id"], list())
		motifList.append(currentMotif)
		
		for position in motif:
			if(position.tag == "match"):
				pass
				# skip for now
			else:
				currentMotif.pswm.add(dict([("A", position.attrib["A"]), ("C", position.attrib["C"]), ('T', position.attrib["T"]), ('G', position.attrib["G"])]))
				
	return motifList

'''
importDNAFromFASTA(fasFile):
fasFile: the input strings in FASTA format
The function reads in the FASTA file and splits it into two lists, one with all the labels, one with all the code, and returns a tuple containing those two lists. The user can then easily locate, for example, string 5, which is labels[4] and code[4].
'''

def importDNAFromFASTA(fasFile):
	sourceDataAsList = fasFile.readlines()
	labels = sourceDataAsList[0::2]
	code = sourceDataAsList[1::2]
	dnalist = list()
	
	if(len(labels) != len(code)):
		raise ValueError("There must be one label for each DNA sequence.")
	
	else:
		for x in xrange(0, len(labels)):
			dnalist.append(models.DNASequence(labels[x], code[x]))
	
	return dnalist


def importFromFIMOtxt(filePointer):
	'''
	importFromFIMOtxt(filePointer):
		Imports from FIMO text output. Format is as follows:
			#pattern name	sequence name	start	stop	strand	score	p-value	q-value	matched sequence
		example:
			CTCTTTTC	3MMU004994	33	40	+	15.7143	1.79e-05	0.0297	CTCTTTTC
		entries seperated by tabs
	'''
	dictOfMatches = dict()
	
	entries = filePointer.readlines()
	
	for x in entries[1:]:
		line = x.split('\t')
		try:
			dictOfMatches[line[0]].append(models.RNA1DMatch(line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8]))
		except (KeyError, AttributeError):
			dictOfMatches[line[0]] = [models.RNA1DMatch(line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8])]
	
	return dictOfMatches