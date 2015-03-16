#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2014 Joseph M. Sleiman

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.

import models

'''
The Looker is in charge of providing some analysis and clarity to the data. It is used as an object so that it can store any and all references to the required data.

The Looker focuses on ALL sequences. The individual sequences with their matches and motifs are Groupings objects.

Looker should contain a structure of Grouping objects. Rather, it should contain a dictionary, where they're tagged via sequence labels (or each Grouping is assigned a position, and they are reffered as such internally?).

{
	"sequence_label_1" : models.Grouping,
	"sequence_label_2" : models.Grouping,
	etc
}

This way the Looker can contain just basic methods to run statistics on each Grouping, but does not perform those statistics itself. In a way, it would act as a sort of intermediary between the Python core and Web UI.

As we can't perform many analytics until we actually try to communicate with the user, this will have to do.

The Looker at this point merely prepares the sequences for analysis.

''' 

class Looker:
	def __init__(self, fasta, motif1dlist, match1d, motif2dlist):
		'''
		__init__(self, fasta, motif1dlist, motif2dlist):
		
		fasta: fasta list
		
		motif2dlist: a list containing all of the 2d motifs (not the matches)
		'''
		# since the motifs now contain references to the matches, this is possible.
		self.fasta = fasta
		self.motif1d = motif1dlist
		self.motif2d = motif2dlist
		
		self.groupMotif1D(match1d)
		
		self.sequenceGroupings = self.groupBySequence()
		
	def groupMotif1D(self, match1d):
		'''
		simple enough to handle:
		since each entry in self.match1d is tagged with the corresponding motif pattern, sort the motif1dlist, and then attach them'''
		 
		tempdict = dict()
		 
		for entry in self.motif1d:
			tempdict[entry.sequence] = entry
		
		for k,v in match1d.items():
			
			for match in v:
				dna = 0
				
				for sequence in self.fasta:
					# ok this is probably really sloppy because we don't know how labels were stripped in fimo
					if(match.sequenceLabel in sequence.label):
						match.sequenceNumber = self.fasta.index(sequence)
						tempdict[k].store(match)
						break
				
		
	def groupBySequence(self):
		'''
		groupBySequence(self):
			This function attempts to group all matches based on the sequences in which they occured.
		'''
		tempgroups = list()
		
		for x in xrange(0, len(self.fasta)):
			tempgroups.append(models.Grouping(self.fasta[x], list(), list(), x))
		
		for motif in self.motif1d:
			# essentially we can find out which of the matches
			for key in motif.matches:
				tempgroups[key].listOf1DMatches.append(motif.matches[key][0])
				tempgroups[key].add1DMotif(motif)
		
		for motif in self.motif2d:
			for key in motif.matches:
				tempgroups[key].listOf2DMatches.append(motif.matches[key][0])
				tempgroups[key].add2DMotif(motif)
		
		for grouping in tempgroups:
			grouping.rank2DByEnergy()
		
		return tempgroups
	