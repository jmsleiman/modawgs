#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2014 Joseph M. Sleiman

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version. 

from abc import ABCMeta, abstractmethod

class RNA1DMotif(object):
	'''
	RNA1DMotif is used to represent the pattern or motif of a 1D Match. Most of the values are taken straight from DREME, and might require some re-structuring for long-term goals.
	'''
	def __init__(self, sequence, evalue, pvalue, nsites, n, p, length, unerased_evalue, idOfMotif, inPSWM):
		self.sequence = sequence
		self.evalue = evalue
		self.pvalue = pvalue
		self.nsites = nsites
		self.n = n
		self.p = p
		self.length = length
		self.unerased_evalue = unerased_evalue
		self.idOfMotif = idOfMotif
		self.number = int(self.idOfMotif[1:])
		self.pswm = PSWM(inPSWM)
		self.matches = dict()
	
	def store(self, match):
		'''
		store(self, match):
		
		Simply attemps to add a particular match to our list of matches. If multiple matches occured in one sequence, we note that as well.
		'''
		
		''' this requires a big rethink. notably, how are we going to tag matches to motifs... and more importantly: why?'''
		try:
			match.matchNumber = len(self.matches[match.sequenceNumber])
			self.matches[match.sequenceNumber].append(match)
		except (KeyError, AttributeError):
			self.matches[match.sequenceNumber] = [match]
			match.matchNumber = 0
		
	def __str__(self):
		return "rna1dmotif"

class RNA1DMatch(object):
	'''
	RNA1DMatch is used to represent the actual matches found by FIMO. Like RNA1DMotif, it's ripped a bit straight from FIMO, so it might require some love.
	'''
	def __init__(self, patternName, sequenceLabel, start, stop, strand, score, pvalue, qvalue, matchedSequence):
		self.patternName = patternName # references the motif pattern that we match with: eg CCNNCC
		self.sequenceLabel = sequenceLabel # references the DNA label we're paired to
		self.sequenceNumber = 0
		self.vm = None
		self.start = int(start)# the position at which we start
		self.stop = int(stop)# the position at which we stop
		self.strand = strand # top or bottom strand
		self.score = score # how good of a match is this
		self.pvalue = pvalue
		self.qvalue = qvalue
		self.matchedSequence = matchedSequence # what did we match with: eg CCAACC
		self.matchNumber = 0
	
	def __len__(self):
		return len(self.matchedSequence)

class PSWM(object):
	'''
	PSWM is just a simple class for holding PSWMs.
	It was kind of convenient to have a PSWM structure that could print out differently.
	'''
	def __init__(self, listOfDicts):
		self.positions = listOfDicts
	
	def __str__(self):
		s = ""
		
		for entry in self.positions:
			s = s + "A: " + entry["A"] + ",\t" + "T: " + entry["T"] + ",\t" + "G: " + entry["G"] + ",\t" "C: " + entry["C"] + "\n"
		
		return s
	
	def add(self, l):
		self.positions.append(l)
	
class RNA2DMotif(object):
	'''
	RNA2DMotif is the object type that stores the motifs, the patterns found.
	'''
	def __init__(self, motifID, genericNucleotideSequence, genericNucleotideStructure):
		self.motifID = motifID
		self.genericNucleotideSequence = genericNucleotideSequence
		self.genericNucleotideStructure = genericNucleotideStructure
		self.matches = {}
		self.matchesByEnergy = {}
	
	def notify(self, match):
		'''
		notify(self, match):
		
		match: RNA2DMatch object that is linked to this motif, by virtue of being its match. Now the matches can identify their parental motif, and the motifs can identify their instances.
		'''
		self.store(match)
		self.estore(match)
	
	def store(self, match):
		'''
		store(self, match):
		
		Simply attemps to add a particular match to our list of matches. If multiple matches occured in one sequence, we note that as well.
		'''
		try:
			match.matchNumber = len(self.matches[match.sequenceNumber])
			self.matches[match.matchLocation].append(match)
		except (KeyError, AttributeError):
			self.matches[match.matchLocation] = [match]
			match.matchNumber = 0
		
		
	def estore(self, match):
		'''
		estore(self, match):
		
		Simply attemps to add a particular match to our list of matches, based on its energy value. In the event of collision, use a list instead of just the element
		'''
		try:
			self.matchesByEnergy[match.matchEnergy].append(match)
		except (KeyError, AttributeError):
			self.matchesByEnergy[match.matchEnergy] = [match]
		
	def __str__(self):
		return "<<\t" + "motifID: " + str(self.motifID) + "\n\tgeneric nucleotide sequence: " + str(self.genericNucleotideSequence) + "\n\tgeneric nucleotide structure: " + str(self.genericNucleotideStructure) + ">>"
	
class RNA2DMatch(object):
	'''
	RNA2DMatch is the object type that stores the individual matches.
	'''
	def __init__(self, genericMotif, matchLocation, matchPosition, matchEnergy, sequenceID, specificStructure):
		self.matchNumber = None
		self.genericMotif = genericMotif
		self.matchLocation = int(matchLocation)
		self.sequenceNumber = self.matchLocation
		self.matchPosition = int(matchPosition)
		self.start = self.matchPosition
		self.matchEnergy = float(matchEnergy)
		self.sequenceID = sequenceID
		self.specificStructure = specificStructure
		self.h5sh3 = None
		self.parseDotBracket()
		self.genericMotif.notify(self)
		self.vm = None
		# print self.h5sh3
	
	def __str__(self):
		return "<<\tgenericMotif: " + str(self.genericMotif) + "\n\tmatchLocation: " + str(self.matchLocation) + "\n\tmatchPosition: " + str(self.matchPosition) + "\n\tmatchEnergy: " + str(self.matchEnergy) + "\n\tsequenceID: " + str(self.sequenceID) + "\n\tspecificStructure: " + str(self.specificStructure) + "\n>>"
	
	#def __repr__(self):
	#	return self.__str__()

	def parseDotBracket(self):
		'''
		parseDotBracket(self):
		
		This function is used to try and parse the dot bracket notation. It will work with any synthetically designed dot-bracket structure, and should not fail on invalid data, but rather throw up a ValueError. It's assumed input from Seed will be clean.
		'''
		dotcount = 0
		popcount = 0
		position = 0
		stack = list()
		
		while(position != len(self.specificStructure)):
			if (self.specificStructure[position] == '('):
				stack.append(self.specificStructure[position])
				position = position + 1
			
			elif (self.specificStructure[position] == '.'):
				dotcount = 0
				while(position != len(self.specificStructure)):
					if (self.specificStructure[position] == '.'):
						dotcount = dotcount + 1
						position = position + 1
					else:
						break
				stack.append(Spacer(dotcount))
				
			elif(self.specificStructure[position] == ')'):
				top = list()
				while(stack[-1] != '('):
					top.append(stack.pop())
				popcount = 0
				#print "top: "
				#print top
				while(position != len(self.specificStructure)):
					if (self.specificStructure[position] == ')' and stack[-1] == '('):
						#print "can into here"
						stack.pop()
						popcount = popcount + 1
					else:
						break
					position = position + 1
				
				five = H5(popcount)
				three = H3(popcount)
				
				five.pair = three
				three.pair = five
				stack.append(five)
				stack.append(list(reversed(top)))
				stack.append(three)
				
			else:
				raise ValueError("models.py: RNA2DMatch.parseDotBracket -> not a valid symbol for dot-bracket: " + self.specificStructure[position])
		
		self.h5sh3 = stack


class Terms(object):
	'''
	Abstract class for defining Terms.
	'''
	__metaclass__ = ABCMeta
	
	def __init__(self):
		self.length = None
		self.nextSegment = None
		self.startPosition = None
		self.pair = None
	
	@abstractmethod
	def __str__(self): pass
	
	@abstractmethod
	def __len__(self): pass
	
class H5(Terms):
	'''
	H5, H3, and S terms are handled very much alike. Some of these variables are currently unused, but could come in use later on.
	'''
	def __init__(self, length):
		# quite possibly tuples of this data?
		self.length = length
		self.nextSegment = None
		self.startPosition = None
		self.pair = None
		self.spacer = False
		self.sym = "h5"
	
	def __str__(self):
		return "H5: length: {0}".format(self.length)
	
	def __len__(self):
		return self.length
	
class Spacer(Terms):
	'''
	H5, H3, and S terms are handled very much alike. Some of these variables are currently unused, but could come in use later on.
	'''
	def __init__(self, length):
		# quite possibly tuples of this data?
		self.length = length
		self.nextSegment = None
		self.startPosition = None
		self.pair = None
		self.spacer = True
		self.sym = "s"
	
	def __str__(self):
		return "S: length: {0}".format(self.length)
	
	def __len__(self):
		return self.length
	
class H3(Terms):
	'''
	H5, H3, and S terms are handled very much alike. Some of these variables are currently unused, but could come in use later on.
	'''
	def __init__(self, length):
		# quite possibly tuples of this data?
		self.length = length
		self.nextSegment = None
		self.startPosition = None
		self.pair = None
		self.spacer = False
		self.sym = "h3"
	
	def __str__(self):
		return "H3: length: {0}".format(self.length)
	
	def __len__(self):
		return self.length

class DNASequence(object):
	'''
	class DNASequence:

	Used to store each sequence and its label.
	'''
	def __init__(self, label, code):
		self.label = label.strip('\n')
		self.code = code.strip('\n')
	
	def __len__(self):
		return len(self.code)

class Grouping(object):
	'''
	Grouping:
	
	Groups DNA sequences and motifs together.
	'''
	def __init__(self, DNAseq, listOf1DMatches, listOf2DMatches, seqNumber):
		self.DNASequenceStructure = DNAseq
		self.number = seqNumber
		self.listOf1DMatches = listOf1DMatches
		self.listOf2DMatches = listOf2DMatches
		
		
		self.energyGraphData = dict()
		self.energyGraphKeys = list()
		
		self.rankedByEnergy = None # a dictionary containing values in the form of {MFE : [listOfMatches]}
		self.rankedOnStatistics = None
		
		self.listOf1DMotifs = list()
		self.listOf2DMotifs = list()
		
		
	def add1DMotif(self, motif):
		if(motif in self.listOf1DMotifs):
			return False
		else:
			self.listOf1DMotifs.append(motif)
	
	def add2DMotif(self, motif):
		if(motif in self.listOf2DMotifs):
			return False
		else:
			self.listOf2DMotifs.append(motif)
	
	
	def rank2DByEnergy(self):
		'''
		rank2DByEnergy(self):
		
		Creates a ranking system to quickly find matches that have an energy value at a certain rank. The contained structure is simply a dictionary:
		
		{MinimumFreeEnergy : ListOfApplicableMatches}
		
		For filtering out unwanted motifs, it is possible to just check each match for its corresponding motif. For now. :)
		
		'''
		tempdict = {}
		for match in self.listOf2DMatches:
			try:
				len(match) # so did we have multiple matches of the same motif
				for submatch in match: # so we better collect them all
					try:
						if(submatch.sequenceNumber == self.number):
							tempdict[submatch.matchEnergy].append(submatch) 
					except (TypeError, KeyError): # so there isn't a list here yet
						if(submatch.sequenceNumber == self.number):
							tempdict[submatch.matchEnergy] = [submatch]
			except TypeError: # if not
				try:
					if(match.sequenceNumber == self.number):
						tempdict[match.matchEnergy].append(match) 
				except (TypeError, KeyError): # so there isn't a list here yet
					if(match.sequenceNumber == self.number):
						tempdict[match.matchEnergy] = [match]
		
		
		self.rankedByEnergy = tempdict
		
		
		for key in self.rankedByEnergy.keys():
			self.energyGraphData[key] = len(self.rankedByEnergy[key])
			self.energyGraphKeys.append(key)
			
		self.energyGraphKeys = sorted(self.energyGraphKeys)
		
