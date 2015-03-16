#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2014 Joseph M. Sleiman

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.

import readers
import models
import looker
import copy

'''
Unfortunately this is some kind of hybrid Controller-ViewModel. It will likely be refactored into just a Controller, with most of these classes going into seperate files.
'''


def flatten(*args):
	for x in args:
		if hasattr(x, '__iter__'):
			for y in flatten(*x):
				yield y
		else:
			yield x

class LookerViewModel(object):
	'''
	LookerViewModel, handles all of the requests that come out of the ModawgsWS, often passes them along to the Grouping that is currently being viewed.
	'''

	def __init__(self, modawgsmodel):
		self.model = modawgsmodel.thelooker
		self.controls = None
		# remember to populate these as VMs
		self.listOfGroupingVM = list()
		for entry in self.model.sequenceGroupings:
			self.listOfGroupingVM.append(GroupingViewModel(entry))
		
		self.sequenceList = self.getSequenceList()
		self.currentGroupingIndex = 0
		self.currentGroupingVM = None
		self.currentEnergyCutoff = 1


	'''
	Cool python tricks: decorators and properties
	'''
	@property
	def currentGroupingIndex(self):
		return self.__currentGroupingIndex
	
	@currentGroupingIndex.setter
	def currentGroupingIndex(self, x):
		self.__currentGroupingIndex = x
		self.currentGroupingVM = self.listOfGroupingVM[x]
	
	def getCurrentEnergyGraph(self):
		currentEnergyGraph = list()
		i = 0
		
		print "\t\tLookerVM: Adjusted MFE: Current cuttoff is {0}".format(self.currentEnergyCutoff)
		
		for x in self.currentGroupingVM.energyGraphKeys:
			#print "\t MFE: {0}, amount of matches: {1}".format(x, self.currentGroupingVM.energyGraphData[x])
			if(self.currentEnergyCutoff == -1):
				currentEnergyGraph.append({"mfe" : x, "frequency" : self.currentGroupingVM.energyGraphData[x], "cut" : False})
			elif(i >= self.currentEnergyCutoff):
				currentEnergyGraph.append({"mfe" : x, "frequency" : self.currentGroupingVM.energyGraphData[x], "cut" : True})
			else:
				currentEnergyGraph.append({"mfe" : x, "frequency" : self.currentGroupingVM.energyGraphData[x], "cut" : False})
			i = i + 1
		return currentEnergyGraph
	
	def getSequenceList(self):
		x = list()
		for grouping in self.listOfGroupingVM:
			x.append(grouping.dnaLabel)
		return x
	
	def getControls(self):
		return self.currentGroupingVM.getControls()
	
	def getMotifInfoBlock(self, motifID):
		return None
	
	def energyChange(self):
		x = self.currentGroupingVM.updateEnergyPermissions(self.currentEnergyCutoff - 1)
		return x
	
	def addMatch(self, obj):
		x = self.currentGroupingVM.addMatch(obj["match-type"], int(obj["motif-number"]), int(obj["match-number"]))
		return x
	
	def removeMatch(self, obj):
		x = self.currentGroupingVM.removeMatch(obj["match-type"], int(obj["motif-number"]), int(obj["match-number"]))
		return (x, obj)
	
	def prepareChordData(self):
		try:
			x = self.currentGroupingVM.chordData
			return x
		except Exception as e:
			print e
			return None
		
class Motif1DViewModel(object):
	def __init__(self, motifModel):
		self.model = motifModel
		self.selected = False
		self.valid = True
		self.number = self.model.number
		
		self.infoBlock = "sequence: {0}\ne-val: {1}\np-val: {2}\nnsites: {3}\nn {4}\np: {5}\nlength: {6}\nunerased_evalue: {7}\nMotif ID: {8}\nPSWM {9}".format(self.model.sequence, self.model.evalue, self.model.pvalue, self.model.nsites, self.model.n, self.model.p, self.model.length, self.model.unerased_evalue, self.model.idOfMotif, self.model.pswm)

class Motif2DViewModel(object):
	def __init__(self, motifModel):
		self.model = motifModel
		self.selected = False
		self.valid = True
		self.number = self.model.motifID
		
		self.infoBlock = "seed motifID: {0}\nGeneric Sequence: {1}\nGeneric Structure: {2}".format(self.model.motifID, self.model.genericNucleotideSequence, self.model.genericNucleotideStructure)
	
	
class Match1DViewModel(object):
	def __init__(self, matchModel, num, parentnum):
		self.model = matchModel
		self.model.vm = self
		self.selected = False
		self.valid = True
		self.number = num
		self.motifnumber = parentnum
		self.hitnumber = 0

class Match2DViewModel(object):
	def __init__(self, matchModel, num, parentnum):
		self.model = matchModel
		self.model.vm = self
		self.selected = False
		self.valid = False
		self.number = num
		self.motifnumber = parentnum
		self.hitnumber = 0
		self.flatlist = None

class GroupingViewModel(object):
	def __init__(self, groupingModel):
		self.model = groupingModel
		self.dnaLabel = copy.deepcopy(self.model.DNASequenceStructure.label)
		self.controls = None
		self.listOf1DMotifVM = list()
		self.listOf2DMotifVM = list()
		self.listOf1DMatchVM = list()
		self.listOf2DMatchVM = list()
		
		self.labels = dict()
		self.connections = list()
		self.groupcounter = 0
		self.occupationGrid = [None]*len(self.model.DNASequenceStructure.code)
		
		self.hitlist = list()
		
		self.chordData = None
		'''
		The following blocks of code spawn MotifViewModels for each of the motifs in the models.
		'''
		
		for entry in self.model.listOf1DMotifs:
			x = Motif1DViewModel(entry)
			self.listOf1DMotifVM.append(x)
			n = 0
			for key in entry.matches:
				if(key == self.model.number):
					try:
						for elem in entry.matches[key]:
							self.listOf1DMatchVM.append(Match1DViewModel(elem, n, x.number))
							n = n + 1
					except Exception:
						self.listOf1DMatchVM.append(Match1DViewModel(entry.matches[key], n, x.number))
					
		
		for entry in self.model.listOf2DMotifs:
			x = Motif2DViewModel(entry)
			self.listOf2DMotifVM.append(x)
			n = 0
			for key in entry.matches:
				if(key == self.model.number):
					try:
						for elem in entry.matches[key]:
							self.listOf2DMatchVM.append(Match2DViewModel(elem, n, x.number))
							n = n + 1
					except Exception:
						self.listOf2DMatchVM.append(Match2DViewModel(entry.matches[key], n, x.number))
					
		
		
		self.energyGraphKeys = self.model.energyGraphKeys
		self.energyGraphData = self.model.energyGraphData
	
	@property
	def chordData(self):
		self.labels = dict()
		self.connections = list()
		self.__chordData = {"labels" : self.labels, "connections" : self.connections}
		
		'''
		so we have a hitlist, which contains all the motifs we want to allow
		the format for the chordData is:
		
		connections = [ [{group : unique_id, value : size}, ...] ,  [{group : unique_id, value : size}, ...] , ...]
			Where each entry in connections represents a new segment
			Where each entry in connections[i] represents a possible chord (if it contains many dict), or a spacer (if it contains only one dict)
			Entries cannnot be duplicated.
		
		labels = {unique_id : "label", unique_id : "label", ...}
			Where each unique_id corresponds with one from the connections list
			Where all unique_id form a total order amongst each other
			Entries cannot be duplicated
		'''
		
		i = 0
		spaceCounter = 0
		labelcount = 0
		layaway = set() # we put Terms here, on the assumption we will go back for them once we've found the pair
		reserved = list() # reserved positions. silly D3.
		
		while(i < len(self.occupationGrid)):
			# order doesn't seem to matter for connections
			# it does matter for labels
			if(self.occupationGrid[i] == None):
				spaceCounter += 1
			else:
				if(spaceCounter > 0):
					self.labels[labelcount] = "(spacer)"
					self.connections.append([{"group" : labelcount, "value" : spaceCounter}])
					spaceCounter = 0
					labelcount += 1
				
				# now to put this ID in a set
				if(not(self.occupationGrid[i].termref in layaway)):
					t = ""
					if(self.occupationGrid[i].termref == self.occupationGrid[i].objref):
						t = "M1.{1}.{2}".format(t, self.occupationGrid[i].objref.motifnumber, self.occupationGrid[i].objref.model.matchNumber)
						layaway.add(self.occupationGrid[i].termref)
					else:
						t = "M2.{1}.{2}".format(t, self.occupationGrid[i].objref.motifnumber, self.occupationGrid[i].objref.model.matchNumber)
						layaway.add(self.occupationGrid[i].termref)
					
					reserved.append([self.occupationGrid[i], t, labelcount])
					labelcount += 1
					# M{0}.{1}.{2} ({3}) -- where {0} comes from the class, 
			i += 1
		if(spaceCounter > 0):
			self.labels[labelcount] = "(spacer)"
			self.connections.append([{"group" : labelcount, "value" : spaceCounter}])
			spaceCounter = 0
			labelcount += 1
		
		pending2D = list()
		pendingRefs = list()
		
		'''
		It's a bit complicated, but we're trying to assemble all the match information now that the spacers are out of the way. There aren't many good datastructures to represent ideas like "range 0, 10 are owned by A" that also allow for accessors that behave like lists (ie range[3] returns the same thing as range[4]).
		
		Then again, that does sound feasible.
		'''
		
		for entry in reserved:
			if(entry[0].linear):
				self.connections.append([{"group" : entry[2], "value" : len(entry[0].objref.model)}])
				self.labels[entry[2]] = entry[1]
			else:
				print entry
				
				if(entry[0].termref.pair in pendingRefs):
					self.labels[entry[2]] = entry[1] + "({0})".format(entry[0].termref.sym)
					x = pendingRefs.index(entry[0].termref.pair)
					pendingRefs.pop(x)
					piece = pending2D.pop(x)
					
					mypiece = {"group" : entry[2], "value" : len(entry[0].termref)}
					
					self.connections.append([piece, mypiece])
					
					
				else:
					self.labels[entry[2]] = entry[1] + "({0})".format(entry[0].termref.sym)
					pending2D.append({"group" : entry[2], "value" : len(entry[0].termref)})
					pendingRefs.append(entry[0].termref)
				
				#self.connections.append([{"group" : entry[0][2], "value" : len(entry[0][3].objref.termref)}])
		
		return self.__chordData
		#return {'connections': [[{'group': 0, 'value': 33}], [{'group': 2, 'value': 11}], [{'group': 1, 'value': 13}]],'labels': {0: '(spacer)', 1: 'M1.1.0', 2: '(spacer)'}}
	
	@chordData.setter
	def chordData(self, x):
		self.__chordData = x
	
	def getControls(self):
		'''
		The controls never change, so they're only generated once.
		'''
		if(self.controls == None):
			self.controls = list()
			self.getMotif1DControls()
			self.getMotif2DControls()
			self.getMatch1DControls()
			self.getMatch2DControls()
		
		return self.controls
	
	def getMotif1DControls(self):
		'''
		Each method for controls is seperated, based on what object we're working with. Different motif types have different behavior.
		'''
		for motifVM in self.listOf1DMotifVM:
			d = {"number" : motifVM.number, "type" : "motif_1d", "valid" : motifVM.valid, "selected" : motifVM.selected}
			self.controls.append(d)
		
	def getMotif2DControls(self):
		for motifVM in self.listOf2DMotifVM:
			d = {"number" : motifVM.number, "type" : "motif_2d", "valid" : motifVM.valid, "selected" : motifVM.selected}
			self.controls.append(d)
	
	def getMatch1DControls(self):
		for matchVM in self.listOf1DMatchVM:
			# need to test if they're compatible at all
			if(matchVM.model.sequenceNumber == self.model.number):
				
				d = {"number" : matchVM.model.matchNumber, "type" : "match_1d", "valid" : matchVM.valid, "selected" : matchVM.selected, "parent" : matchVM.motifnumber}
				self.controls.append(d)
	
	def getMatch2DControls(self):
		for matchVM in self.listOf2DMatchVM:
			if(matchVM.model.sequenceNumber == self.model.number):
				
				d = {"number" : matchVM.model.matchNumber, "type" : "match_2d", "valid" : matchVM.valid, "selected" : matchVM.selected, "parent" : matchVM.motifnumber}
			self.controls.append(d)
	
	def updateEnergyPermissions(self, cutoff):
		kl = self.model.energyGraphKeys
		responses = list()
		
		for key in self.model.rankedByEnergy.keys():
			if(key <= kl[cutoff]):
				# if you are within the acceptable rank
				for m in self.model.rankedByEnergy[key]:
					# let's check if you were once removed
					if(m.vm.valid == False):
						# well you are now!
						m.vm.valid = True
						responses.append({"match-number" : m.vm.number, "motif-number" : m.vm.motifnumber, "enabled" : True})
			else:
				for m in self.model.rankedByEnergy[key]:
					# now you must be disabled
					if(m.vm.valid == True):
						# no need sending a False->False
						m.vm.valid = False
						responses.append({"match-number" : m.vm.number, "motif-number" : m.vm.motifnumber, "enabled" : False})
		
		return responses
	
	def addMatch(self, matchtype, motifnumber, matchnumber):
		'''
		This is what happens when you click a checkbox. It's a bit tangly. I don't know why. It is a bit complicated to try and figure out collisions.
		
		In the future, I'd like a way to score compatible and incompatible matches during creation -- but I fear that would be overly complicated. Even scoring it each time a change occurs means n-1 calls per selected matches -- and a user won't just pick one match. It's not terrible, but it's just hard to work out what makes a pair compatible or incompatible. It's true, right now I just check whether they overlap or not, but there are also other concerns.
		
		One thing I didn't like about this was that I'm sure I could have had a compare feature, to compare two motifs and see their compatibility. Feature creep.
		'''
		if(matchtype == "1d"):
			for hit in self.listOf1DMatchVM:
				if(hit.motifnumber == motifnumber):
					for x in xrange(hit.model.start, hit.model.stop):
						if(self.occupationGrid[x] == None):
							hit.selected = True
							self.hitlist.append(hit)
							print "hitlist inbound"
							print self.hitlist
							
							for m in xrange(hit.model.start, hit.model.stop):
								self.occupationGrid[m] = Block(True, hit, hit)
					
							
							return {"selected" : True}
						else:
							hit.selected = False
							
							t = ""
							if(isinstance(self.occupationGrid[x].objref, Match1DViewModel)):
								t = "1d"
							else:
								t = "2d"
							
							return {"selected" : False, "conflicting" : {"match-number" : self.occupationGrid[x].objref.model.matchNumber, "motif-number" : self.occupationGrid[x].objref.motifnumber, "motif-type" : t}, "objective" : {"match-number" : matchnumber, "motif-number" : motifnumber, "motif-type" : matchtype}}
					
					
		elif(matchtype == "2d"):
			for hit in self.listOf2DMatchVM:
				if(hit.motifnumber == motifnumber):
					# test if it's permissible
					start = hit.model.matchPosition
					
					hit.flatlist = list(flatten(hit.model.h5sh3))
					
					# check and see if we can even start here
					if(self.occupationGrid[start] == None):
						# ok, now start looking at the terms
						w = 0
						for term in hit.flatlist:
							if(not term.spacer):
								# if our term isn't just a spacer
								x = 0
								while(x < len(term)):
									# as long as there are spots in the term
									# is the place we're looking at not-empty or a not-spacer?
									if(self.occupationGrid[start+x+w] != None):
										hit.selected = False
										# make sure to deal with complications
										t = ""
										if(isinstance(self.occupationGrid[start+x+w].objref, Match1DViewModel)):
											t = "1d"
										else:
											t = "2d"
										
										return {"selected" : False, "conflicting" : {"match-number" : self.occupationGrid[start+w+x].objref.model.matchNumber, "motif-number" : self.occupationGrid[start+w+x].objref.motifnumber, "motif-type" : t}, "objective" : {"match-number" : matchnumber, "motif-number" : motifnumber, "motif-type" : matchtype}}
									x = x + 1
								w = w + x
							else:
								w = w + len(term)
							
						# at this point, none of the blocking terms have been excluded. there might be an embedded motif but that is OK
						hit.selected = True
						
						# a flatter list means we don't have to worry about infinitely-nested depth traversal. yuck at that. we don't really need to consder it right now. we're just looking at collisions.
						w = 0
						for term in hit.flatlist:
							if(not term.spacer):
								# so if it's essentially we reserve this spot
								x = 0
								while(x < len(term)):
									self.occupationGrid[start+x+w] = Block(True, term, hit)
									x = x + 1
								w = w + x
							else:
								w = w + len(term)
							
						self.hitlist.append(hit)
						return {"selected" : True}
					else:
						return {"selected" : False, "conflicting" : {"match-number" : self.occupationGrid[start].objref.model.matchNumber, "motif-number" : self.occupationGrid[start].objref.motifnumber, "motif-type" : "2d"}, "objective" : {"match-number" : matchnumber, "motif-number" : motifnumber, "motif-type" : matchtype}}
		else:
			return None
	
	def removeMatch(self, matchtype, motifnumber, matchnumber):
		'''
		Returns True when a match can be removed
		Returns False when a match could not be removed
		'''
		try:
			x = 0
			hit = None
			listToWreck = list()
			
			while(x < len(self.hitlist)):
				hit = self.hitlist[x]
				if(hit.motifnumber == motifnumber):
					hit.selected = False
					for b in xrange(0, len(self.occupationGrid)):
						if(self.occupationGrid[b] == None):
							pass
						elif(self.occupationGrid[b].objref == hit):
							self.occupationGrid[b].blocking = False
							self.occupationGrid[b].termref = None
							self.occupationGrid[b].objref = None
							listToWreck.append(b)
							
					break
				x = x+1
			self.hitlist.remove(hit)
			print "did indeed remove that match"
			
			for target in listToWreck:
				self.occupationGrid[target] = None
			
			return True
		except Exception as e:
			print "Exception occured"
			print e
			return False
		

class Block(object):
	def __init__(self, blocking, referenceToTerm, objectReference):
		self.blocking = blocking
		self.termref = referenceToTerm
		self.objref = objectReference
		self.marked = False
		self.linear = False
		if(objectReference == referenceToTerm):
			self.linear = True

class ModawgsModel(object):
	def __init__(self, seedSource, dremeSource, fastaSource, fimoSource):
		seedResults = readers.importFromSeedXML(seedSource)
		dremeResults = readers.importFromDREMEXML(dremeSource)
		fastaResults = readers.importDNAFromFASTA(fastaSource)
		fimoResults = readers.importFromFIMOtxt(fimoSource)
		
		self.thelooker = looker.Looker(fastaResults, dremeResults, fimoResults, seedResults[0])
	
	def spawnVM(self):
		return LookerViewModel(self)