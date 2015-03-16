'''
The MIT License (MIT)

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
'''
import signal, sys, ssl, logging
import traceback
from SimpleWebSocketServer import WebSocket, SimpleWebSocketServer, SimpleSSLWebSocketServer
from optparse import OptionParser
import json

import coordinator

logging.basicConfig(format='%(asctime)s %(message)s', level=logging.DEBUG)

class ModawgsWS(WebSocket):
	'''
	structure of messages:
	{
	action: "string"
	data: ["entry1", "entry 2", etc]
	}
	
	But what kind of data are we sending to and fro?
	
	probably just a json array of values
	
	'''
	
	'''
	now the biggest problem is: how do i communicate back and forth with a backend core?
	
	every time we spawn a connection, spawn a new view model, but retain all the old models. that way 1 server = 1 data set, but 1 window = 1 view
	
	'''
	
	def __init__(self, server, sock, address, modawgsmodel):
		self.vm = modawgsmodel.spawnVM()
		super(ModawgsWS, self).__init__(server, sock, address)
	
	def handleMessage(self):
		
		
		if self.data is None:
			self.data = ''
		
		try:
			
			msg = json.loads(self.data.decode("utf-8"))
			print "\t\t", self.data.decode("utf-8")
			
			print "message accepted: {0}".format(msg)
			if(msg["action"] == "request-change-sequence"):
				self.vm.currentGroupingIndex = int(msg["data"])
				self.sendEnergyGraph()
				self.sendUpdatedControls()
				self.sendEnergyMatchChanges()
				
			elif(msg["action"] == "request-energy-graph"):
				self.sendEnergyGraph()
			
			elif(msg["action"] == "request-update-controls"):
				self.sendUpdatedControls()
			
			elif(msg["action"] == "request-set-energy-cutoff"):
				self.vm.currentEnergyCutoff = int(msg["data"])
				self.sendEnergyMatchChanges()
				self.sendEnergyGraph()
			
			elif(msg["action"] == "request-add-match"):
				x = self.vm.addMatch(msg["data"])
				# not exactly a useful message
				msg = {"action" : "response-add-match", "data" : x}
				print msg
				
				# {"selected" : False, "conflicting" : {"match-number" : self.occupationGrid[x].model.matchNumber, "motif-number" : self.occupationGrid[x].model.motifnumber, "motif-type" : "2d"}}
				
				self.sendMessage(json.dumps(msg))
				
			elif(msg["action"] == "request-remove-match"):
				x = self.vm.removeMatch(msg["data"])
				# not exactly a useful message
				msg = {"action" : "response-remove-match", "data" : {"enabled" : x[0], "entry" : x[1]}}
				self.sendMessage(json.dumps(msg))
			
			elif(msg["action"] == "request-chord-data"):
				print "preparing to send the new chord"
				self.sendChordData()
			
			else:
				print "error: received unknown message: {0}".format(msg)
		
		except Exception as n:
			print "Exception occured: {0}, {1}".format(n, sys.exc_traceback.tb_lineno)
	
	# ====== data tx/rx ================================================================== #
	def sendEnergyMatchChanges(self):
		x = self.vm.energyChange()
		msg = {"action" : "response-update-matches", "data" : x}
		self.sendMessage(json.dumps(msg))
	
	def sendChordData(self):
		x = self.vm.prepareChordData()
		msg = {"action" : "response-chord-data", "data" : x}
		print x
		self.sendMessage(json.dumps(msg))
	
	def sendEnergyGraph(self):
		print "ModawgsWS -> sendEnergyGraph()"
		x = self.vm.getCurrentEnergyGraph()
		msg = {"action" : "response-energy-graph", "confirmation" : "true", "data" : x}
		self.sendMessage(json.dumps(msg))
	
	def sendUpdatedControls(self):
		print "ModawgsWS -> sendUpdatedControls()"
		try:
			x = self.vm.getControls()
			msg = {"action" : "response-initial-controls", "data" : x}
		except Exception as e:
			print "oh man wtf"
			print e
		self.sendMessage(json.dumps(msg))
	'''
	Connection handlers
	'''
	def handleConnected(self):
		print self.address, 'connected'
		sl = self.vm.getSequenceList()
		msg = {"action" : "response-sequence-list", "confirmation" : "true", "data" : sl}
		self.sendMessage(json.dumps(msg))
		self.vm.currentGroupingIndex = 0
		self.sendEnergyGraph()
		self.sendUpdatedControls()
	
	def handleClose(self):
		self.x += 1
		print "value of x was: ", self.x
		print self.address, 'closed'
	
	
if __name__ == "__main__":
	port = 8000
	host = ""
	
	print """
#     #         ######     #    #     #  #####   #####
##   ##   ####  #     #   # #   #  #  # #     # #     #
# # # #  #    # #     #  #   #  #  #  # #       #
#  #  #  #    # #     # #     # #  #  # #  ####  #####
#     #  #    # #     # ####### #  #  # #     #       #
#     #  #    # #     # #     # #  #  # #     # #     #
#     #   ####  ######  #     #  ## ##   #####   #####"""
	print "Modawgs v1.0b"
	print "At the moment, you'll need to run your own FASTA file in Seed, DREME, and FIMO seperately. Sorry!"

	dremeSource = None
	seedSource = None
	fastaSource = None
	fimoSource = None
	
	seedSource = file("seed.xml", "r")
	dremeSource = file("dreme.xml", "r")
	fastaSource = file("data.fas", "r")
	fimoSource = file("fimo.txt", "r")
	
	#while(True):
		#print "Please provide a file path for Seed's XML output: "
		#seed = raw_input("filename: " )
		#try:
			#seedSource = file(seed, "r")
			#break
		#except Exception as e:
			#print "Sorry, {0} came up. Try again.".format(e)
	
	#while(True):
		#print "Please provide a file path for Dreme's XML output: "
		#dreme = raw_input("filename: " )
		#try:
			#dremeSource = file(dreme, "r")
			#break
		#except Exception as e:
			#print "Sorry, {0} came up. Try again.".format(e)
	
	#while(True):
		#print "Please provide a file path for FIMO's txt output: "
		#fimo = raw_input("filename: " )
		#try:
			#fimoSource = file(fimo, "r")
			#break
		#except Exception as e:
			#print "Sorry, {0} came up. Try again.".format(e)
	
	#while(True):
		#print "Please provide a file path for your FASTA output: "
		#fasta = raw_input("filename: " )
		#try:
			#fastaSource = file(fasta, "r")
			#break
		#except Exception as e:
			#print "Sorry, {0} came up. Try again.".format(e)
	
	print "Thanks! When you're ready, open the index.html file in the Viewer folder to autoconnect to the viewer."
	print "If you've already loaded it, just refresh the page or hit 'connect'."
	cls = ModawgsWS
	
	mainmodel = coordinator.ModawgsModel(seedSource, dremeSource, fastaSource, fimoSource)
	server = SimpleWebSocketServer(host, port, cls)
	server.modawgs = mainmodel
	

	def close_sig_handler(signal, frame):
		server.close()
		sys.exit()

	signal.signal(signal.SIGINT, close_sig_handler)
	server.serveforever()
