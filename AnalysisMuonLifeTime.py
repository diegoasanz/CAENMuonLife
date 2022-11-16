#!/usr/bin/env python
import os, glob
import shutil
import struct
import subprocess as subp
import sys
import time
from ConfigParser import ConfigParser
from optparse import OptionParser

from collections import OrderedDict

import ROOT as ro
import numpy as np
import cPickle as pickle
import ipdb

from Channel_Caen import Channel_Caen
from Settings_Caen import Settings_Caen
from Utils import *
# from memory_profiler import profile

class AnalysisMuonLifeTime:
	def __init__(self, directory='.', config='CAENAnalysisConfig.cfg'):
		print 'Starting Muon Life Time Analysis ...'

		self.config = config
		self.inDir = directory
		self.inputFile = ''
		self.in_root_file = None
		self.in_tree_name = ''
		self.in_root_tree = None
		self.settings = None
		self.signals_enabled = {}
		self.trigger_ch = None
		self.veto_ch = None
		self.start_event = 0
		self.tot_events = 0
		self.max_events = 0
		self.num_signals = 1
		self.off_time = 4
		self.delay_time = 50
		self.stop_window = 30
		self.veto_window = 50
		self.veto_sigmas = 30
		self.veto_threshold = -10000
		self.trigger_threshold = -10000
		self.signals_sigmas = {}
		self.signals_threshold = {}
		self.max_mem = 8.  # in GB
		self.volt_res = 0.2  # mV

		self.ptsWave, self.event, self.events, self.max_events = 10240, np.zeros(1, 'uint32'), 0, 0
		self.eventVect = np.empty(0, 'uint32')
		self.timeVect, self.signalWaveVects, self.triggerWaveVect, self.vetoWaveVect = np.empty(0, 'int16'), {}, np.empty(0, 'int16'), np.empty(0, 'int16')
		self.vetoedEvents = np.empty(0, 'uint32')
		self.vetoedEventsTrig = np.empty(0, 'uint32')
		self.is_event_vetoed = np.empty(0, '?')
		self.is_event_vetoed2 = np.empty(0, '?')
		self.chStopVect, self.chDecayVect, self.chDecayVect2 = np.empty(0, 'int8'), np.empty(0, 'int8'), np.empty(0, 'int8')
		self.timeStopVect, self.timeDecayVect, self.timeDecayVect2 = np.empty(0, 'int16'), np.empty(0, 'int16'), np.empty(0, 'int16')
		self.lifeTimeVect, self.lifeTimeVect2 = np.empty(0, 'int16'), np.empty(0, 'int16')
		self.timeStop = np.zeros(1, 'int16')
		self.timeDecay = np.zeros(1, 'int16')
		self.timeDecay2 = np.zeros(1, 'int16')
		self.lifeTime = np.zeros(1, 'int16')
		self.lifeTime2 = np.zeros(1, 'int16')
		self.chStop = np.zeros(1, 'int8')
		self.chDecay = np.zeros(1, 'int8')
		self.chDecay2 = np.zeros(1, 'int8')

		self.peds = {}
		self.sigmas = {}

		self.cut0 = ro.TCut('cut0', '')

		self.branchesAll = ['vetoedEvent', 'stopChannel', 'decayChannel', 'timeStop', 'timeDecay', 'lifeTime']
		self.hasBranch = {}

		self.utils = Utils()

		self.outDir = ''

		self.canvas = {}
		self.profile = {}
		self.histo = {}
		self.fit = {}
		self.line = {}

		self.Load_Config_File()
		self.outDir = self.inDir
		self.CreateCut0()
		self.LoadInputTree()

		self.analysisFile = None
		self.analysisTree = None
		self.analysisTreeName = self.in_tree_name + '.analysis'
		self.analysisTreeExisted = False

	def Load_Config_File(self):
		'''
		Loads configuration file given as a parameter. It looks for the different sections "RUN", "VETO" "SIGNALS" and "SIGNALSx" where x
		is the specific signal channel
		:return:
		'''
		parser = ConfigParser()
		if os.path.isfile(self.config):
			print 'Reading configuration file:', self.config, '...', ; sys.stdout.flush()
			parser.read(self.config)

			if parser.has_section('RUN'):
				if parser.has_option('RUN', 'file_name'):
					self.inputFile = parser.get('RUN', 'file_name')
				if parser.has_option('RUN', 'max_mem'):
					self.max_mem = abs(parser.getfloat('RUN', 'max_mem'))

			if parser.has_section('VETO'):
				if parser.has_option('VETO', 'window'):
					self.veto_window = parser.getint('VETO', 'window')
				if parser.has_option('VETO', 'sigmas'):
					self.veto_sigmas = parser.getfloat('VETO', 'sigmas')
				if parser.has_option('VETO', 'threshold') and self.veto_sigmas == 0:
					self.veto_threshold = parser.getfloat('VETO', 'threshold')

			if parser.has_section('SIGNALS'):
				if parser.has_option('SIGNALS', 'number'):
					self.num_signals = parser.getint('SIGNALS', 'number')

			for sig in xrange(self.num_signals):
				if parser.has_section('SIGNAL' + str(sig)):
					self.signals_enabled[sig] = False
					if parser.has_option('SIGNAL' + str(sig), 'enabled'):
						self.signals_enabled[sig] = parser.getboolean('SIGNAL' + str(sig), 'enabled')
					if parser.has_option('SIGNAL' + str(sig), 'sigmas'):
						self.signals_sigmas[sig] = parser.getfloat('SIGNAL' + str(sig), 'sigmas')
					self.signals_threshold[sig] = -10000
					if parser.has_option('SIGNAL' + str(sig), 'threshold'):
						self.signals_threshold[sig] = parser.getfloat('SIGNAL' + str(sig), 'threshold')

			if parser.has_section('CUTS'):
				if parser.has_option('CUTS', 'off_time'):
					self.off_time = parser.getint('CUTS', 'off_time')
				if parser.has_option('CUTS', 'delay'):
					self.delay_time = parser.getint('CUTS', 'delay')
				if parser.has_option('CUTS', 'start_event'):
					self.start_event = parser.getint('CUTS', 'start_event')
				if parser.has_option('CUTS', 'max_events'):
					self.max_events = parser.getint('CUTS', 'max_events')
				if parser.has_option('CUTS', 'stop_window'):
					self.stop_window = parser.getint('CUTS', 'stop_window')

			print 'Done'

	def LoadInputTree(self):
		'''
		Loads the input raw tree from the root file. From the branch "time", it extracts the number of points per saved wave
		:return:
		'''
		if not os.path.isdir(self.inDir):
			ExitMessage('The directory {d} does not exist! Exiting...'.format(d=self.inDir), os.EX_DATAERR)
		if os.path.isfile('{d}/{f}'.format(d=self.inDir, f=self.inputFile)):
			self.in_tree_name = self.inputFile.split('.root')[0].split('.raw')[0]
			self.in_root_file = ro.TFile('{d}/{f}'.format(d=self.inDir, f=self.inputFile))
			self.in_root_tree = self.in_root_file.Get(self.in_tree_name)
			self.ptsWave = int(self.in_root_tree.GetLeaf('time').GetLen())
			self.tot_events = self.in_root_tree.GetEntries()
			self.max_events = self.max_events if self.max_events != 0 else self.tot_events - self.start_event
			# self.AddVetoedFriend()
		else:
			ExitMessage('The file {f} is not inside {d}. Exiting...'.format(d=self.inDir, f=self.in_tree_name), os.EX_DATAERR)

	def AnalysisWaves(self, recreate=False):
		'''
		Method that does all the analysis steps.
		:param recreate: if recreate is True, it will recalculate the vetoed events, rewrite the veto event tree, and reanalyse all the
		events. It will also rewrite the analysis tree. If recreate is False, it will try to load the Veto tree and the analysis tree.
		:return:
		'''
		if not recreate:
			self.OpenAnalysisROOTFile('READ')
		if not self.analysisFile:
			self.OpenAnalysisROOTFile('RECREATE')

		if not np.array(self.hasBranch.values(), '?').all():
			self.CloseAnalysisROOTFile()
			self.DetermineThresholds()
			# ipdb.set_trace()
			self.AddVetoedFriend(overwr=recreate)
			# if not self.in_root_tree.GetFriend('vetoedEventsTree'):
			# 	self.FindVetoedEvents()
			# 	self.AddVetoedFriend()
			self.AnalyseEvents()
			self.OpenAnalysisROOTFile('RECREATE')
			self.FillAnalysisTree()
			self.CloseAnalysisROOTFile()
			self.OpenAnalysisROOTFile('READ')

	def OpenAnalysisROOTFile(self, mode='READ'):
		'''
		Method that opens the Analysis Root file in the specified mode
		:param mode: typical opening modes for root files. E.g. READ, RECREATE
		:return:
		'''
		if not os.path.isdir(self.outDir):
			ExitMessage('The directory {d} does not exist. Exiting...'.format(d=self.outDir), os.EX_DATAERR)

		if self.analysisFile:
			if self.analysisFile.IsOpen():
				if self.analysisFile.GetOption().lower() != mode.lower():
					if self.analysisFile.ReOpen(mode) == -1:
						ExitMessage('Could not reopen file {f}.root in mode: {m}. Exiting...'.format(f=self.analysisTreeName, m=mode))
				return
			else:
				self.analysisFile = None
				self.analysisFile = ro.TFile('{o}/{f}.root'.format(o=self.outDir, f=self.analysisTreeName), mode)
				self.LoadAnalysisTree()
		else:
			if mode.lower() in ['new', 'create', 'update', 'recreate']:
				mode2 = mode if os.path.isfile('{d}/{f}.root'.format(d=self.outDir, f=self.analysisTreeName)) else 'RECREATE'
				self.analysisFile = ro.TFile('{o}/{f}.root'.format(o=self.outDir, f=self.analysisTreeName), mode2)
				self.LoadAnalysisTree()
			else:
				if os.path.isfile('{d}/{f}.root'.format(d=self.outDir, f=self.analysisTreeName)):
					self.analysisFile = ro.TFile('{o}/{f}.root'.format(o=self.outDir, f=self.analysisTreeName), mode)
					self.LoadAnalysisTree()
				else:
					print 'Can\'t open the file {f}.root in {m} mode because it does not exist!!! Returning...'.format(f=self.analysisTreeName, m=mode)

	def LoadAnalysisTree(self):
		'''
		Loads the analysis tree from the analysis root file
		:return:
		'''
		self.analysisTreeExisted = True
		if not self.analysisFile.IsOpen():
			print 'Analysis file is closed. Opening it in update mode...'
			self.OpenAnalysisROOTFile('UPDATE')
		self.analysisFile.cd()
		self.analysisTree = self.analysisFile.Get(self.analysisTreeName)
		if not self.analysisTree:
			self.analysisTreeExisted = False
			self.analysisTree = ro.TTree(self.analysisTreeName, self.analysisTreeName)
		# else:
		# 	if not self.analysisTree.GetFriend(self.in_tree_name):
		# 		self.analysisTree.AddFriend(self.in_tree_name, '{d}/{f}.root'.format(d=self.inDir, f=self.in_tree_name))
		self.hasBranch = {branch: self.TreeHasBranch(branch) for branch in self.branchesAll}

	def TreeHasBranch(self, branch):
		'''
		Checks if the input root tree, or the analysis tree has the specified "branch"
		:param branch: name of the branch that might be in the input root tree or in the analysis root tree
		:return: True if the given branch was found in the input root tree or in the analysis root tree. False if it was not found
		'''
		if self.in_root_tree.GetBranch(branch) or self.analysisTree.GetBranch(branch):
			return True
		return False

	def CreateCut0(self):
		'''
		Creates the initial cut, based on the event limits given in the analysis configuration file
		:return:
		'''
		if self.cut0.GetTitle() != '':
			self.cut0.SetTitle('')
		if self.max_events == 0:
			self.cut0 += ro.TCut('eventsCut', '({s}<=event)'.format(s=self.start_event))
		else:
			self.cut0 += ro.TCut('eventsCut', '(({s}<=event)&&(event<{m}))'.format(s=self.start_event, m=self.start_event + self.max_events))

	def ResetCut0(self):
		'''
		Resets the parameter cut0
		:return:
		'''
		self.cut0.Clear()
		self.cut0 = ro.TCut('cut0', '')

	def DetermineThresholds(self):
		'''
		Determines the threshold for the different signals and the veto scintillator, if they were not given in the analysis configuration file.
		The analysis configuration file could have given the threshold as a number of sigmas below the pedestal, or as a fixed value in mV.
		:return:
		'''
		print 'Determining thresholds ...'
		# create 1D histogram where the pedestal distributions will be saved
		minv, maxv = -1900, 100
		bins = int(float(maxv - minv) / self.volt_res)
		temph = ro.TH1F('temph', 'temph', bins, minv, maxv)
		print '\nVeto Scintillator'
		# if the veto threshold is unrealistic
		if self.veto_threshold < -2000:
			# get the voltage signals for 1000 waveforms in mV for the veto scintillator, 500ns before t=0 (where it should be flat)
			self.in_root_tree.Draw('voltageVeto/10>>temph', '(-500>time)', 'goff', 1000, self.start_event)
			# update the threshold as a number of sigmas below the mean of the pedestal for the veto scintillator
			self.veto_threshold = temph.GetMean() - self.veto_sigmas * temph.GetRMS()
			print 'Veto scintillator, pedestal: {m} mV, RMS: {si} mV'.format(m=temph.GetMean(), si=temph.GetRMS())
		print 'Veto scintillator threshold:', self.veto_threshold, 'mV'
		print '\nTrigger Scintillator'
		# if the trigger threshold is unrealistic
		if self.trigger_threshold < -2000:
			# get the voltage signals for 1000 waveforms in mV for the trigger scintillator, 500ns before t=0 (where it should be flat)
			self.in_root_tree.Draw('voltageTrigger/10>>temph', '(-500>time)', 'goff', 1000, self.start_event)
			# update the threshold as a number of sigmas below the mean of the pedestal for the trigger scintillator
			self.trigger_threshold = temph.GetMean() - self.veto_sigmas * temph.GetRMS()
			print 'Trigger scintillator, pedestal: {m} mV, RMS: {si} mV'.format(m=temph.GetMean(), si=temph.GetRMS())
		print 'Trigger scintillator threshold:', self.trigger_threshold, 'mV'
		# loop over all the signal channels
		for sig in xrange(self.num_signals):
			print '\nSignal', sig
			# check if the channel is enabled for analysis
			if self.signals_enabled[sig]:
				print 'Signal enabled'
				# if the signal's trheshold is unrealistic (below -2V)
				if self.signals_threshold[sig] < -2000:
					# get the voltages for 1000 waveforms in mV for the signal channel, 500ns before t=0 (where it should be flat)
					self.in_root_tree.Draw('voltageSignal{s}/10>>temph'.format(s=sig), '(-500>time)||(time>500)', 'goff', 1000, self.start_event)
					# update the threshold as a number of sigmas below the mean of the pedestal for the signal channel
					self.signals_threshold[sig] = temph.GetMean() - self.signals_sigmas[sig] * temph.GetRMS()
					print 'Signal {s}, pedestal: {m} mV, RMS: {si} mV'.format(s=sig, m=temph.GetMean(), si=temph.GetRMS())
				print 'Signal threshold:', self.signals_threshold[sig], 'mV'
			else:
				print 'Signal disabled'

	def FindVetoedEvents(self):
		'''
		Find the events that have a signal in the veto channel in the time window where the trigger happened.
		This method is used before creating the "vetoed" tree friend.
		:return:
		'''
		print 'Finding vetoed events ...'
		# get list of events where the veto voltage is larger than the threshold in the time region where the trigger happens
		templen = self.in_root_tree.Draw('>>listVeto', '({to}<=time)&&(time<={to}+{ti})&&(voltageVeto/10<{t})'.format(t=self.veto_threshold, ti=self.veto_window, to=self.off_time), '', self.max_events, self.start_event)
		while templen > self.in_root_tree.GetEstimate():
			self.in_root_tree.SetEstimate(templen)
			templen = self.in_root_tree.Draw('>>listVeto', '({to}<=time)&&(time<={to}+{ti})&&(voltageVeto/10<{t})'.format(t=self.veto_threshold, ti=self.veto_window, to=self.off_time), '', self.max_events, self.start_event)
		listVeto = ro.gDirectory.Get('listVeto')
		# make array with the events that were vetoed
		self.vetoedEvents = np.array([ev for ev in range(self.start_event, self.tot_events) if listVeto.Contains(ev)], 'uint32')
		# initialize arrays that account for vetoed events with all entries in "False" for all the events
		self.is_event_vetoed = np.zeros(self.tot_events, '?')
		self.is_event_vetoed2 = np.zeros(self.tot_events, '?')
		# exclude the events before the "start event"
		if self.start_event > 0:
			print 'Excluding first', self.start_event, 'events (tagged as vetoed)'
			# mark the first events before "start event" as vetoed
			self.is_event_vetoed.put(range(self.start_event), True)
		# exclude vetoed events
		if self.vetoedEvents.size > 0:
			print 'Found', self.vetoedEvents.size, 'veto events'
			# mark the identified vetoed events as True in the array
			self.is_event_vetoed.put(self.vetoedEvents, True)

		print 'Finding signals signals passing through the trigger and veto scintillators after the original triggering signal'
		# get list of events where the veto voltage is larger than the threshold in the time region where the decay is supposed to happen
		templen = self.in_root_tree.Draw('>>listVeto2', '({to}+{td}<=time)&&(voltageVeto/10<{t})'.format(t=self.veto_threshold, td=self.delay_time, to=self.off_time), '', self.max_events, self.start_event)
		while templen > self.in_root_tree.GetEstimate():
			self.in_root_tree.SetEstimate(templen)
			templen = self.in_root_tree.Draw('>>listVeto2', '({to}+{td}<=time)&&(voltageVeto/10<{t})'.format(t=self.veto_threshold, td=self.delay_time, to=self.off_time), '', self.max_events, self.start_event)
		listVeto2 = ro.gDirectory.Get('listVeto2')
		# make array with the events that have veto signal in the time region where the decay is supposed to occur
		vetoedEvents2 = np.array([ev for ev in range(self.start_event, self.tot_events) if listVeto2.Contains(ev)], 'uint32')
		# get list of events where the trigger voltage is larger than the threshold in the time region where the decay is supposed to happen
		templen = self.in_root_tree.Draw('>>listTrig2', '({to}+{td}<=time)&&(voltageTrigger/10<{t})'.format(t=self.trigger_threshold, td=self.delay_time, to=self.off_time), '', self.max_events, self.start_event)
		while templen > self.in_root_tree.GetEstimate():
			self.in_root_tree.SetEstimate(templen)
			templen = self.in_root_tree.Draw('>>listTrig2', '({to}+{td}<=time)&&(voltageTrigger/10<{t})'.format(t=self.trigger_threshold, td=self.delay_time, to=self.off_time), '', self.max_events, self.start_event)
		listTrig2 = ro.gDirectory.Get('listTrig2')
		# make array with the events that have trigger signal in the time region where the decay is supposed to occur
		trigEvents2 = np.array([ev for ev in range(self.start_event, self.tot_events) if listTrig2.Contains(ev)], 'uint32')
		# update the array with True for the identified events that had either a veto or trigger signal in the time region where the decay occurs
		if vetoedEvents2.size > 0:
			print 'Found', vetoedEvents2.size, 'events with veto signals in the decay window'
			self.is_event_vetoed2.put(vetoedEvents2, True)
		if trigEvents2.size > 0:
			print 'Found', trigEvents2.size, 'events with trigger signals in the decay window'
			self.is_event_vetoed2.put(trigEvents2, True)

	def AddVetoedFriend(self, overwr=False):
		'''
		Adds the vetoed event tree as a friend of the input root Tree
		:param overwr: if True, it will find the vetoed events and (re)create the root tree.
		If false, it will try to load the vetoed event tree by looking for the file that should contain in. If it does not exist, it will create it
		:return:
		'''
		print 'Adding vetoed friend ...', ; sys.stdout.flush()
		if not self.in_root_tree.GetFriend('vetoedEventsTree'):
			if os.path.isfile('{d}/{f}.vetoed.root'.format(d=self.outDir, f=self.in_tree_name)) and not overwr:
				self.in_root_tree.AddFriend('vetoedEventsTree', '{d}/{f}.vetoed.root'.format(d=self.outDir, f=self.in_tree_name))
				self.hasBranch['vetoedEvent'] = True
				self.cut0 += ro.TCut('vetoedCut', '(!vetoedEvent)')
			elif self.veto_threshold >= -2000 and np.all([self.signals_threshold[sig] >= -2000 for sig in self.signals_enabled.keys() if self.signals_enabled[sig]]):
				self.FindVetoedEvents()
				self.CreateVetoedFriend()
				self.AddVetoedFriend()
			else:
				print 'Can\'t add vetoed friend. Must calculate first the thresholds'
				return
		else:
			if not os.path.isfile('{d}/{f}.vetoed.root'.format(d=self.outDir, f=self.in_tree_name)):
				self.in_root_tree.RemoveFriend('vetoedEventsTree')
				self.AddVetoedFriend(overwr)
		print 'Done'

	def CreateVetoedFriend(self):
		'''
		Creates the vetoedEventsTree that contains the branch vetoedEvent, which is 0 if the event was not vetoed by the veto scintillator,
		or 1 if the event was vetoed. The other variable "vetoedEventDecay" is 1 if the trigger or veto scintillator detect a particle after
		the trigger window, where the decay event is expected, and is 0 if no particle is triggered during the decay time.
		:return:
		'''
		# create root file
		vetoedFile = ro.TFile('{d}/{f}.vetoed.root'.format(d=self.outDir, f=self.in_tree_name), 'RECREATE')
		# create root tree
		vetoedTree = ro.TTree('vetoedEventsTree', 'vetoedEventsTree')
		# create branches
		vetevt = np.zeros(1, '?')
		vetevt2 = np.zeros(1, '?')
		vetoedTree.Branch('vetoedEvent', vetevt, 'vetoedEvent/O')
		vetoedTree.Branch('vetoedEventDecay', vetevt2, 'vetoedEventDecay/O')
		self.utils.CreateProgressBar(self.tot_events + 1)
		self.utils.bar.start()
		# loop over all events
		for ev in xrange(self.tot_events):
			vetevt.itemset(False)
			vetevt2.itemset(False)
			if self.is_event_vetoed[ev]:
				vetevt.itemset(True)
			if self.is_event_vetoed2[ev]:
				vetevt2.itemset(True)
			vetoedTree.Fill()
			self.utils.bar.update(ev + 1)
		vetoedFile.Write()
		vetoedFile.Close()
		self.utils.bar.finish()
		print 'Finished creating vetoedEventsTree'

	def AnalyseEvents(self):
		'''
		Loops over events that are not vetoed (vetoedEvent==0) and finds "stop" events from the muons, and then "decay" events from the muons. The respective
		vectors are updated with the found parameters. If more than one "stop" signal is found it will choose the "earliest" one. The same
		happens with the "decay" signals. Two "decay" variables are set. The first one (without the number at the end 2) only saves the information
		for "decay" events that happened in the same channel where the muon was "stopped". The second set of "decay" variables
		(the ones with a '2' at the end) consider a "decay" event in all the channels.
		:return:
		'''
		print 'Analysing events'
		# get list of events that are not vetoed (vetoedEvent==0)
		leng = self.in_root_tree.Draw('event', self.cut0.GetTitle(), 'goff', self.max_events, self.start_event)
		while leng > self.in_root_tree.GetEstimate():
			self.in_root_tree.SetEstimate(leng)
			leng = self.in_root_tree.Draw('event', self.cut0.GetTitle(), 'goff', self.max_events, self.start_event)
		print 'Getting analysable events'
		tempEvts = self.in_root_tree.GetVal(0)
		self.eventVect = np.array([tempEvts[i] for i in xrange(leng)], 'uint32')
		# initialize arrays with the default values for each event
		self.chStopVect = np.full(self.tot_events, -1, 'int8')
		self.chDecayVect = np.full(self.tot_events, -1, 'int8')
		self.chDecayVect2 = np.full(self.tot_events, -1, 'int8')
		self.timeStopVect = np.full(self.tot_events, -10000, 'int16')
		self.timeDecayVect = np.full(self.tot_events, -10000, 'int16')
		self.timeDecayVect2 = np.full(self.tot_events, -10000, 'int16')
		self.lifeTimeVect = np.full(self.tot_events, -10000, 'int16')
		self.lifeTimeVect2 = np.full(self.tot_events, -10000, 'int16')
		# initialize arrays that will be used to access the data for each event
		self.timeVect = np.zeros(self.ptsWave, 'int16')
		self.signalWaveVects = {sig: np.zeros(self.ptsWave, 'int16') for sig in xrange(self.num_signals)}
		# Enable reading of all variables in the input tree
		self.in_root_tree.SetBranchStatus('*', 1)
		# Disable the access to the branches of the trigger and veto scintillators to speed up the access of the entries
		self.in_root_tree.SetBranchStatus('voltageTrigger', 0)
		self.in_root_tree.SetBranchStatus('voltageVeto', 0)

		onesVect = np.ones(self.ptsWave, '?')
		maxpbar = self.eventVect.size
		itbar = 0
		print 'Looping over analysable events'
		self.utils.CreateProgressBar(maxpbar)
		self.utils.bar.start()
		# loop over events with vetoedEvent==0
		for ev in self.eventVect:
			# temporal variable initialization
			timeStop = -10000
			chStop = -1
			tempTimeStop = -10000
			timeDecay = -10000
			chDecay = -1
			timeDecay2 = -10000
			chDecay2 = -1
			# Acquire the data for the current event
			self.in_root_tree.GetEntry(ev)
			np.putmask(self.timeVect, onesVect, self.in_root_tree.time)
			for sig in xrange(self.num_signals):
				if self.signals_enabled[sig]:
					np.putmask(self.signalWaveVects[sig], onesVect, self.in_root_tree.voltageSignal0 if sig == 0 else self.in_root_tree.voltageSignal1 if sig == 1 else self.in_root_tree.voltageSignal2 if sig == 2 else self.in_root_tree.voltageSignal3)
			# condition that defines the window in time for the stop events
			cond_stop_window = np.bitwise_and(self.off_time < self.timeVect, self.timeVect <= self.off_time + self.stop_window)
			# get the indices for the time array where the stop window is
			timeStopInd = np.argwhere(cond_stop_window).flatten()
			# initialize lists with the information of the signal channels and the times where stop signals are found
			chsStop = []
			chsTimeStop = []
			# find the indices for the falling edge transitions that cross the corresponding threshold for each signal
			fallingEdgePos = {sig: self.GetFallingEdgeTriggers(self.signalWaveVects[sig].astype('int32'), self.signals_threshold[sig]) for sig in xrange(self.num_signals)}
			## Look for stop events
			# loop over each of the signal channels
			for sig in xrange(self.num_signals):
				foundStop = False
				posStop = -1
				# if there is at least one falling edge transition that crosses the threshold
				if fallingEdgePos[sig].size > 0:
					# loop over all the identified falling edge transitions
					for posFE in xrange(fallingEdgePos[sig].size):
						# if one of the transitions is in the stop window
						if fallingEdgePos[sig][posFE] in timeStopInd:
							foundStop = True
							posStop = fallingEdgePos[sig][posFE]
							break
				if foundStop:
					# append to the lists the signal number and the time where the stop was detected
					chsStop.append(sig)
					tempTimeStop = self.timeVect[posStop]
					chsTimeStop.append(tempTimeStop)
			if len(chsStop) > 0:
				# select the earliest stop signal information from the recorded ones
				pos = np.argmin(chsTimeStop)
				timeStop = chsTimeStop[pos]
				chStop = chsStop[pos]
				self.timeStopVect.itemset(ev, timeStop)
				self.chStopVect.itemset(ev, chStop)
				## Look for Decay events
				# initialize lists for decay events
				chsDecay = []
				chsDecay2 = []
				chsTimeDecay = []
				chsTimeDecay2 = []
				# condition that defines the window in time for the decay events
				cond_decay_window = (timeStop + self.delay_time) < self.timeVect
				# get the indices for the time array where the decay window is
				timeDecayInd = np.argwhere(cond_decay_window).flatten()
				# loop over signals
				for sig in xrange(self.num_signals):
					# decay in a channel different from the saved stop event
					# if True:
					if sig != chStop:
						foundDecay = False
						posDecay = -1
						# if there is any falling edge transition that crosses the threshold
						if fallingEdgePos[sig].size > 0:
							# loop over all the identified transitions
							for posFE in xrange(fallingEdgePos[sig].size):
								# check if the transition is in the decay window region
								if fallingEdgePos[sig][posFE] in timeDecayInd:
									foundDecay = True
									posDecay = fallingEdgePos[sig][posFE]
									break
						if foundDecay:
							# append to the lists the signal number and the time where the decay was detected
							chsDecay2.append(sig)
							chsTimeDecay2.append(self.timeVect[posDecay])
				# check decay in the channel where the stop event was found (chStop)
				foundDecay = False
				posDecay = -1
				# check if there is more than one falling edge transition (one was for the identified stop event)
				if fallingEdgePos[chStop].size > 1:
					for posFE in xrange(fallingEdgePos[chStop].size):
						# check if the transition is in the decay window region
						if fallingEdgePos[chStop][posFE] in timeDecayInd:
							foundDecay = True
							posDecay = fallingEdgePos[chStop][posFE]
							break
				if foundDecay:
					# append to the lists the signal number and the time where the decay was detected
					chsDecay.append(chStop)
					chsDecay2.append(chStop)
					chsTimeDecay.append(self.timeVect[posDecay])
					chsTimeDecay2.append(self.timeVect[posDecay])
				# get minimum decay time in array for decay in same stop channel
				if len(chsDecay) > 0:
					pos = np.argmin(chsTimeDecay)
					timeDecay = chsTimeDecay[pos]
					chDecay = chsDecay[pos]
					# update the arrays with the decay information for the current event (for the case where the decay and stop happened in the same signal channel)
					self.timeDecayVect.itemset(ev, timeDecay)
					self.chDecayVect.itemset(ev, chDecay)
					self.lifeTimeVect.itemset(ev, timeDecay - timeStop)
				# get minimum decay time in array for any decay
				if len(chsDecay2) > 0:
					pos = np.argmin(chsTimeDecay2)
					timeDecay2 = chsTimeDecay2[pos]
					chDecay2 = chsDecay2[pos]
					# update the arrays with the decay information for the current event (for any case)
					self.timeDecayVect2.itemset(ev, timeDecay2)
					self.chDecayVect2.itemset(ev, chDecay2)
					self.lifeTimeVect2.itemset(ev, timeDecay2 - timeStop)

			itbar += 1
			self.utils.bar.update(itbar)
		print 'Finished analysing events'

	def GetFallingEdgeTriggers(self, signal, threshold):
		'''
		Finds the indices of events that have a falling edge that passes through a given threshold
		:param signal: array containing the voltage signal. Each unit is in 100uV.
		:param threshold: given threshold to evaluate the falling edge position. it is given in mV
		:return: returns a numpy array with the position of the events that have a falling edge that passes through the given threshold.
		If the signal never crosses the threshold as a falling edge, the returned array is empty.
		'''
		return np.flatnonzero((signal[:-1] >= 10. * threshold) & (signal[1:] <= 10. * threshold)) +1

	def CloseAnalysisROOTFile(self):
		'''
		Closes the Analysis root file and resets the variables for its tree and file
		:return:
		'''
		if self.analysisFile:
			if self.analysisFile.IsOpen():
				self.analysisFile.Close()
			if self.analysisTree:
				del self.analysisTree
		self.analysisTree = None
		self.analysisFile = None

	def CloseInputROOTFiles(self):
		'''
		Closes the input root file and resets the variables for its tree and file
		:return:
		'''
		if self.in_root_file:
			if self.in_root_file.IsOpen():
				self.in_root_file.Close()
			if self.in_root_tree:
				del self.in_root_tree
		self.in_root_file = None
		self.in_root_tree = None

	def FillAnalysisTree(self):
		'''
		Fills the analysis tree
		:return:
		'''
		print 'Filling analysis tree ...'
		# Create branches in analysis tree
		self.analysisTree.Branch('stopChannel', self.chStop, 'stopChannel/B')
		self.analysisTree.Branch('decayChannel', self.chDecay, 'decayChannel/B')
		self.analysisTree.Branch('decayChannel2', self.chDecay, 'decayChannel2/B')
		self.analysisTree.Branch('timeStop', self.timeStop, 'timeStop/S')
		self.analysisTree.Branch('timeDecay', self.timeDecay, 'timeDecay/S')
		self.analysisTree.Branch('timeDecay2', self.timeDecay2, 'timeDecay2/S')
		self.analysisTree.Branch('lifeTime', self.lifeTime, 'lifeTime/S')
		self.analysisTree.Branch('lifeTime2', self.lifeTime2, 'lifeTime2/S')
		self.utils.CreateProgressBar(self.tot_events)
		self.utils.bar.start()
		# loop over all events
		for ev in xrange(self.tot_events):
			self.chStop.itemset(self.chStopVect[ev])
			self.chDecay.itemset(self.chDecayVect[ev])
			self.chDecay2.itemset(self.chDecayVect2[ev])
			self.timeStop.itemset(self.timeStopVect[ev])
			self.timeDecay.itemset(self.timeDecayVect[ev])
			self.timeDecay2.itemset(self.timeDecayVect2[ev])
			self.lifeTime.itemset(self.lifeTimeVect[ev])
			self.lifeTime2.itemset(self.lifeTimeVect2[ev])
			self.analysisTree.Fill()
			self.utils.bar.update(ev + 1)
		# check if the input raw tree is a friend of the analysis tree. If not, make it a friend
		if not self.analysisTree.GetFriend(self.in_tree_name):
			self.analysisTree.AddFriend(self.in_tree_name, '{d}/{f}'.format(d=os.path.abspath(self.inDir), f=self.inputFile))
		# check if the vetoed event tree is a friend of the analysis tree. If not, make it a friend
		if not self.analysisTree.GetFriend('vetoedEventsTree'):
			self.analysisTree.AddFriend('vetoedEventsTree', '{d}/{f}.vetoed.root'.format(d=self.outDir, f=self.in_tree_name))
		# write the filled analysis tree
		self.analysisTree.Write()
		self.utils.bar.finish()

	def SaveCanvasInlist(self, lista):
		'''
		Save the canvas whose keys are given in lista
		:param lista: list of the keys that identify the canvases in the dictionary 'canvas'
		:return:
		'''
		if not os.path.isdir('{d}'.format(d=self.inDir)):
			ExitMessage('The directory does not exist!!!!', os.EX_UNAVAILABLE)
		for canv in lista:
			if self.canvas.has_key(canv):
				self.canvas[canv].SaveAs('{d}/{c}.png'.format(d=self.inDir, c=canv))
				self.canvas[canv].SaveAs('{d}/{c}.pdf'.format(d=self.inDir, c=canv))
				self.canvas[canv].SaveAs('{d}/{c}.root'.format(d=self.inDir, c=canv))

	def SaveAllCanvas(self):
		'''
		Saves all the canvas present in the dictionary 'canvas'
		:return:
		'''
		self.SaveCanvasInlist(self.canvas.keys())

def main():
	parser = OptionParser()
	parser.add_option('-d', '--inDir', dest='inDir', default='.', type='string', help='Directory containing the run files')
	parser.add_option('-c', '--configFile', dest='config', default='CAENAnalysisConfig.cfg', type='string', help='Path to file containing Analysis configuration file')
	parser.add_option('-a', '--automatic', dest='auto', default=False, help='Toggles automatic basic analysis', action='store_true')
	parser.add_option('-r', '--recreate', default=False, help='Recreates analysis file if it was already created', action='store_true')

	(options, args) = parser.parse_args()
	directory = str(options.inDir)
	config = str(options.config)
	autom = bool(options.auto)
	recreate = bool(options.recreate)

	ana = AnalysisMuonLifeTime(directory, config)

	# ana.LoadAnalysisTree()
	# ana.LoadPickles()
	if autom:
		ana.AnalysisWaves(recreate)
		# ana.SaveAllCanvas()
	return ana

	# if auto:
	# 	ccd.StartHVControl()
	# 	ccd.GetBaseLines()
	# 	ccd.SavePickles()
	# 	written_events = ccd.GetData()
	# 	ccd.settings.num_events = written_events
	# 	ccd.SavePickles()  # update pickles with the real amount of written events
	# 	ccd.settings.MoveBinaryFiles()
	# 	ccd.settings.RenameDigitiserSettings()
	# 	ccd.CloseHVClient()
	# 	if not ccd.settings.simultaneous_conversion:
	# 		ccd.CreateRootFile(files_moved=True)
	# 		while ccd.pconv.poll() is None:
	# 			continue
	# 		ccd.CloseSubprocess('converter', stdin=False, stdout=False)
	#
	# print 'Finished :)'
	# sys.stdout.write('\a\a\a')
	# sys.stdout.flush()

if __name__ == '__main__':
	ana = main()
