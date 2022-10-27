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
from Langaus import LanGaus
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
		self.veto_sigmas = 0
		self.veto_threshold = -10000
		self.trigger_threshold = -10000
		self.signals_sigmas = {}
		self.signals_threshold = {}
		self.max_mem = 8.  # in GB
		self.volt_res = 0.25  # mV

		self.ptsWave, self.event, self.events, self.max_events = 10240, np.zeros(1, 'uint32'), 0, 0
		self.eventVect = np.empty(0, 'uint32')
		self.timeVect, self.signalWaveVects, self.triggerWaveVect, self.vetoWaveVect = np.empty(0, 'int16'), {}, np.empty(0, 'int16'), np.empty(0, 'int16')
		self.vetoedEvents = np.empty(0, 'uint32')
		self.vetoedEventsTrig = np.empty(0, 'uint32')
		self.is_event_vetoed = np.empty(0, '?')
		self.is_event_vetoed_trig = np.empty(0, '?')
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

	def AnalysisWaves(self):
		self.OpenAnalysisROOTFile('READ')
		if not self.analysisFile:
			self.OpenAnalysisROOTFile('RECREATE')

		if not np.array(self.hasBranch.values(), '?').all():
			self.CloseAnalysisROOTFile()
			self.DetermineThresholds()
			# ipdb.set_trace()
			self.AddVetoedFriend()
			# if not self.in_root_tree.GetFriend('vetoedEventsTree'):
			# 	self.FindVetoedEvents()
			# 	self.AddVetoedFriend()
			self.AnalyseEvents()
			self.OpenAnalysisROOTFile('RECREATE')
			self.FillAnalysisTree()
			self.CloseAnalysisROOTFile()
			self.OpenAnalysisROOTFile('READ')

	def OpenAnalysisROOTFile(self, mode='READ'):
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
		if self.in_root_tree.GetBranch(branch) or self.analysisTree.GetBranch(branch):
			return True
		return False

	def CreateCut0(self):
		if self.cut0.GetTitle() != '':
			self.cut0.SetTitle('')
		if self.max_events == 0:
			self.cut0 += ro.TCut('eventsCut', '({s}<=event)'.format(s=self.start_event))
		else:
			self.cut0 += ro.TCut('eventsCut', '(({s}<=event)&&(event<{m}))'.format(s=self.start_event, m=self.start_event + self.max_events))

	def ResetCut0(self):
		self.cut0.Clear()
		self.cut0 = ro.TCut('cut0', '')

	def DetermineThresholds(self):
		print 'Determining thresholds ...'
		minv, maxv = -1900, 100
		bins = int(float(maxv - minv) / self.volt_res)
		temph = ro.TH1F('temph', 'temph', bins, minv, maxv)
		print '\nVeto Scintillator'
		if self.veto_threshold < -2000:
			self.in_root_tree.Draw('voltageVeto/10>>temph', '(-500>time)||(time>500)', 'goff', 10, self.start_event)
			self.veto_threshold = temph.GetMean() - self.veto_sigmas * temph.GetRMS()
			print 'Veto scintillator, pedestal: {m} mV, RMS: {si} mV'.format(m=temph.GetMean(), si=temph.GetRMS())
		print 'Veto scintillator threshold:', self.veto_threshold, 'mV'
		print '\nTrigger Scintillator'
		if self.trigger_threshold < -2000:
			self.in_root_tree.Draw('voltageTrigger/10>>temph', '(-500>time)||(time>500)', 'goff', 10, self.start_event)
			self.trigger_threshold = temph.GetMean() - self.veto_sigmas * temph.GetRMS()
			print 'Trigger scintillator, pedestal: {m} mV, RMS: {si} mV'.format(m=temph.GetMean(), si=temph.GetRMS())
		print 'Trigger scintillator threshold:', self.trigger_threshold, 'mV'
		for sig in xrange(self.num_signals):
			print '\nSignal', sig
			if self.signals_enabled[sig]:
				print 'Signal enabled'
				if self.signals_threshold[sig] < -2000:
					self.in_root_tree.Draw('voltageSignal{s}/10>>temph'.format(s=sig), '(-500>time)||(time>500)', 'goff', 10, self.start_event)
					self.signals_threshold[sig] = temph.GetMean() - self.signals_sigmas[sig] * temph.GetRMS()
					print 'Signal {s}, pedestal: {m} mV, RMS: {si} mV'.format(s=sig, m=temph.GetMean(), si=temph.GetRMS())
				print 'Signal threshold:', self.signals_threshold[sig], 'mV'
			else:
				print 'Signal disabled'

	def FindVetoedEvents(self):
		print 'Finding vetoed events ...'
		templen = self.in_root_tree.Draw('>>listVeto', '({to}<=time)&&(time<={to}+{ti})&&(voltageVeto/10<{t})'.format(t=self.veto_threshold, ti=self.veto_window, to=self.off_time), '', self.max_events, self.start_event)
		# templen = self.in_root_tree.Draw('>>listVeto', '(0<=time)&&(time<={ti})&&(voltageVeto/10<{t})'.format(t=self.veto_threshold, ti=self.off_time + self.stop_window), '', self.max_events, self.start_event)
		while templen > self.in_root_tree.GetEstimate():
			self.in_root_tree.SetEstimate(templen)
			templen = self.in_root_tree.Draw('>>listVeto', '({to}<=time)&&(time<={to}+{ti})&&(voltageVeto/10<{t})'.format(t=self.veto_threshold, ti=self.veto_window, to=self.off_time), '', self.max_events, self.start_event)
		# templen = self.in_root_tree.Draw('>>listVeto', '(0<=time)&&(time<={ti})&&(voltageVeto/10<{t})'.format(t=self.veto_threshold, ti=self.off_time + self.stop_window), '', self.max_events, self.start_event)
		listVeto = ro.gDirectory.Get('listVeto')
		self.vetoedEvents = np.array([ev for ev in range(self.start_event, self.tot_events) if listVeto.Contains(ev)], 'uint32')
		self.is_event_vetoed = np.zeros(self.tot_events, '?')
		if self.start_event > 0:
			print 'Excluding first', self.start_event, 'events (tagged as vetoed)'
			self.is_event_vetoed.put(range(self.start_event), True)
		if self.vetoedEvents.size > 0:
			print 'Found', self.vetoedEvents.size, 'veto events'
			self.is_event_vetoed.put(self.vetoedEvents, True)

		# print 'Finding second triggers ...'
		# templen = self.in_root_tree.Draw('>>listVetoTrig', '(10<=time)&&(time<={ti})&&(voltageVeto/10<{t})'.format(t=self.veto_threshold, ti=self.veto_window), '', self.max_events, self.start_event)
		# while templen > self.in_root_tree.GetEstimate():
		# 	self.in_root_tree.SetEstimate(templen)
		# 	templen = self.in_root_tree.Draw('>>listVetoTrig', '(10<=time)&&(time<={ti})&&(voltageVeto/10<{t})'.format(t=self.veto_threshold, ti=self.veto_window), '', self.max_events, self.start_event)
		# listVetoTrig = ro.gDirectory.Get('listVetoTrig')
		# self.vetoedEventsTrig = np.array([ev for ev in range(self.start_event, self.tot_events) if listVetoTrig.Contains(ev)], 'uint32')
		# self.is_event_vetoed_trig = np.zeros(self.tot_events, '?')
		# if self.start_event > 0:
		# 	print 'Excluding first', self.start_event, 'events (tagged as vetoed trigger)'
		# 	self.is_event_vetoed_trig.put(range(self.start_event), True)
		# if self.vetoedEventsTrig.size > 0:
		# 	print 'Found', self.vetoedEventsTrig.size, 'veto trigger events'
		# 	self.is_event_vetoed_trig.put(self.vetoedEventsTrig, True)

	def AddVetoedFriend(self, overwr=False):
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
		vetoedFile = ro.TFile('{d}/{f}.vetoed.root'.format(d=self.outDir, f=self.in_tree_name), 'RECREATE')
		vetoedTree = ro.TTree('vetoedEventsTree', 'vetoedEventsTree')
		vetevt = np.zeros(1, '?')
		# vetevtTrig = np.zeros(1, '?')
		vetoedTree.Branch('vetoedEvent', vetevt, 'vetoedEvent/O')
		# vetoedTree.Branch('vetoedEventTrig', vetevtTrig, 'vetoedEventTrig/O')
		self.utils.CreateProgressBar(self.tot_events + 1)
		self.utils.bar.start()
		for ev in xrange(self.tot_events):
			vetevt.itemset(False)
			# vetevtTrig.itemset(False)
			if self.is_event_vetoed[ev]:
				vetevt.itemset(True)
			# if self.is_event_vetoed_trig[ev]:
			# 	vetevtTrig.itemset(True)
			vetoedTree.Fill()
			self.utils.bar.update(ev + 1)
		vetoedFile.Write()
		vetoedFile.Close()
		self.utils.bar.finish()
		print 'Finished creating vetoedEventsTree'

	def AnalyseEvents(self):
		print 'Analysing events'
		leng = self.in_root_tree.Draw('event', self.cut0.GetTitle(), 'goff', self.max_events, self.start_event)
		while leng > self.in_root_tree.GetEstimate():
			self.in_root_tree.SetEstimate(leng)
			leng = self.in_root_tree.Draw('event', self.cut0.GetTitle(), 'goff', self.max_events, self.start_event)
		print 'Getting analysable events'
		tempEvts = self.in_root_tree.GetVal(0)
		self.eventVect = np.array([tempEvts[i] for i in xrange(leng)], 'uint32')
		self.chStopVect = np.full(self.tot_events, -1, 'int8')
		self.chDecayVect = np.full(self.tot_events, -1, 'int8')
		self.chDecayVect2 = np.full(self.tot_events, -1, 'int8')
		self.timeStopVect = np.full(self.tot_events, -10000, 'int16')
		self.timeDecayVect = np.full(self.tot_events, -10000, 'int16')
		self.timeDecayVect2 = np.full(self.tot_events, -10000, 'int16')
		self.in_root_tree.SetBranchStatus('*', 1)
		self.in_root_tree.SetBranchStatus('voltageTrigger', 0)
		self.in_root_tree.SetBranchStatus('voltageVeto', 0)
		self.timeVect = np.zeros(self.ptsWave, 'int16')
		self.signalWaveVects = {sig: np.zeros(self.ptsWave, 'int16') for sig in xrange(self.num_signals)}
		onesVect = np.ones(self.ptsWave, '?')
		maxpbar = self.eventVect.size
		itbar = 0
		print 'Looping over analysable events'
		self.utils.CreateProgressBar(maxpbar)
		self.utils.bar.start()
		cont = 0
		for ev in self.eventVect:
			self.in_root_tree.GetEntry(ev)
			np.putmask(self.timeVect, onesVect, self.in_root_tree.time)
			for sig in xrange(self.num_signals):
				if self.signals_enabled[sig]:
					np.putmask(self.signalWaveVects[sig], onesVect, self.in_root_tree.voltageSignal0 if sig == 0 else self.in_root_tree.voltageSignal1 if sig == 1 else self.in_root_tree.voltageSignal2 if sig == 2 else self.in_root_tree.voltageSignal3)
			cond_stop_window = np.bitwise_and(self.off_time < self.timeVect, self.timeVect <= self.off_time + self.stop_window)
			timeStopInd = np.argwhere(cond_stop_window).flatten()
			chsStop = []
			chsTimeStop = []
			# chsTimeStopEnd = []
			cond_sig_threshold = {sig: self.signalWaveVects[sig] < 10 * self.signals_threshold[sig] for sig in xrange(self.num_signals)}
			for sig in xrange(self.num_signals):
				cond_sig_threshold_i = cond_sig_threshold[sig]
				foundStop = cond_sig_threshold_i[timeStopInd].any()
				if foundStop:
					chsStop.append(sig)
					cond_time_stop = np.bitwise_and(cond_stop_window, cond_sig_threshold_i)
					tempTimeStop = self.timeVect[cond_time_stop.argmax()]
					chsTimeStop.append(tempTimeStop)
					# cond_stop_end_window = self.timeVect > tempTimeStop
					# cond_time_stop_end = np.bitwise_and(cond_stop_end_window, np.bitwise_not(cond_sig_threshold_i))
					# tempTimeStopEnd = self.timeVect[cond_time_stop_end.argmax()]
					# chsTimeStopEnd.append(tempTimeStopEnd)
			if len(chsStop) > 0:
				pos = np.argmin(chsTimeStop)
				timeStop = chsTimeStop[pos]
				chStop = chsStop[pos]
				# timeStopEnd = chsTimeStopEnd[pos]
				self.timeStopVect.itemset(ev, timeStop)
				self.chStopVect.itemset(ev, chStop)
				chsDecay = []
				chsDecay2 = []
				chsTimeDecay = []
				chsTimeDecay2 = []
				cond_decay_window = (timeStop + self.delay_time) < self.timeVect
				timeDecayInd = np.argwhere(cond_decay_window).flatten()
				for sig in xrange(self.num_signals):
					# decay in another channel
					# if True:
					if sig != chStop:
						cond_sig_threshold_i = cond_sig_threshold[sig]
						foundDecay = cond_sig_threshold_i[timeDecayInd].any()
						if foundDecay:
							chsDecay.append(sig)
							chsDecay2.append(sig)
							cond_time_decay = np.bitwise_and(cond_decay_window, cond_sig_threshold_i)
							chsTimeDecay.append(self.timeVect[cond_time_decay.argmax()])
							chsTimeDecay2.append(self.timeVect[cond_time_decay.argmax()])
				# check stop signal
				# cond_decay_window = (timeStopEnd + self.delay_time) < self.timeVect
				cond_decay_window = (timeStop + self.delay_time) < self.timeVect
				timeDecayInd = np.argwhere(cond_decay_window).flatten()
				cond_sig_threshold_i = cond_sig_threshold[chStop]
				foundDecay = cond_sig_threshold_i[timeDecayInd].any()
				if foundDecay:
					chsDecay2.append(chStop)
					cond_time_decay = np.bitwise_and(cond_decay_window, cond_sig_threshold_i)
					chsTimeDecay2.append(self.timeVect[cond_time_decay.argmax()])
				if len(chsDecay) > 0:
					pos = np.argmin(chsTimeDecay)
					timeDecay = chsTimeDecay[pos]
					chDecay = chsDecay[pos]
					self.timeDecayVect.itemset(ev, timeDecay)
					self.chDecayVect.itemset(ev, chDecay)
				if len(chsDecay2) > 0:
					pos = np.argmin(chsTimeDecay2)
					if len(chsDecay2) > 1:
						print 'Decays: ', chsTimeDecay2
						cont += 1
					timeDecay2 = chsTimeDecay2[pos]
					chDecay2 = chsDecay2[pos]
					self.timeDecayVect2.itemset(ev, timeDecay2)
					self.chDecayVect2.itemset(ev, chDecay2)

			itbar += 1
			self.utils.bar.update(itbar)
		print 'Cont', cont
		self.lifeTimeVect = np.subtract(self.timeDecayVect, self.timeStopVect, dtype='int16')
		self.lifeTimeVect2 = np.subtract(self.timeDecayVect2, self.timeStopVect, dtype='int16')
		print 'Finished analysing events'

	def CloseAnalysisROOTFile(self):
		if self.analysisFile:
			if self.analysisFile.IsOpen():
				self.analysisFile.Close()
			if self.analysisTree:
				del self.analysisTree
		self.analysisTree = None
		self.analysisFile = None

	def CloseInputROOTFiles(self):
		if self.in_root_file:
			if self.in_root_file.IsOpen():
				self.in_root_file.Close()
			if self.in_root_tree:
				del self.in_root_tree
		self.in_root_file = None
		self.in_root_tree = None

	def FillAnalysisTree(self):
		print 'Filling analysis tree ...'
		self.analysisTree.Branch('stopChannel', self.chStop, 'stopChannel/B')
		self.analysisTree.Branch('decayChannel', self.chDecay, 'decayChannel/B')
		self.analysisTree.Branch('timeStop', self.timeStop, 'timeStop/S')
		self.analysisTree.Branch('timeDecay', self.timeDecay, 'timeDecay/S')
		self.analysisTree.Branch('timeDecay2', self.timeDecay2, 'timeDecay2/S')
		self.analysisTree.Branch('lifeTime', self.lifeTime, 'lifeTime/S')
		self.analysisTree.Branch('lifeTime2', self.lifeTime2, 'lifeTime2/S')
		self.utils.CreateProgressBar(self.tot_events)
		self.utils.bar.start()
		for ev in xrange(self.tot_events):
			self.chStop.itemset(self.chStopVect[ev])
			self.chDecay.itemset(self.chDecayVect[ev])
			self.timeStop.itemset(self.timeStopVect[ev])
			self.timeDecay.itemset(self.timeDecayVect[ev])
			self.timeDecay2.itemset(self.timeDecayVect2[ev])
			self.lifeTime.itemset(self.lifeTimeVect[ev])
			self.lifeTime2.itemset(self.lifeTimeVect2[ev])
			self.analysisTree.Fill()
			self.utils.bar.update(ev + 1)
		if not self.analysisTree.GetFriend(self.in_tree_name):
			self.analysisTree.AddFriend(self.in_tree_name, '{d}/{f}'.format(d=os.path.abspath(self.inDir), f=self.inputFile))
		if not self.analysisTree.GetFriend('vetoedEventsTree'):
			self.analysisTree.AddFriend('vetoedEventsTree', '{d}/{f}.vetoed.root'.format(d=self.outDir, f=self.in_tree_name))
		# self.analysisTree.Write('', ro.TObject.kOverwrite)
		self.analysisTree.Write()
		self.utils.bar.finish()

	def SaveCanvasInlist(self, lista):
		if not os.path.isdir('{d}'.format(d=self.inDir)):
			ExitMessage('The directory does not exist!!!!', os.EX_UNAVAILABLE)
		for canv in lista:
			if self.canvas.has_key(canv):
				self.canvas[canv].SaveAs('{d}/{c}.png'.format(d=self.inDir, c=canv))
				self.canvas[canv].SaveAs('{d}/{c}.root'.format(d=self.inDir, c=canv))

	def SaveAllCanvas(self):
		self.SaveCanvasInlist(self.canvas.keys())

def main():
	parser = OptionParser()
	parser.add_option('-d', '--inDir', dest='inDir', default='.', type='string', help='Directory containing the run files')
	parser.add_option('-c', '--configFile', dest='config', default='CAENAnalysisConfig.cfg', type='string', help='Path to file containing Analysis configuration file')
	parser.add_option('-a', '--automatic', dest='auto', default=False, help='Toggles automatic basic analysis', action='store_true')

	(options, args) = parser.parse_args()
	directory = str(options.inDir)
	config = str(options.config)
	autom = bool(options.auto)

	ana = AnalysisMuonLifeTime(directory, config)

	# ana.LoadAnalysisTree()
	# ana.LoadPickles()
	if autom:
		ana.AnalysisWaves()
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
