#!/usr/bin/env python
import visa
import csv
import numpy as np
from struct import unpack
import time, os, sys
from optparse import OptionParser
import progressbar
import ipdb
from pykeyboard import PyKeyboard
from ConfigParser import ConfigParser
import subprocess as subp
import struct
import ROOT as ro
import shutil
from Utils import *
from collections import OrderedDict
# from DataAcquisition import DataAcquisition


class CCD_Analysis:
	def __init__(self, config='CAENAnalysisConfig.cfg', infile='', bias=0.0, verbose=False):
		self.config, self.outDir, self.inputFile, self.bias, self.verb = config, '/'.join(infile.split('/')[:-1]), infile.split('/')[-1], bias, verbose
		self.treeName = '.'.join(self.inputFile.split('.')[:-1])
		self.fileRaw, self.treeRaw = None, None
		self.ptsWave, self.event, self.events, self.max_events = 0, np.zeros(1, 'I'), 0, 0
		self.eventVect = np.empty(0, 'f8')
		self.timeVect, self.signalWaveVect, self.triggerWaveVect, self.vetoWaveVect = np.empty(0, 'f8'), np.empty(0, 'f8'), np.empty(0, 'f8'), np.empty(0, 'f8')

		self.pedVect, self.pedSigmaVect, self.sigAndPedVect, self.sigAndPedSigmaVect, self.sigVect = None, None, None, None, None
		self.ped, self.pedSigma, self.sigAndPed, self.sigAndPedSigma, self.sig = np.zeros(1, 'f'), np.zeros(1, 'f'), np.zeros(1, 'f'), np.zeros(1, 'f'), np.zeros(1, 'f')
		self.vetoedEvent, self.badShape, self.badPedestal = np.empty(0, '?'), np.empty(0, np.dtype('int8')), np.empty(0, '?')
		self.voltageHV, self.currentHV = np.empty(0, 'f8'), np.empty(0, 'f8')
		self.timeHV = np.empty(0, 'f8')

		self.signalWaveMeanVect, self.signalWaveSigmaVect = None, None

		self.pedestalIntegrationTime = 0.4e-6
		self.pedestalTEndPos = -20e-9
		self.peakTime = 2.121e-6 if self.bias >= 0 else 2.120e-6
		self.doPeakPos = True
		self.peakForward = self.pedestalIntegrationTime/2.0
		self.peakBackward = self.pedestalIntegrationTime/2.0
		self.doBadPedestalCut = True
		self.badShapeCut = 2
		self.doVetoedEventCut = True
		self.peakPosCut = 1e-6
		self.currentCut = 10e-9

		self.Load_Config_File()

		self.cut0 = ro.TCut('cut0', '')

		self.branches1DTotal = ['event', 'vetoedEvent', 'badShape', 'badPedestal', 'voltageHV', 'currentHV', 'timeHV', 'peakPosition', 'pedestal', 'pedestalSigma', 'signalAndPedestal', 'signalAndPedestalSigma', 'signal']
		self.branches1DType = {'event': 'uint32', 'vetoedEvent': 'bool', 'badShape': 'int8', 'badPedestal': 'bool', 'voltageHV': 'float32', 'currentHV': 'float32', 'timeHV': 'float64', 'peakPosition': 'float32', 'pedestal': 'float32', 'pedestalSigma': 'float32', 'signalAndPedestal': 'float32', 'signalAndPedestalSigma': 'float32', 'signal': 'float32'}
		self.branchesWavesTotal = ['time', 'voltageSignal', 'voltageTrigger', 'voltageVeto']
		self.branchesWavesType = {'time': 'float64', 'voltageSignal': 'float64', 'voltageTrigger': 'float64', 'voltageVeto': 'float64'}
		self.branchesAll = self.branches1DTotal + self.branchesWavesTotal

		self.branches1DLoad = ['event', 'voltageHV', 'currentHV', 'timeHV', 'peakPosition']
		self.branchesWavesLoad = ['time', 'voltageSignal']

		self.analysisScalarsBranches = ['pedestal', 'pedestalSigma', 'signalAndPedestal', 'signalAndPedestalSigma', 'signal']

		self.dicBraVect1D = OrderedDict()
		self.dicBraVectWaves = OrderedDict()

		self.hasBranch = {}

		self.pedestalTimeIndices = None  # has the indices for the pedestal for each event
		self.peak_positions = None  # has the peak position in time for the peak of the signal for each event
		self.peak_position = np.zeros(1, 'f')
		self.signalTimeIndices = None  # has the indices for the integral of the signal for each event

		self.utils = Utils()

	def Load_Config_File(self):
		parser = ConfigParser()
		if os.path.isfile(self.config):
			print 'Reading configuration file:', self.config, '...'
			parser.read(self.config)

			if parser.has_section('ANALYSIS'):
				if parser.has_option('ANALYSIS', 'bias'):
					self.bias = parser.getfloat('ANALYSIS', 'bias')
				if parser.has_option('ANALYSIS', 'input_file') and self.inputFile == '':
					self.inputFile = parser.get('ANALYSIS', 'input_file')
				if parser.has_option('ANALYSIS', 'out_dir') and self.outDir == '':
					self.outDir = parser.get('ANALYSIS', 'out_dir')
				if parser.has_option('ANALYSIS', 'max_events'):
					self.max_events = parser.getint('ANALYSIS', 'max_events')
				if parser.has_option('ANALYSIS', 'peak_time'):
					self.peakTime = parser.getfloat('ANALYSIS', 'peak_time') * 1e-6
				if parser.has_option('ANALYSIS', 'integration_time'):
					self.pedestalIntegrationTime = parser.getfloat('ANALYSIS', 'integration_time') * 1e-6
				if parser.has_option('ANALYSIS', 'do_peak_positioning'):
					self.doPeakPos = parser.getboolean('ANALYSIS', 'do_peak_positioning')
				if parser.has_option('ANALYSIS', 'transition_time'):
					self.pedestalTEndPos = parser.getfloat('ANALYSIS', 'transition_time') * 1e-9
				if parser.has_option('ANALYSIS', 'forward_bakward_ratio'):
					self.peakForward = self.pedestalIntegrationTime / (1.0 + 1.0 / parser.getfloat('ANALYSIS', 'forward_bakward_ratio'))
					self.peakBackward = self.pedestalIntegrationTime / (1.0 + parser.getfloat('ANALYSIS', 'forward_bakward_ratio'))

			if parser.has_section('CUTS'):
				if parser.has_option('CUTS', 'bad_pedestal'):
					self.doBadPedestalCut = parser.getboolean('CUTS', 'bad_pedestal')
				if parser.has_option('CUTS', 'bad_shape'):
					self.badShapeCut = parser.getint('CUTS', 'bad_shape')
				if parser.has_option('CUTS', 'vetoed_events'):
					self.doVetoedEventCut = parser.getboolean('CUTS', 'vetoed_events')
				if parser.has_option('CUTS', 'peak_position'):
					self.peakPosCut = parser.getfloat('CUTS', 'peak_position') * 1e-6
				if parser.has_option('CUTS', 'current_cut'):
					self.currentCut = parser.getfloat('CUTS', 'current_cut') * 1e-9

	def AnalysisWaves(self, doCuts0=True):
		self.OpenROOTFile('UPDATE')
		self.LoadTree()
		if doCuts0: self.CreateCut0()
		if not self.hasBranch['peakPosition']:
			self.LoadVectorsFromTree()
			self.ExplicitVectorsFromDictionary()
			if self.doPeakPos:
				self.FindRealPeakPosition()
			else:
				self.peak_positions = np.full(self.peakTime)
			self.FillTreePeakPositions()
			self.CloseROOTFile()
			self.OpenROOTFile('UPDATE')
			self.Reset_Braches_Lists_And_Dictionaries()
			self.LoadTree()
		self.AddPeakPositionCut()
		if not np.array([self.hasBranch[key0] for key0 in self.analysisScalarsBranches]).all():
			self.LoadVectorsFromTree()
			self.ExplicitVectorsFromDictionary()
			self.FindPedestalPosition()
			self.FindSignalPositions(self.peakBackward, self.peakForward)
			self.CalculatePedestalsAndSignals()
			self.FillPedestalsAndSignals()
			self.CloseROOTFile()
		self.OpenROOTFile('READ')
		self.Reset_Braches_Lists_And_Dictionaries()
		self.LoadTree()

	def OpenROOTFile(self, mode='READ'):
		if not os.path.isdir(self.outDir):
			print 'Directory:', self.outDir, '; does not exist. Exiting!'
			exit()
		if not os.path.isfile('{o}/{f}'.format(o=self.outDir, f=self.inputFile)):
			print 'File:', self.inputFile, '; does not exist in:', self.outDir, '. Exiting!'
			exit()
		if self.fileRaw:
			if self.fileRaw.IsOpen():
				if self.fileRaw.GetOption().lower() != mode.lower():
					if self.fileRaw.ReOpen(mode) == -1:
						print 'Could not reopen file:', self.fileRaw, 'in mode:', mode, '. Exiting!'
						exit()
				return
			else:
				self.fileRaw = None
		self.fileRaw = ro.TFile('{o}/{f}'.format(o=self.outDir, f=self.inputFile), mode)

	def LoadTree(self):
		if not self.fileRaw.IsOpen():
			print 'First OpenROOTFile before Loading Tree'
			return
		self.treeRaw = self.fileRaw.Get(self.treeName)
		self.hasBranch = {branch: self.TreeHasBranch(branch) for branch in self.branchesAll}
		self.UpdateBranchesLists()
		self.IsTimeHVaTimeStamp()

	def TreeHasBranch(self, branch):
		if self.treeRaw.GetLeaf(branch):
			return True
		else:
			return False

	def UpdateBranchesLists(self):
		for branch in self.branches1DTotal[:]:
			if not self.hasBranch[branch]:
				self.branches1DTotal.remove(branch)
		for branch in self.branches1DLoad[:]:
			if not self.hasBranch[branch]:
				self.branches1DLoad.remove(branch)

		for branch in self.branchesWavesTotal[:]:
			if not self.hasBranch[branch]:
				self.branchesWavesTotal.remove(branch)
		for branch in self.branchesWavesLoad[:]:
			if not self.hasBranch[branch]:
				self.branchesWavesLoad.remove(branch)

	def IsTimeHVaTimeStamp(self):
		if self.hasBranch['timeHV']:
			if self.treeRaw.GetLeaf('timeHV').GetTypeName() != 'TDatime':
				# self.branches1DTotal = ['timeHV.AsDouble()' if branch == 'timeHV' else branch for branch in self.branches1DTotal]
				self.branches1DLoad = ['timeHV.AsDouble()' if branch == 'timeHV' else branch for branch in self.branches1DLoad]
				del self.branches1DType['timeHV']
				self.branches1DType['timeHV.AsDouble()'] = 'float64'
				# self.branches1DType = {'event': 'uint32', 'vetoedEvent': 'bool', 'badShape': 'int8', 'badPedestal': 'bool', 'voltageHV': 'float32', 'currentHV': 'float32', 'timeHV.AsDouble()': 'float64', 'peakPosition': 'float64'}
			else:
				# self.branches1DTotal = ['timeHV.Convert()' if branch == 'timeHV' else branch for branch in self.branches1DTotal]
				self.branches1DLoad = ['timeHV.Convert()' if branch == 'timeHV' else branch for branch in self.branches1DLoad]
				del self.branches1DType['timeHV']
				self.branches1DType['timeHV.Convert()'] = 'uint32'
				# self.branches1DType = {'event': 'uint32', 'vetoedEvent': 'bool', 'badShape': 'int8', 'badPedestal': 'bool', 'voltageHV': 'float32', 'currentHV': 'float32', 'timeHV.Convert()': 'uint32', 'peakPosition': 'float64'}

	def CreateCut0(self):
		if self.cut0.GetTitle() != '':
			self.cut0.SetTitle('')
		if self.doBadPedestalCut and 'badPedestal' in self.branches1DTotal:
			self.cut0 += ro.TCut('badPedCut', 'badPedestal==0')
		if self.doVetoedEventCut and 'vetoedEvent' in self.branches1DTotal:
			self.cut0 += ro.TCut('vetoedEventCut', 'vetoedEvent==0')
		if self.badShapeCut == 1 and 'badShape' in self.branches1DTotal:
			self.cut0 += ro.TCut('badShapeCut', 'badShape!=1')
		elif self.badShapeCut == 2 and 'badShape' in self.branches1DTotal:
			self.cut0 += ro.TCut('badShapeCut', 'badShape==0')
		if 'currentHV' in self.branches1DTotal:
			self.cut0 += ro.TCut('currentCut', 'abs(currentHV)<{cc}'.format(cc=self.currentCut))

	def LoadVectorsFromTree(self):
		self.max_events = self.treeRaw.GetEntries() if self.max_events == 0 else self.max_events
		self.ptsWave = self.treeRaw.GetLeaf('time').GetLen()

		leng = self.treeRaw.Draw(':'.join(self.branches1DLoad), self.cut0, 'goff para', self.max_events)
		if leng == -1:
			print 'Error, could not load the branches: {b}. Try again :('.format(b=':'.join(self.branches1DLoad))
			return
		while leng > self.treeRaw.GetEstimate():
			self.treeRaw.SetEstimate(leng)
			leng = self.treeRaw.Draw(':'.join(self.branches1DLoad), self.cut0, 'goff para', self.max_events)
		self.events = leng
		for pos, branch in enumerate(self.branches1DLoad):
			if self.verb: print 'Vectorising branch:', branch, '...', ; sys.stdout.flush()
			temp = self.treeRaw.GetVal(pos)
			self.dicBraVect1D[branch] = np.array([temp[ev] for ev in xrange(self.events)], dtype=np.dtype(self.branches1DType[branch]))
			if self.verb: print 'Done'

		leng = self.treeRaw.Draw(':'.join(self.branchesWavesLoad), self.cut0, 'goff para', self.max_events)
		if leng == -1:
			print 'Error, could not load the branches {b}. Try again :('.format(b=':'.join(self.branchesWavesLoad))
			return
		while leng > self.treeRaw.GetEstimate():
			self.treeRaw.SetEstimate(leng)
			leng = self.treeRaw.Draw(':'.join(self.branchesWavesLoad), self.cut0, 'goff para', self.max_events)
		for pos, branch in enumerate(self.branchesWavesLoad):
			if self.verb: print 'Vectorising branch:', branch, '...', ; sys.stdout.flush()
			temp = self.treeRaw.GetVal(pos)
			self.dicBraVectWaves[branch] = np.array([[temp[ev * self.ptsWave + pt] for pt in xrange(self.ptsWave)] for ev in xrange(self.events)], dtype=np.dtype(self.branchesWavesType[branch]))
			if self.verb: print 'Done'

	def ExplicitVectorsFromDictionary(self):
		if self.hasBranch['voltageSignal']:
			self.signalWaveVect = self.dicBraVectWaves['voltageSignal']
		# if self.hasBranch['voltageTrigger']:
		# 	self.triggerWaveVect = self.dicBraVectWaves['voltageTrigger']
		# if self.hasBranch['voltageVeto']:
		# 	self.vetoWaveVect = self.dicBraVectWaves['voltageVeto']
		if self.hasBranch['time']:
			self.timeVect = self.dicBraVectWaves['time']
		if self.hasBranch['event']:
			self.eventVect = self.dicBraVect1D['event']
		# if self.hasBranch['vetoedEvent']:
		# 	self.vetoedEvent = self.dicBraVectWaves['vetoedEvent']
		# if self.hasBranch['badShape']:
		# 	self.badShape = self.dicBraVectWaves['badShape']
		# if self.hasBranch['badPedestal']:
		# 	self.badPedestal = self.dicBraVectWaves['badPedestal']
		if self.hasBranch['voltageHV']:
			self.voltageHV = self.dicBraVect1D['voltageHV']
		if self.hasBranch['currentHV']:
			self.currentHV = self.dicBraVect1D['currentHV']
		if self.hasBranch['timeHV']:
			key = 'timeHV.Convert()' if 'timeHV.Convert()' in self.branches1DLoad else 'timeHV.AsDouble()'
			self.timeHV = self.dicBraVect1D[key]
		if self.hasBranch['peakPosition']:
			self.peak_positions = self.dicBraVect1D['peakPosition']

	def FindRealPeakPosition(self):
		print 'Getting real peak positions...', ;sys.stdout.flush()
		mpos = self.signalWaveVect.argmin(axis=1) if self.bias >= 0 else self.signalWaveVect.argmax(axis=1)
		time_mpos = self.timeVect[:, mpos].diagonal()
		xmin, xmax = time_mpos - self.pedestalIntegrationTime / 2.0, time_mpos + self.pedestalIntegrationTime / 2.0
		fit = [ro.TGraph(len(timei), timei, self.signalWaveVect[it]).Fit('pol2', 'QMN0FS', '', xmin[it], xmax[it]) for it, timei in enumerate(self.timeVect)]
		b, a = np.array([fiti.Parameter(1) for fiti in fit]), np.array([fiti.Parameter(2) for fiti in fit])
		self.peak_positions = np.divide(-b, 2 * a)
		print 'Done'

	def FillTreePeakPositions(self):
		print 'Filling tree with peak positions...'
		peakPosBra = self.treeRaw.Branch('peakPosition', self.peak_position, 'peakPosition/F')
		self.utils.CreateProgressBar(self.treeRaw.GetEntries())
		self.utils.bar.start()
		for ev in xrange(self.treeRaw.GetEntries()):
			self.treeRaw.GetEntry(ev)
			if ev in self.eventVect:
				try:
					self.peak_position.itemset(self.peak_positions[np.argwhere(ev == self.eventVect).flatten()])
				except ValueError:
					print 'could not fill event', ev, 'it should have a peak position of:', self.peak_positions[np.argwhere(ev == self.eventVect).flatten()]
					exit()
			else:
				self.peak_position.itemset(0)
			peakPosBra.Fill()
			self.utils.bar.update(ev + 1)
		self.treeRaw.Write()
		self.utils.bar.finish()

	def CloseROOTFile(self):
		self.treeRaw.Delete()
		self.fileRaw.Close()

	def Reset_Braches_Lists_And_Dictionaries(self):
		self.branches1DTotal = ['event', 'vetoedEvent', 'badShape', 'badPedestal', 'voltageHV', 'currentHV', 'timeHV', 'peakPosition', 'pedestal', 'pedestalSigma', 'signalAndPedestal', 'signalAndPedestalSigma', 'signal']
		self.branches1DType = {'event': 'uint32', 'vetoedEvent': 'bool', 'badShape': 'int8', 'badPedestal': 'bool', 'voltageHV': 'float32', 'currentHV': 'float32', 'timeHV': 'float64', 'peakPosition': 'float32', 'pedestal': 'float32', 'pedestalSigma': 'float32', 'signalAndPedestal': 'float32', 'signalAndPedestalSigma': 'float32', 'signal': 'float32'}
		self.branchesWavesTotal = ['time', 'voltageSignal', 'voltageTrigger', 'voltageVeto']
		self.branchesWavesType = {'time': 'float64', 'voltageSignal': 'float64', 'voltageTrigger': 'float64', 'voltageVeto': 'float64'}
		self.branchesAll = self.branches1DTotal + self.branchesWavesTotal
		self.branches1DLoad = ['event', 'voltageHV', 'currentHV', 'timeHV', 'peakPosition']
		self.branchesWavesLoad = ['time', 'voltageSignal']
		self.analysisScalarsBranches = ['pedestal', 'pedestalSigma', 'signalAndPedestal', 'signalAndPedestalSigma', 'signal']
		self.dicBraVect1D = OrderedDict()
		self.dicBraVectWaves = OrderedDict()
		self.hasBranch = {}

	def AddPeakPositionCut(self):
		self.cut0 += ro.TCut('peakTimeCut', 'abs(peakPosition-{pp})<={ppc}'.format(pp=self.peakTime, ppc=self.peakPosCut))

	def FindPedestalPosition(self):
		print 'Calculating position of pedestals...', ;sys.stdout.flush()
		self.pedestalTimeIndices = [np.argwhere(np.bitwise_and(self.pedestalTEndPos - self.pedestalIntegrationTime <= timeVectEvi, timeVectEvi <= self.pedestalTEndPos)).flatten() for timeVectEvi in self.timeVect]
		print 'Done'

	def FindSignalPositions(self, backward, forward):
		print 'Calculating position of signals...', ;sys.stdout.flush()
		self.signalTimeIndices = [np.argwhere(abs(timeVectEvi - self.peak_positions[it] + (forward - backward)/2.0) <= (forward + backward)/2.0).flatten() for it, timeVectEvi in enumerate(self.timeVect)]
		print 'Done'

	def CalculatePedestalsAndSignals(self):
		print 'Calculating pedestals and signals...', ;sys.stdout.flush()
		self.pedVect = np.array([self.signalWaveVect[ev, pedTimeIndxs].mean() if pedTimeIndxs.size > 0 else -10 for ev, pedTimeIndxs in enumerate(self.pedestalTimeIndices)])
		self.pedSigmaVect = np.array([self.signalWaveVect[ev, pedTimeIndxs].std() if pedTimeIndxs.size > 1 else -10 for ev, pedTimeIndxs in enumerate(self.pedestalTimeIndices)])
		self.sigAndPedVect = np.array([self.signalWaveVect[ev, sigTimeIndxs].mean() if sigTimeIndxs.size > 0 else -10 for ev, sigTimeIndxs in enumerate(self.signalTimeIndices)])
		self.sigAndPedSigmaVect = np.array([self.signalWaveVect[ev, sigTimeIndxs].std() if sigTimeIndxs.size > 1 else -10 for ev, sigTimeIndxs in enumerate(self.signalTimeIndices)])
		self.sigVect = np.subtract(self.sigAndPedVect, self.pedVect)
		print 'Done'

	def FillPedestalsAndSignals(self):
		print 'Filling tree with peak positions...'
		pedBra = self.treeRaw.Branch('pedestal', self.ped, 'pedestal/F')
		pedSigmaBra = self.treeRaw.Branch('pedestalSigma', self.pedSigma, 'pedestalSigma/F')
		pedSignalBra = self.treeRaw.Branch('signalAndPedestal', self.sigAndPed, 'signalAndPedestal/F')
		pedSignalSigmaBra = self.treeRaw.Branch('signalAndPedestalSigma', self.sigAndPed, 'signalAndPedestalSigma/F')
		sigBra = self.treeRaw.Branch('signal', self.sig, 'signal/F')
		self.utils.CreateProgressBar(self.max_events)
		self.utils.bar.start()
		for ev in xrange(self.max_events):
			self.treeRaw.GetEntry(ev)
			if ev in self.eventVect:
				try:
					argum = np.argwhere(ev == self.eventVect).flatten()
					self.ped.itemset(self.pedVect[argum])
					self.pedSigma.itemset(self.pedSigmaVect[argum])
					self.sigAndPed.itemset(self.sigAndPedVect[argum])
					self.sigAndPedSigma.itemset(self.sigAndPedSigmaVect[argum])
					self.sig.itemset(self.sigVect[argum])
				except ValueError:
					print 'Could not fill event', ev
					exit()
			else:
				self.ped.itemset(0)
				self.pedSigma.itemset(0)
				self.sigAndPed.itemset(0)
				self.sigAndPedSigma.itemset(0)
				self.sig.itemset(0)
			pedBra.Fill()
			pedSigmaBra.Fill()
			pedSignalBra.Fill()
			pedSignalSigmaBra.Fill()
			sigBra.Fill()
			self.utils.bar.update(ev + 1)
		self.treeRaw.Write()
		self.utils.bar.finish()

	def ExtractMeanOfWaveforms(self):
		self.signalWaveMeanVect = self.signalWaveVect.mean(axis=0)
		self.signalWaveSigmaVect = self.signalWaveVect.std(axis=0)

	def ResetTreeToOriginal(self, keepBranches=['event','time','voltageSignal','voltageTrigger','voltageVeto','vetoedEvent','badShape','badPedestal','voltageHV','currentHV','timeHV']):
		print 'Restoring tree with the following branches:', keepBranches, '...'
		raw_input('Press a key and Enter to continue: ')
		self.OpenROOTFile('READ')
		self.LoadTree()
		self.treeRaw.SetBranchStatus('*', 0)
		for branch in keepBranches:
			if self.TreeHasBranch(branch):
				self.treeRaw.SetBranchStatus(branch, 1)
		newFile = ro.TFile('{o}/temp.root'.format(o=self.outDir), 'recreate')
		newTree = self.treeRaw.CloneTree()
		newTree.Print()
		newFile.Write()
		del self.fileRaw
		del newFile
		self.fileRaw = None
		checkFile = ro.TFile('{o}/temp.root'.format(o=self.outDir), 'READ')
		checkTree = checkFile.Get(self.treeName)
		doMoveFile = True
		if checkTree:
			for branch in keepBranches:
				if not checkTree.GetLeaf(branch):
					doMoveFile = False
					break
			if doMoveFile:
				print 'The file was cloned successfully :)'
				checkFile.Close()
				del checkFile
				shutil.move('{o}/temp.root'.format(o=self.outDir), '{o}/{f}'.format(o=self.outDir, f=self.inputFile))
				return
		print 'The file was not cloned successfully :S. Check original tree and "temp.root"'

if __name__ == '__main__':

	parser = OptionParser()
	parser.add_option('-c', '--configFile', dest='config', default='CAENAnalysisConfig.cfg', type='string', help='Path to file containing Analysis configuration file')
	parser.add_option('-i', '--inputFile', dest='infile', default='', type='string', help='Path to root file to be analysed')
	parser.add_option('-b', '--bias', dest='bias', default=0, type='float', help='Bias voltage used')
	parser.add_option('-v', '--verbose', dest='verb', default=False, action='store_true', help='Toggles verbose')
	(options, args) = parser.parse_args()
	config = str(options.config)
	infile = str(options.infile)
	bias = float(options.bias)
	verb = bool(options.verb)

	z = CCD_Analysis(config, infile, bias, verb)
