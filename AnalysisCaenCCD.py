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

from Channel_Caen import Channel_Caen
from Settings_Caen import Settings_Caen
from HV_Control import HV_Control
from Utils import *
# from memory_profiler import profile

trig_rand_time = 0.2
wait_time_hv = 7
BRANCHES1DTOTAL = ['event', 'vetoedEvent', 'badShape', 'badPedestal', 'voltageHV', 'currentHV', 'timeHV', 'peakPosition', 'pedestal', 'pedestalSigma', 'signalAndPedestal', 'signalAndPedestalSigma', 'signal']
BRANCHES1DTYPE = {'event': 'uint32', 'vetoedEvent': 'bool', 'badShape': 'int8', 'badPedestal': 'bool', 'voltageHV': 'float32', 'currentHV': 'float32', 'timeHV.AsDouble()': 'float64', 'timeHV.Convert()': 'uint32', 'peakPosition': 'float32', 'pedestal': 'float32', 'pedestalSigma': 'float32', 'signalAndPedestal': 'float32','signalAndPedestalSigma': 'float32', 'signal': 'float32'}
BRANCHESWAVESTOTAL = ['time', 'voltageSignal', 'voltageTrigger', 'voltageVeto']
BRANCHESWAVESTYPE = {'time': 'float64', 'voltageSignal': 'float64', 'voltageTrigger': 'float64', 'voltageVeto': 'float64'}
BRANCHES1DLOAD = ['event', 'voltageHV', 'currentHV', 'timeHV.Convert()', 'timeHV.AsDouble()', 'peakPosition']
BRANCHESWAVESLOAD = ['time', 'voltageSignal']
ANALYSISSCALARBRANCHES = ['pedestal', 'pedestalSigma', 'signalAndPedestal', 'signalAndPedestalSigma', 'signal']

class AnalysisCaenCCD:
	def __init__(self, directory='.', config='CAENAnalysisConfig.cfg', infile='', bias=0.0, verbose=False):
		print 'Starting CCD Analysis ...'

		self.config = config
		self.verb = verbose
		self.inputFile = ''
		self.inDir = directory
		self.in_root_file = None
		self.in_tree_name = ''
		self.in_root_tree = None
		self.settings = None
		self.signal_ch = None
		self.trigger_ch = None
		self.veto_ch = None
		self.max_events = 0
		self.bias = bias
		self.pedestalIntegrationTime = 0.4e-6
		self.pedestalTEndPos = -20e-9
		self.peakTime = 2.121e-6 if self.bias >= 0 else 2.120e-6
		self.doPeakPos = True
		self.peakForward = self.pedestalIntegrationTime / 2.0
		self.peakBackward = self.pedestalIntegrationTime / 2.0
		self.doBadPedestalCut = True
		self.badShapeCut = 2
		self.doVetoedEventCut = True
		self.peakPosCut = 1e-6
		self.currentCut = 10e-9

		self.analysisTreeExisted = False

		self.ptsWave, self.event, self.events, self.max_events = 0, np.zeros(1, 'I'), 0, 0
		self.eventVect = np.empty(0, 'f8')
		self.timeVect, self.signalWaveVect, self.triggerWaveVect, self.vetoWaveVect = np.empty(0, 'f8'), np.empty(0, 'f8'), np.empty(0, 'f8'), np.empty(0, 'f8')

		self.pedVect, self.pedSigmaVect, self.sigAndPedVect, self.sigAndPedSigmaVect, self.sigVect = None, None, None, None, None
		self.ped, self.pedSigma, self.sigAndPed, self.sigAndPedSigma, self.sig = np.zeros(1, 'f'), np.zeros(1, 'f'), np.zeros(1, 'f'), np.zeros(1, 'f'), np.zeros(1, 'f')
		self.vetoedEvent, self.badShape, self.badPedestal = np.empty(0, '?'), np.empty(0, np.dtype('int8')), np.empty(0, '?')
		self.voltageHV, self.currentHV = np.empty(0, 'f8'), np.empty(0, 'f8')
		self.timeHV = np.empty(0, 'f8')

		self.signalWaveMeanVect, self.signalWaveSigmaVect = None, None

		self.cut0 = ro.TCut('cut0', '')

		self.branches1DTotal = BRANCHES1DTOTAL[:]
		self.branches1DType = BRANCHES1DTYPE.copy()
		self.branchesWavesTotal = BRANCHESWAVESTOTAL[:]
		self.branchesWavesType = BRANCHESWAVESTYPE.copy()
		self.branchesAll = self.branches1DTotal + self.branchesWavesTotal

		self.branches1DLoad = BRANCHES1DLOAD[:]
		self.branchesWavesLoad = BRANCHESWAVESLOAD[:]
		self.dic1DVectLoaded = {}
		self.dicWavesVectLoaded = {}

		self.analysisScalarsBranches = ANALYSISSCALARBRANCHES[:]

		self.dicBraVect1D = OrderedDict()
		self.dicBraVectWaves = OrderedDict()

		self.hasBranch = {}

		self.pedestalTimeIndices = None  # has the indices for the pedestal for each event
		self.peak_positions = None  # has the peak position in time for the peak of the signal for each event
		self.peak_position = np.zeros(1, 'f')
		self.signalTimeIndices = None  # has the indices for the integral of the signal for each event

		self.cut0 = ro.TCut('cut0', '')

		self.utils = Utils()

		self.outDir = ''

		self.canvas = {}
		self.profile = {}
		self.histo = {}


		if infile == '' and directory != '.':
			print 'Is analysis of data after 06/18...'
			self.LoadInputTree()
			self.LoadPickles()
			self.outDir = self.inDir
			self.inputFile = self.in_tree_name + '.root'

		elif infile != '' and directory == '.':
			print 'Is analysis of data before 06/18...'
			self.outDir, self.inputFile = '/'.join(infile.split('/')[:-1]), infile.split('/')[-1]
			self.in_tree_name = '.'.join(self.inputFile.split('.')[:-1])
			self.inDir = self.outDir
			self.LoadInputTree()
			self.LoadPickles()
		else:
			ExitMessage('I don\'t know what to do. If you want to run the analysis for old files, give inputFile, configFile, bias. If you want to run the analysis for new files, just give the directory where the pickles, data files and root files are.', os.EX_CONFIG)

		self.analysisFile = None
		self.analysisTree = None
		self.analysisTreeName = self.in_tree_name + '.analysis'

		self.Load_Config_File()


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
		self.OpenAnalysisROOTFile('UPDATE')
		if doCuts0: self.CreateCut0()

		if not self.hasBranch['peakPosition']:
			self.LoadVectorsFromTree()
			self.ExplicitVectorsFromDictionary()
			if self.doPeakPos:
				self.FindRealPeakPosition()
			else:
				self.peak_positions = np.full(self.events, self.peakTime)
			self.FillTreePeakPositions()
			self.CloseAnalysisROOTFile()
			self.CloseInputROOTFiles()
			self.Reset_Braches_Lists_And_Dictionaries()
			self.LoadInputTree()
			self.OpenAnalysisROOTFile('UPDATE')
		# self.AddPeakPositionCut()
		if not np.array([self.hasBranch[key0] for key0 in self.analysisScalarsBranches]).all():
			self.LoadVectorsFromTree()
			self.ExplicitVectorsFromDictionary()
			self.FindPedestalPosition()
			self.FindSignalPositions(self.peakBackward, self.peakForward)
			self.CalculatePedestalsAndSignals()
			self.FillPedestalsAndSignals()
		self.CloseAnalysisROOTFile()
		self.CloseInputROOTFiles()
		self.Reset_Braches_Lists_And_Dictionaries()
		self.LoadInputTree()
		self.OpenAnalysisROOTFile('READ')
		self.PlotPeakPositionDistributions()

	def LoadInputTree(self):
		if not os.path.isdir(self.inDir):
			ExitMessage('The directory {d} does not exist! Exiting...'.format(d=self.inDir), os.EX_DATAERR)
		if self.in_tree_name == '':
			root_files = glob.glob('{d}/*.root'.format(d=self.inDir))
			if len(root_files) == 1:
				self.in_tree_name = root_files[0].split('/')[-1].split('.root')[0]
			elif len(root_files) == 0:
				ExitMessage('There is no root file inside the directory {d}. Exiting...'.format(d=self.inDir), os.EX_DATAERR)
			else:
				print 'The following files were encountered:'
				for it in root_files:
					print it
				file_to_open = raw_input('Copy and paste the one that you want to open: ')
				if file_to_open in root_files:
					self.in_tree_name = file_to_open.split('/')[-1].split('.root')[0]
				else:
					ExitMessage('The given file: {f} is not part of the given list. Exiting...'.format(f=file_to_open), os.EX_DATAERR)
		if os.path.isfile('{d}/{f}.root'.format(d=self.inDir, f=self.in_tree_name)):
			self.in_root_file = ro.TFile('{d}/{f}.root'.format(d=self.inDir, f=self.in_tree_name))
			self.in_root_tree = self.in_root_file.Get(self.in_tree_name)
		else:
			ExitMessage('The file {f} is not inside {d}. Exiting...'.format(d=self.inDir, f=self.in_tree_name), os.EX_DATAERR)

	def LoadPickles(self):
		if self.in_tree_name != '':
			if os.path.isdir(self.inDir):
				if os.path.isfile('{d}/{f}.settings'.format(d=self.inDir, f=self.in_tree_name)):
					with open('{d}/{f}.settings'.format(d=self.inDir, f=self.in_tree_name), 'rb') as pklsets:
						self.settings = pickle.load(pklsets)
				if os.path.isfile('{d}/{f}.signal_ch'.format(d=self.inDir, f=self.in_tree_name)):
					with open('{d}/{f}.signal_ch'.format(d=self.inDir, f=self.in_tree_name)) as pklsign:
						self.signal_ch = pickle.load(pklsign)
				if os.path.isfile('{d}/{f}.trigger_ch'.format(d=self.inDir, f=self.in_tree_name)):
					with open('{d}/{f}.trigger_ch'.format(d=self.inDir, f=self.in_tree_name)) as pkltrig:
						self.trigger_ch = pickle.load(pkltrig)
				if os.path.isfile('{d}/{f}.veto_ch'.format(d=self.inDir, f=self.in_tree_name)):
					with open('{d}/{f}.veto_ch'.format(d=self.inDir, f=self.in_tree_name)) as pklveto:
						self.veto_ch = pickle.load(pklveto)
				elif os.path.isfile('{d}/{f}.veto'.format(d=self.inDir, f=self.in_tree_name)):
					with open('{d}/{f}.veto'.format(d=self.inDir, f=self.in_tree_name)) as pklveto:
						self.veto_ch = pickle.load(pklveto)

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
					ExitMessage('Can\'t open the file {f}.root in {m} mode because it does not exist!!! Exiting...'.format(f=self.analysisTreeName, m=mode), os.EX_DATAERR)

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
		self.IsTimeHVaTimeStamp()
		self.UpdateBranchesLists()

	def TreeHasBranch(self, branch):
		if self.in_root_tree.GetBranch(branch) or self.analysisTree.GetBranch(branch):
			return True
		return False

	def IsTimeHVaTimeStamp(self):
		if self.hasBranch['timeHV']:
			if self.in_root_tree.GetLeaf('timeHV').GetTypeName() != 'TDatime':
				# self.branches1DLoad = ['timeHV.AsDouble()' if branch == 'timeHV' else branch for branch in self.branches1DLoad]
				if 'timeHV.Convert()' in self.branches1DLoad: self.branches1DLoad.remove('timeHV.Convert()')
				# del self.branches1DType['timeHV']
				# self.branches1DType['timeHV.AsDouble()'] = 'float64'
				if self.branches1DType.has_key('timeHV.Convert()'): del self.branches1DType['timeHV.Convert()']
			else:
				# self.branches1DLoad = ['timeHV.Convert()' if branch == 'timeHV' else branch for branch in self.branches1DLoad]
				if 'timeHV.AsDouble()' in self.branches1DLoad: self.branches1DLoad.remove('timeHV.AsDouble()')
				# del self.branches1DType['timeHV']
				# self.branches1DType['timeHV.Convert()'] = 'uint32'
				if self.branches1DType.has_key('timeHV.AsDouble()'): del self.branches1DType['timeHV.AsDouble()']

	def UpdateBranchesLists(self):
		for branch in self.branches1DTotal[:]:
			if not self.hasBranch[branch]:
				self.branches1DTotal.remove(branch)
		for branch in self.branches1DLoad[:]:
			if branch.startswith('timeHV'):
				if not self.hasBranch['timeHV']:
					self.branches1DLoad.remove(branch)
				else:
					if self.dicBraVect1D.has_key(branch):
						self.dic1DVectLoaded[branch] = True if self.dicBraVect1D[branch] else False
					else:
						self.dic1DVectLoaded[branch] = False
			elif not self.hasBranch[branch]:
				self.branches1DLoad.remove(branch)
			else:
				if self.dicBraVect1D.has_key(branch):
					self.dic1DVectLoaded[branch] = True if self.dicBraVect1D[branch] else False
				else:
					self.dic1DVectLoaded[branch] = False

		for branch in self.branchesWavesTotal[:]:
			if not self.hasBranch[branch]:
				self.branchesWavesTotal.remove(branch)
		for branch in self.branchesWavesLoad[:]:
			if not self.hasBranch[branch]:
				self.branchesWavesLoad.remove(branch)
			else:
				if self.dicBraVectWaves.has_key(branch):
					self.dicWavesVectLoaded[branch] = True if self.dicBraVectWaves[branch] else False
				else:
					self.dicWavesVectLoaded[branch] = False

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

	def ResetCut0(self):
		self.cut0.Clear()
		self.cut0 = ro.TCut('cut0', '')

	def LoadVectorsFromTree(self):
		self.max_events = self.in_root_tree.GetEntries() if self.max_events == 0 else self.max_events
		self.ptsWave = self.in_root_tree.GetLeaf('time').GetLen()
		working_tree = self.analysisTree if self.analysisTreeExisted else self.in_root_tree
		branches_to_load_1D = [branch for branch in self.branches1DLoad if not self.dic1DVectLoaded[branch]]
		leng = working_tree.Draw(':'.join(branches_to_load_1D), self.cut0, 'goff para', self.max_events)
		if leng == -1:
			print 'Error, could not load the branches: {b}. Try again :('.format(b=':'.join(branches_to_load_1D))
			return
		while leng > working_tree.GetEstimate():
			working_tree.SetEstimate(leng)
			leng = working_tree.Draw(':'.join(branches_to_load_1D), self.cut0, 'goff para', self.max_events)
		self.events = leng
		for pos, branch in enumerate(branches_to_load_1D):
			if self.verb: print 'Vectorising branch:', branch, '...', ; sys.stdout.flush()
			temp = working_tree.GetVal(pos)
			self.dicBraVect1D[branch] = np.array([temp[ev] for ev in xrange(self.events)], dtype=np.dtype(self.branches1DType[branch]))
			self.dic1DVectLoaded[branch] = True
			if self.verb: print 'Done'

		branches_to_load_waves = [branch for branch in self.branchesWavesLoad if not self.dicWavesVectLoaded[branch]]
		leng = working_tree.Draw(':'.join(branches_to_load_waves), self.cut0, 'goff para', self.max_events)
		if leng == -1:
			print 'Error, could not load the branches {b}. Try again :('.format(b=':'.join(branches_to_load_waves))
			return
		while leng > working_tree.GetEstimate():
			working_tree.SetEstimate(leng)
			leng = working_tree.Draw(':'.join(branches_to_load_waves), self.cut0, 'goff para', self.max_events)
		for pos, branch in enumerate(branches_to_load_waves):
			if self.verb: print 'Vectorising branch:', branch, '...', ; sys.stdout.flush()
			temp = working_tree.GetVal(pos)
			self.dicBraVectWaves[branch] = np.array([[temp[ev * self.ptsWave + pt] for pt in xrange(self.ptsWave)] for ev in xrange(self.events)], dtype=np.dtype(self.branchesWavesType[branch]))
			self.dicWavesVectLoaded[branch] = True
			if self.verb: print 'Done'

	def ExplicitVectorsFromDictionary(self):
		if self.hasBranch['voltageSignal']:
			if self.dicWavesVectLoaded['voltageSignal']:
				self.signalWaveVect = self.dicBraVectWaves['voltageSignal']
		if self.hasBranch['time']:
			if self.dicWavesVectLoaded['time']:
				self.timeVect = self.dicBraVectWaves['time']
		if self.hasBranch['event']:
			if self.dic1DVectLoaded['event']:
				self.eventVect = self.dicBraVect1D['event']
		if self.hasBranch['voltageHV']:
			if self.dic1DVectLoaded['voltageHV']:
				self.voltageHV = self.dicBraVect1D['voltageHV']
		if self.hasBranch['currentHV']:
			if self.dic1DVectLoaded['currentHV']:
				self.currentHV = self.dicBraVect1D['currentHV']
		if self.hasBranch['timeHV']:
			key = 'timeHV.Convert()' if 'timeHV.Convert()' in self.branches1DLoad else 'timeHV.AsDouble()'
			if self.dic1DVectLoaded[key]:
				self.timeHV = self.dicBraVect1D[key]
		if self.hasBranch['peakPosition']:
			if self.dic1DVectLoaded['peakPosition']:
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
		peakPosBra = self.analysisTree.Branch('peakPosition', self.peak_position, 'peakPosition/F')
		self.utils.CreateProgressBar(self.in_root_tree.GetEntries())
		self.utils.bar.start()
		for ev in xrange(self.in_root_tree.GetEntries()):
			# self.in_root_tree.GetEntry(ev)
			if ev in self.eventVect:
				try:
					self.peak_position.itemset(self.peak_positions[np.argwhere(ev == self.eventVect).flatten()])
				except ValueError:
					ExitMessage('Could not fill event {ev}; it should have a peak position of: {v}. Exiting'.format(ev=ev, v=self.peak_positions[np.argwhere(ev == self.eventVect).flatten()]), os.EX_DATAERR)
			else:
				self.peak_position.itemset(0)
			# peakPosBra.Fill()
			self.analysisTree.Fill()
			self.utils.bar.update(ev + 1)
		self.analysisFile.cd()
		if not self.analysisTree.GetFriend(self.in_tree_name):
			self.analysisTree.AddFriend(self.in_tree_name, '{d}/{f}.root'.format(d=self.inDir, f=self.in_tree_name))
		self.analysisTree.Write()
		self.utils.bar.finish()

	def CloseAnalysisROOTFile(self):
		if self.analysisFile.IsOpen():
			self.analysisFile.Close()
		if self.analysisTree:
			del self.analysisTree
		self.analysisTree = None
		self.analysisFile = None

	def CloseInputROOTFiles(self):
		if self.in_root_file.IsOpen():
			self.in_root_file.Close()
		if self.in_root_tree:
			del self.in_root_tree
		self.in_root_file = None
		self.in_root_tree = None

	def Reset_Braches_Lists_And_Dictionaries(self):
		self.branches1DTotal = BRANCHES1DTOTAL[:]
		self.branches1DType = BRANCHES1DTYPE.copy()
		self.branchesWavesTotal = BRANCHESWAVESTOTAL[:]
		self.branchesWavesType = BRANCHESWAVESTYPE.copy()
		self.branchesAll = self.branches1DTotal + self.branchesWavesTotal
		self.branches1DLoad = BRANCHES1DLOAD[:]
		self.branchesWavesLoad = BRANCHESWAVESLOAD[:]
		self.analysisScalarsBranches = ANALYSISSCALARBRANCHES[:]
		self.dicBraVect1D = OrderedDict()
		self.dicBraVectWaves = OrderedDict()
		self.hasBranch = {}

	def AddPeakPositionCut(self):
		self.cut0 += ro.TCut('peakPosCut', 'abs(peakPosition-{pp})<={ppc}'.format(pp=self.peakTime, ppc=self.peakPosCut))

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
		print 'Filling tree with scalars...'
		pedBra = self.analysisTree.Branch('pedestal', self.ped, 'pedestal/F')
		pedSigmaBra = self.analysisTree.Branch('pedestalSigma', self.pedSigma, 'pedestalSigma/F')
		pedSignalBra = self.analysisTree.Branch('signalAndPedestal', self.sigAndPed, 'signalAndPedestal/F')
		pedSignalSigmaBra = self.analysisTree.Branch('signalAndPedestalSigma', self.sigAndPed, 'signalAndPedestalSigma/F')
		sigBra = self.analysisTree.Branch('signal', self.sig, 'signal/F')
		self.utils.CreateProgressBar(self.max_events)
		self.utils.bar.start()
		for ev in xrange(self.in_root_tree.GetEntries()):
			# self.in_root_tree.GetEntry(ev)
			self.analysisTree.GetEntry(ev)
			if ev in self.eventVect:
				try:
					argum = np.argwhere(ev == self.eventVect).flatten()
					self.ped.itemset(self.pedVect[argum])
					self.pedSigma.itemset(self.pedSigmaVect[argum])
					self.sigAndPed.itemset(self.sigAndPedVect[argum])
					self.sigAndPedSigma.itemset(self.sigAndPedSigmaVect[argum])
					self.sig.itemset(self.sigVect[argum])
				except ValueError:
					ExitMessage('Could not fill event {ev}. Exiting...'.format(ev=ev), os.EX_DATAERR)
			else:
				self.ped.itemset(0)
				self.pedSigma.itemset(0)
				self.sigAndPed.itemset(0)
				self.sigAndPedSigma.itemset(0)
				self.sig.itemset(0)
			self.analysisTree.Fill()
			# pedBra.Fill()
			# pedSigmaBra.Fill()
			# pedSignalBra.Fill()
			# pedSignalSigmaBra.Fill()
			# sigBra.Fill()
			self.utils.bar.update(ev + 1)
		self.analysisFile.cd()
		if not self.analysisTree.GetFriend(self.in_tree_name):
			self.analysisTree.AddFriend(self.in_tree_name, '{d}/{f}.root'.format(d=self.inDir, f=self.in_tree_name))
		self.analysisTree.Write()
		self.utils.bar.finish()

	def ExtractMeanOfWaveforms(self):
		self.signalWaveMeanVect = self.signalWaveVect.mean(axis=0)
		self.signalWaveSigmaVect = self.signalWaveVect.std(axis=0)

	def DrawHisto(self, name, xmin, xmax, deltax, var, varname, option='e'):
		ro.TFormula.SetMaxima(100000)
		if self.histo[name]:
			self.histo[name].Delete()
			del self.histo[name]
		self.histo[name] = ro.TH1F('h_' + name, 'h_' + name, int(np.floor((xmax - xmin) / deltax + 0.5)), xmin, xmax)
		self.histo[name].GetXaxis().SetTitle(varname)
		self.histo[name].GetYaxis().SetTitle('entries')
		if 'goff' not in option:
			if self.canvas[name]:
				self.canvas[name].Close()
				del self.canvas[name]
			self.canvas[name] = ro.TCanvas('c_' + name, 'c_' + name, 1)
			self.canvas[name].cd()
		self.analysisTree.Draw('{v}>>h_{n}'.format(v=var, n=name), self.cut0, option)
		if 'goff' not in option:
			self.canvas[name].SetGridx()
			self.canvas[name].SetGridy()
			self.canvas[name].SetTicky()
			ro.gPad.Update()
			SetDefault1DStats(self.histo[name])
		ro.TFormula.SetMaxima(1000)

	def DrawHisto2D(self, name, varx, xmin, xmax, deltax, xname, vary, ymin, ymax, deltay, yname, cuts='', option='colz'):
		ro.TFormula.SetMaxima(100000)
		if self.histo[name]:
			self.histo[name].Delete()
			del self.histo[name]
		self.histo[name] = ro.TH2F('h_' + name, 'h_' + name, int(np.floor((xmax - xmin) / deltax + 0.5) + 2), xmin - deltax, xmax + deltax, int(np.floor((ymax - ymin) / deltay + 0.5) + 2), ymin - deltay, ymax + deltay)
		self.histo[name].GetXaxis().SetTitle(xname)
		self.histo[name].GetYaxis().SetTitle(yname)
		self.histo[name].GetZaxis().SetTitle('entries')
		if 'goff' not in option:
			if self.canvas[name]:
				self.canvas[name].Close()
				del self.canvas[name]
			self.canvas[name] = ro.TCanvas('c_' + name, 'c_' + name, 1)
			self.canvas[name].cd()
		cuts0 = self.cut0.GetTitle() if cuts == '' else cuts
		self.analysisTree.Draw('{y}:{x}>>h_{n}'.format(y=vary, x=varx, n=name), cuts0, option)
		if 'goff' not in option:
			self.canvas[name].SetGridx()
			self.canvas[name].SetGridy()
			self.canvas[name].SetTickx()
			self.canvas[name].SetTicky()
			ro.gPad.Update()
			SetDefault2DStats(self.histo[name])
		ro.TFormula.SetMaxima(1000)

	def PlotPeakPositionDistributions(self, name='peakPosDist', low_t=1e-9, up_t=5001e-9, nbins=500):
		self.DrawHisto(name, low_t, up_t, (up_t - low_t) / nbins, 'peakPosition', 'Peak Position [s]', 'e')
		self.histo[name].GetXaxis().SetRangeUser(self.histo[name].GetMean() - 5 * self.histo[name].GetRMS(), self.histo[name].GetMean() + 5 * self.histo[name].GetRMS())
		func = ro.TF1('fit_' + name, 'gaus', low_t, up_t)
		func.SetNpx(1000)
		params = np.array((self.histo[name].GetEntries(), self.histo[name].GetMean(), self.histo[name].GetRMS()), 'float64')
		func.SetParameters(params)
		fit = self.histo[name].Fit('fit_' + name, 'QEMS', '', self.histo[name].GetMean() - self.histo[name].GetRMS(), self.histo[name].GetMean() + self.histo[name].GetRMS())
		SetDefaultFitStats(self.histo[name], func)
		self.peakTime = fit.Parameter(1)
		self.AddPeakPositionCut()

	def PlotSignal

	def PlotWaveforms(self, name='SignalWaveform', type='signal', cuts='', option='colz'):
		var = 'voltageSignal' if 'signal' in type.lower() else 'voltageTrigger' if 'trig' in type.lower() else 'voltageVeto' if 'veto' in type.lower() else ''
		if var == '':
			print 'type should be "signal", "trigger" or "veto"'
			return
		# cuts0 = self.cut0.GetTitle() if cuts == '' else cuts
		vname = ('signal' if var == 'voltageSignal' else 'trigger' if var == 'voltageTrigger' else 'veto' if var == 'voltageVeto' else '') + ' [V]'
		tmin, tmax, deltat = self.analysisTree.GetMinimum('time'), self.analysisTree.GetMaximum('time'), self.settings.time_res
		vmin, vmax, deltav = self.analysisTree.GetMinimum(var), self.analysisTree.GetMaximum(var), (self.signal_ch.offseted_adc_to_volts_cal['p1'] if var == 'voltageSignal' else self.settings.sigRes)
		self.DrawHisto2D(name, 'time', tmin - deltat/2.0, tmax + deltat/2.0, deltat, 'time[s]', var, vmin, vmax, deltav, vname, cuts, option)


	def ResetTreeToOriginal(self, keepBranches=['event','time','voltageSignal','voltageTrigger','voltageVeto','vetoedEvent','badShape','badPedestal','voltageHV','currentHV','timeHV']):
		print 'Restoring tree with the following branches:', keepBranches, '...'
		raw_input('Press a key and Enter to continue: ')
		self.OpenAnalysisROOTFile('READ')
		self.LoadAnalysisTree()
		self.in_root_tree.SetBranchStatus('*', 0)
		for branch in keepBranches:
			if self.TreeHasBranch(branch):
				self.in_root_tree.SetBranchStatus(branch, 1)
		newFile = ro.TFile('{o}/temp.root'.format(o=self.outDir), 'recreate')
		newTree = self.in_root_tree.CloneTree()
		newTree.Print()
		newFile.Write()
		del self.in_root_file
		del newFile
		self.in_root_file = None
		checkFile = ro.TFile('{o}/temp.root'.format(o=self.outDir), 'READ')
		checkTree = checkFile.Get(self.in_tree_name)
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


	def PrintPlotLimits(self, ti=-5.12e-7, tf=4.606e-6, vmin=-0.7, vmax=0.05):
		print np.double([(tf-ti)/float(self.settings.time_res) +1, ti-self.settings.time_res/2.0,
		                 tf+self.settings.time_res/2.0, (vmax-vmin)/self.settings.sigRes, vmin, vmax])

def main():
	parser = OptionParser()
	parser.add_option('-d', '--inDir', dest='inDir', default='.', type='string', help='Directory containing the run files')
	parser.add_option('-c', '--configFile', dest='config', default='CAENAnalysisConfig.cfg', type='string', help='Path to file containing Analysis configuration file')
	parser.add_option('-i', '--inputFile', dest='infile', default='', type='string', help='Path to root file to be analysed')
	parser.add_option('-b', '--bias', dest='bias', default=0, type='float', help='Bias voltage used')
	parser.add_option('-v', '--verbose', dest='verb', default=False, help='Toggles verbose', action='store_true')
	# parser.add_option('-a', '--automatic', dest='auto', default=False, help='Toggles automatic conversion and analysis afterwards', action='store_true')

	(options, args) = parser.parse_args()
	directory = str(options.inDir)
	config = str(options.config)
	infile = str(options.infile)
	bias = float(options.bias)
	verb = bool(options.verb)

	ana = AnalysisCaenCCD(directory, config, infile, bias, verb)

	# ana.LoadAnalysisTree()
	# ana.LoadPickles()
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
