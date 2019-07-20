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
import pickle
import shutil
from Settings_Caen import Settings_Caen
from Channel_Caen import Channel_Caen
from copy import deepcopy
from Utils import *


# from DataAcquisition import DataAcquisition

class Converter_Caen:
	def __init__(self, settings_object_path='', data_path=''):
		self.settings_object_path = settings_object_path
		self.settings = Settings_Caen()
		self.settings_full_path = os.path.abspath(settings_object_path)
		self.output_dir = '/'.join(self.settings_full_path.split('/')[:-1])
		self.raw_dir = self.output_dir if data_path == '' else data_path
		self.filename = self.settings_full_path.split('/')[-1].split('.settings')[0]
		# load pickles
		self.settings = pickle.load(open('{d}/{f}.settings'.format(d=self.output_dir, f=self.filename), 'rb')) if not self.output_dir == '/1bla1' else Settings_Caen()
		self.signal_ch = {sigi: pickle.load(open('{d}/{f}.signal_ch{s}'.format(s=sigi, d=self.output_dir, f=self.filename), 'rb')) if not self.output_dir == '/1bla1' else Channel_Caen() for sigi in xrange(self.settings.num_signals)}
		self.trigger_ch = pickle.load(open('{d}/{f}.trigger_ch'.format(d=self.output_dir, f=self.filename), 'rb')) if not self.output_dir == '/1bla1' else Channel_Caen()
		self.veto_ch = pickle.load(open('{d}/{f}.veto'.format(d=self.output_dir, f=self.filename), 'rb')) if not self.output_dir == '/1bla1' else Channel_Caen()
		# define paths for raw data
		self.signal_path = {sigi: data_path + '/raw_wave{chs}.dat'.format(chs=self.settings.sigCh[sigi]) if self.settings.simultaneous_conversion else data_path + '/' + self.filename + '_signal{s}.dat'.format(s=sigi) for sigi in xrange(self.settings.num_signals)}
		self.trigger_path = data_path + '/raw_wave{cht}.dat'.format(cht=self.settings.trigCh) if self.settings.simultaneous_conversion else data_path + '/' + self.filename + '_trigger.dat'
		self.veto_path = data_path + '/raw_wave{cha}.dat'.format(cha=self.settings.acCh) if self.settings.simultaneous_conversion else data_path + '/' + self.filename + '_veto.dat'
		# get settings variables
		self.points = self.settings.points
		self.num_events = self.settings.num_events
		self.struct_len = self.settings.struct_len
		self.struct_fmt = self.settings.struct_fmt
		self.adc_res = self.settings.sigRes
		self.adc_offset = 0
		self.sig_offset = {sigi: self.signal_ch[sigi].dc_offset_percent for sigi in xrange(self.settings.num_signals)}
		self.trig_offset = self.trigger_ch.dc_offset_percent
		self.anti_co_offset = self.veto_ch.dc_offset_percent
		self.time_res = self.settings.time_res
		self.post_trig_percent = self.settings.post_trig_percent
		self.trig_value = self.trigger_ch.thr_volts
		# self.veto_value = self.veto_ch.thr_volts
		self.veto_value = self.veto_ch.Volts_to_ADC(self.veto_ch.thr_volts)
		self.dig_bits = self.settings.dig_bits
		self.simultaneous_conversion = self.settings.simultaneous_conversion
		self.time_recal = self.settings.time_calib
		self.sig_polarity = {sigi: self.signal_ch[sigi].polarity for sigi in xrange(self.settings.num_signals)}
		# define trigger and veto windows
		self.trigger_search_window = 0.1e-6
		self.veto_window_around_trigg = 0.1e-6
		self.doVeto = True

		self.array_points = np.arange(self.points, dtype=np.dtype('int32'))

		self.raw_file = None
		self.raw_tree = None
		self.sigBra = {sigi: None for sigi in xrange(self.settings.num_signals)}
		self.eventBra = self.trigBra = self.vetoBra = self.timeBra = self.vetoedBra = self.badShapeBra = self.badPedBra = None
		# self.hourBra = self.minuteBra = self.secondBra = None
		self.hourMinSecBra = None

		self.t0 = time.time()

		self.signal_written_events = self.trigger_written_events = self.anti_co_written_events = None
		self.fs = self.ft = self.fa = None
		self.fs = {sigi: None for sigi in xrange(self.settings.num_signals)}
		self.wait_for_data = None

		self.datas = self.datat = self.dataa = None
		self.datas = {sigi: None for sigi in xrange(self.settings.num_signals)}
		self.sigADC = self.trigADC = self.vetoADC = None
		self.sigADC = {sigi: None for sigi in xrange(self.settings.num_signals)}
		self.sigVolts = self.trigVolts = self.vetoVolts = None
		self.sigVolts = {sigi: None for sigi in xrange(self.settings.num_signals)}
		self.trigPos = None
		self.timeVect = None
		self.vetoed_event = None
		self.condition_base_line = None
		self.condition_peak_pos = None
		self.hour_event = self.minute_event = self.second_event = None
		# self.time_break = int(np.ceil(self.time_recal + 30))
		self.time_break = int(1800)
		# self.hour_min_sec_event = None
		# self.currentTime = None

		self.struct_s = self.struct_t = self.struct_ac = None
		self.struct_s = {sigi: None for sigi in xrange(self.settings.num_signals)}

		self.bar = None

	def SetupRootFile(self):
		if self.simultaneous_conversion:
			print 'Start creating root file simultaneously with data taking'
		else:
			print 'Start creating root file'
		self.raw_file = ro.TFile('{wd}/{r}.root'.format(wd=self.output_dir, r=self.filename), 'RECREATE')
		self.raw_tree = ro.TTree(self.filename, self.filename)
		self.raw_tree.SetAutoFlush(100)
		self.raw_tree.SetAutoSave(-10485760)
		if self.doVeto:
			self.vetoBra = np.zeros(self.points, 'f8')
			self.vetoedBra = np.zeros(1, '?')
		self.eventBra = np.zeros(1, 'I')
		self.sigBra = {sigi: np.zeros(self.points, 'f8') for sigi in xrange(self.settings.num_signals)}
		self.trigBra = np.zeros(self.points, 'f8')
		self.timeBra = np.zeros(self.points, 'f8')
		# self.badShapeBra = np.zeros(1, dtype=np.dtype('int8'))  # signed char
		# self.badPedBra = np.zeros(1, '?')
		self.raw_tree.Branch('event', self.eventBra, 'event/i')
		self.raw_tree.Branch('time', self.timeBra, 'time[{s}]/D'.format(s=self.points))
		for sigi in xrange(self.settings.num_signals):
			self.raw_tree.Branch('voltageSignal{s}'.format(s=sigi), self.sigBra[sigi], 'voltageSignal{si}[{s}]/D'.format(si=sigi, s=self.points))
		self.raw_tree.Branch('voltageTrigger', self.trigBra, 'voltageTrigger[{s}]/D'.format(s=self.points))
		self.raw_tree.Branch('voltageVeto', self.vetoBra, 'voltageVeto[{s}]/D'.format(s=self.points))
		self.raw_tree.Branch('vetoedEvent', self.vetoedBra, 'vetoedEvent/O')
		# self.raw_tree.Branch('badShape', self.badShapeBra, 'badShape/B')  # signed char
		# self.raw_tree.Branch('badPedestal', self.badPedBra, 'badPedestal/O')

	def GetBinariesNumberWrittenEvents(self):
		self.signal_written_events = int(round(os.path.getsize(self.signal_path[0]) / self.struct_len)) if os.path.isfile(self.signal_path[0]) else 0
		self.trigger_written_events = int(round(os.path.getsize(self.trigger_path) / self.struct_len)) if os.path.isfile(self.trigger_path) else 0
		self.anti_co_written_events = int(round(os.path.getsize(self.veto_path) / self.struct_len)) if os.path.isfile(self.veto_path) else 0

	def OpenRawBinaries(self):
		self.fs = {sigi: open(self.signal_path[sigi], 'rb') for sigi in xrange(self.settings.num_signals)}
		self.ft = open(self.trigger_path, 'rb')
		self.fa = open(self.veto_path, 'rb')

	def CreateProgressBar(self, maxVal=1):
		# widgets = [
		# 	'Processed: ', progressbar.Counter(),
		# 	' out of {mv} '.format(mv=maxVal), progressbar.Percentage(),
		# 	' ', progressbar.Bar(marker='>'),
		# 	' ', progressbar.Timer(),
		# 	' ', progressbar.ETA()
		# 	# ' ', progressbar.AdaptativeETA(),
		# 	#  ' ', progressbar.AdaptativeTransferSpeed()
		# ]
		# self.bar = progressbar.ProgressBar(widgets=widgets, maxval=maxVal)
		self.bar = CreateProgressBarUtils(maxVal)

	def ConvertEvents(self):
		self.bar.start()
		for ev in xrange(self.num_events):
			self.CheckFilesSizes(ev)
			self.WaitForData(ev)
			self.ReadData(ev)
			self.CheckData()
			self.struct_s = {sigi: struct.Struct(self.struct_fmt).unpack_from(datasi) for sigi, datasi in self.datas.iteritems()}
			self.sigADC = {sigi: np.array(struct_si, 'H') for sigi, struct_si in self.struct_s.iteritems()}
			self.sigVolts = {sigi: self.ADC_to_Volts('signal', sigi) for sigi in xrange(self.settings.num_signals)}
			self.struct_t = struct.Struct(self.struct_fmt).unpack_from(self.datat)
			self.trigADC = np.array(self.struct_t, 'H')
			self.trigVolts = self.ADC_to_Volts('trigger')
			self.LookForTime0()
			self.timeVect = np.linspace(-self.trigPos * self.time_res, self.time_res * (self.points - 1 - self.trigPos), self.points, dtype='f8')
			self.struct_ac = struct.Struct(self.struct_fmt).unpack_from(self.dataa)
			self.vetoADC = np.array(self.struct_ac, 'H')
			self.vetoVolts = self.ADC_to_Volts('veto')
			self.vetoed_event = self.IsEventVetoed()
			# self.DefineSignalBaseLineAndPeakPosition()
			# self.bad_shape_event = self.IsEventBadShape()
			# self.bad_pedstal_event = self.IsPedestalBad()
			self.FillBranches(ev)
			if ev == 10:
				self.raw_tree.OptimizeBaskets()
			self.raw_tree.Fill()
			self.bar.update(ev + 1)

	def CheckFilesSizes(self, ev):
		self.wait_for_data = (self.signal_written_events <= ev) or (self.trigger_written_events <= ev) or (self.anti_co_written_events <= ev)

	def WaitForData(self, ev):
		t1 = time.time()
		while self.wait_for_data:
			if self.simultaneous_conversion:
				if time.time() - t1 > self.time_break:
					ExitMessage('No data has been saved in file for event {ev} in the past {t} seconds... exiting!'.format(ev=ev, t=self.time_break), os.EX_NOINPUT)
				for fsi in self.fs.itervalues():
					if not fsi.closed:
						fsi.close()
				if not self.ft.closed:
					self.ft.close()
				if not self.fa.closed:
					self.fa.close()
				self.GetBinariesNumberWrittenEvents()
				self.CheckFilesSizes(ev)
				if not self.wait_for_data:
					self.OpenRawBinaries()
			else:
				ExitMessage('The data is corrupted... exiting', os.EX_DATAERR)

	def ReadData(self, ev):
		for sigi, fsi in self.fs.iteritems():
			fsi.seek(ev * self.struct_len, 0)
			self.datas[sigi] = fsi.read(self.struct_len)
		self.ft.seek(ev * self.struct_len, 0)
		self.datat = self.ft.read(self.struct_len)
		self.fa.seek(ev * self.struct_len, 0)
		self.dataa = self.fa.read(self.struct_len)

	def CheckData(self):
		if not self.datas[0] or not self.datat or not self.dataa:
			ExitMessage('No event in signal or trigger files... exiting', os.EX_DATAERR)

	def LookForTime0(self):
		guess_pos = int(round(self.points * (100.0 - self.post_trig_percent)/100.0))
		condition_trigg = np.array(np.abs(self.array_points - guess_pos) <= int(round(self.trigger_search_window/self.time_res)), dtype='?')
		condition_no_trigg = np.array(1 - condition_trigg, dtype='?')
		# mean = np.extract(condition_no_trigg, self.trigVolts).mean()
		# sigma = np.extract(condition_no_trigg, self.trigVolts).std()
		temp_trig_volts = np.copy(self.trigVolts)
		np.putmask(temp_trig_volts, condition_no_trigg, 100)
		volt_min_pos = temp_trig_volts.argmin()
		condition_trigg = np.bitwise_and(condition_trigg, np.array(self.array_points <= volt_min_pos))
		np.putmask(temp_trig_volts, np.bitwise_not(condition_trigg), 100)
		self.trigPos = np.abs(temp_trig_volts - self.trig_value).argmin()
		# self.trigPos = guess_pos
		del guess_pos, condition_trigg, condition_no_trigg, temp_trig_volts, volt_min_pos

	def IsEventVetoed(self):
		condition_veto_base_line = np.array(np.abs(self.array_points - self.trigPos) > int(round(self.veto_window_around_trigg / float(self.time_res))), dtype='?')
		condition_search = np.array(1 - condition_veto_base_line, dtype='?')
		# meanbl = np.extract(condition_veto_base_line, self.vetoADC).mean()
		# veto_event = bool((np.extract(condition_search, self.vetoADC) - meanbl + self.veto_value).min() <= 0)
		veto_event = bool((np.extract(condition_search, self.vetoADC) - self.veto_value).min() <= 0)
		del condition_search, condition_veto_base_line #, meanbl
		return veto_event

	def DefineSignalBaseLineAndPeakPosition(self):
		self.condition_base_line = np.array(self.array_points <= self.trigPos, dtype='?')
		# # values 2.2 and 0.5 us come from shape of signal
		# self.condition_peak_pos = np.array(np.abs(self.array_points - (self.peak_pos_estimate/float(self.time_res) + self.trigPos)) <= self.peak_pos_window/float(self.time_res), dtype='?')

	def IsEventBadShape(self):
		# mean = np.extract(self.condition_base_line, self.sigADC).mean()
		sigma = np.extract(self.condition_base_line, self.sigADC).std()
		# lim_inf = self.condition_peak_pos.argmax()
		# lim_sup = self.points - self.condition_peak_pos[::-1].argmax() - 1
		# peak_pos = self.sigADC.argmin() if self.polarity == 1 else self.sigADC.argmax()
		# if lim_inf < peak_pos < lim_sup:
		# 	# The event has a good shape
		# 	return 0
		# else:
		# 	modified_adc = self.sigADC - sigma if self.polarity == 1 else self.sigADC + sigma
		# 	modified_adc[lim_inf] += 2*sigma if self.polarity == 1 else -2*sigma
		# 	modified_adc[lim_sup] += 2*sigma if self.polarity == 1 else -2*sigma
		# 	peak_pos = modified_adc.argmin() if self.polarity == 1 else modified_adc.argmax()
		# 	if lim_inf < peak_pos < lim_sup:
		# 		# Can't tell if the event has a bad shape
		# 		return -1
		# 	else:
		# 		# Event has bad shape
		# 		return 1

	def IsPedestalBad(self):
		sigma = np.extract(self.condition_base_line, self.sigADC).std()
		mean = np.extract(self.condition_base_line, self.sigADC).mean()
		self.adc_res = self.signal_ch.offseted_adc_to_volts_cal['p1']
		sigma_volts = sigma * self.adc_res
		mean_volts = mean * self.adc_res
		diff_volts = abs(int(self.sigADC[0]) - int(self.sigADC[self.trigPos])) * self.adc_res
		if sigma_volts >= 2e-3 or diff_volts >= 15e-3:  # if a signal is not flat enough due to previous unstable states
			return True
		else:
			last_sig = self.sigADC[-1] * self.adc_res
			if abs(last_sig - mean_volts) > 75e-3:  # if a signal does not return to base line due to multiple effects
				return True
			else:  # when the signal pedestal is well behaved and the signal returns to baseline
				return False

	def FillBranches(self, ev):
		self.eventBra.fill(ev)
		np.putmask(self.timeBra, np.ones(self.points, '?'), self.timeVect)
		for sigi, sigvolti in self.sigVolts.iteritems():
			np.putmask(self.sigBra[sigi], np.ones(self.points, '?'), sigvolti)
		# np.putmask(self.sigBra, np.bitwise_not(np.zeros(self.points, '?')), self.sigVolts)
		np.putmask(self.trigBra, np.ones(self.points, '?'), self.trigVolts)
		np.putmask(self.vetoBra, np.ones(self.points, '?'), self.vetoVolts)
		self.vetoedBra.fill(self.vetoed_event)
		# self.badShapeBra.fill(self.bad_shape_event)
		# self.badPedBra.fill(self.bad_pedstal_event)

	def CloseAll(self):
		self.bar.finish()
		self.raw_file.Write()
		self.raw_file.Close()
		for fsi in self.fs.itervalues():
			fsi.close()
			del fsi
		self.ft.close()
		del self.ft
		if self.doVeto:
			self.fa.close()
			del self.fa
		self.t0 = time.time() - self.t0
		print 'Time creating root tree:', self.t0, 'seconds'
		exit()

	def IsPedestalNotFlat(self, signalADC, points, trigPos, time_res):
		array_points = np.arange(points, dtype=np.dtype('int32'))
		condition_base_line = np.array(array_points - trigPos <= 0, dtype='?')

	def ADC_to_Volts(self, sig_type, signum=0):
		adcs = 0
		if sig_type == 'signal':
			adcs = self.sigADC[signum]
			result = self.signal_ch[signum].ADC_to_Volts(adcs)
		elif sig_type == 'trigger':
			adcs = self.trigADC
			result = self.trigger_ch.ADC_to_Volts(adcs)
		elif sig_type == 'veto':
			adcs = self.vetoADC
			result = self.veto_ch.ADC_to_Volts(adcs)
		else:
			ExitMessage('Wrong type. Exiting', os.EX_SOFTWARE)
		return result

if __name__ == '__main__':
	settings_object_path = str(sys.argv[1])  # settings pickle path
	if len(sys.argv) > 2:
		data_path = str(sys.argv[2])  # path where the binary data in adcs is
	else:
		data_path = ''
	converter = Converter_Caen(settings_object_path=settings_object_path, data_path=data_path)

	converter.SetupRootFile()
	converter.GetBinariesNumberWrittenEvents()
	converter.OpenRawBinaries()
	converter.CreateProgressBar(converter.num_events)
	converter.ConvertEvents()
	converter.CloseAll()




