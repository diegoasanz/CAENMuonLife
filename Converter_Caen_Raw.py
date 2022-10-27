#!/usr/bin/env python
import visa
import csv
import numpy as np
from struct import unpack
import time, os, sys
from optparse import OptionParser
import progressbar
import ipdb
# from pykeyboard import PyKeyboard
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

class Converter_Caen_Raw:
	def __init__(self, settings_object_path='', data_path='', simultaneous_data_conv=True):
		self.settings_object_path = settings_object_path
		self.settings = Settings_Caen()
		self.settings_full_path = os.path.abspath(settings_object_path)
		self.output_dir = '/'.join(self.settings_full_path.split('/')[:-1])
		self.raw_dir = self.output_dir if data_path == '' else data_path
		self.filename = self.settings_full_path.split('/')[-1].split('.settings')[0]
		# load pickles
		self.settings = pickle.load(open('{d}/{f}.settings'.format(d=self.output_dir, f=self.filename), 'rb')) if not self.output_dir == '/1bla1' else Settings_Caen()
		self.settings.simultaneous_conversion = simultaneous_data_conv  # overrides the flag used while taking data, if it is converted offline
		self.signal_ch = {sigi: pickle.load(open('{d}/{f}.signal_ch{s}'.format(s=sigi, d=self.output_dir, f=self.filename), 'rb')) if not self.output_dir == '/1bla1' else Channel_Caen() for sigi in xrange(self.settings.num_signals)}
		self.trigger_ch = pickle.load(open('{d}/{f}.trigger_ch'.format(d=self.output_dir, f=self.filename), 'rb')) if not self.output_dir == '/1bla1' else Channel_Caen()
		self.veto_ch = pickle.load(open('{d}/{f}.veto'.format(d=self.output_dir, f=self.filename), 'rb')) if not self.output_dir == '/1bla1' else Channel_Caen()
		# define paths for raw data
		self.signal_path = {sigi: data_path + '/raw_wave{chs}.dat'.format(chs=self.settings.sigCh[sigi]) if self.settings.simultaneous_conversion else data_path + '/' + self.filename + '_signal{s}.dat'.format(s=sigi) for sigi in xrange(self.settings.num_signals)}
		self.trigger_path = data_path + '/raw_wave{cht}.dat'.format(cht=self.settings.trigCh) if self.settings.simultaneous_conversion else data_path + '/' + self.filename + '_trigger.dat'
		self.veto_path = data_path + '/raw_wave{cha}.dat'.format(cha=self.settings.vetoCh) if self.settings.simultaneous_conversion else data_path + '/' + self.filename + '_veto.dat'
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
		self.trigger_search_window = 0.5e-6
		self.veto_window_around_trigg = 0.5e-6
		self.sig_window_around_trigg = 0.5e-6
		self.doVeto = True

		self.array_points = np.arange(self.points, dtype=np.dtype('int32'))

		self.raw_file = None
		self.raw_tree = None
		self.sigBra = {sigi: None for sigi in xrange(self.settings.num_signals)}
		# self.sigADCBra = {sigi: None for sigi in xrange(self.settins.num_signals)}
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

		self.sigm_stop_threshold = 4
		self.sigm_decay_threshold = 4

		self.cond_sig_ped_pos = None
		self.cond_sig_decay_window = None
		self.cond_sig_stop_window = None

		self.bar = None

	def SetupRootFile(self):
		if self.simultaneous_conversion:
			print 'Start creating root file simultaneously with data taking'
		else:
			print 'Start creating root file'
		self.raw_file = ro.TFile('{wd}/{r}.raw.root'.format(wd=self.output_dir, r=self.filename), 'RECREATE')
		self.raw_tree = ro.TTree(self.filename, self.filename)
		self.raw_tree.SetAutoFlush(100)
		self.raw_tree.SetAutoSave(-10485760)
		if self.doVeto:
			self.vetoBra = np.zeros(self.points, 'int16')
			self.vetoedBra = np.zeros(1, '?')
		self.eventBra = np.zeros(1, 'I')
		self.sigBra = {sigi: np.zeros(self.points, 'int16') for sigi in xrange(self.settings.num_signals)}
		self.trigBra = np.zeros(self.points, 'int16')
		self.timeBra = np.zeros(self.points, 'int16')

		self.raw_tree.Branch('event', self.eventBra, 'event/i')
		self.raw_tree.Branch('time', self.timeBra, 'time[{s}]/S'.format(s=self.points))
		for sigi in xrange(self.settings.num_signals):
			self.raw_tree.Branch('voltageSignal{s}'.format(s=sigi), self.sigBra[sigi], 'voltageSignal{si}[{s}]/S'.format(si=sigi, s=self.points))
		self.raw_tree.Branch('voltageTrigger', self.trigBra, 'voltageTrigger[{s}]/S'.format(s=self.points))
		self.raw_tree.Branch('voltageVeto', self.vetoBra, 'voltageVeto[{s}]/S'.format(s=self.points))

	def GetBinariesNumberWrittenEvents(self):
		self.signal_written_events = int(round(os.path.getsize(self.signal_path[0]) / self.struct_len)) if os.path.isfile(self.signal_path[0]) else 0
		self.trigger_written_events = int(round(os.path.getsize(self.trigger_path) / self.struct_len)) if os.path.isfile(self.trigger_path) else 0
		self.anti_co_written_events = int(round(os.path.getsize(self.veto_path) / self.struct_len)) if os.path.isfile(self.veto_path) else 0

	def OpenRawBinaries(self):
		self.fs = {sigi: open(self.signal_path[sigi], 'rb') for sigi in xrange(self.settings.num_signals)}
		self.ft = open(self.trigger_path, 'rb')
		self.fa = open(self.veto_path, 'rb')

	def CreateProgressBar(self, maxVal=1):
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
			self.sigVolts = {sigi: np.floor(0.5+(self.ADC_to_Volts('signal', sigi) / 0.0001)).astype('int16') for sigi in xrange(self.settings.num_signals)}
			self.struct_t = struct.Struct(self.struct_fmt).unpack_from(self.datat)
			self.trigADC = np.array(self.struct_t, 'H')
			self.trigVolts = np.floor(0.5 + (self.ADC_to_Volts('trigger') / 0.0001)).astype('int16')
			# self.LookForTime0()

			self.points = 10240
			self.settings.post_trig_percent = 95
			post_trig_percent = self.settings.post_trig_percent

			self.timeVect = np.linspace(-self.time_res * 1e9 * self.points * (100 - post_trig_percent)/100., self.time_res * 1e9 * (self.points - 1 - (self.points * (100 - post_trig_percent)/100.)), self.points, dtype='f8')
			self.timeVect = np.floor(0.5 + self.timeVect).astype('int16')
			self.struct_ac = struct.Struct(self.struct_fmt).unpack_from(self.dataa)
			self.vetoADC = np.array(self.struct_ac, 'H')
			self.vetoVolts = np.floor(0.5 + self.ADC_to_Volts('veto') / 0.0001).astype('int16')
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

	def GetInfoSignal(self):
		self.cond_sig_ped_pos = np.array((self.array_points - (self.trigPos - self.sig_window_around_trigg / float(self.time_res))).astype('float64') < 0, dtype='?')
		self.cond_sig_stop_window = np.array(np.abs(self.array_points - self.trigPos) <= int(round(self.sig_window_around_trigg / float(self.time_res))), dtype='?')
		self.cond_sig_decay_window = np.array(self.array_points - self.trigPos - self.sig_window_around_trigg / float(self.time_res) > 0, dtype='?')
		sigPedVolts = {sigi: np.extract(self.cond_sig_ped_pos, self.sigVolts[sigi]) for sigi in xrange(self.settings.num_signals)}
		pedMean = {sigi: sigPedVolts[sigi].mean() for sigi in xrange(self.settings.num_signals)}
		pedSigma = {sigi: sigPedVolts[sigi].std() for sigi in xrange(self.settings.num_signals)}
		return {'mean': pedMean, 'sigma': pedSigma}

	def FillBranches(self, ev):
		self.eventBra.fill(ev)
		np.putmask(self.timeBra, np.ones(self.points, '?'), self.timeVect)
		for sigi, sigvolti in self.sigVolts.iteritems():
			np.putmask(self.sigBra[sigi], np.ones(self.points, '?'), sigvolti)
		np.putmask(self.trigBra, np.ones(self.points, '?'), self.trigVolts)
		np.putmask(self.vetoBra, np.ones(self.points, '?'), self.vetoVolts)

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
	# first argument is the path to the settings pickle file
	# second argument is the path of the directory that contains the raw data.
	# By default, it assumes simultaneous data conversion. If the conversion is done offline (aka. not simultaneous), then the 3rd parameter has to be given and should be '0'

	if len(sys.argv) < 2:
		print 'Usage is: Converter_Caen.py <settings_pickle_path> <dir_with_raw_data> 0 for offline conversion)'
		exit()
	settings_object = str(sys.argv[1])  # settings pickle path
	if settings_object in ['-h', '--help']:
		print 'Usage is: Converter_Caen.py <settings_pickle_path> <dir_with_raw_data> 0 for offline conversion)'
		exit()
	print 'settings object', settings_object
	if len(sys.argv) > 2:
		data_path = str(sys.argv[2])  # path where the binary data in adcs is. It is a directory path containing the raw files.
		print 'data_path', data_path
	else:
		data_path = ''
		print 'data_path empty ""'
	is_simultaneous_data_conv = True
	if len(sys.argv) > 3:
		if IsInt(str(sys.argv[3])):
			is_simultaneous_data_conv = bool(int(str(sys.argv[3])))
			print 'simultaneous is now', is_simultaneous_data_conv

	converter = Converter_Caen_Raw(settings_object_path=settings_object, data_path=data_path, simultaneous_data_conv=is_simultaneous_data_conv)

	converter.SetupRootFile()
	converter.GetBinariesNumberWrittenEvents()
	converter.OpenRawBinaries()
	converter.CreateProgressBar(converter.num_events)
	converter.ConvertEvents()
	converter.CloseAll()




